import os
import re
import time
import glob
import shutil
import subprocess
from typing import Tuple, Optional, Dict

import pandas as pd
import requests
import certifi

# =========================
# PARAMÈTRES (dossier projet)
# =========================
ROOT = "blast_local"

PROTEOMES_DIR = os.path.join(ROOT, "proteomes")
DB_DIR = os.path.join(ROOT, "db")
QUERIES_FASTA = os.path.join(ROOT, "queries", "nos_prots.fasta")
RESULTS_DIR = os.path.join(ROOT, "results")

TABLE_PROT = os.path.join(ROOT, "queries", "proteines_connues.xlsx")
UNIPROT_COL = "ID Uniprot"
NAME_COL = "Nom de la protéine"
ORG_COL = "Organisme"  # optionnel si présent

# =========================
# BLAST params
# =========================
EVALUE_MAX = 1e-5  # défaut (vertébrés)

# >>> Invertébrés = 1e-3 (selon TES dossiers)
EVALUE_BY_SPECIES = {
    "daphnia_magna": 1e-3,
    "dreissena_polymorpha": 1e-3,
}

MAX_TARGET_SEQS = 200
NUM_THREADS = 8

FILTER_MEMBRANE = True
ADD_AF_SCORE = True
AF_SLEEP_SEC = 0.2

FORCE_REBUILD_DB = False


# =========================
# UTILITAIRES
# =========================
def ensure_dir(path: str) -> None:
    """Crée un dossier et tous ses parents si besoin. Ne plante pas si le 
    dossier existe déjà."""
    os.makedirs(path, exist_ok=True)


def file_nonempty(path: str) -> bool:
    """Retourn true si le fichier fasta existe et qu'il n'est pas vide."""
    return os.path.isfile(path) and os.path.getsize(path) > 0


def species_key_from_filename(path: str) -> str:
    """Transforme un nom de fichier de proteome en cle d'espece propre.
    Exemple : 'blast_local/proteomes/danio_rerio.fasta' --> 'danio_rerio'"""
    base = os.path.basename(path)
    for ext in [".fasta", ".fa", ".faa", ".fna", ".fas"]:
        if base.endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


def check_tool_exists(tool: str) -> None:
    """Fonction qui verifie que blastp et makeblastdb sont accessibles dans le 
    PATH."""
    if shutil.which(tool) is None:
        raise RuntimeError(
            f"Commande '{tool}' introuvable. "
            f"Tu es bien dans ton environnement pixi ? (pixi shell)"
        )


def evalue_for_species(species_key: str) -> float:
    """Retourne la evalue spécifique de l'orga, sinon evalue par defaut."""
    return EVALUE_BY_SPECIES.get(species_key, EVALUE_MAX)


# =========================
# PARSING UniProt / BLAST
# =========================
def parse_uniprot_acc(id_str: str) -> str:
    """Extrait l'accession Uniprot propre a partir de la premiere ligne 
    fasta."""
    s = str(id_str)
    parts = s.split("|")
    if len(parts) >= 2:
        return parts[1].split(".")[0]
    return s


def parse_subject_desc(desc: str) -> Tuple[str, Optional[str], Optional[str]]:
    """A partir de stitle de BLAST, extraire nom de la proteine, gene et
    organisme."""
    if desc is None:
        return "", None, None
    text = str(desc)
    prot = text.split(" OS=", 1)[0].strip() if " OS=" in text else text.strip()
    m_gn = re.search(r"\bGN=([^ =]+)", text)
    gene = m_gn.group(1) if m_gn else None
    m_os = re.search(r"\bOS=([^=]+?)(?: OX=| GN=| PE=| SV=|$)", text)
    organism = m_os.group(1).strip() if m_os else None
    return prot, gene, organism


# =========================
# UniProt JSON + filtre membrane (cache)
# =========================
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{accession}.json"
_uniprot_json_cache: Dict[str, Optional[dict]] = {}
_membrane_cache: Dict[str, bool] = {}


def get_uniprot_json(accession: str) -> Optional[dict]:
    accession = str(accession).strip()
    if not accession:
        return None

    if accession in _uniprot_json_cache:
        return _uniprot_json_cache[accession]

    url = UNIPROT_JSON.format(accession=accession)
    try:
        r = requests.get(url, timeout=15, verify=certifi.where())
    except Exception:
        _uniprot_json_cache[accession] = None
        return None

    if r.status_code != 200:
        _uniprot_json_cache[accession] = None
        return None

    try:
        data = r.json()
    except Exception:
        _uniprot_json_cache[accession] = None
        return None

    _uniprot_json_cache[accession] = data
    return data


def is_membrane_uniprot(accession: str) -> bool:
    accession = str(accession).strip()
    if not accession:
        return False

    if accession in _membrane_cache:
        return _membrane_cache[accession]

    data = get_uniprot_json(accession)
    if not data:
        _membrane_cache[accession] = False
        return False

    def has_cell_membrane_term(text: str) -> bool:
        t = (text or "").lower()
        if "membran" not in t:
            return False
        strong = [
            "plasma membrane",
            "cell membrane",
            "cell surface",
            "basolateral plasma membrane",
            "apical plasma membrane",
            "basal plasma membrane",
        ]
        if any(k in t for k in strong):
            return True
        return "cell" in t

    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            for p in ref.get("properties", []):
                if p.get("key") == "GoTerm":
                    if has_cell_membrane_term(p.get("value", "")):
                        _membrane_cache[accession] = True
                        return True

    for com in data.get("comments", []):
        if com.get("commentType") == "SUBCELLULAR_LOCATION":
            for t in com.get("texts", []):
                if has_cell_membrane_term(t.get("value", "")):
                    _membrane_cache[accession] = True
                    return True
            for loc in com.get("subcellularLocations", []):
                val = loc.get("location", {}).get("value", "")
                if has_cell_membrane_term(val):
                    _membrane_cache[accession] = True
                    return True

    _membrane_cache[accession] = False
    return False


# =========================
# AlphaFold score_AF (cache)
# =========================
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
_score_af_cache: Dict[str, Optional[float]] = {}


def normalize_to_uniprot(raw: str) -> str:
    s = str(raw).strip().replace(" ", "")
    m = re.search(r"(?:AF-)?([A-Z0-9]{6,10})(?:-F\d+)?(?:-model_v\d+)?", s)
    return m.group(1) if m else s


def get_af_pdb_url(uniprot_id: str) -> Optional[str]:
    url = AF_API.format(uniprot_id=uniprot_id)
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return None
    data = r.json()
    if isinstance(data, list) and data:
        return data[0].get("pdbUrl")
    return None


def mean_plddt_from_pdb(pdb_text: str) -> Optional[float]:
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                vals.append(float(line[60:66]))  # B-factor = pLDDT
            except ValueError:
                pass
    if not vals:
        return None
    return round(sum(vals) / len(vals), 2)


def get_mean_plddt(uniprot_id: str, sleep_sec: float = 0.2) -> Optional[float]:
    uniprot_id = normalize_to_uniprot(uniprot_id)
    if not uniprot_id:
        return None
    if uniprot_id in _score_af_cache:
        return _score_af_cache[uniprot_id]

    pdb_url = get_af_pdb_url(uniprot_id)
    if not pdb_url:
        _score_af_cache[uniprot_id] = None
        return None

    r = requests.get(pdb_url, timeout=60)
    if r.status_code != 200:
        _score_af_cache[uniprot_id] = None
        return None

    score = mean_plddt_from_pdb(r.text)
    _score_af_cache[uniprot_id] = score
    time.sleep(sleep_sec)
    return score


# =========================
# BLAST DB + run
# =========================
def build_blast_db(proteome_fasta: str, db_prefix: str) -> None:
    expected = [db_prefix + ext for ext in [".pin", ".psq", ".phr"]]
    db_exists = all(os.path.exists(p) for p in expected)

    if db_exists and not FORCE_REBUILD_DB:
        print(f"[INFO] DB déjà présente: {db_prefix}")
        return

    print(f"[INFO] Construction DB: {db_prefix}")
    cmd = ["makeblastdb", "-in", proteome_fasta, "-dbtype",
           "prot", "-parse_seqids", "-out", db_prefix]
    subprocess.run(cmd, check=True)


def run_blastp(query_fasta: str, db_prefix: str, out_tsv: str, species_key: str) -> float:
    cutoff = evalue_for_species(species_key)
    cmd = [
        "blastp",
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", out_tsv,
        "-evalue", str(cutoff),
        "-max_target_seqs", str(MAX_TARGET_SEQS),
        "-num_threads", str(NUM_THREADS),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle ssciname",
    ]
    subprocess.run(cmd, check=True)
    return cutoff


def load_query_metadata(table_xlsx: str) -> Tuple[Dict[str, str], Dict[str, Optional[str]]]:
    df_prots = pd.read_excel(table_xlsx)
    df_prots[UNIPROT_COL] = df_prots[UNIPROT_COL].astype(str).str.strip()

    name_map = dict(zip(df_prots[UNIPROT_COL], df_prots.get(
        NAME_COL, df_prots[UNIPROT_COL])))

    org_map = {}
    if ORG_COL in df_prots.columns:
        org_map = dict(zip(df_prots[UNIPROT_COL], df_prots[ORG_COL]))
    return name_map, org_map


def postprocess_to_excels(
    tsv_path: str,
    out_dir: str,
    species_key: str,
    name_map: Dict[str, str],
    org_map: Dict[str, Optional[str]],
    evalue_cutoff: float,
) -> None:
    cols = ["query", "subject", "identity", "length", "evalue",
            "bitscore", "subject_raw", "subject_ssciname"]
    df = pd.read_csv(tsv_path, sep="\t", names=cols)

    df["db"] = df["subject"].astype(str).str.slice(0, 2)
    df["reviewed"] = df["db"].eq("sp")
    df["query_uniprot"] = df["query"].apply(parse_uniprot_acc)
    df["subject_uniprot"] = df["subject"].apply(parse_uniprot_acc)

    parsed = df["subject_raw"].apply(parse_subject_desc)
    df["subject_name"] = parsed.apply(lambda x: x[0])
    df["subject_gene"] = parsed.apply(lambda x: x[1])
    df["subject_organism_os"] = parsed.apply(lambda x: x[2])

    df["subject_organism"] = df["subject_ssciname"].astype("object")
    mask_empty_org = df["subject_organism"].isna() | (
        df["subject_organism"] == "")
    df.loc[mask_empty_org, "subject_organism"] = df.loc[mask_empty_org,
                                                        "subject_organism_os"]

    df["query_name"] = df["query_uniprot"].map(
        name_map).fillna(df["query_uniprot"])
    df["query_organism"] = df["query_uniprot"].map(
        org_map) if org_map else None

    df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df = df[(df["evalue"] >= 0) & (df["evalue"] <=
                                   evalue_cutoff)].reset_index(drop=True)

    print(f"[INFO] {species_key}: hits après filtre e-value "
          f"(<= {evalue_cutoff:g}): {len(df)}")

    if FILTER_MEMBRANE:
        print(f"[INFO] {species_key}: filtre membrane via UniProt…")
        df["is_membrane"] = df["subject_uniprot"].astype(
            str).apply(is_membrane_uniprot)
        df = df[df["is_membrane"]].reset_index(drop=True)
        print(f"[INFO] {species_key}: hits après filtre membrane: {len(df)}")

    if ADD_AF_SCORE:
        unique_subjects = df["subject_uniprot"].dropna().astype(
            str).drop_duplicates().tolist()
        print(f"[INFO] {species_key}: AlphaFold score sur "
              f"{len(unique_subjects)} sujets uniques…")
        score_map = {}
        for i, acc in enumerate(unique_subjects, start=1):
            print(f"  [AF] {species_key} ({i}/{len(unique_subjects)}) {acc}")
            score_map[acc] = get_mean_plddt(acc, sleep_sec=AF_SLEEP_SEC)
        df["score_AF"] = df["subject_uniprot"].astype(str).map(score_map)
    else:
        df["score_AF"] = None

    df["review_status"] = df["reviewed"].map(
        {True: "Swiss-Prot (reviewed)", False: "TrEMBL (unreviewed)"}
    )

    ensure_dir(out_dir)
    out_xlsx = os.path.join(out_dir, f"{species_key}_blast_detailed.xlsx")
    out_unique = os.path.join(out_dir, f"{species_key}_blast_unique.xlsx")

    cols_final = [
        "query_name", "query_uniprot", "query_organism",
        "subject_uniprot", "subject_name", "subject_gene", "subject_organism",
        "identity", "length", "evalue", "bitscore",
        "review_status",
        "score_AF",
    ]
    df[cols_final].to_excel(out_xlsx, index=False)

    unique_df = df.drop_duplicates(subset="subject_uniprot").copy()
    cols_unique = ["subject_uniprot", "subject_name", "subject_gene",
                   "subject_organism", "review_status", "score_AF"]
    unique_df[cols_unique].to_excel(out_unique, index=False)

    print(f"[OK] {species_key}: {len(unique_df)} sujets uniques exportés.")
    print(f"     - Detailed: {out_xlsx}")
    print(f"     - Unique  : {out_unique}")


# =========================
# MAIN LOOP
# =========================
def main():
    check_tool_exists("makeblastdb")
    check_tool_exists("blastp")

    if not file_nonempty(QUERIES_FASTA):
        raise FileNotFoundError(f"Queries introuvables: {QUERIES_FASTA}")
    if not file_nonempty(TABLE_PROT):
        raise FileNotFoundError(f"Excel queries introuvable: {TABLE_PROT}")

    ensure_dir(DB_DIR)
    ensure_dir(RESULTS_DIR)

    name_map, org_map = load_query_metadata(TABLE_PROT)

    proteome_files = sorted(
        glob.glob(os.path.join(PROTEOMES_DIR, "*.fasta"))
        + glob.glob(os.path.join(PROTEOMES_DIR, "*.faa"))
        + glob.glob(os.path.join(PROTEOMES_DIR, "*.fa"))
    )
    if not proteome_files:
        raise FileNotFoundError(f"Aucun proteome trouvé dans {PROTEOMES_DIR}")

    print(f"[INFO] {len(proteome_files)} protéomes détectés.")

    for proteome in proteome_files:
        sp = species_key_from_filename(proteome)
        print("\n" + "=" * 80)
        print(f"[INFO] Traitement: {sp}")
        print(f"[INFO] Proteome: {proteome}")
        print(f"[INFO] E-value utilisée: {evalue_for_species(sp):g}")

        db_prefix = os.path.join(DB_DIR, sp)
        build_blast_db(proteome, db_prefix)

        species_out_dir = os.path.join(RESULTS_DIR, sp)
        ensure_dir(species_out_dir)
        tsv_out = os.path.join(species_out_dir, f"{sp}_blast.tsv")

        cutoff_used = run_blastp(QUERIES_FASTA, db_prefix, tsv_out, sp)

        postprocess_to_excels(
            tsv_path=tsv_out,
            out_dir=species_out_dir,
            species_key=sp,
            name_map=name_map,
            org_map=org_map,
            evalue_cutoff=cutoff_used,
        )

    print("\n[OK] Terminé : tous les protéomes ont été traités.")


if __name__ == "__main__":
    main()
