import os
import re
import gzip
import requests
import pandas as pd

EXCEL_FILE = "1_proteine_test.xlsx"
SHEET_NAME = 0

ORG_COLS = {
    "Poisson zebre": "zebrafish",
    "Xenope": "xenope",
    "Rat": "rattus",
    "Dreissna": "moule",
    "Daphnia": "daphnia",
}

# On tente d'abord l'API, puis des URLs directes.
AF_DIRECT_PATTERNS = [
    # EBI files (souvent le plus stable)
    "https://alphafold.ebi.ac.uk/files/AF-{u}-F1-model_v{v}.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{u}-F1-model_v{v}.pdb.gz",
]

AF_VERSIONS_TO_TRY = [6, 5, 4, 3, 2]


def get_uniprot(cell):
    if pd.isna(cell):
        return None
    txt = str(cell).strip()
    m = re.search(r"\b([A-NR-Z0-9][A-Z0-9]{5}([A-Z0-9]{4})?)\b", txt)
    return m.group(1) if m else None


def uniprot_accessions(acc, session):
    """primary + secondary (pour les A0A... qui pointent parfois vers une primaire différente)"""
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    try:
        r = session.get(url, timeout=30)
        if r.status_code != 200:
            return [acc]
        js = r.json()
        primary = js.get("primaryAccession")
        secondary = js.get("secondaryAccessions") or []
        out = []
        if acc:
            out.append(acc)
        if primary and primary not in out:
            out.append(primary)
        for s in secondary:
            if s and s not in out:
                out.append(s)
        return out
    except Exception:
        return [acc]


def alphafold_api_pdb_url(acc, session):
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
    r = session.get(api, timeout=30)
    if r.status_code != 200:
        return None
    data = r.json()
    if not data:
        return None
    # Quand dispo, c'est l'URL officielle du PDB
    return data[0].get("pdbUrl")


def download_pdb_from_url(url, out_pdb, session):
    r = session.get(url, timeout=120)
    if r.status_code != 200 or not r.content or len(r.content) < 1000:
        return False

    if url.endswith(".gz"):
        try:
            pdb_bytes = gzip.decompress(r.content)
        except Exception:
            return False
        if len(pdb_bytes) < 1000:
            return False
        with open(out_pdb, "wb") as f:
            f.write(pdb_bytes)
        return True

    with open(out_pdb, "wb") as f:
        f.write(r.content)
    return True


def try_direct_files(acc, out_pdb, session):
    for v in AF_VERSIONS_TO_TRY:
        for pat in AF_DIRECT_PATTERNS:
            url = pat.format(u=acc, v=v)
            if download_pdb_from_url(url, out_pdb, session):
                return url
    return None


df = pd.read_excel(EXCEL_FILE, sheet_name=SHEET_NAME,
                   engine="openpyxl").fillna("")

sess = requests.Session()
sess.headers.update({"User-Agent": "alphafold-downloader/2.0"})

for _, row in df.iterrows():
    gene = str(row.get("Gene", "")).strip().lower()
    if not gene:
        continue

    for col, org in ORG_COLS.items():
        acc0 = get_uniprot(row.get(col, ""))
        if not acc0:
            continue

        out_pdb = f"{gene}_{org}.pdb"
        if os.path.exists(out_pdb):
            print(f"SKIP  {out_pdb}")
            continue

        # candidats: accession brute, sans suffixe isoforme, + remap UniProt
        candidates = []
        candidates.extend([acc0, acc0.split("-")[0]]
                          if "-" in acc0 else [acc0])

        expanded = []
        for a in candidates:
            expanded.extend(uniprot_accessions(a, sess))

        # dédoublonner
        seen = set()
        candidates = []
        for a in expanded:
            if a and a not in seen:
                candidates.append(a)
                seen.add(a)

        got = False

        for acc in candidates:
            # 1) API
            api_url = alphafold_api_pdb_url(acc, sess)
            if api_url and download_pdb_from_url(api_url, out_pdb, sess):
                print(f"OK    {out_pdb}  (API via {acc})")
                got = True
                break

            # 2) fallback direct v6..v2 (+ .gz)
            direct_url = try_direct_files(acc, out_pdb, sess)
            if direct_url:
                print(f"OK    {out_pdb}  (direct {
                      direct_url.split('/')[-1]} via {acc})")
                got = True
                break

        if not got:
            print(f"FAIL  {gene} {
                  org} ({acc0}) -> aucun PDB AlphaFold trouvé (API + direct v6..v2)")
