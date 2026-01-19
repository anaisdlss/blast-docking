#!/usr/bin/env python3
import csv
from pathlib import Path

# =====================
# CONFIG
# =====================
ROOT = Path(".").resolve()
WORK = ROOT / "work"
RESULTS = ROOT / "results"

IN_CSV = RESULTS / "summary_refine_filtered.csv"
OUT_CSV = RESULTS / "summary_refine_filtered_with_pocket_props.csv"

CUTOFF_A = 4.0

# Où chercher les .sdf (le script essaie plusieurs dossiers)
SDF_SEARCH_DIRS = [
    WORK / "ligands" / "in",
    WORK / "ligands" / "3d",
    WORK / "ligands",
    ROOT / "liguands",   # ton dossier existe (orthographe)
    ROOT / "ligands",    # au cas où
    ROOT,
]

# Résidus
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR_UNCHARGED = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
POS = {"LYS", "ARG"}
NEG = {"ASP", "GLU"}

POCKET_COLS = {
    "pocket_residues_n",
    "hydrophobic_pct",
    "polar_pct",
    "pos_count",
    "neg_count",
    "net_charge",
}

HBOND_COLS = {
    "hbond_count",
}

# atomes accepteurs/donneurs (proxy simple)
HBOND_ACCEPTORS = {"O", "N", "S"}
HBOND_DONORS = {"O", "N"}  # sans H explicites, on reste simple
HBOND_DIST_A = 3.5

# Colonnes ligands qu'on veut remplir de façon robuste
LIGAND_COLS = {
    "lig_cyclic",
    "lig_len_aa",
    "lig_polarity_0_1",
    "lig_charge_ph74",
    "lig_flex_0_1",
    # + infos “brutes” utiles pour justifier les scores
    "lig_mw",
    "lig_tpsa",
    "lig_rot_bonds",
    "lig_ring_count",
    "lig_formal_charge",
    "lig_charge_kind",
}

MANUAL_LIG_CYCLIC = {
    "microcystineLR": "cyclique",
    "cyanopeptolinA": "cyclique",
    "aerucyclamideA": "cyclique",
    "microcyclamide7806B": "cyclique",
    "aeruginosin98A": "lineaire",
    "microgininFR1": "lineaire",
}

# =====================
# RDKit (robuste)
# =====================
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
except Exception as e:
    raise SystemExit(
        "\n[ERROR] RDKit n'est pas disponible dans ton environnement.\n"
        "Installe-le (conda-forge recommandé), puis relance.\n"
        "Ex: conda install -c conda-forge rdkit\n"
    )

# =====================
# Parsers pocket
# =====================


def parse_pdb_atoms(pdb_path: Path):
    atoms = []
    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip().upper()
            resn = line[17:20].strip().upper()
            chain = line[21:22].strip()
            resi = line[22:26].strip()

            # élément : colonne 77-78 si propre, sinon fallback sur atom name
            elem = line[76:78].strip().upper()
            if not elem:
                atom_name = line[12:16].strip().upper()
                elem = "".join([c for c in atom_name if c.isalpha()])[:1]

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atoms.append({
                "chain": chain, "resi": resi, "resn": resn,
                "x": x, "y": y, "z": z,
                "elem": elem,
                "atom_name": atom_name,
            })
    return atoms


def parse_pdbqt_atoms(pdbqt_path: Path):
    atoms = []
    with open(pdbqt_path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            # atom type pdbqt est souvent en fin de ligne (ex: OA, NA, HD...)
            parts = line.split()
            atype = parts[-1] if parts else ""
            atype = atype.upper()

            # déduire élément simple
            elem = atype[:1]  # O, N, S, C...
            if elem not in {"C", "N", "O", "S", "P", "F", "I", "B", "H", "CL", "BR"}:
                elem = line[12:16].strip().upper()[:1]

            atoms.append({"x": x, "y": y, "z": z, "elem": elem})
    return atoms


BACKBONE_NAMES = {"N", "O", "C", "CA", "OXT"}


def count_hbond_like_pocket_residues(protein_atoms, ligand_atoms, pocket_set, cutoff_a=3.2):
    if not protein_atoms or not ligand_atoms or not pocket_set:
        return 0
    cutoff2 = cutoff_a * cutoff_a

    lig = [a for a in ligand_atoms if a.get("elem") in {"O", "N"}]
    if not lig:
        return 0

    pocket_keys = set(pocket_set)

    prot = []
    for a in protein_atoms:
        key = (a["chain"], a["resi"], a["resn"])
        if key not in pocket_keys:
            continue
        if a.get("elem") not in {"O", "N"}:
            continue
        if a.get("atom_name", "") in BACKBONE_NAMES:
            continue
        prot.append(a)

    hit = set()
    for pa in prot:
        pcoord = (pa["x"], pa["y"], pa["z"])
        pres = (pa["chain"], pa["resi"], pa["resn"])
        for la in lig:
            if dist2(pcoord, (la["x"], la["y"], la["z"])) <= cutoff2:
                hit.add(pres)
                break
    return len(hit)


def count_hbonds_residue_level(protein_atoms, ligand_atoms, cutoff_a=3.5):
    if not protein_atoms or not ligand_atoms:
        return 0

    cutoff2 = cutoff_a * cutoff_a
    hbond_residues = set()

    prot = [
        a for a in protein_atoms
        if a.get("elem") in {"O", "N"}
    ]
    lig = [
        a for a in ligand_atoms
        if a.get("elem") in {"O", "N"}
    ]

    for pa in prot:
        pcoord = (pa["x"], pa["y"], pa["z"])
        pres = (pa["chain"], pa["resi"], pa["resn"])

        for la in lig:
            if dist2(pcoord, (la["x"], la["y"], la["z"])) <= cutoff2:
                hbond_residues.add(pres)
                break   # ⛔ un résidu = 1 max

    return len(hbond_residues)


def dist2(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz


def pocket_residues(protein_atoms, ligand_coords, cutoff_a: float):
    if not protein_atoms or not ligand_coords:
        return set()
    cutoff2 = cutoff_a * cutoff_a
    pocket = set()
    for at in protein_atoms:
        acoord = (at["x"], at["y"], at["z"])
        for lc in ligand_coords:
            if dist2(acoord, lc) <= cutoff2:
                pocket.add((at["chain"], at["resi"], at["resn"]))
                break
    return pocket


def safe_float(x):
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", "None", "nan", "NaN"):
        return None
    try:
        return float(s)
    except ValueError:
        return None

# =====================
# Ligand helpers (robuste)
# =====================


def ligand_stem_from_row(row: dict) -> str:
    """
    Exemple:
    - microcystineLR_2D.sdf -> microcystineLR
    - microcyclamide7806B_3D_rdkit_3d.sdf -> microcyclamide7806B
    """
    s = Path(row.get("ligand", "")).stem
    s = s.replace("_rdkit_3d", "")
    s = s.replace("_2D", "").replace("_3D", "")
    return s


def find_ligand_sdf(row: dict) -> Path | None:
    lig = row.get("ligand", "")
    if not isinstance(lig, str) or not lig:
        return None

    p = Path(lig)
    if p.exists():
        return p

    filename = p.name  # ex: microcystineLR_2D.sdf
    stem_raw = Path(filename).stem  # microcystineLR_2D
    stem = stem_raw.replace("_2D", "").replace("_3D", "")  # microcystineLR

    # toutes les variantes qu’on accepte chez toi
    candidates = [
        filename,                 # microcystineLR_2D.sdf
        f"{stem}.sdf",            # microcystineLR.sdf
        f"{stem}_2D.sdf",
        f"{stem}_3D.sdf",
        f"{stem}_rdkit_3d.sdf",           # microcystineLR_rdkit_3d.sdf
        f"{stem}_3D_rdkit_3d.sdf",        # microcyclamide7806B_3D_rdkit_3d.sdf
    ]

    for d in SDF_SEARCH_DIRS:
        if not d.exists():
            continue

        # recherche directe
        for name in candidates:
            cand = d / name
            if cand.exists():
                return cand

        # recherche récursive (si tu as des sous-dossiers)
        for fp in d.rglob("*.sdf"):
            if fp.name in candidates:
                return fp

    return None


AMIDE = Chem.MolFromSmarts("C(=O)N")  # proxy peptide-like


def compute_ligand_props_from_sdf(sdf_path: Path) -> dict:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = None
    for m in suppl:
        if m is not None:
            mol = m
            break

    if mol is None:
        return {}

    mw = Descriptors.MolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    ring = rdMolDescriptors.CalcNumRings(mol)
    formal_charge = Chem.GetFormalCharge(mol)   # ✅ ICI

    # ⚠️ RDKit "ring>0" = présence d'au moins un cycle, ça ne veut PAS dire cyclopeptide.
    # On force une vérité biochimique avec un mapping manuel.
    stem = sdf_path.name
    stem = stem.replace(".sdf", "")
    stem = stem.replace("_rdkit_3d", "")
    stem = stem.replace("_2D", "").replace("_3D", "")
    stem = stem.replace("_3D_rdkit_3d", "")
    lig_cyclic = MANUAL_LIG_CYCLIC.get(stem, "")

    lig_polarity_0_1 = min(tpsa / 200.0, 1.0)
    lig_flex_0_1 = rot / (rot + ring + 1.0)

    amide_bonds = len(mol.GetSubstructMatches(AMIDE))
    lig_len_aa = amide_bonds + 1 if amide_bonds > 0 else ""

    return {
        "lig_cyclic": lig_cyclic,
        "lig_len_aa": lig_len_aa,
        "lig_polarity_0_1": round(lig_polarity_0_1, 3),
        "lig_flex_0_1": round(lig_flex_0_1, 3),
        "lig_mw": round(mw, 2),
        "lig_tpsa": round(tpsa, 2),
        "lig_rot_bonds": int(rot),
        "lig_ring_count": int(ring),
        "lig_formal_charge": int(formal_charge),
        "lig_charge_ph74": int(formal_charge),
        "lig_charge_kind": "formal_charge_rdkit",
    }


def ligand_already_done(row: dict) -> bool:
    # considéré “fait” si toutes colonnes ligand existent ET lig_tpsa non vide (stable)
    if not LIGAND_COLS.issubset(row.keys()):
        return False
    return str(row.get("lig_tpsa", "")).strip() != ""


def pocket_already_done(row: dict) -> bool:
    if not POCKET_COLS.issubset(row.keys()):
        return False
    return str(row.get("pocket_residues_n", "")).strip() != ""

# =====================
# MAIN
# =====================


def main():
    print("[INFO] IN_CSV  =", IN_CSV.resolve())
    print("[INFO] OUT_CSV =", OUT_CSV.resolve())

    if not IN_CSV.exists():
        raise SystemExit(f"[ERROR] Input CSV not found: {IN_CSV}")

    # 1) rows d'entrée
    with open(IN_CSV, newline="") as f:
        in_rows = list(csv.DictReader(f))

    # 2) reprendre OUT si existe
    existing_by_tag = {}
    if OUT_CSV.exists():
        print(f"[INFO] Reprise depuis fichier existant: {OUT_CSV}")
        with open(OUT_CSV, newline="") as f:
            for r in csv.DictReader(f):
                t = r.get("tag", "")
                if t:
                    existing_by_tag[t] = r

    total = len(in_rows)
    print(f"[INFO] Lignes à traiter: {total}")

    out_rows = []
    missing_files = 0
    missing_sdf = 0

    computed_pockets = 0
    computed_ligands = 0
    computed_hbonds = 0

    for i, base in enumerate(in_rows, start=1):
        tag = base.get("tag", "")
        row = dict(existing_by_tag.get(tag, base))

        if i == 1 or i % 10 == 0:
            print(f"[PROGRESS] {i}/{total} | {tag}")

        # ============================
        # LIGAND PROPS
        # ============================
        sdf_path = find_ligand_sdf(row)

        if sdf_path is None:
            missing_sdf += 1
            for k in LIGAND_COLS:
                row.setdefault(k, "")
        else:
            if not ligand_already_done(row):
                props = compute_ligand_props_from_sdf(sdf_path)
                if props:
                    computed_ligands += 1
                    for k in LIGAND_COLS:
                        if str(row.get(k, "")).strip() == "":
                            row[k] = props.get(k, "")
                else:
                    missing_sdf += 1
                    for k in LIGAND_COLS:
                        row.setdefault(k, "")

            stem_row = ligand_stem_from_row(row)
            row["lig_cyclic"] = MANUAL_LIG_CYCLIC.get(
                stem_row, row.get("lig_cyclic", ""))

        # ============================
        # FICHIERS DOCKING
        # ============================
        results_dir = Path(row.get("results_dir", ""))
        rec_pdb = results_dir / "receptor_nowat.pdb"
        lig_pdbqt = results_dir / "out.pdbqt"

        if (not rec_pdb.exists()) or (not lig_pdbqt.exists()):
            missing_files += 1
            for k in POCKET_COLS:
                row.setdefault(k, "")
            row["hbond_count"] = ""
            out_rows.append(row)
            continue

        protein_atoms = parse_pdb_atoms(rec_pdb)
        ligand_atoms = parse_pdbqt_atoms(lig_pdbqt)
        ligand_coords = [(a["x"], a["y"], a["z"]) for a in ligand_atoms]

        # ============================
        # POCKET (toujours calculé ici car utile pour hbonds strict)
        # ============================
        pocket = pocket_residues(protein_atoms, ligand_coords, CUTOFF_A)

        # ============================
        # POCKET PROPS (si manquant)
        # ============================
        if not pocket_already_done(row):
            residues = [r[2] for r in pocket]
            n_total = len(residues)

            if n_total == 0:
                hyd_pct = pol_pct = 0.0
                pos_n = neg_n = net = 0
            else:
                hyd_n = sum(1 for rn in residues if rn in HYDROPHOBIC)
                pol_n = sum(1 for rn in residues if rn in POLAR_UNCHARGED)
                pos_n = sum(1 for rn in residues if rn in POS)
                neg_n = sum(1 for rn in residues if rn in NEG)
                net = pos_n - neg_n
                hyd_pct = round(100.0 * hyd_n / n_total, 1)
                pol_pct = round(100.0 * pol_n / n_total, 1)

            row["pocket_residues_n"] = str(n_total)
            row["hydrophobic_pct"] = str(hyd_pct)
            row["polar_pct"] = str(pol_pct)
            row["pos_count"] = str(pos_n)
            row["neg_count"] = str(neg_n)
            row["net_charge"] = str(net)
            computed_pockets += 1
        else:
            for k in POCKET_COLS:
                row.setdefault(k, "")

        # ============================
        # HBONDS STRICT (poche only) ✅
        # ============================
        row["hbond_count"] = str(
            count_hbond_like_pocket_residues(
                protein_atoms, ligand_atoms, pocket, cutoff_a=3.2
            )
        )
        computed_hbonds += 1

        # ✅ ✅ ✅ OUBLI MAJEUR : ajouter la ligne
        out_rows.append(row)

    # tri affinité
    def aff_key(r):
        a = safe_float(r.get("best_affinity_kcal_mol"))
        return a if a is not None else 9999.0

    out_rows.sort(key=aff_key)

    all_keys = set()
    for r in out_rows:
        all_keys.update(r.keys())

    preferred = [
        "status", "protein", "ligand", "type", "tag",
        "best_affinity_kcal_mol", "min_rmsd_lb_modes2plus",
        "results_dir", "pymol_open",
        "lig_cyclic", "lig_len_aa", "lig_polarity_0_1", "lig_flex_0_1",
        "lig_charge_ph74", "lig_charge_kind", "lig_formal_charge",
        "lig_mw", "lig_tpsa", "lig_rot_bonds", "lig_ring_count",
        "pocket_residues_n", "hydrophobic_pct", "polar_pct",
        "pos_count", "neg_count", "net_charge",
        "hbond_count",
    ]
    fieldnames = [k for k in preferred if k in all_keys] + \
        sorted([k for k in all_keys if k not in preferred])

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    print(f"[OK] Wrote: {OUT_CSV}")
    print(f"[INFO] Ligands calculés cette fois: {computed_ligands}")
    print(f"[INFO] Pockets calculées cette fois: {computed_pockets}")
    print(f"[INFO] Hbonds calculées cette fois: {computed_hbonds}")

    if missing_sdf:
        print(
            f"[WARN] {missing_sdf} lignes: SDF ligand introuvable/illisible.")
    if missing_files:
        print(f"[WARN] "
              f"{missing_files} lignes: fichiers docking manquants (receptor_nowat.pdb ou out.pdbqt).")


if __name__ == "__main__":
    main()
