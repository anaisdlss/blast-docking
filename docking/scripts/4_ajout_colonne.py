#!/usr/bin/env python3
import csv
from pathlib import Path
import subprocess
import shutil
from plip.structure.preparation import PDBComplex

# =====================
# CONFIG
# =====================
ROOT = Path(".").resolve()
WORK = ROOT / "work"
RESULTS = ROOT / "results"

IN_CSV = RESULTS / "summary_refine_filtered.csv"
OUT_CSV = RESULTS / "summary_refine_filtered_with_pocket_props.csv"

CUTOFF_A = 4.0

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

# =====================
# RDKit (robuste)
# =====================
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
except Exception:
    raise SystemExit(
        "\n[ERROR] RDKit n'est pas disponible.\n"
        "Installe-le puis relance.\n"
        "Ex: conda install -c conda-forge rdkit\n"
    )

# =====================
# LIGANDS (RDKit)
# =====================
SDF_SEARCH_DIRS = [
    WORK / "ligands" / "in",
    WORK / "ligands" / "3d",
    WORK / "ligands",
    ROOT / "liguands",
    ROOT / "ligands",
    ROOT,
]

LIGAND_COLS = {
    "lig_cyclic",
    "lig_len_aa",
    "lig_polarity_0_1",
    "lig_charge_ph74",
    "lig_flex_0_1",
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

AMIDE = Chem.MolFromSmarts("C(=O)N")

# =====================
# HELPERS
# =====================


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


def run_cmd(cmd, quiet=False):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        if not quiet:
            print("[CMD ERROR]", " ".join(cmd))
            print("STDOUT:", r.stdout[:500])
            print("STDERR:", r.stderr[:500])
        raise RuntimeError("Command failed")
    return r.stdout


def find_obabel_cmd():
    """
    Retourne une liste de commande base:
    - ["obabel"] si dispo dans PATH
    - sinon ["pixi","run","obabel"] (utile si obabel n'existe que dans pixi env)
    """
    if shutil.which("obabel"):
        return ["obabel"]
    # fallback: pixi run obabel
    if shutil.which("pixi"):
        return ["pixi", "run", "obabel"]
    raise RuntimeError("OpenBabel (obabel) introuvable (ni PATH ni pixi).")


def pdbqt_to_pdb(pdbqt_path: Path, out_pdb_path: Path) -> bool:
    try:
        obabel = find_obabel_cmd()
        run_cmd(obabel + ["-ipdbqt", str(pdbqt_path),
                "-opdb", "-O", str(out_pdb_path)], quiet=True)
        return out_pdb_path.exists() and out_pdb_path.stat().st_size > 0
    except Exception:
        return False


def make_complex_pdb(receptor_pdb: Path, ligand_pdb: Path, out_complex_pdb: Path) -> None:
    """
    Concatène receptor + ligand dans un seul PDB.
    Force le resname du ligand à LIG + HETATM.
    """
    with open(out_complex_pdb, "w") as w:
        # Protéine
        with open(receptor_pdb, "r", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    w.write(line)
        w.write("TER\n")

        # Ligand
        with open(ligand_pdb, "r", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                line = "HETATM" + line[6:]                 # force HETATM
                line = line[:17] + "LIG" + line[20:]       # force resname LIG
                w.write(line)

        w.write("END\n")


def plip_hbond_count(complex_pdb: Path) -> int:
    """
    Compte hbonds PLIP (hbonds_pdon + hbonds_ldon).
    """
    mol = PDBComplex()
    mol.load_pdb(str(complex_pdb))
    mol.analyze()

    total = 0
    for _bsid, s in mol.interaction_sets.items():
        total += len(getattr(s, "hbonds_pdon", []))
        total += len(getattr(s, "hbonds_ldon", []))
    return int(total)

# =====================
# POCKET (simple)
# =====================


def parse_pdb_atoms(pdb_path: Path):
    atoms = []
    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            resn = line[17:20].strip().upper()
            chain = line[21:22].strip()
            resi = line[22:26].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atoms.append({"chain": chain, "resi": resi,
                         "resn": resn, "x": x, "y": y, "z": z})
    return atoms


def parse_pdbqt_coords(pdbqt_path: Path):
    coords = []
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
            coords.append((x, y, z))
    return coords


def dist2(a, b):
    dx = a[0]-b[0]
    dy = a[1]-b[1]
    dz = a[2]-b[2]
    return dx*dx + dy*dy + dz*dz


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

# =====================
# LIGAND RDKit (identique à ton approche)
# =====================


def ligand_stem_from_row(row: dict) -> str:
    s = Path(row.get("ligand", "")).stem
    s = s.replace("_rdkit_3d", "").replace("_2D", "").replace("_3D", "")
    return s


def find_ligand_sdf(row: dict) -> Path | None:
    lig = row.get("ligand", "")
    if not isinstance(lig, str) or not lig:
        return None
    p = Path(lig)
    if p.exists():
        return p

    filename = p.name
    stem_raw = Path(filename).stem
    stem = stem_raw.replace("_2D", "").replace("_3D", "")

    candidates = [
        filename,
        f"{stem}.sdf",
        f"{stem}_2D.sdf",
        f"{stem}_3D.sdf",
        f"{stem}_rdkit_3d.sdf",
        f"{stem}_3D_rdkit_3d.sdf",
    ]

    for d in SDF_SEARCH_DIRS:
        if not d.exists():
            continue
        for name in candidates:
            cand = d / name
            if cand.exists():
                return cand
        for fp in d.rglob("*.sdf"):
            if fp.name in candidates:
                return fp
    return None


def compute_ligand_props_from_sdf(sdf_path: Path) -> dict:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return {}

    mw = Descriptors.MolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    ring = rdMolDescriptors.CalcNumRings(mol)
    formal_charge = Chem.GetFormalCharge(mol)

    stem = sdf_path.name.replace(".sdf", "")
    stem = stem.replace("_rdkit_3d", "").replace(
        "_2D", "").replace("_3D", "").replace("_3D_rdkit_3d", "")
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

    with open(IN_CSV, newline="") as f:
        in_rows = list(csv.DictReader(f))

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

        print(f"[PROGRESS] {i}/{total} | {tag}")

        # --- Ligand props (RDKit) ---
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
            stem_row = ligand_stem_from_row(row)
            row["lig_cyclic"] = MANUAL_LIG_CYCLIC.get(
                stem_row, row.get("lig_cyclic", ""))

        # --- Docking files ---
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

        # --- Pocket props ---
        protein_atoms = parse_pdb_atoms(rec_pdb)
        ligand_coords = parse_pdbqt_coords(lig_pdbqt)
        pocket = pocket_residues(protein_atoms, ligand_coords, CUTOFF_A)

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

        # --- HBonds via PLIP ---
        lig_pdb = results_dir / "ligand_plip.pdb"
        complex_pdb = results_dir / "complex_plip.pdb"

        if pdbqt_to_pdb(lig_pdbqt, lig_pdb):
            make_complex_pdb(rec_pdb, lig_pdb, complex_pdb)
            try:
                row["hbond_count"] = str(plip_hbond_count(complex_pdb))
                computed_hbonds += 1
            except Exception as e:
                # utile pour debug
                print(f"[WARN] PLIP failed for {tag}: {e}")
                row["hbond_count"] = ""
        else:
            print(f"[WARN] obabel conversion failed for {tag}")
            row["hbond_count"] = ""

        out_rows.append(row)

    # sort by affinity if present
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
    print(f"[INFO] Ligands calculés: {computed_ligands}")
    print(f"[INFO] Pockets calculées: {computed_pockets}")
    print(f"[INFO] Hbonds PLIP calculées: {computed_hbonds}")
    if missing_sdf:
        print(f"[WARN] {missing_sdf} lignes: SDF introuvable/illisible.")
    if missing_files:
        print(f"[WARN] {missing_files} lignes: fichiers docking manquants.")


if __name__ == "__main__":
    main()
