#!/usr/bin/env python3
import re
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(".").resolve()

PROTEINS_DIR = ROOT / "proteins" / "known"
LIGANDS_DIR = ROOT / "liguands"

WORK = ROOT / "work"
WORK_REC_CLEAN = WORK / "receptors" / "clean"
WORK_REC_PDBQT = WORK / "receptors" / "pdbqt"
WORK_LIG_IN = WORK / "ligands" / "in"
WORK_LIG_3D = WORK / "ligands" / "3d"
WORK_LIG_PDBQT = WORK / "ligands" / "pdbqt"

RESULTS = ROOT / "results"
LOGS = ROOT / "logs"

PAIRS = [
    # ACE + microginine
    ("ace_homo.pdb", "microgininFR1_2D.sdf", "ACE__microgininFR1"),

    # Transporteurs + microcystine
    ("oatp1b1_homo.pdb", "microcystineLR_2D.sdf",
     "OATP1B1_homo__microcystineLR"),
    ("oatp1d1_truitearc.pdb", "microcystineLR_2D.sdf",
     "OATP1D1_truite__microcystineLR"),
    ("oatp1d1_zebrafish.pdb", "microcystineLR_2D.sdf",
     "OATP1D1_zfish__microcystineLR"),
    ("oatp1f2_zebrafish.pdb", "microcystineLR_2D.sdf",
     "OATP1F2_zfish__microcystineLR"),

    # MRP2 (ABCC2) + microcystine
    ("mrp2_homo.pdb",         "microcystineLR_2D.sdf", "MRP2_homo__microcystineLR"),
]

DEFAULT_EXH = 16
DEFAULT_MODES = 50

BOX_ACE = (24, 24, 24)
BOX_OATP = (30, 30, 30)


def ensure_dirs():
    for d in [
        WORK_REC_CLEAN, WORK_REC_PDBQT, WORK_LIG_IN, WORK_LIG_3D, WORK_LIG_PDBQT,
        RESULTS, LOGS, LOGS/"prep_receptors", LOGS/"prep_ligands"
    ]:
        d.mkdir(parents=True, exist_ok=True)


def which_or_fail(exe: str) -> Path:
    p = shutil.which(exe)
    if not p:
        raise FileNotFoundError(f"Executable not found in PATH: {exe}")
    return Path(p)


def pixi_bin(script_name: str) -> Path:
    py = which_or_fail("python")
    bin_dir = py.parent
    script = bin_dir / script_name
    if not script.exists():
        raise FileNotFoundError(f"Cannot find {script_name} in "
                                f"{bin_dir}. Are you inside `pixi shell`?")
    return script


def run(cmd, log_path: Path | None = None):
    cmd = [str(x) for x in cmd]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = (p.stdout or "") + ("\n" if p.stdout else "") + (p.stderr or "")
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(out)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed (exit={p.returncode}): "
                           f"{' '.join(cmd)}\n\n{out}")
    return out


def clean_pdb_atom_only(pdb_in: Path, pdb_out: Path):
    """
    Nettoyage robuste pour docking Vina/Meeko:
    - garde uniquement ATOM/TER/END (donc supprime tous les HETATM : glycans, ions, ligands, etc.)
    - supprime implicitement les eaux (souvent en HETATM de toute façon)
    """
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "TER", "END")):
                fout.write(line)


def receptor_center_from_pdbqt(pdbqt: Path):
    xs, ys, zs = [], [], []
    with open(pdbqt, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    if not xs:
        raise ValueError(f"No atoms read from {pdbqt}")
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))


def parse_vina_best_affinity(log_text: str):
    # prend l'affinité de la ligne du mode 1
    for ln in log_text.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+(\.\d+)?)\s+", ln)
        if m:
            return float(m.group(1))
    return None


def prepare_receptor(pdb_name: str) -> Path:
    pdb_in = PROTEINS_DIR / pdb_name
    if not pdb_in.exists():
        raise FileNotFoundError(pdb_in)

    raw = WORK_REC_CLEAN / f"{pdb_in.stem}_raw.pdb"
    nowat = WORK_REC_CLEAN / f"{pdb_in.stem}_nowat.pdb"
    shutil.copy2(pdb_in, raw)
    clean_pdb_atom_only(raw, nowat)

    mk_prepare_receptor = pixi_bin("mk_prepare_receptor.py")
    out_base = WORK_REC_PDBQT / pdb_in.stem

    # IMPORTANT : pas de -v ici (box), sinon Meeko exige un centre/taille.
    log_path = LOGS / "prep_receptors" / f"{pdb_in.stem}.log"
    run(["python", mk_prepare_receptor, "-i", nowat,
        "-o", out_base, "-p", "-j", "-a"], log_path)

    out_pdbqt = Path(str(out_base) + ".pdbqt")
    if not out_pdbqt.exists():
        raise RuntimeError(f"Receptor PDBQT not produced: "
                           f"{out_pdbqt}\nSee {log_path}")
    return out_pdbqt


def rdkit_make_3d(sdf_in: Path, sdf_out: Path, log_path: Path):
    code = f"""
from rdkit import Chem
from rdkit.Chem import AllChem
inp = r"{sdf_in}"
out = r"{sdf_out}"
suppl = Chem.SDMolSupplier(inp, removeHs=False)
mol = suppl[0] if suppl and len(suppl)>0 else None
if mol is None:
    raise SystemExit("RDKit cannot read SDF")
mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.useSmallRingTorsions = True
params.useMacrocycleTorsions = True
ok = AllChem.EmbedMolecule(mol, params)
if ok != 0:
    ok = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    if ok != 0:
        raise SystemExit("3D embedding failed")
if AllChem.MMFFHasAllMoleculeParams(mol):
    AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s', maxIters=2000)
else:
    AllChem.UFFOptimizeMolecule(mol, maxIters=4000)
w = Chem.SDWriter(out)
w.write(mol); w.close()
print("OK", out)
"""
    run(["python", "-c", code], log_path)


def prepare_ligand(lig_sdf_name: str) -> Path:
    sdf_in = LIGANDS_DIR / lig_sdf_name
    if not sdf_in.exists():
        raise FileNotFoundError(sdf_in)

    WORK_LIG_IN.mkdir(parents=True, exist_ok=True)
    sdf_work = WORK_LIG_IN / lig_sdf_name
    shutil.copy2(sdf_in, sdf_work)

    WORK_LIG_3D.mkdir(parents=True, exist_ok=True)
    sdf_3d = WORK_LIG_3D / \
        (Path(lig_sdf_name).stem.replace("_2D", "") + "_rdkit_3d.sdf")

    rdkit_make_3d(
        sdf_work, sdf_3d,
        LOGS / "prep_ligands" / f"{Path(lig_sdf_name).stem}_rdkit.log"
    )

    mk_prepare_ligand = pixi_bin("mk_prepare_ligand.py")
    WORK_LIG_PDBQT.mkdir(parents=True, exist_ok=True)
    out_pdbqt = WORK_LIG_PDBQT / \
        (Path(lig_sdf_name).stem.replace("_2D", "") + ".pdbqt")

    run(
        ["python", mk_prepare_ligand, "-i", sdf_3d, "-o", out_pdbqt],
        LOGS / "prep_ligands" / f"{Path(lig_sdf_name).stem}_pdbqt.log"
    )

    if not out_pdbqt.exists():
        raise RuntimeError(f"Ligand PDBQT not produced: {out_pdbqt}")
    return out_pdbqt


def write_vina_conf(path: Path, receptor: Path, ligand: Path, center, size_xyz, exh, modes):
    cx, cy, cz = center
    sx, sy, sz = size_xyz
    txt = f"""receptor = {receptor}
ligand   = {ligand}

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = {sx}
size_y = {sy}
size_z = {sz}

exhaustiveness = {exh}
num_modes = {modes}
"""
    path.write_text(txt)


def obabel_convert(inp: Path, out: Path):
    obabel = which_or_fail("obabel")
    run([obabel, inp, "-O", out])


def dock_one(tag: str, rec_pdbqt: Path, rec_nowat_pdb: Path, lig_pdbqt: Path, size_xyz):
    vina = which_or_fail("vina")
    out_dir = RESULTS / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    center = receptor_center_from_pdbqt(rec_pdbqt)
    conf = out_dir / "vina.conf"
    write_vina_conf(conf, rec_pdbqt, lig_pdbqt, center,
                    size_xyz, DEFAULT_EXH, DEFAULT_MODES)

    out_pdbqt = out_dir / "out.pdbqt"
    log_path = out_dir / "vina.log"

    # Vina 1.2.x : pas de --log, on capture stdout/stderr dans vina.log
    vina_log = run([vina, "--config", conf, "--out",
                   out_pdbqt, "--verbosity", "1"], log_path)

    out_pdb = out_dir / "out.pdb"
    obabel_convert(out_pdbqt, out_pdb)

    shutil.copy2(rec_nowat_pdb, out_dir / "receptor_nowat.pdb")

    best = parse_vina_best_affinity(vina_log)
    return {
        "tag": tag,
        "receptor": rec_pdbqt.name,
        "ligand": lig_pdbqt.name,
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "size_x": size_xyz[0],
        "size_y": size_xyz[1],
        "size_z": size_xyz[2],
        "exhaustiveness": DEFAULT_EXH,
        "num_modes": DEFAULT_MODES,
        "best_affinity_kcal_mol": best,
        "results_dir": str(out_dir),
        "pymol_open": f"pymol {out_dir/'receptor_nowat.pdb'} {out_dir/'out.pdb'}"
    }


def main():
    ensure_dirs()

    rec_cache = {}
    lig_cache = {}
    summary = []

    for protein_pdb, ligand_sdf, tag in PAIRS:
        if protein_pdb not in rec_cache:
            rec_pdbqt = prepare_receptor(protein_pdb)
            rec_nowat_pdb = WORK_REC_CLEAN / \
                f"{Path(protein_pdb).stem}_nowat.pdb"
            rec_cache[protein_pdb] = (rec_pdbqt, rec_nowat_pdb)

        if ligand_sdf not in lig_cache:
            lig_pdbqt = prepare_ligand(ligand_sdf)
            lig_cache[ligand_sdf] = lig_pdbqt

        rec_pdbqt, rec_nowat_pdb = rec_cache[protein_pdb]
        lig_pdbqt = lig_cache[ligand_sdf]

        size_xyz = BOX_ACE if protein_pdb.startswith("ace_") else BOX_OATP

        row = dock_one(tag, rec_pdbqt, rec_nowat_pdb, lig_pdbqt, size_xyz)
        summary.append(row)
        print(f"✅ Done: {tag} | best={row['best_affinity_kcal_mol']} kcal/mol")

    out_csv = RESULTS / "summary_scores.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    print("\nCSV written:", out_csv)
    for r in summary:
        print(f"{r['tag']}: best="
              f"{r['best_affinity_kcal_mol']} | {r['results_dir']}")
        print("Open in PyMOL:", r["pymol_open"])


if __name__ == "__main__":
    main()
