#!/usr/bin/env python3
import csv
import re
import shutil
import subprocess
from pathlib import Path

# =====================
# PATHS
# =====================
ROOT = Path(".").resolve()

WORK = ROOT / "work"
WORK_REC_CLEAN = WORK / "receptors" / "clean"
WORK_REC_PDBQT = WORK / "receptors" / "pdbqt"
WORK_LIG_PDBQT = WORK / "ligands" / "pdbqt"

RESULTS = ROOT / "results"
FAST_CSV = RESULTS / "summary_fast.csv"

# =====================
# REFINE PARAMETERS
# =====================
REF_EXH = 16
REF_MODES = 50
BOX_SIZE = (30, 30, 30)

# =====================
# FILTERS
# =====================
TRANSPORTER_PREFIXES = (
    "abca", "abcb", "abcc", "abcg", "cftr", "slco", "oatp", "mrp"
)

ENZYMES = ("tmprss2", "corin", "ace", "ppp5c")

AFF_TRANSPORTER = -6.5
RMSD_TRANSPORTER = 3.0

AFF_ENZYME = -7.5
RMSD_ENZYME = 2.0


# =====================
# UTILS
# =====================
def which_or_fail(exe: str) -> Path:
    p = shutil.which(exe)
    if not p:
        raise RuntimeError("Executable not found in PATH: " + exe)
    return Path(p)


def receptor_center_from_pdbqt(pdbqt: Path):
    xs, ys, zs = [], [], []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))


def parse_vina_table(text: str):
    best = None
    rmsds = []
    for ln in text.splitlines():
        m = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", ln)
        if not m:
            continue
        mode = int(m.group(1))
        aff = float(m.group(2))
        rmsd = float(m.group(3))
        if mode == 1:
            best = aff
        else:
            rmsds.append(rmsd)
    return best, min(rmsds) if rmsds else None


def write_vina_conf(path, receptor, ligand, center):
    cx, cy, cz = center
    sx, sy, sz = BOX_SIZE
    path.write_text(
        f"receptor = {receptor}\n"
        f"ligand = {ligand}\n\n"
        f"center_x = {cx}\n"
        f"center_y = {cy}\n"
        f"center_z = {cz}\n\n"
        f"size_x = {sx}\n"
        f"size_y = {sy}\n"
        f"size_z = {sz}\n\n"
        f"exhaustiveness = {REF_EXH}\n"
        f"num_modes = {REF_MODES}\n"
    )


# =====================
# MAIN
# =====================
def main():

    if not FAST_CSV.exists():
        raise RuntimeError("summary_fast.csv not found")

    vina = which_or_fail("vina")

    with open(FAST_CSV) as f:
        rows = list(csv.DictReader(f))

    hits = []

    for r in rows:
        prot = r["protein"].lower()
        aff = float(r["best_affinity_kcal_mol"])
        rmsd = float(r["min_rmsd_lb_modes2plus"])

        is_transporter = prot.startswith(TRANSPORTER_PREFIXES)
        is_enzyme = any(e in prot for e in ENZYMES)

        if is_transporter:
            if aff <= AFF_TRANSPORTER and rmsd <= RMSD_TRANSPORTER:
                hits.append(r)

        elif is_enzyme:
            if aff <= AFF_ENZYME and rmsd <= RMSD_ENZYME:
                hits.append(r)

    print(f"[INFO] Hits retenus pour REFINE : {len(hits)}")

    for r in hits:

        tag = r["tag"]
        protein = Path(r["protein"]).stem
        ligand = Path(r["ligand"]).stem.replace("_2D", "")

        rec_pdbqt = WORK_REC_PDBQT / f"{protein}.pdbqt"
        rec_pdb = WORK_REC_CLEAN / f"{protein}_nowat.pdb"
        lig_pdbqt = WORK_LIG_PDBQT / f"{ligand}.pdbqt"

        out_dir = RESULTS / "REFINE" / tag
        out_dir.mkdir(parents=True, exist_ok=True)

        out_pdbqt = out_dir / "out.pdbqt"
        log_path = out_dir / "vina.log"

        if out_pdbqt.exists() and log_path.exists():
            print(f"[SKIP] {tag}")
            continue

        center = receptor_center_from_pdbqt(rec_pdbqt)
        conf = out_dir / "vina.conf"
        write_vina_conf(conf, rec_pdbqt, lig_pdbqt, center)

        print(f"[REFINE] {tag}")

        p = subprocess.run(
            [vina, "--config", conf, "--out", out_pdbqt],
            capture_output=True,
            text=True
        )

        log_path.write_text((p.stdout or "") + (p.stderr or ""))

        if p.returncode != 0:
            print(f"[ERROR] Vina failed for {tag}")
            continue

    print("\n[OK] REFINE terminé")


if __name__ == "__main__":
    main()
