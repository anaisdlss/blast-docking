#!/usr/bin/env python3
import csv
import re
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

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

# Parallélisme (Mac M3, tu peux mettre 4 ou 6 si tu veux)
MAX_WORKERS = 4
CPU_PER_DOCK = 1  # IMPORTANT : 1 sinon ça explose

# Conversion optionnelle
CONVERT_PDB = False

# =====================
# FILTERS (FINAL strict) ✅ inchangés
# =====================
TRANSPORTER_PREFIXES = ("abca", "abcb", "abcc", "abcg",
                        "cftr", "slco", "oatp", "mrp")
ENZYMES = ("tmprss2", "corin", "ace", "ppp5c")

AFF_TRANSPORTER = -7.0
RMSD_TRANSPORTER = 2.5
AFF_ENZYME = -7.5
RMSD_ENZYME = 2.0

# =====================
# FAST PREFILTER (buffer) ✅ plus permissif
# =====================
FAST_AFF_TRANSPORTER = -6.0
FAST_RMSD_TRANSPORTER = 4.0
FAST_AFF_ENZYME = -7.0
FAST_RMSD_ENZYME = 2.5


# =====================
# UTILS
# =====================
def which_or_fail(exe: str) -> str:
    p = shutil.which(exe)
    if not p:
        raise RuntimeError("Executable not found in PATH: " + exe)
    return p


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


def classify_protein(protein_filename: str) -> str:
    s = Path(protein_filename).stem.lower()
    if s.startswith(TRANSPORTER_PREFIXES):
        return "transporteur"
    if any(e in s for e in ENZYMES):
        return "enzyme"
    return "autre"


def pass_thresholds_final(protein_filename: str, aff: float, rmsd_lb: float) -> bool:
    ptype = classify_protein(protein_filename)
    if ptype == "transporteur":
        return (aff <= AFF_TRANSPORTER) and (rmsd_lb <= RMSD_TRANSPORTER)
    if ptype == "enzyme":
        return (aff <= AFF_ENZYME) and (rmsd_lb <= RMSD_ENZYME)
    return False


def pass_thresholds_fast_prefilter(protein_filename: str, aff: float, rmsd_lb: float) -> bool:
    ptype = classify_protein(protein_filename)
    if ptype == "transporteur":
        return (aff <= FAST_AFF_TRANSPORTER) and (rmsd_lb <= FAST_RMSD_TRANSPORTER)
    if ptype == "enzyme":
        return (aff <= FAST_AFF_ENZYME) and (rmsd_lb <= FAST_RMSD_ENZYME)
    return False


def receptor_center_from_pdbqt(pdbqt: Path):
    xs, ys, zs = [], [], []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    if not xs:
        raise RuntimeError("No ATOM lines read from: " + str(pdbqt))
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def parse_vina_table(text: str):
    best = None
    rmsd_lb_list = []
    for ln in text.splitlines():
        m = re.match(
            r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)", ln)
        if not m:
            continue
        mode = int(m.group(1))
        aff = float(m.group(2))
        rmsd_lb = float(m.group(3))
        if mode == 1:
            best = aff
        elif mode >= 2:
            rmsd_lb_list.append(rmsd_lb)
    return best, (min(rmsd_lb_list) if rmsd_lb_list else None)


def obabel_convert(obabel: str, inp: Path, out: Path) -> bool:
    p = subprocess.run([obabel, str(inp), "-O", str(out)],
                       capture_output=True, text=True)
    return p.returncode == 0


def run_one_refine(job: dict) -> dict:
    """
    Worker SAFE : ne crashe jamais sans renvoyer un dict.
    """
    try:
        vina = job["vina"]
        obabel = job.get("obabel", None)

        tag = job["tag"]
        protein = job["protein"]
        ligand = job["ligand"]

        rec_pdbqt = Path(job["rec_pdbqt"])
        rec_nowat = Path(job["rec_nowat_pdb"])
        lig_pdbqt = Path(job["lig_pdbqt"])
        cx, cy, cz = job["center"]
        out_dir = Path(job["out_dir"])

        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / "out.pdbqt"
        log_path = out_dir / "vina.log"

        # copie receptor pour PyMOL
        rec_copy = out_dir / "receptor_nowat.pdb"
        if rec_nowat.exists() and (not rec_copy.exists()):
            try:
                shutil.copy2(rec_nowat, rec_copy)
            except Exception:
                pass

        # SKIP si déjà fait
        if out_pdbqt.exists() and log_path.exists():
            txt = log_path.read_text(errors="ignore")
            best, min_rmsd = parse_vina_table(txt)
            return {
                "status": "SKIP",
                "protein": protein,
                "ligand": ligand,
                "type": classify_protein(protein),
                "tag": tag,
                "best_affinity_kcal_mol": best,
                "min_rmsd_lb_modes2plus": min_rmsd,
                "results_dir": str(out_dir),
                "pymol_open": f"pymol {out_dir/'receptor_nowat.pdb'} {out_dir/'out.pdb'}",
            }

        sx, sy, sz = BOX_SIZE
        cmd = [
            vina,
            "--receptor", str(rec_pdbqt),
            "--ligand", str(lig_pdbqt),
            "--center_x", f"{cx:.3f}",
            "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", str(sx),
            "--size_y", str(sy),
            "--size_z", str(sz),
            "--exhaustiveness", str(REF_EXH),
            "--num_modes", str(REF_MODES),
            "--cpu", str(CPU_PER_DOCK),
            "--out", str(out_pdbqt),
            "--verbosity", "1",
        ]

        p = subprocess.run(cmd, capture_output=True, text=True)
        log_text = (p.stdout or "") + \
            ("\n" if p.stdout else "") + (p.stderr or "")
        log_path.write_text(log_text)

        if p.returncode != 0:
            return {
                "status": "ERROR",
                "protein": protein,
                "ligand": ligand,
                "type": classify_protein(protein),
                "tag": tag,
                "best_affinity_kcal_mol": None,
                "min_rmsd_lb_modes2plus": None,
                "results_dir": str(out_dir),
                "pymol_open": "",
            }

        best, min_rmsd = parse_vina_table(log_text)

        # conversion optionnelle
        if CONVERT_PDB and obabel:
            out_pdb = out_dir / "out.pdb"
            if not out_pdb.exists():
                obabel_convert(obabel, out_pdbqt, out_pdb)

        return {
            "status": "OK",
            "protein": protein,
            "ligand": ligand,
            "type": classify_protein(protein),
            "tag": tag,
            "best_affinity_kcal_mol": best,
            "min_rmsd_lb_modes2plus": min_rmsd,
            "results_dir": str(out_dir),
            "pymol_open": f"pymol {out_dir/'receptor_nowat.pdb'} {out_dir/'out.pdb'}",
        }

    except Exception as e:
        # Worker fail-safe
        return {
            "status": "ERROR",
            "protein": job.get("protein", ""),
            "ligand": job.get("ligand", ""),
            "type": classify_protein(job.get("protein", "")) if job.get("protein") else "autre",
            "tag": job.get("tag", ""),
            "best_affinity_kcal_mol": None,
            "min_rmsd_lb_modes2plus": None,
            "results_dir": job.get("out_dir", ""),
            "pymol_open": "",
        }


def main():
    if not FAST_CSV.exists():
        raise RuntimeError("summary_fast.csv not found: " + str(FAST_CSV))

    vina = which_or_fail("vina")
    obabel = shutil.which("obabel") if CONVERT_PDB else None

    with open(FAST_CSV) as f:
        rows = list(csv.DictReader(f))

    # 1) FAST prefilter (buffer)
    hits = []
    for r in rows:
        aff = safe_float(r.get("best_affinity_kcal_mol"))
        rmsd = safe_float(r.get("min_rmsd_lb_modes2plus"))
        if aff is None or rmsd is None:
            continue
        if classify_protein(r["protein"]) == "autre":
            continue
        if pass_thresholds_fast_prefilter(r["protein"], aff, rmsd):
            hits.append(r)

    print(f"[INFO] Hits retenus pour REFINE (buffer FAST) : {len(hits)}")

    # 2) cache centers
    centers = {}
    for r in hits:
        protein_stem = Path(r["protein"]).stem
        rec_pdbqt = WORK_REC_PDBQT / f"{protein_stem}.pdbqt"
        if rec_pdbqt not in centers and rec_pdbqt.exists():
            centers[rec_pdbqt] = receptor_center_from_pdbqt(rec_pdbqt)

    # 3) build jobs
    jobs = []
    for r in hits:
        tag = r["tag"]
        protein_stem = Path(r["protein"]).stem
        ligand_stem = Path(r["ligand"]).stem.replace("_2D", "")

        rec_pdbqt = WORK_REC_PDBQT / f"{protein_stem}.pdbqt"
        rec_nowat = WORK_REC_CLEAN / f"{protein_stem}_nowat.pdb"
        lig_pdbqt = WORK_LIG_PDBQT / f"{ligand_stem}.pdbqt"

        if not (rec_pdbqt.exists() and rec_nowat.exists() and lig_pdbqt.exists()):
            continue
        if rec_pdbqt not in centers:
            continue

        out_dir = RESULTS / "REFINE" / tag

        jobs.append({
            "vina": vina,
            "obabel": obabel,
            "tag": tag,
            "protein": r["protein"],
            "ligand": r["ligand"],
            "rec_pdbqt": str(rec_pdbqt),
            "rec_nowat_pdb": str(rec_nowat),
            "lig_pdbqt": str(lig_pdbqt),
            "center": centers[rec_pdbqt],
            "out_dir": str(out_dir),
        })

    print(f"[INFO] Dockings REFINE à exécuter / skipper : {len(jobs)}")
    print(f"[INFO] Parallèle = {MAX_WORKERS} jobs | cpu/job={CPU_PER_DOCK}")

    refine_rows_all = []
    refine_rows_filtered = []

    # 4) run parallel
    done = 0
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(run_one_refine, j) for j in jobs]
        for fut in as_completed(futures):
            row = fut.result()
            done += 1
            status = row.get("status", "OK")
            tag = row.get("tag", "")
            best = row.get("best_affinity_kcal_mol", None)
            mr = row.get("min_rmsd_lb_modes2plus", None)
            print(
                f"[{done}/{len(jobs)}] {status} {tag} | best={best} | minRMSDlb={mr}")

            refine_rows_all.append(row)
            if best is not None and mr is not None and pass_thresholds_final(row["protein"], float(best), float(mr)):
                refine_rows_filtered.append(row)

    # 5) write CSV (extrasaction ignore = impossible de planter)
    RESULTS.mkdir(parents=True, exist_ok=True)

    keys = [
        "status", "protein", "ligand", "type", "tag",
        "best_affinity_kcal_mol", "min_rmsd_lb_modes2plus",
        "results_dir", "pymol_open"
    ]

    if refine_rows_all:
        out_all = RESULTS / "summary_refine_all.csv"
        with open(out_all, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(refine_rows_all)

        print(f"[OK] CSV REFINE (tous) : {out_all}")
        print(f"[INFO] Nombre de dockings REFINE : {len(refine_rows_all)}")

    # Tri par affinité : du plus négatif (meilleur) au moins négatif
    refine_rows_filtered.sort(
        key=lambda x: x["best_affinity_kcal_mol"]
    )

    if refine_rows_filtered:
        out_filt = RESULTS / "summary_refine_filtered.csv"
        with open(out_filt, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(refine_rows_filtered)

        print(f"[OK] CSV REFINE filtré : {out_filt}")
        print(f"[INFO] Après filtre strict : {len(refine_rows_filtered)}")
    else:
        print("[INFO] Aucun résultat REFINE ne passe les seuils stricts.")

    print("\n[OK] REFINE terminé.")


if __name__ == "__main__":
    main()
