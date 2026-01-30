#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build summary_fast.csv from docking FAST results directories.

Expected layout:
  ./docking/results/FAST/<tag>/
      out.pdbqt or out.pdb (vina remarks preferred)
      receptor_nowat.pdb
      (optional) *.sdf (ligand copy)
      (optional) *.pdb (original receptor input)
      (optional) *.log / *.txt

Output:
  ./results/summary_fast.csv

No filtering: we include ALL result folders in the CSV.
We still compute passed_affinity / passed_rmsd booleans for convenience.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


# Defaults inferred from your example
DEFAULT_AFFINITY_THRESHOLD = -6.5  # kcal/mol (more negative is better)
DEFAULT_RMSD_THRESHOLD = 2.0       # Angstrom (lower is better)


@dataclass
class VinaPose:
    affinity: float
    rmsd_lb: float
    rmsd_ub: float


VINA_REMARK_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)"
)

# Parses the table form sometimes found in logs:
# mode | affinity | rmsd lb | rmsd ub
VINA_TABLE_LINE_RE = re.compile(
    r"^\s*(\d+)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*$"
)


def read_text_quiet(p: Path) -> Optional[str]:
    try:
        return p.read_text(errors="ignore")
    except Exception:
        return None


def parse_vina_poses_from_text(text: str) -> List[VinaPose]:
    poses: List[VinaPose] = []

    # 1) Prefer REMARK VINA RESULT lines (most robust, mode order already correct)
    for m in VINA_REMARK_RE.finditer(text):
        aff = float(m.group(1))
        lb = float(m.group(2))
        ub = float(m.group(3))
        poses.append(VinaPose(affinity=aff, rmsd_lb=lb, rmsd_ub=ub))

    if poses:
        return poses

    # 2) Otherwise, try log table lines
    # We accept any block of lines that looks like: "1  -7.2  0.0  0.0"
    for line in text.splitlines():
        m = VINA_TABLE_LINE_RE.match(line)
        if m:
            # mode = int(m.group(1))  # not needed except for counting
            aff = float(m.group(2))
            lb = float(m.group(3))
            ub = float(m.group(4))
            poses.append(VinaPose(affinity=aff, rmsd_lb=lb, rmsd_ub=ub))

    return poses


def parse_vina_poses_in_dir(d: Path) -> List[VinaPose]:
    """
    Try common files in a FAST result directory to extract vina poses.
    Priority:
      out.pdbqt -> out.pdb -> any *.pdbqt/*.pdb -> any *.log/*.txt
    """
    candidates: List[Path] = []

    for name in ["out.pdbqt", "out.pdb"]:
        p = d / name
        if p.exists():
            candidates.append(p)

    # Other pdbqt/pdb
    candidates += sorted([p for p in d.glob("*.pdbqt")
                         if p.name not in {"out.pdbqt"}])
    candidates += sorted([p for p in d.glob("*.pdb")
                         if p.name not in {"out.pdb", "receptor_nowat.pdb"}])

    # Logs/texts
    candidates += sorted(list(d.glob("*.log")))
    candidates += sorted(list(d.glob("*.txt")))

    for p in candidates:
        txt = read_text_quiet(p)
        if not txt:
            continue
        poses = parse_vina_poses_from_text(txt)
        if poses:
            return poses

    return []


def pick_protein_name(tag: str, d: Path) -> str:
    """
    Prefer a protein pdb file present in dir (excluding out/receptor_nowat).
    Otherwise fallback to first part of tag + .pdb
    """
    pdbs = [p for p in d.glob(
        "*.pdb") if p.name not in {"out.pdb", "receptor_nowat.pdb"}]
    if len(pdbs) == 1:
        return pdbs[0].name

    # fallback from tag
    protein_part = tag.split("__", 1)[0] if "__" in tag else tag
    if not protein_part.lower().endswith(".pdb"):
        protein_part += ".pdb"
    return protein_part


def pick_ligand_name(tag: str, d: Path) -> str:
    """
    Prefer an sdf present in dir.
    Otherwise fallback to second part of tag + .sdf
    """
    sdfs = list(d.glob("*.sdf"))
    if len(sdfs) == 1:
        return sdfs[0].name
    if len(sdfs) > 1:
        # If multiple, pick the biggest (often the real ligand, not a tiny stub)
        sdfs_sorted = sorted(sdfs, key=lambda p: p.stat(
        ).st_size if p.exists() else 0, reverse=True)
        return sdfs_sorted[0].name

    # fallback from tag
    ligand_part = tag.split("__", 1)[1] if "__" in tag else "ligand"
    if not ligand_part.lower().endswith(".sdf"):
        ligand_part += ".sdf"
    return ligand_part


def safe_float_or_nan(x: Optional[float]) -> float:
    return float("nan") if x is None else float(x)


def main() -> int:
    # Paths relative to current working directory
    root = Path(__file__).resolve().parents[1]
    results_base = root / "docking" / "results" / "FAST"
    out_dir = root / "results"
    out_csv = out_dir / "summary_fast.csv"

    # Optional CLI args:
    #   python make_summary_fast.py [affinity_threshold] [rmsd_threshold]
    affinity_thr = DEFAULT_AFFINITY_THRESHOLD
    rmsd_thr = DEFAULT_RMSD_THRESHOLD
    if len(sys.argv) >= 2:
        affinity_thr = float(sys.argv[1])
    if len(sys.argv) >= 3:
        rmsd_thr = float(sys.argv[2])

    if not results_base.exists():
        print(f"ERROR: results folder not found: {
              results_base}", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for d in sorted([p for p in results_base.iterdir() if p.is_dir()]):
        tag = d.name

        protein = pick_protein_name(tag, d)
        ligand = pick_ligand_name(tag, d)

        poses = parse_vina_poses_in_dir(d)

        best_aff: Optional[float] = None
        min_rmsd_lb_modes2plus: Optional[float] = None

        if poses:
            # "best" is the most negative affinity.
            best_aff = min(p.affinity for p in poses)

            # modes2plus: exclude the first pose (mode 1)
            if len(poses) >= 2:
                min_rmsd_lb_modes2plus = min(p.rmsd_lb for p in poses[1:])

        passed_affinity = (best_aff is not None) and (best_aff <= affinity_thr)
        passed_rmsd = (min_rmsd_lb_modes2plus is not None) and (
            min_rmsd_lb_modes2plus <= rmsd_thr)

        # absolute paths like your example
        results_dir_abs = str(d.resolve())

        receptor_nowat = d / "receptor_nowat.pdb"
        out_pdb = d / "out.pdb"
        # If out.pdb doesn't exist but out.pdbqt exists, still create a command
        if not out_pdb.exists() and (d / "out.pdbqt").exists():
            out_pdb = d / "out.pdbqt"

        pymol_open = f"pymol {receptor_nowat.resolve()} {out_pdb.resolve()}"

        rows.append({
            "protein": protein,
            "ligand": ligand,
            "tag": tag,
            "best_affinity_kcal_mol": safe_float_or_nan(best_aff),
            "min_rmsd_lb_modes2plus": safe_float_or_nan(min_rmsd_lb_modes2plus),
            "passed_affinity": bool(passed_affinity),
            "passed_rmsd": bool(passed_rmsd),
            "results_dir": results_dir_abs,
            "pymol_open": pymol_open,
        })

    # Write CSV
    fieldnames = [
        "protein",
        "ligand",
        "tag",
        "best_affinity_kcal_mol",
        "min_rmsd_lb_modes2plus",
        "passed_affinity",
        "passed_rmsd",
        "results_dir",
        "pymol_open",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            # Make NaN look like blank (optional). If you prefer "nan", remove this block.
            for k in ["best_affinity_kcal_mol", "min_rmsd_lb_modes2plus"]:
                v = r[k]
                if isinstance(v, float) and math.isnan(v):
                    r[k] = ""
            w.writerow(r)

    print(f"OK: wrote {len(rows)} rows -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
