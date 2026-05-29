#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract target peaks from mzML files for a given polarity.

Usage:
    python script/extract_targets_pra.py \
        --mode pos|neg \
        --db db/ms1_metabolome_db.csv \
        --work_dir ./pos \
        --mz_tolerance 10 \
        --abs_int_threshold 500 \
        --rel_int_threshold 0.5 \
        --n_workers 4
"""

import os
import argparse
import pandas as pd
from pyopenms import MSExperiment, MzMLFile
import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

# Polarity specific adducts
ADDUCTS_POS = {
    "[M+H]+": {"charge": 1, "mass_diff": 1.007276, "priority": 1},
    "[M+Na]+": {"charge": 1, "mass_diff": 22.989218, "priority": 2},
    "[M+NH4]+": {"charge": 1, "mass_diff": 18.033823, "priority": 3},
    "[M+K]+": {"charge": 1, "mass_diff": 38.963158, "priority": 4},
    "[M+H-H2O]+": {"charge": 1, "mass_diff": -17.003288, "priority": 5},
}

ADDUCTS_NEG = {
    "[M-H]-": {"charge": -1, "mass_diff": -1.007276, "priority": 1},
    "[M+Cl]-": {"charge": -1, "mass_diff": 34.969402, "priority": 2},
    "[M+FA-H]-": {"charge": -1, "mass_diff": 44.998201, "priority": 3},
}


def load_database(db_path, mode):
    """Load the metabolite database and expand with adducts."""
    df = pd.read_csv(db_path)
    adducts = ADDUCTS_POS if mode == "pos" else ADDUCTS_NEG
    expanded = []
    for _, r in df.iterrows():
        base = r["MonoisotopicMass"]
        for name, info in adducts.items():
            mz = (base + info["mass_diff"]) / abs(info["charge"])  # simplified
            expanded.append({
                "compound_name": r["compound_name"],
                "mz": mz,
                "adduct": name,
                "priority": info["priority"],
            })
    return pd.DataFrame(expanded)


def process_file(file, rt_range, abs_thr, rel_thr):
    exp = MSExperiment()
    try:
        MzMLFile().load(file, exp)
    except Exception:
        return []
    data = []
    for spec in exp:
        rt = spec.getRT() / 60  # minutes
        if rt < rt_range[0] or rt > rt_range[1]:
            continue
        mzs, ints = spec.get_peaks()
        mask = ints >= abs_thr
        if rel_thr > 0 and mask.any():
            max_int = ints[mask].max()
            mask &= ints >= rel_thr * max_int
        idx = np.where(mask)[0]
        for i in idx:
            data.append({
                "sample": os.path.splitext(os.path.basename(file))[0],
                "mz": mzs[i],
                "rt": rt,
                "intensity": ints[i],
            })
    return data


def worker(file, rt_range, abs_thr, rel_thr):
    return process_file(file, rt_range, abs_thr, rel_thr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["pos", "neg"])
    parser.add_argument("--db", required=True)
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--mz_tolerance", type=float, default=10.0)
    parser.add_argument("--abs_int_threshold", type=float, default=500.0)
    parser.add_argument("--rel_int_threshold", type=float, default=0.5)
    parser.add_argument("--rt_min", type=float, default=0.0)
    parser.add_argument("--rt_max", type=float, default=10.0)
    parser.add_argument("--n_workers", type=int, default=mp.cpu_count() - 1)
    args = parser.parse_args()

    db = load_database(args.db, args.mode)

    # Find mzML files
    mzml_files = [os.path.join(args.work_dir, f) for f in os.listdir(args.work_dir) if f.lower().endswith(".mzml")]
    if not mzml_files:
        print("No mzML files found.")
        return

    rt_range = (args.rt_min, args.rt_max)
    workers = args.n_workers
    pool = mp.Pool(workers)
    all_data = []
    for res in tqdm(pool.imap(partial(worker, rt_range=rt_range, abs_thr=args.abs_int_threshold, rel_thr=args.rel_int_threshold), mzml_files), total=len(mzml_files)):
        all_data.extend(res)
    pool.close()
    pool.join()

    if not all_data:
        print("No peaks extracted.")
        return

    peak_df = pd.DataFrame(all_data)
    # Merge with database by mz tolerance
    # For simplicity, perform a nearest neighbor merge on mz within tolerance
    mz_tol = args.mz_tolerance
    merged = []
    for _, r in peak_df.iterrows():
        mz = r["mz"]
        # find closest db entry
        diffs = (db["mz"] - mz).abs()
        idx = diffs.idxmin()
        if abs(db.loc[idx, "mz"] - mz) <= mz_tol:
            row = db.loc[idx].to_dict()
            row.update({"sample": r["sample"], "intensity": r["intensity"], "rt": r["rt"], "mz": mz})
            merged.append(row)
    out_df = pd.DataFrame(merged)
    out_path = os.path.join(args.work_dir, f"targets_{args.mode}_filtered.csv")
    out_df.to_csv(out_path, index=False)
    print(f"Output written to {out_path}")

if __name__ == "__main__":
    main()
"