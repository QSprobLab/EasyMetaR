# LC‑MS/MS Metabolomics Data Processing Workflow

This repository contains a set of Python scripts that automate the analysis of metabolomics data acquired with liquid chromatography tandem mass spectrometry (LC‑MS/MS). The workflow is designed to work on **Windows, Linux** and **macOS** and can be run directly from the command line or integrated into a GitHub Actions CI pipeline.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Folder Structure](#folder-structure)
- [Scripts & Usage](#scripts--usage)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

The workflow performs the following steps:
1. **Convert raw data** (e.g. Thermo `.d` → `.mzML`) using MSConvert.
2. **Extract target peaks** for each polarity (`pos`/`neg`).
3. **Merge polarity data** to produce a consolidated feature table.
4. **Impute missing values** and normalize to TIC = 1000.

The scripts are located in `script/` and rely on a simple CSV database of target metabolites.

---

## Prerequisites

| Component | Version | Notes |
|---|---|---|
| Python | 3.8+ | The scripts use only the standard library and the packages listed in `requirements.txt` |
| MSConvert | [ProteoWizard](https://proteowizard.sourceforge.io/download.html) | Converts Thermo `.d` to `.mzML` |
| Conda (optional) | 4.9+ | Helpful for installing `pyopenms` on Windows |

---

## Installation

```bash
# Create a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# Upgrade pip and install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

On Windows you might prefer to use Conda for `pyopenms`:

```bash
conda install -c conda-forge pyopenms
```

---

## Folder Structure

```
project-root/
├── pos/                      # positive‑ion data
│   ├── BLK/                  # Blank samples
│   ├── QC/                   # QC samples
│   ├── MSMS/                 # MS2 QC files
│   └── <group‑folders>/      # Experimental groups (e.g. GFFBif_end)
├── neg/                      # negative‑ion data (same layout as pos)
├── script/                   # Python scripts
├── db/                       # Database files
│   └── ms1_metabolome_db.csv
└── metadata.csv              # Optional sample metadata
```

Rename each raw file to `groupname_originalfilename.mzML`.  For example:

| Sample | Type | Batch |
|---|---|---|
| BLK_1.mzML | Blank | Batch3 |
| QC_1.mzML | QC | Batch1 |
| GFFBif_end_6.mzML | Sample | Batch2 |

---

## Scripts & Usage

### 1. `extract_targets_pra.py`

Extracts target peaks for a single polarity.

```bash
python script/extract_targets_pra.py \
    --mode <pos|neg> \
    --db db/ms1_metabolome_db.csv \
    --work_dir ./<polarity‑folder> \
    --relaxed_filtering \
    --mz_tolerance 10 \
    --abs_int_threshold 500 \
    --skip_ctrl_filter \
    --n_workers 10
```

**Parameters**
| Flag | Default | Description |
|---|---|---|
| `--mode` | **required** | `pos` or `neg` |
| `--db` | **required** | Path to the CSV database |
| `--work_dir` | **required** | Directory containing `.mzML` files |
| `--mz_tolerance` | 10 | ppm tolerance for m/z matching |
| `--rt_tolerance` | 0.5 | minutes tolerance for retention time |
| `--abs_int_threshold` | 500 | Absolute intensity threshold |
| `--n_workers` | CPU‑1 | Number of parallel workers |

The script outputs a CSV file named `targets_<mode>_<polarity>_filtered.csv` in the same directory.

### 2. `postprocess_qc_rsc.py`

Merges positive and negative polarity feature tables and generates a final feature table.

```bash
python script/postprocess_qc_rsc.py \
    --pos targets_pos_filtered.csv \
    --neg targets_neg_filtered.csv \
    --output-tab feature_table_raw.csv \
    --output-info feature_info.csv \
    --db db/ms1_metabolome_db.csv
```

The outputs are:
- `feature_table_raw.csv` – rows are metabolites, columns are samples.
- `feature_info.csv` – metadata per feature (m/z, RT, adduct, DB match, …).

### 3. `fill_and_tic.py`

Imputes zeros with the 1/10 of the minimal non‑zero value per row and normalises each column to a TIC of 1000.

```bash
python script/fill_and_tic.py
```

It expects `feature_table_raw.csv` in the current directory and writes `feature_table.csv`.

---

## Full Workflow Example

```bash
# 1. Set working directory
cd /path/to/project-root

# 2. Verify files are in the correct folders (pos/, neg/, script/, db/)

# 3. Extract positive polarity peaks
python script/extract_targets_pra.py \
    --mode pos \
    --db db/ms1_metabolome_db.csv \
    --work_dir pos \
    --relaxed_filtering \
    --mz_tolerance 10 \
    --abs_int_threshold 500 \
    --skip_ctrl_filter \
    --n_workers 10

# 4. Extract negative polarity peaks
python script/extract_targets_pra.py \
    --mode neg \
    --db db/ms1_metabolome_db.csv \
    --work_dir neg \
    --relaxed_filtering \
    --mz_tolerance 10 \
    --abs_int_threshold 500 \
    --skip_ctrl_filter \
    --n_workers 10

# 5. Merge polarities
python script/postprocess_qc_rsc.py \
    --pos pos/targets_pos_filtered.csv \
    --neg neg/targets_neg_filtered.csv \
    --output-tab feature_table_raw.csv \
    --output-info feature_info.csv \
    --db db/ms1_metabolome_db.csv

# 6. Impute and normalise
python script/fill_and_tic.py
```

---

## Troubleshooting

| Issue | Possible Fix |
|---|---|
| `pyopenms` install fails on Windows | Use Conda or install the pre‑compiled wheel from `https://github.com/pyopenms/pyopenms/releases` |
| File encoding errors (e.g. GBK) | The scripts automatically try UTF‑8, GBK, GB2312. Ensure your CSV files are UTF‑8 if possible |
| Out‑of‑memory | Reduce `--n_workers` or split data into batches |
| No QC samples | The scripts skip RT calibration if no QC is detected |

---

## License

This project is released under the MIT license.  See the `LICENSE` file for details.

---

## Citation

If you use this workflow in a publication, please cite:
- *Metabolomics data analysis using LC‑MS/MS and Python*.
- `MSConvert` (ProteoWizard), `pyOpenMS`, and the scripts in this repository.
