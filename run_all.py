"""
run_all.py
Runs the full mitotyping pipeline end to end:
    1. 01_run_pipeline.py    - pseudobulk -> MitoPathway scores -> mitoPPS
    2. 02_age_regression.py  - mitoPPS ~ age + sex OLS regression
    3. 03_outlier_qc.py      - outlier flagging + sensitivity analysis + QC plots
    4. 04_visualization.py   - all publication-style figures

Each script also runs standalone (they just read/write CSVs in OUTDIR),
so re-run any single stage on its own once its inputs exist.

Usage:
    python run_all.py
"""

import runpy

STAGES = [
    "01_run_pipeline.py",
    "02_age_regression.py",
    "03_outlier_qc.py",
    "04_visualization.py",
]

if __name__ == "__main__":
    for stage in STAGES:
        print(f"\n{'='*70}\nRunning {stage}\n{'='*70}")
        runpy.run_path(stage, run_name="__main__")
