#!/usr/bin/env python3
"""
Master runner script for ESKAPE AMR analysis pipeline.
Runs all analysis scripts in the correct order.
"""
import subprocess
import sys
import os

SCRIPTS = [
    "01_build_amr_matrix.py",
    "02_mdr_analysis.py",
    "03_core_accessory_amr.py",
    "04_final_summary.py",
    "05_heatmaps.py",
]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

print("=" * 60)
print(" ESKAPE Analysis Pipeline Runner")
print("=" * 60)

for script in SCRIPTS:
    path = os.path.join(BASE_DIR, script)
    if not os.path.exists(path):
        print(f"\n❌ Script not found: {path}")
        sys.exit(1)

    print(f"\n{'─' * 60}")
    print(f"▶ Running {script}...")
    print(f"{'─' * 60}")

    result = subprocess.run([sys.executable, path])
    if result.returncode != 0:
        print(f"\n❌ {script} failed with exit code {result.returncode}")
        sys.exit(1)

print("\n" + "=" * 60)
print(" ✅ ALL SCRIPTS COMPLETED SUCCESSFULLY!")
print("=" * 60)
