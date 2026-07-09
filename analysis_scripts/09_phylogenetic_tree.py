#!/usr/bin/env python3
"""
09_phylogenetic_tree.py
=======================
Build standard maximum-likelihood phylogenies from Roary core-genome
alignments using FastTree (Price et al. 2009, 2010).

Because IQ-TREE 2 was not available in the configured conda channels,
FastTree 2.1.10 is used as the published, model-based ML fallback:
  * General Time Reversible (GTR) nucleotide substitution model
  * Gamma20 rate heterogeneity (-gamma)
  * 1,000 local bootstrap replicates (-boot 1000)

Inputs
------
roary_out/<species>/core_gene_alignment.aln   (FASTA)

Outputs
-------
analysis_results/09_phylogenetic_tree/<species>_core_tree.nwk
analysis_results/09_phylogenetic_tree/<species>_fasttree.log
analysis_results/09_phylogenetic_tree/tree_summary.csv

Note
----
Run this script inside the `eskape_phylogeny` conda environment, or with
FastTree (capital F) on the system PATH.
"""

import os
import subprocess
import shutil
import sys
import time
import multiprocessing
import numpy as np
import pandas as pd
from Bio import AlignIO

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "09_phylogenetic_tree")
os.makedirs(OUT_DIR, exist_ok=True)

# FastTree is multi-threaded via OpenMP. Cap per-species threads to avoid
# oversubscription when building all six trees in parallel.
FASTTREE_THREADS = int(os.environ.get("FASTTREE_THREADS", "8"))


def check_fasttree():
    """Ensure FastTree is available and report its version."""
    fasttree = shutil.which("FastTree")
    if fasttree is None:
        print(
            "\n❌ FastTree executable not found on PATH.\n"
            "   Please activate the 'eskape_phylogeny' conda environment:\n"
            "      conda activate eskape_phylogeny\n"
            "   or install FastTree and ensure it is on PATH.\n",
            flush=True,
        )
        sys.exit(1)

    try:
        version = subprocess.check_output(
            [fasttree, "-help"], stderr=subprocess.STDOUT, text=True, timeout=10
        ).splitlines()[0]
    except Exception:
        version = "FastTree (version unknown)"
    return fasttree, version


def alignment_stats(alignment):
    """
    Return alignment length and number of polymorphic (SNP) sites.

    A site is counted as a SNP if, after excluding gaps ('-') and ambiguous
    'N'/'n' characters, it contains at least two different nucleotides.
    The implementation is vectorized with NumPy for large core-genome
    alignments.
    """
    n_strains = len(alignment)
    aln_len = alignment.get_alignment_length()

    # Convert sequences to a 2D uint8 array of ASCII byte values.
    seqs = [str(record.seq) for record in alignment]
    arr_bytes = np.array([list(s.encode()) for s in seqs], dtype=np.uint8)

    # Mask invalid characters (gaps and Ns).
    valid = (
        (arr_bytes != ord("-"))
        & (arr_bytes != ord("N"))
        & (arr_bytes != ord("n"))
    )

    valid_count = valid.sum(axis=0)
    has_valid = valid.any(axis=0)

    # Min / max valid ASCII value per column; columns with no valid chars are
    # set to identical sentinel values so they are not counted as SNPs.
    min_valid = np.where(
        has_valid,
        np.min(np.where(valid, arr_bytes, 255), axis=0),
        0,
    )
    max_valid = np.where(
        has_valid,
        np.max(np.where(valid, arr_bytes, 0), axis=0),
        0,
    )

    polymorphic = (valid_count >= 2) & (min_valid != max_valid)
    n_snps = int(polymorphic.sum())

    return n_strains, aln_len, n_snps


def build_tree_for_species(args):
    """Worker function: build one FastTree ML tree for a single species."""
    species, fasttree_bin = args
    species_cap = species.capitalize()
    aln_path = os.path.join(ROARY_DIR, species, "core_gene_alignment.aln")
    tree_out = os.path.join(OUT_DIR, f"{species}_core_tree.nwk")
    log_out = os.path.join(OUT_DIR, f"{species}_fasttree.log")

    if not os.path.exists(aln_path):
        return {
            "status": "skipped",
            "Species": species_cap,
            "message": f"alignment not found: {aln_path}",
        }

    print(f"\n🌲 {species_cap}: loading alignment...", flush=True)
    try:
        alignment = AlignIO.read(aln_path, "fasta")
    except Exception as e:
        return {
            "status": "error",
            "Species": species_cap,
            "message": f"Failed to parse alignment: {e}",
        }

    print(f"   {species_cap}: computing alignment statistics...", flush=True)
    t0 = time.time()
    n_strains, aln_len, n_snps = alignment_stats(alignment)
    stats_time = time.time() - t0
    print(
        f"   {species_cap}: {aln_len:,} bp × {n_strains} strains "
        f"({n_snps:,} SNP sites; stats computed in {stats_time:.1f}s)",
        flush=True,
    )

    names = [record.id for record in alignment]
    if len(names) != len(set(names)):
        print(
            f"   ⚠️  {species_cap}: duplicate sequence IDs detected; "
            "FastTree may fail.",
            flush=True,
        )

    # FastTree command: nucleotide, GTR, Gamma, 1000 local bootstrap replicates
    cmd = [
        fasttree_bin,
        "-nt",
        "-gtr",
        "-gamma",
        "-boot", "1000",
        "-nopr",
        aln_path,
    ]

    print(
        f"   {species_cap}: running FastTree ML inference "
        f"(up to {FASTTREE_THREADS} threads)...",
        flush=True,
    )
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(FASTTREE_THREADS)
    try:
        with open(tree_out, "w") as tree_fh, open(log_out, "w") as log_fh:
            subprocess.run(
                cmd,
                stdout=tree_fh,
                stderr=log_fh,
                check=True,
                text=True,
                env=env,
            )
    except subprocess.CalledProcessError as e:
        return {
            "status": "error",
            "Species": species_cap,
            "message": f"FastTree failed (exit code {e.returncode}); see {log_out}",
        }

    elapsed = time.time() - t0
    print(
        f"   ✅ {species_cap}: tree saved ({elapsed:.1f}s): {tree_out}",
        flush=True,
    )

    return {
        "status": "success",
        "Species": species_cap,
        "Strains": n_strains,
        "Alignment_Bp": aln_len,
        "SNP_Sites": n_snps,
        "Substitution_Model": "GTR+Gamma",
        "Bootstrap_Replicates": 1000,
        "Bootstrap_Type": "FastTree local bootstrap",
        "Runtime_Sec": round(elapsed, 1),
        "Tree": tree_out,
        "Log": log_out,
    }


# ===============================
# BUILD ML TREES WITH FASTTREE
# ===============================
fasttree_bin, fasttree_version = check_fasttree()

species_list = sorted([
    d for d in os.listdir(ROARY_DIR)
    if os.path.isdir(os.path.join(ROARY_DIR, d))
])

print("=" * 70, flush=True)
print("Building Core Genome Maximum-Likelihood Trees with FastTree", flush=True)
print(f"Executable: {fasttree_bin}", flush=True)
print(f"Version:    {fasttree_version}", flush=True)
print("Model:      GTR+Gamma, 1000 local bootstrap replicates", flush=True)
print(f"Parallel species: {len(species_list)} (max {FASTTREE_THREADS} threads each)", flush=True)
print("=" * 70, flush=True)

start_all = time.time()

# Build all species trees in parallel using one process per species.
with multiprocessing.Pool(processes=len(species_list)) as pool:
    results = pool.map(build_tree_for_species, [(s, fasttree_bin) for s in species_list])

tree_results = []
errors = []
skipped = []
for res in results:
    if res["status"] == "success":
        tree_results.append({
            k: v for k, v in res.items() if k != "status"
        })
    elif res["status"] == "skipped":
        skipped.append(res)
    else:
        errors.append(res)

# ===============================
# SUMMARY
# ===============================
print("\n" + "=" * 70, flush=True)
print("SUMMARY", flush=True)
print("=" * 70, flush=True)

if tree_results:
    summary_df = pd.DataFrame(tree_results)
    summary_path = os.path.join(OUT_DIR, "tree_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print("\n📋 Tree Summary:", flush=True)
    print(summary_df.to_string(index=False), flush=True)
    print(f"\n✅ All trees and logs saved to: {OUT_DIR}", flush=True)

if skipped:
    print("\n⚠️  Skipped species:", flush=True)
    for item in skipped:
        print(f"   - {item['Species']}: {item['message']}", flush=True)

if errors:
    print("\n❌ Errors:", flush=True)
    for item in errors:
        print(f"   - {item['Species']}: {item['message']}", flush=True)

if not tree_results:
    print("\n⚠️  No trees were generated.", flush=True)

print(
    f"\n⏱️  Total wall-clock time: {time.time() - start_all:.1f}s",
    flush=True,
)
print("\n🎉 Phylogenetic tree construction complete!", flush=True)
