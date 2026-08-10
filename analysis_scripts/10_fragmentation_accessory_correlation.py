#!/usr/bin/env python3
"""
Script 10: Fragmentation vs accessory genome size correlation analysis

Assesses whether draft assembly fragmentation (total contigs, contig N50)
correlates with per-strain accessory genome size. This addresses Reviewer #3
comment 20 in DECISION_2.pdf.

Inputs:
    - analysis_results/04_final_summary/checkm2_summary.csv
    - roary_out/{species}/gene_presence_absence.csv

Outputs:
    - analysis_results/10_fragmentation_accessory_correlation/
        * fragmentation_accessory_correlation.csv   (summary statistics)
        * per_strain_accessory_fragmentation.csv    (raw per-strain data)
        * fragmentation_accessory_correlation.png   (Supplementary Figure S9)

Dependencies: pandas, numpy, scipy, matplotlib, seaborn
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Correlate assembly fragmentation with accessory genome size."
    )
    parser.add_argument(
        "--base-dir",
        default="/Users/admin/Desktop/ESKAPE",
        help="Project root directory.",
    )
    parser.add_argument(
        "--core-threshold",
        type=float,
        default=0.99,
        help="Presence threshold for defining core genes (default: 0.99).",
    )
    return parser.parse_args()


def load_checkm2(path):
    df = pd.read_csv(path)
    df["Species"] = df["Species"].str.lower().str.strip()
    df["Name"] = df["Name"].astype(str).str.strip()
    return df


def get_strain_columns(roary_df):
    """Return strain sample columns (Roary columns from index 14 onward)."""
    return roary_df.columns[14:]


def compute_accessory_counts(roary_path, core_threshold):
    roary = pd.read_csv(roary_path, low_memory=False)
    strain_cols = get_strain_columns(roary)
    n_strains = len(strain_cols)

    presence = roary[strain_cols].notnull().sum(axis=1)
    is_core = presence >= core_threshold * n_strains
    accessory_df = roary.loc[~is_core, strain_cols].notnull().sum(axis=0).reset_index()
    accessory_df.columns = ["Name", "Accessory_Count"]
    accessory_df["Name"] = accessory_df["Name"].astype(str).str.strip()
    return accessory_df, n_strains


def correlate(df, metric):
    x = df["Accessory_Count"].values
    y = df[metric].values
    if len(x) > 2:
        r_pearson, p_pearson = pearsonr(x, y)
        r_spearman, p_spearman = spearmanr(x, y)
    else:
        r_pearson = p_pearson = r_spearman = p_spearman = np.nan
    return r_pearson, p_pearson, r_spearman, p_spearman


def main():
    args = parse_args()
    base_dir = args.base_dir
    core_threshold = args.core_threshold

    roary_dir = os.path.join(base_dir, "roary_out")
    checkm2_file = os.path.join(
        base_dir, "analysis_results", "04_final_summary", "checkm2_summary.csv"
    )
    out_dir = os.path.join(
        base_dir, "analysis_results", "10_fragmentation_accessory_correlation"
    )
    os.makedirs(out_dir, exist_ok=True)

    species_list = [
        "acinetobacter",
        "enterobacter",
        "enterococcus",
        "klebsiella",
        "pseudomonas",
        "staphylococcus",
    ]

    ch = load_checkm2(checkm2_file)
    quality_cols = ["Name", "Total_Contigs", "Contig_N50"]

    results = []
    all_data = []

    for species in species_list:
        roary_path = os.path.join(roary_dir, species, "gene_presence_absence.csv")
        if not os.path.exists(roary_path):
            logging.warning("Missing Roary file for %s: %s", species, roary_path)
            continue

        logging.info("Processing %s ...", species)
        accessory_df, n_strains = compute_accessory_counts(roary_path, core_threshold)
        accessory_df["Species"] = species

        # Merge CheckM2 quality metrics with Roary strain names.
        # Strain names in checkm2_summary.csv must match those in the Roary matrix.
        ch_sp = ch[ch["Species"] == species][quality_cols].copy()
        merged = accessory_df.merge(ch_sp, on="Name", how="inner")
        if merged.empty:
            logging.warning("No matching strains for %s after merging with CheckM2.", species)
            continue
        if len(merged) < n_strains:
            logging.warning(
                "%s: matched %d of %d strains to CheckM2 quality metrics.",
                species.capitalize(),
                len(merged),
                n_strains,
            )

        all_data.append(merged)

        for metric in ["Total_Contigs", "Contig_N50"]:
            r_p, p_p, r_s, p_s = correlate(merged, metric)
            results.append({
                "Species": species.capitalize(),
                "Metric": metric,
                "N": len(merged),
                "Pearson_r": round(r_p, 3) if not np.isnan(r_p) else np.nan,
                "Pearson_p": p_p,
                "Spearman_rho": round(r_s, 3) if not np.isnan(r_s) else np.nan,
                "Spearman_p": p_s,
            })

    results_df = pd.DataFrame(results)
    results_path = os.path.join(out_dir, "fragmentation_accessory_correlation.csv")
    results_df.to_csv(results_path, index=False)
    logging.info("Wrote summary statistics to %s", results_path)

    all_df = pd.concat(all_data, ignore_index=True)
    # Keep only the columns needed for publication
    all_df["Species"] = all_df["Species"].str.capitalize()
    all_df = all_df[["Name", "Species", "Accessory_Count", "Total_Contigs", "Contig_N50"]]
    data_path = os.path.join(out_dir, "per_strain_accessory_fragmentation.csv")
    all_df.to_csv(data_path, index=False)
    logging.info("Wrote per-strain data to %s", data_path)

    # Build figure
    fig, axes = plt.subplots(nrows=6, ncols=2, figsize=(12, 18))
    fig.suptitle(
        "Relationship between assembly fragmentation and accessory genome size per strain",
        fontsize=14,
        weight="bold",
    )

    for i, species in enumerate(species_list):
        label = species.capitalize()
        df_sp = all_df[all_df["Species"] == label]

        ax1 = axes[i, 0]
        ax1.scatter(
            df_sp["Total_Contigs"],
            df_sp["Accessory_Count"],
            color="steelblue",
            edgecolor="black",
            s=60,
        )
        ax1.set_xlabel("Total contigs", fontsize=10)
        ax1.set_ylabel("Accessory gene count", fontsize=10)
        ax1.set_title(f"{label}: Accessory vs. Total contigs", fontsize=11)

        row = results_df[
            (results_df["Species"] == label) & (results_df["Metric"] == "Total_Contigs")
        ].iloc[0]
        ax1.text(
            0.05,
            0.95,
            f"Pearson r = {row['Pearson_r']}\np = {row['Pearson_p']:.3f}",
            transform=ax1.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        ax2 = axes[i, 1]
        ax2.scatter(
            df_sp["Contig_N50"],
            df_sp["Accessory_Count"],
            color="darkorange",
            edgecolor="black",
            s=60,
        )
        ax2.set_xlabel("Contig N50 (bp)", fontsize=10)
        ax2.set_ylabel("Accessory gene count", fontsize=10)
        ax2.set_title(f"{label}: Accessory vs. Contig N50", fontsize=11)
        ax2.ticklabel_format(style="scientific", axis="x", scilimits=(0, 0))

        row2 = results_df[
            (results_df["Species"] == label) & (results_df["Metric"] == "Contig_N50")
        ].iloc[0]
        ax2.text(
            0.05,
            0.95,
            f"Pearson r = {row2['Pearson_r']}\np = {row2['Pearson_p']:.3f}",
            transform=ax2.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    fig_path = os.path.join(out_dir, "fragmentation_accessory_correlation.png")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info("Wrote figure to %s", fig_path)

    print("\n" + results_df.to_string(index=False))


if __name__ == "__main__":
    main()
