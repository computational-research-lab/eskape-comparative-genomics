import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
AMR_DIR = os.path.join(BASE_DIR, "analysis_results", "03_core_accessory_amr")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "07_fisher_core_accessory")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# HELPERS
# ===============================
def fisher_exact_r(a, b, c, d):
    """
    Compute exact odds ratio, p-value, and 95% CI using R's fisher.test.
    This ensures the CI and p-value are mutually consistent (both exact).
    Table layout:
        [[a, b],
         [c, d]]
    """
    r_code = (
        f'result <- fisher.test(matrix(c({a},{c},{b},{d}), nrow=2)); '
        f'cat(result$estimate, "\\n", result$p.value, "\\n", '
        f'result$conf.int[1], "\\n", result$conf.int[2], "\\n", sep="")'
    )
    try:
        result = subprocess.run(
            ['R', '--vanilla', '--slave', '-e', r_code],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        lines = result.stdout.strip().split('\n')
        or_val = float(lines[0])
        p_val = float(lines[1])
        ci_lo = float(lines[2])
        ci_hi = float(lines[3])
        return or_val, p_val, ci_lo, ci_hi
    except Exception as e:
        # Fallback to scipy if R is unavailable
        print(f"  ⚠️  R fisher.test failed ({e}), falling back to scipy")
        table = [[a, b], [c, d]]
        or_val, p_val = fisher_exact(table)
        return or_val, p_val, np.nan, np.nan


def fisher_with_ci(a, b, c, d, alpha=0.05):
    """
    Run Fisher's exact test on 2x2 table.
    Returns: odds_ratio, p_value, ci_lower, ci_upper
    Uses R's exact conditional method for CI to guarantee agreement with p-value.
    """
    return fisher_exact_r(a, b, c, d)


def benjamini_hochberg(p_values):
    """
    Compute Benjamini-Hochberg FDR q-values.
    Returns array of q-values in the same order as input p_values.
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    ranks = np.arange(1, n + 1)
    # BH-adjusted p-values (non-monotonic)
    raw_q = sorted_p * n / ranks
    # Enforce monotonicity from largest to smallest
    monotonic_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    # Clip at 1.0
    monotonic_q = np.clip(monotonic_q, 0.0, 1.0)
    # Restore original order
    q_values = np.empty(n)
    q_values[order] = monotonic_q
    return q_values


# ===============================
# MAIN
# ===============================
species_list = sorted([d for d in os.listdir(ROARY_DIR)
                       if os.path.isdir(os.path.join(ROARY_DIR, d))])

results = []
raw_p_values = []

print("=" * 70)
print("Fisher's Exact Test: AMR Enrichment in Core vs Accessory Genome")
print("=" * 70)

for species in species_list:
    species_cap = species.capitalize()

    # Load AMR core/accessory data
    amr_path = os.path.join(AMR_DIR, species, "amr_core_accessory.csv")
    if not os.path.exists(amr_path):
        print(f"⚠️  Missing AMR data for {species_cap}")
        continue

    amr_df = pd.read_csv(amr_path)

    # Count UNIQUE AMR genes per category
    unique_amr = amr_df.groupby("Category")["AMR_gene"].nunique()
    core_amr = unique_amr.get("Core", 0)
    acc_amr = unique_amr.get("Accessory", 0)
    unk_amr = unique_amr.get("Unknown", 0)

    # Total pangenome genes from Roary
    roary = pd.read_csv(os.path.join(ROARY_DIR, species, "gene_presence_absence.csv"), low_memory=False)
    strain_cols = roary.columns[14:]
    total_strains = len(strain_cols)
    roary["Presence_Count"] = roary[strain_cols].notnull().sum(axis=1)
    roary["Category"] = roary["Presence_Count"].apply(
        lambda x: "Core" if x >= 0.99 * total_strains else "Accessory"
    )

    total_core = (roary["Category"] == "Core").sum()
    total_acc = (roary["Category"] == "Accessory").sum()

    # Non-AMR background genes
    nonamr_core = total_core - core_amr
    nonamr_acc = total_acc - acc_amr

    # 2x2 table:
    #           Core    Accessory
    # AMR       a       b
    # Non-AMR   c       d
    a, b = core_amr, acc_amr
    c, d = nonamr_core, nonamr_acc

    or_val, p_val, ci_lo, ci_hi = fisher_with_ci(a, b, c, d)
    raw_p_values.append(p_val)

    # Direction (based on uncorrected p < 0.05 for backwards compatibility)
    if or_val > 1 and p_val < 0.05:
        direction = "Enriched in CORE"
    elif or_val < 1 and p_val < 0.05:
        direction = "Enriched in ACCESSORY"
    else:
        direction = "No significant enrichment"

    results.append({
        "Species": species_cap,
        "Core_AMR_Genes": core_amr,
        "Accessory_AMR_Genes": acc_amr,
        "Unknown_AMR_Genes": unk_amr,
        "Total_Core_Genes": total_core,
        "Total_Accessory_Genes": total_acc,
        "OR": round(or_val, 3),
        "CI_Lower": round(ci_lo, 3),
        "CI_Upper": round(ci_hi, 3),
        "p_value": p_val,
        "Significant": "Yes" if p_val < 0.05 else "No",
        "Direction": direction
    })

    print(f"\n{species_cap}")
    print(f"  AMR genes: Core={core_amr}, Accessory={acc_amr}, Unknown={unk_amr}")
    print(f"  Total genes: Core={total_core}, Accessory={total_acc}")
    print(f"  2x2 table: [{a},{b} / {c},{d}]")
    print(f"  OR = {or_val:.3f} (95% CI: {ci_lo:.3f}–{ci_hi:.3f})")
    print(f"  p = {p_val:.4f} → {direction}")

# ===============================
# APPLY BENJAMINI-HOCHBERG FDR CORRECTION
# ===============================
if results:
    fdr_q_values = benjamini_hochberg(raw_p_values)
    for i, q_val in enumerate(fdr_q_values):
        results[i]["FDR_q"] = round(q_val, 6)
        results[i]["Significant_after_FDR"] = "Yes" if q_val < 0.05 else "No"

# ===============================
# SAVE SUMMARY
# ===============================
if results:
    res_df = pd.DataFrame(results)
    out_csv = os.path.join(OUT_DIR, "fisher_core_accessory_enrichment.csv")
    res_df.to_csv(out_csv, index=False)
    print(f"\n✅ Summary saved: {out_csv}")
    print("\n📋 Final Summary:")
    print(res_df[["Species", "Core_AMR_Genes", "Accessory_AMR_Genes", "OR", "CI_Lower", "CI_Upper", "p_value", "FDR_q", "Significant_after_FDR"]].to_string(index=False))

    # ===============================
    # PLOT: Odds Ratios
    # ===============================
    plot_df = res_df.copy()
    plot_df["Significant_FDR"] = plot_df["Significant_after_FDR"] == "Yes"
    plot_df["Color"] = plot_df.apply(
        lambda row: "#F44336" if (row["OR"] < 1 and row["Significant_FDR"])
        else "#4CAF50" if (row["OR"] > 1 and row["Significant_FDR"])
        else "#9E9E9E",
        axis=1
    )
    plot_df["MarkerFaceColor"] = plot_df.apply(
        lambda row: row["Color"] if row["Significant_FDR"] else "white",
        axis=1
    )

    plt.figure(figsize=(11, 6))

    for idx, row in plot_df.iterrows():
        color = row["Color"]
        facecolor = row["MarkerFaceColor"]
        or_val = row["OR"]
        ci_lo = row["CI_Lower"]
        ci_hi = row["CI_Upper"]

        # Handle edge cases where OR == 0 or CI bounds cross
        if or_val == 0:
            # Plot a zero OR at a small positive value so it is visible on the log scale.
            # The reported OR and CI are unchanged; this is purely for visualization.
            or_plot = 0.001
            err_lo = or_plot * 0.999  # lower cap near zero; clipped at axis floor
            err_hi = ci_hi - or_plot
        else:
            or_plot = or_val
            err_lo = or_plot - ci_lo
            if err_lo < 0:
                err_lo = or_plot * 0.01
            err_hi = ci_hi - or_plot

        if err_hi < 0:
            err_hi = or_plot * 0.01

        plt.errorbar(
            x=idx, y=or_plot,
            yerr=[[err_lo], [err_hi]],
            fmt='o', color=color, ecolor=color, capsize=5, capthick=2, markersize=10,
            markerfacecolor=facecolor,
            markeredgecolor='black', markeredgewidth=1.2
        )

    plt.xticks(range(len(plot_df)), plot_df["Species"], rotation=30, ha='right')
    plt.axhline(y=1, color='black', linestyle='--', linewidth=1)
    plt.ylabel("Odds Ratio (AMR in Core vs Accessory)", fontsize=12)
    plt.xlabel("Species", fontsize=12)
    plt.title("AMR Gene Enrichment: Core vs Accessory Genome\n(OR > 1 = Core-enriched; OR < 1 = Accessory-enriched; FDR q < 0.05 = filled)",
              fontsize=13, weight='bold')
    plt.yscale('log')
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    # Custom legend
    from matplotlib.patches import Patch
    import matplotlib.lines as mlines
    legend_elements = [
        Patch(facecolor="#F44336", edgecolor='black', label="Accessory-enriched (FDR q<0.05)"),
        Patch(facecolor="#4CAF50", edgecolor='black', label="Core-enriched (FDR q<0.05)"),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white',
                      markeredgecolor='black', markersize=10, label="Not significant (FDR q≥0.05)")
    ]
    plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, "fisher_core_accessory_oddsratios.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n✅ Plot saved: {plot_path}")

print("\n🎉 Fisher's exact test analysis complete!")
