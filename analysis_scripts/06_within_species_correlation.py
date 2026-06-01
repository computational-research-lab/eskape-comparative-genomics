import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
MDR_PATH = os.path.join(BASE_DIR, "analysis_results", "02_mdr_analysis", "mdr_from_amrfinder_class.csv")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "06_within_species_correlation")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# LOAD MDR DATA
# ===============================
mdr_df = pd.read_csv(MDR_PATH)
mdr_df["Species"] = mdr_df["Strain"].apply(lambda x: x.split("_")[0])
mdr_df["Strain_Base"] = mdr_df["Strain"].apply(lambda x: "_".join(x.split("_")[1:]))

# ===============================
# PROCESS EACH SPECIES
# ===============================
species_list = sorted([d for d in os.listdir(ROARY_DIR)
                       if os.path.isdir(os.path.join(ROARY_DIR, d))])

results = []
plot_data = []

for species in species_list:
    species_cap = species.capitalize()
    roary_path = os.path.join(ROARY_DIR, species, "gene_presence_absence.csv")

    if not os.path.exists(roary_path):
        print(f"⚠️  Roary file not found for {species}")
        continue

    roary = pd.read_csv(roary_path, low_memory=False)
    strain_cols = roary.columns[14:]
    total_strains = len(strain_cols)

    # Classify core vs accessory
    roary["Presence_Count"] = roary[strain_cols].notnull().sum(axis=1)
    roary["Category"] = roary["Presence_Count"].apply(
        lambda x: "Core" if x >= 0.99 * total_strains else "Accessory"
    )

    # Count accessory genes per strain
    accessory_counts = {}
    for strain in strain_cols:
        accessory_counts[strain] = roary[roary["Category"] == "Accessory"][strain].notnull().sum()

    acc_df = pd.DataFrame(list(accessory_counts.items()), columns=["Strain_Base", "Accessory_Genes"])

    # Merge with MDR data for this species
    species_mdr = mdr_df[mdr_df["Species"] == species_cap][["Strain_Base", "Total_Classes", "MDR"]].copy()
    merged = pd.merge(acc_df, species_mdr, on="Strain_Base", how="inner")

    if len(merged) < 3:
        print(f"⚠️  {species}: insufficient matched strains ({len(merged)})")
        continue

    # Correlations
    r_pearson, p_pearson = pearsonr(merged["Accessory_Genes"], merged["Total_Classes"])
    r_spear, p_spear = spearmanr(merged["Accessory_Genes"], merged["Total_Classes"])

    results.append({
        "Species": species_cap,
        "N_Strains": len(merged),
        "Pearson_r": round(r_pearson, 3),
        "Pearson_p": round(p_pearson, 4),
        "Spearman_rho": round(r_spear, 3),
        "Spearman_p": round(p_spear, 4),
        "Mean_Accessory": round(merged["Accessory_Genes"].mean(), 1),
        "Mean_AMR_Classes": round(merged["Total_Classes"].mean(), 2)
    })

    # Add species label for combined plot
    merged["Species"] = species_cap
    plot_data.append(merged)

    print(f"✅ {species_cap:15s} | n={len(merged):2d} | Pearson r={r_pearson:+.3f} (p={p_pearson:.3f}) | "
          f"Spearman ρ={r_spear:+.3f} (p={p_spear:.3f})")

# ===============================
# SAVE SUMMARY TABLE
# ===============================
if results:
    summary_df = pd.DataFrame(results)
    summary_path = os.path.join(OUT_DIR, "within_species_correlation_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\n✅ Summary saved: {summary_path}")
    print("\n📋 Within-Species Correlation Summary:")
    print(summary_df.to_string(index=False))

# ===============================
# COMBINED FACETED PLOT
# ===============================
if plot_data:
    combined = pd.concat(plot_data, ignore_index=True)

    g = sns.FacetGrid(combined, col="Species", col_wrap=3, sharex=False, sharey=False,
                      height=3.5, aspect=1.1)
    g.map_dataframe(sns.scatterplot, x="Accessory_Genes", y="Total_Classes",
                    hue="MDR", palette={"Yes": "#F44336", "No": "#4CAF50"},
                    s=80, edgecolor="black", linewidth=0.5)

    # Add regression line per facet
    def regline(data, color, **kwargs):
        if len(data) > 2:
            x = data["Accessory_Genes"]
            y = data["Total_Classes"]
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            plt.plot(x.sort_values(), p(x.sort_values()), linestyle="--", color="gray", alpha=0.7)

    g.map_dataframe(regline)

    g.set_axis_labels("Accessory Genome Size (genes)", "AMR Classes")
    g.set_titles(col_template="{col_name}", size=11, weight='bold')
    g.add_legend(title="MDR", adjust_subtitles=True)
    g.tight_layout()

    plot_path = os.path.join(OUT_DIR, "within_species_accessory_vs_amr.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n✅ Faceted plot saved: {plot_path}")

print("\n🎉 Within-species correlation analysis complete!")
