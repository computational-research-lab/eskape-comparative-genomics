import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_results")
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
CHECKM2_DIR = os.path.join(BASE_DIR, "checkm2_results")
OUT_DIR = os.path.join(ANALYSIS_DIR, "04_final_summary")
os.makedirs(OUT_DIR, exist_ok=True)


def parse_roary_stats(species):
    """Read Roary summary_statistics.txt and return key numbers."""
    path = os.path.join(ROARY_DIR, species, "summary_statistics.txt")
    stats = {}
    if not os.path.exists(path):
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                key = parts[0].strip().lower().replace(" ", "_")
                try:
                    stats[key] = int(parts[-1].strip())
                except ValueError:
                    pass
    return stats


# ===============================
# STEP 1: Load MDR data
# ===============================
mdr_path = os.path.join(ANALYSIS_DIR, "02_mdr_analysis", "mdr_from_amrfinder_class.csv")
mdr_df = pd.read_csv(mdr_path)
mdr_df["Species"] = mdr_df["Strain"].apply(lambda x: x.split("_")[0].lower())

mdr_summary = mdr_df.groupby("Species")["MDR"].value_counts().unstack(fill_value=0)
mdr_summary["Total"] = mdr_summary.sum(axis=1)
mdr_summary["MDR_%"] = (mdr_summary.get("Yes", 0) / mdr_summary["Total"]) * 100
avg_classes = mdr_df.groupby("Species")["Total_Classes"].mean()

# ===============================
# STEP 2: Parse pangenome stats dynamically from Roary
# ===============================
species_list = sorted([d for d in os.listdir(ROARY_DIR)
                       if os.path.isdir(os.path.join(ROARY_DIR, d))])

pangenome_data = []
for sp in species_list:
    stats = parse_roary_stats(sp)
    pangenome_data.append({
        "Species": sp,
        "Core_Genes": stats.get("core_genes", 0),
        "Soft_Core_Genes": stats.get("soft_core_genes", 0),
        "Shell_Genes": stats.get("shell_genes", 0),
        "Cloud_Genes": stats.get("cloud_genes", 0),
        "Total_Genes": stats.get("total_genes", 0)
    })

pangenome = pd.DataFrame(pangenome_data)
# Accessory = everything not strict core (soft core + shell + cloud)
pangenome["Accessory_Genes"] = (
    pangenome["Soft_Core_Genes"] + pangenome["Shell_Genes"] + pangenome["Cloud_Genes"]
)
pangenome = pangenome.set_index("Species")

# ===============================
# STEP 3: Load core/accessory AMR summary (optional)
# ===============================
ca_path = os.path.join(ANALYSIS_DIR, "03_core_accessory_amr", "final_summary_core_accessory.csv")
if os.path.exists(ca_path):
    ca_df = pd.read_csv(ca_path)
    ca_df["Species"] = ca_df["Species"].str.lower()
    ca_df = ca_df.set_index("Species")
else:
    ca_df = pd.DataFrame()

# ===============================
# STEP 4: Build final merged summary
# ===============================
final = pd.concat([mdr_summary["MDR_%"], avg_classes, pangenome], axis=1)
final.columns = [
    "MDR_%", "Avg_AMR_Classes", "Core_Genes", "Soft_Core_Genes",
    "Shell_Genes", "Cloud_Genes", "Total_Genes", "Accessory_Genes"
]

if not ca_df.empty:
    final = final.join(ca_df[["Core", "Accessory"]], rsuffix="_AMR")
    final = final.rename(columns={"Core": "Core_AMR_Hits", "Accessory": "Accessory_AMR_Hits"})

final = final.reset_index()
out_summary = os.path.join(OUT_DIR, "final_species_summary.csv")
final.to_csv(out_summary, index=False)

print("✅ Final species summary saved!")
print(final.to_string(index=False))

# ===============================
# STEP 5: Statistical correlation
# ===============================
if len(final) > 2:
    r_pearson, p_pearson = pearsonr(final["Accessory_Genes"], final["MDR_%"])
    r_spear, p_spear = spearmanr(final["Accessory_Genes"], final["MDR_%"])
    print(f"\n📊 Pearson correlation (Accessory vs MDR%):  r = {r_pearson:.3f}, p = {p_pearson:.3f}")
    print(f"📊 Spearman correlation (Accessory vs MDR%): ρ = {r_spear:.3f}, p = {p_spear:.3f}")
else:
    print("\n⚠️  Not enough species for correlation test (need >2).")

# ===============================
# STEP 6: Scatter plot (Accessory vs MDR%)
# ===============================
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=final,
    x="Accessory_Genes",
    y="MDR_%",
    hue="Species",
    s=150,
    edgecolor="black",
    linewidth=1.2
)

# Add trend line if enough points
if len(final) > 2:
    z = np.polyfit(final["Accessory_Genes"], final["MDR_%"], 1)
    p = np.poly1d(z)
    plt.plot(final["Accessory_Genes"], p(final["Accessory_Genes"]),
             linestyle="--", color="gray", alpha=0.7, label="Trend")

plt.xlabel("Accessory Genome Size (Genes)", fontsize=12)
plt.ylabel("MDR Prevalence (%)", fontsize=12)
plt.title("Relationship between Genome Plasticity and MDR", fontsize=14, weight='bold')
plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

plot_path = os.path.join(OUT_DIR, "final_mdr_vs_accessory_clean.png")
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"\n✅ Scatter plot saved: {plot_path}")

# ===============================
# STEP 7: CheckM2 quality summary
# ===============================
checkm2_data = []
for sp in species_list:
    path = os.path.join(CHECKM2_DIR, sp, "quality_report.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "Species", sp.capitalize())
        checkm2_data.append(df)

if checkm2_data:
    checkm2_combined = pd.concat(checkm2_data, ignore_index=True)
    checkm2_combined.to_csv(os.path.join(OUT_DIR, "checkm2_summary.csv"), index=False)
    print(f"✅ CheckM2 per-strain summary saved.")

    # Species-level averages
    numeric_cols = ["Completeness", "Contamination", "Genome_Size", "Contig_N50", "Total_Contigs"]
    available_cols = [c for c in numeric_cols if c in checkm2_combined.columns]
    if available_cols:
        checkm2_species = checkm2_combined.groupby("Species")[available_cols].mean().round(2)
        checkm2_species.to_csv(os.path.join(OUT_DIR, "checkm2_species_summary.csv"))
        print(f"✅ CheckM2 species summary saved.")
        print("\n📋 CheckM2 Species Averages:")
        print(checkm2_species.to_string())
else:
    print("⚠️  No CheckM2 reports found.")

print("\n🎉 Final summary complete!")
