import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "05_heatmaps")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# LOAD MATRIX
# ===============================
df = pd.read_csv(os.path.join(BASE_DIR, "analysis_results", "01_build_amr_matrix", "amr_matrix.csv"))
df["Species"] = df["Strain"].apply(lambda x: x.split("_")[0])

# Binary colormap: light yellow = absent, purple = present
cmap = ListedColormap(["#ffffcc", "#6a0dad"])

# ===============================
# GENERATE HEATMAPS PER SPECIES
# ===============================
for species in sorted(df["Species"].unique()):
    print(f"\nProcessing {species}...")

    sub = df[df["Species"] == species].copy()
    sub = sub.drop(columns=["Species"])
    sub = sub.set_index("Strain")

    # Remove rare genes (present in <= 2 strains)
    sub = sub.loc[:, sub.sum() > 2]

    if sub.shape[1] == 0:
        print("  ⚠️  No genes after filtering, skipping")
        continue

    # Keep top 30 most frequent genes for readability
    top_genes = sub.sum().sort_values(ascending=False).head(30).index
    sub = sub[top_genes]

    # Dynamic figure sizing
    fig_width = max(12, sub.shape[1] * 0.4)
    fig_height = max(10, sub.shape[0] * 0.3)

    g = sns.clustermap(
        sub,
        cmap=cmap,
        figsize=(fig_width, fig_height),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"ticks": [0, 1]},
        linewidths=0.5,
        linecolor="gray"
    )

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=9)
    g.ax_cbar.set_yticklabels(["Absent", "Present"])
    g.ax_cbar.set_title("AMR", fontsize=10)
    plt.title(f"{species} AMR Heatmap", y=1.05)

    out_path = os.path.join(OUT_DIR, f"{species}_heatmap.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ Saved: {out_path} ({sub.shape[0]} strains × {sub.shape[1]} genes)")

print("\n🎉 All heatmaps generated!")
