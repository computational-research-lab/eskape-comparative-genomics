import os
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AMR_DIR = os.path.join(BASE_DIR, "amrfinder_results")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "02_mdr_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# STEP 1: Extract drug classes per strain
# ===============================
species_list = sorted([d for d in os.listdir(AMR_DIR)
                       if os.path.isdir(os.path.join(AMR_DIR, d))])

data = []

for species in species_list:
    species_path = os.path.join(AMR_DIR, species)
    species_cap = species.capitalize()
    files = [f for f in os.listdir(species_path) if f.endswith("_amr.tsv")]

    for f in files:
        base = f.replace("_amr.tsv", "")
        strain = f"{species_cap}_{base}"
        filepath = os.path.join(species_path, f)

        try:
            df = pd.read_csv(filepath, sep="\t")
        except Exception as e:
            print(f"⚠️  Error reading {filepath}: {e}")
            continue

        # Detect class column (handles different AMRFinder versions)
        class_col = None
        for col in ["Class", "Drug Class", "class", "Drug class"]:
            if col in df.columns:
                class_col = col
                break

        if class_col is None:
            print(f"⚠️  No class column in {filepath}")
            continue

        classes = df[class_col].dropna().unique()
        row = {"Strain": strain}
        for c in classes:
            row[c] = 1
        data.append(row)

# ===============================
# STEP 2: Build class dataframe
# ===============================
class_df = pd.DataFrame(data).fillna(0)
class_df = class_df.set_index("Strain")
class_df = class_df.astype(int)

# Count total classes per strain
class_df["Total_Classes"] = class_df.sum(axis=1)

# MDR definition: resistance to >= 3 antimicrobial classes
class_df["MDR"] = class_df["Total_Classes"].apply(lambda x: "Yes" if x >= 3 else "No")

# ===============================
# STEP 3: Save
# ===============================
out_csv = os.path.join(OUT_DIR, "mdr_from_amrfinder_class.csv")
class_df.to_csv(out_csv)
print(f"✅ MDR data saved: {out_csv}")
print(f"   Total strains: {len(class_df)}")
print(f"   Avg classes per strain: {class_df['Total_Classes'].mean():.2f}")
print("\n🔍 MDR Summary:")
print(class_df["MDR"].value_counts())

# ===============================
# STEP 4: Plot MDR distribution
# ===============================
df_plot = class_df.reset_index()
df_plot["Species"] = df_plot["Strain"].apply(lambda x: x.split("_")[0])

summary = df_plot.groupby(["Species", "MDR"]).size().unstack(fill_value=0)
summary.columns = ["Non-MDR", "MDR"]

plt.figure(figsize=(10, 6))
ax = summary.plot(
    kind="bar",
    stacked=True,
    color=["#4CAF50", "#F44336"],
    edgecolor="black",
    width=0.7
)

plt.title("MDR Distribution Across Species", fontsize=14, weight='bold')
plt.xlabel("Species", fontsize=12)
plt.ylabel("Number of Strains", fontsize=12)
plt.xticks(rotation=30, ha='right')
plt.legend(title="MDR Status", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.5)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
out_png = os.path.join(OUT_DIR, "mdr_amrfinder_plot.png")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ MDR plot saved: {out_png}")
