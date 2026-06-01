import os
import pandas as pd

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AMR_DIR = os.path.join(BASE_DIR, "amrfinder_results")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "01_build_amr_matrix")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# STEP 1: Collect AMR genes per strain
# ===============================
species_list = sorted([d for d in os.listdir(AMR_DIR)
                       if os.path.isdir(os.path.join(AMR_DIR, d))])

all_genes = set()
data = {}  # strain_name -> set of genes

for species in species_list:
    species_path = os.path.join(AMR_DIR, species)
    species_cap = species.capitalize()

    files = sorted([f for f in os.listdir(species_path) if f.endswith("_amr.tsv")])

    for f in files:
        base = f.replace("_amr.tsv", "")
        strain = f"{species_cap}_{base}"
        filepath = os.path.join(species_path, f)

        try:
            df = pd.read_csv(filepath, sep="\t")
        except Exception as e:
            print(f"⚠️  Skipping {filepath}: {e}")
            continue

        if "Element symbol" not in df.columns:
            print(f"⚠️  No 'Element symbol' column in {filepath}")
            continue

        genes = set(df["Element symbol"].dropna())
        data[strain] = genes
        all_genes.update(genes)

# ===============================
# STEP 2: Build binary matrix
# ===============================
all_genes = sorted(all_genes)
matrix = []

for strain, genes in data.items():
    row = [1 if gene in genes else 0 for gene in all_genes]
    matrix.append(row)

df_matrix = pd.DataFrame(matrix, columns=all_genes)
df_matrix.insert(0, "Strain", list(data.keys()))

# Save full matrix
out_matrix = os.path.join(OUT_DIR, "amr_matrix.csv")
df_matrix.to_csv(out_matrix, index=False)
print(f"✅ AMR matrix saved: {out_matrix}")
print(f"   Strains: {len(data)} | Genes: {len(all_genes)}")

# ===============================
# STEP 3: Save "cleaned" matrix (remove all-zero gene columns)
# ===============================
gene_cols = [c for c in df_matrix.columns if c != "Strain"]
nonzero_cols = [c for c in gene_cols if df_matrix[c].sum() > 0]
df_clean = df_matrix[["Strain"] + nonzero_cols]
clean_path = os.path.join(OUT_DIR, "cleaned_amr.csv")
df_clean.to_csv(clean_path, index=False)
print(f"✅ Cleaned AMR matrix saved: {clean_path}")

# ===============================
# STEP 4: Split by species
# ===============================
for species in species_list:
    species_cap = species.capitalize()
    sub = df_matrix[df_matrix["Strain"].str.startswith(species_cap)]
    if not sub.empty:
        out_path = os.path.join(OUT_DIR, f"amr_{species}.csv")
        sub.to_csv(out_path, index=False)
        print(f"   📁 {species}: {len(sub)} strains → {out_path}")

print("\n🎉 Done building AMR matrices!")
