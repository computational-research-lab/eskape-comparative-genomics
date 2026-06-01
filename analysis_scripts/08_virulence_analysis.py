import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from matplotlib.colors import ListedColormap

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROKKA_DIR = os.path.join(BASE_DIR, "prokka_out")
DB_PATH = os.path.join(BASE_DIR, "databases", "vfdb.dmnd")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "08_virulence_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# Filtering thresholds
MIN_PIDENT = 80.0
MIN_COV = 70.0
MAX_EVAL = 1e-10

# Full path to DIAMOND executable
DIAMOND = "/Users/admin/miniconda3/envs/eskaPE_env/bin/diamond"

# ===============================
# VERIFY DATABASE
# ===============================
if not os.path.exists(DB_PATH):
    print(f"❌ DIAMOND database not found: {DB_PATH}")
    print("   Please run: diamond makedb --in VFDB_setB_pro.fas --db vfdb.dmnd")
    exit(1)

if not os.path.exists(DIAMOND):
    print(f"❌ DIAMOND executable not found: {DIAMOND}")
    exit(1)

# ===============================
# STEP 1: Run DIAMOND BLASTP for all strains
# ===============================
species_list = sorted([d for d in os.listdir(PROKKA_DIR)
                       if os.path.isdir(os.path.join(PROKKA_DIR, d))])

all_hits = {}
all_vfs = set()

print("=" * 60)
print("Running DIAMOND BLASTP against VFDB")
print("=" * 60)

for species in species_list:
    species_cap = species.capitalize()
    faa_pattern = os.path.join(PROKKA_DIR, species, "*", "*.faa")
    faa_files = glob(faa_pattern)
    print(f"\n{species_cap}: {len(faa_files)} genomes")

    for faa in sorted(faa_files):
        strain = os.path.basename(faa).replace(".faa", "")
        strain_full = f"{species_cap}_{strain}"
        out_tsv = os.path.join(OUT_DIR, f"{strain_full}_vfdb.tsv")

        if not os.path.exists(out_tsv):
            cmd = [
                DIAMOND, "blastp",
                "-d", DB_PATH,
                "-q", faa,
                "-o", out_tsv,
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle", "qlen", "slen",
                "--max-target-seqs", "1",
                "--evalue", str(MAX_EVAL),
                "--threads", "4"
            ]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                print(f"  ⚠️  DIAMOND failed for {strain_full}: {e}")
                continue

        # Parse results
        if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
            all_hits[strain_full] = set()
            continue

        try:
            cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle", "qlen", "slen"]
            df = pd.read_csv(out_tsv, sep="\t", names=cols, header=None)
        except Exception as e:
            print(f"  ⚠️  Parse error for {strain_full}: {e}")
            all_hits[strain_full] = set()
            continue

        # Apply filters
        df["cov"] = (df["length"] / df["slen"]) * 100
        filtered = df[(df["pident"] >= MIN_PIDENT) & (df["cov"] >= MIN_COV)]

        # Extract VF gene name from stitle (e.g., "VFG026115(gb|WP_006484930) (wcbP) ...")
        # Use the gene symbol in parentheses if available, otherwise sseqid
        vf_genes = set()
        for _, row in filtered.iterrows():
            stitle = str(row["stitle"])
            # Try to extract gene name like (wcbP) or (blaTEM-1)
            import re
            m = re.search(r'\(([^)]+)\)', stitle)
            if m:
                gene = m.group(1).strip()
                # Skip accession numbers that look like WP_123 or gb|...
                if not re.match(r'^(WP_|gb\|)', gene):
                    vf_genes.add(gene)
                else:
                    vf_genes.add(row["sseqid"].split("(")[0])
            else:
                vf_genes.add(row["sseqid"].split("(")[0])

        all_hits[strain_full] = vf_genes
        all_vfs.update(vf_genes)

# ===============================
# STEP 2: Build binary matrix
# ===============================
all_vfs = sorted(all_vfs)
matrix = []
for strain, vfs in all_hits.items():
    row = [1 if vf in vfs else 0 for vf in all_vfs]
    matrix.append(row)

vf_matrix = pd.DataFrame(matrix, columns=all_vfs)
vf_matrix.insert(0, "Strain", list(all_hits.keys()))
vf_matrix.to_csv(os.path.join(OUT_DIR, "virulence_matrix.csv"), index=False)
print(f"\n✅ Virulence matrix saved: {len(all_hits)} strains × {len(all_vfs)} VFs")

# ===============================
# STEP 3: Species summary
# ===============================
vf_matrix["Species"] = vf_matrix["Strain"].apply(lambda x: x.split("_")[0])
species_summary = vf_matrix.groupby("Species")[all_vfs].sum().T
species_summary.to_csv(os.path.join(OUT_DIR, "virulence_species_summary.csv"))

# Strain-level summary
strain_counts = vf_matrix[all_vfs].sum(axis=1)
summary_df = pd.DataFrame({
    "Strain": vf_matrix["Strain"],
    "Species": vf_matrix["Species"],
    "VF_Count": strain_counts
})
summary_df.to_csv(os.path.join(OUT_DIR, "virulence_strain_summary.csv"), index=False)

print("\n📋 Average VF genes per strain:")
print(summary_df.groupby("Species")["VF_Count"].agg(["mean", "std", "min", "max"]).round(1))

# ===============================
# STEP 4: Plot species-wise VF counts
# ===============================
plt.figure(figsize=(10, 6))
species_order = sorted(summary_df["Species"].unique())
sns.boxplot(data=summary_df, x="Species", y="VF_Count", order=species_order, palette="Set2")
sns.stripplot(data=summary_df, x="Species", y="VF_Count", order=species_order,
              color="black", alpha=0.5, size=4)
plt.title("Virulence Factor Gene Count per Strain", fontsize=14, weight='bold')
plt.xlabel("Species", fontsize=12)
plt.ylabel("Number of VF Genes", fontsize=12)
plt.xticks(rotation=30, ha='right')
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "virulence_count_boxplot.png"), dpi=300, bbox_inches='tight')
plt.close()
print(f"✅ Boxplot saved")

# ===============================
# STEP 5: Heatmap (top VFs overall)
# ===============================
vf_totals = vf_matrix[all_vfs].sum().sort_values(ascending=False)
top_vfs = vf_totals.head(40).index.tolist()

sub = vf_matrix[["Strain"] + top_vfs].copy()
sub = sub.set_index("Strain")

fig_width = max(14, len(top_vfs) * 0.35)
fig_height = max(10, len(sub) * 0.25)

cmap = ListedColormap(["#ffffcc", "#d62728"])
g = sns.clustermap(
    sub,
    cmap=cmap,
    figsize=(fig_width, fig_height),
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"ticks": [0, 1]},
    linewidths=0.3,
    linecolor="gray",
    row_cluster=True,
    col_cluster=True
)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=7)
g.ax_cbar.set_yticklabels(["Absent", "Present"])
g.ax_cbar.set_title("VF", fontsize=10)
plt.title("Top 40 Virulence Factors Across Strains", y=1.02)
plt.savefig(os.path.join(OUT_DIR, "virulence_heatmap_top40.png"), dpi=300, bbox_inches='tight')
plt.close()
print(f"✅ Heatmap saved")

print("\n🎉 Virulence analysis complete!")
