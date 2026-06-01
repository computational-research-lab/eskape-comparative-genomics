import os
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
AMR_DIR = os.path.join(BASE_DIR, "amrfinder_results")
PROKKA_DIR = os.path.join(BASE_DIR, "prokka_out")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "03_core_accessory_amr")
os.makedirs(OUT_DIR, exist_ok=True)

# ===============================
# HELPERS
# ===============================
def load_gff(gff_file):
    """Parse Prokka GFF and return CDS records with locus tags."""
    records = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            contig, _, feature, start, end, _, strand, _, attributes = parts
            if feature != "CDS":
                continue

            attr_dict = {}
            for item in attributes.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v

            locus = attr_dict.get("locus_tag")
            if locus:
                records.append({
                    "contig": contig,
                    "start": int(start),
                    "end": int(end),
                    "locus_tag": locus
                })
    return pd.DataFrame(records)


def parse_roary_summary(species):
    """Parse Roary summary_statistics.txt into a dict."""
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
# MAIN PROCESSING FUNCTION
# ===============================
def process_species(species):
    print(f"\n{'='*50}\nProcessing: {species.capitalize()}\n{'='*50}")

    roary_path = os.path.join(ROARY_DIR, species, "gene_presence_absence.csv")
    if not os.path.exists(roary_path):
        print(f"❌ Roary file not found: {roary_path}")
        return None

    roary = pd.read_csv(roary_path, low_memory=False)
    strain_cols = roary.columns[14:]
    total_strains = len(strain_cols)

    # Define core (>=99% strains) vs accessory
    roary["Presence_Count"] = roary[strain_cols].notnull().sum(axis=1)
    roary["Category"] = roary["Presence_Count"].apply(
        lambda x: "Core" if x >= 0.99 * total_strains else "Accessory"
    )

    # Build locus -> category map
    locus_map = {}
    for _, row in roary.iterrows():
        category = row["Category"]
        for strain in strain_cols:
            val = row[strain]
            if pd.notnull(val):
                loci = str(val).split("\t")
                for locus in loci:
                    locus_map[locus.strip()] = category

    print(f"✔ Roary genes: {len(roary)} | Locus map: {len(locus_map)} | Strains: {total_strains}")

    # Load GFF files for this species
    gff_pattern = os.path.join(PROKKA_DIR, species, "*", "*.gff")
    gff_files = glob(gff_pattern)
    print(f"✔ GFF files found: {len(gff_files)}")

    if not gff_files:
        print("❌ No GFF files found — skipping")
        return None

    gff_list = []
    for g in gff_files:
        df = load_gff(g)
        if not df.empty:
            gff_list.append(df)

    if not gff_list:
        print("❌ GFF parsing failed — no CDS entries")
        return None

    gff_df = pd.concat(gff_list, ignore_index=True)
    print(f"✔ Total CDS records: {len(gff_df)}")

    # Load AMR files
    amr_pattern = os.path.join(AMR_DIR, species, "*_amr.tsv")
    amr_files = glob(amr_pattern)
    print(f"✔ AMR files found: {len(amr_files)}")

    if not amr_files:
        print("❌ No AMR files found — skipping")
        return None

    results = []
    species_cap = species.capitalize()

    for file in amr_files:
        base = os.path.basename(file).replace("_amr.tsv", "")
        strain = f"{species_cap}_{base}"

        try:
            df = pd.read_csv(file, sep="\t")
        except Exception as e:
            print(f"⚠️  Error reading {file}: {e}")
            continue

        for _, row in df.iterrows():
            try:
                contig = row["Contig id"]
                start = int(row["Start"])
                end = int(row["Stop"])
                gene = row["Element symbol"]
            except Exception:
                continue

            # Match AMR coordinates to CDS in GFF
            match = gff_df[
                (gff_df["contig"] == contig) &
                (gff_df["start"] <= start) &
                (gff_df["end"] >= end)
            ]

            if not match.empty:
                locus = match.iloc[0]["locus_tag"]
                category = locus_map.get(locus, "Unknown")
            else:
                locus = "NA"
                category = "Unknown"

            results.append([strain, gene, locus, category])

    res_df = pd.DataFrame(results, columns=["Strain", "AMR_gene", "Locus_tag", "Category"])
    print("\n📊 Category counts:")
    print(res_df["Category"].value_counts().to_string())

    # Save per-species output
    out_dir = os.path.join(OUT_DIR, species)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "amr_core_accessory.csv")
    res_df.to_csv(out_path, index=False)
    print(f"✔ Saved: {out_path}")

    return res_df


# ===============================
# RUN FOR ALL SPECIES
# ===============================
species_list = sorted([d for d in os.listdir(ROARY_DIR)
                       if os.path.isdir(os.path.join(ROARY_DIR, d))])

all_results = []
for sp in species_list:
    df = process_species(sp)
    if df is not None:
        all_results.append(df)

# ===============================
# COMBINE & SUMMARIZE
# ===============================
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined_path = os.path.join(OUT_DIR, "amr_core_accessory_all.csv")
    combined.to_csv(combined_path, index=False)
    print(f"\n✅ Combined output saved: {combined_path}")

    # Per-strain core/accessory counts
    strain_counts = combined.groupby("Strain")["Category"].value_counts().unstack(fill_value=0)
    strain_counts = strain_counts.reindex(columns=["Core", "Accessory", "Unknown"], fill_value=0)
    strain_counts["Species"] = strain_counts.index.str.split("_").str[0]

    # Per-species summary
    species_summary = strain_counts.groupby("Species")[["Core", "Accessory", "Unknown"]].sum().reset_index()
    summary_path = os.path.join(OUT_DIR, "final_summary_core_accessory.csv")
    species_summary.to_csv(summary_path, index=False)
    print(f"✅ Species summary saved: {summary_path}")
    print("\n📋 Core/Accessory Summary:")
    print(species_summary.to_string(index=False))

    # ===============================
    # PLOT: Core vs Accessory AMR
    # ===============================
    plt.figure(figsize=(10, 6))
    melted = species_summary.melt(id_vars="Species", var_name="Category", value_name="Count")
    palette = {"Core": "#4CAF50", "Accessory": "#F44336", "Unknown": "#9E9E9E"}

    sns.barplot(data=melted, x="Species", y="Count", hue="Category", palette=palette)
    plt.title("Core vs Accessory AMR Genes Across Species", fontsize=14, weight='bold')
    plt.xlabel("Species", fontsize=12)
    plt.ylabel("Number of AMR Gene Hits", fontsize=12)
    plt.xticks(rotation=30, ha='right')
    plt.legend(title="Genome Category", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()

    plot_path = os.path.join(OUT_DIR, "core_vs_accessory_amr_genes.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Plot saved: {plot_path}")

print("\n🎉 Core/Accessory AMR analysis complete!")
