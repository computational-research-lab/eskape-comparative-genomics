import os
import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

# ===============================
# CONFIG
# ===============================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROARY_DIR = os.path.join(BASE_DIR, "roary_out")
OUT_DIR = os.path.join(BASE_DIR, "analysis_results", "09_phylogenetic_tree")
os.makedirs(OUT_DIR, exist_ok=True)


def alignment_to_array(alignment):
    """Convert BioPython alignment to numpy char array."""
    seqs = []
    names = []
    for record in alignment:
        names.append(record.id)
        seqs.append(list(str(record.seq)))
    return np.array(seqs), names


def fast_distance_matrix(arr):
    """Compute pairwise Hamming distances on SNP sites only."""
    # Find polymorphic columns (SNPs)
    # Ignore gaps ('-') and 'N' by treating them as matching anything
    # Actually, for core genome alignment from Roary, gaps should be minimal
    polymorphic = np.any(arr != arr[0, :], axis=0)
    snp_arr = arr[:, polymorphic]

    n = snp_arr.shape[0]
    dm = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            # Count mismatches, ignoring gaps in either sequence
            mask = (snp_arr[i, :] != '-') & (snp_arr[j, :] != '-')
            if mask.sum() == 0:
                dist = 0.0
            else:
                mismatches = np.sum(snp_arr[i, mask] != snp_arr[j, mask])
                dist = mismatches / mask.sum()
            dm[i, j] = dist
            dm[j, i] = dist

    return dm, snp_arr.shape[1]


# ===============================
# BUILD NJ TREES FROM ROARY ALIGNMENTS
# ===============================
species_list = sorted([d for d in os.listdir(ROARY_DIR)
                       if os.path.isdir(os.path.join(ROARY_DIR, d))])

print("=" * 60)
print("Building Core Genome Neighbor-Joining Trees")
print("=" * 60)

tree_results = []

for species in species_list:
    species_cap = species.capitalize()
    aln_path = os.path.join(ROARY_DIR, species, "core_gene_alignment.aln")

    if not os.path.exists(aln_path):
        print(f"⚠️  {species_cap}: alignment not found, skipping")
        continue

    tree_out = os.path.join(OUT_DIR, f"{species}_core_tree.nwk")

    print(f"\n🌲 {species_cap}: loading alignment...")
    try:
        alignment = AlignIO.read(aln_path, "fasta")
    except Exception as e:
        print(f"   ❌ Failed to parse alignment: {e}")
        continue

    n_strains = len(alignment)
    aln_len = alignment.get_alignment_length()
    print(f"   Alignment: {aln_len:,} bp × {n_strains} strains")

    # Convert to array and compute distances
    print(f"   Extracting SNPs and computing distances...")
    arr, names = alignment_to_array(alignment)
    dm_array, n_snps = fast_distance_matrix(arr)
    print(f"   SNP sites: {n_snps:,} | Avg pairwise distance: {dm_array[np.triu_indices_from(dm_array, k=1)].mean():.4f}")

    # Build BioPython DistanceMatrix
    dm_list = []
    for i in range(n_strains):
        row = [dm_array[i, j] for j in range(i + 1)]
        dm_list.append(row)

    dm = DistanceMatrix(names, dm_list)

    # Build NJ tree
    print(f"   Building NJ tree...")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Save tree
    Phylo.write(tree, tree_out, "newick")
    print(f"   ✅ Tree saved: {tree_out}")

    tree_results.append({
        "Species": species_cap,
        "Strains": n_strains,
        "Alignment_Bp": aln_len,
        "SNP_Sites": n_snps,
        "Tree": tree_out
    })

# ===============================
# SUMMARY
# ===============================
if tree_results:
    summary_df = pd.DataFrame(tree_results)
    summary_df.to_csv(os.path.join(OUT_DIR, "tree_summary.csv"), index=False)
    print("\n📋 Tree Summary:")
    print(summary_df.to_string(index=False))
    print(f"\n✅ All trees saved to: {OUT_DIR}")
else:
    print("\n⚠️  No trees were generated.")

print("\n🎉 Phylogenetic tree construction complete!")
