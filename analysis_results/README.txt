ANALYSIS RESULTS DIRECTORY STRUCTURE
=====================================

This folder is organized by analysis script. Each subfolder contains ALL outputs
produced by the corresponding script, making it easy to trace where each result
came from.

------------------------------------------------------------------------------
01_build_amr_matrix/
------------------------------------------------------------------------------
Script: 01_build_amr_matrix.py
Input: amrfinder_results/*/*_amr.tsv
Description: Builds binary presence/absence matrix of AMR genes across all 120 strains

Files:
  - amr_matrix.csv              Full binary matrix (120 strains x 253 genes)
  - cleaned_amr.csv             Same as above with all-zero columns removed
  - amr_acinetobacter.csv       Species subset (20 strains)
  - amr_enterobacter.csv        Species subset (20 strains)
  - amr_enterococcus.csv        Species subset (20 strains)
  - amr_klebsiella.csv          Species subset (20 strains)
  - amr_pseudomonas.csv         Species subset (20 strains)
  - amr_staphylococcus.csv      Species subset (20 strains)

------------------------------------------------------------------------------
02_mdr_analysis/
------------------------------------------------------------------------------
Script: 02_mdr_analysis.py
Input: amrfinder_results/*/*_amr.tsv
Description: Extracts drug classes per strain, defines MDR (>=3 classes), plots distribution

Files:
  - mdr_from_amrfinder_class.csv   Table with drug classes per strain + MDR status
  - mdr_amrfinder_plot.png         Stacked bar chart of MDR vs non-MDR per species

------------------------------------------------------------------------------
03_core_accessory_amr/
------------------------------------------------------------------------------
Script: 03_core_accessory_amr.py
Inputs: amrfinder_results/*/*_amr.tsv + prokka_out/*/*/*.gff + roary_out/*/gene_presence_absence.csv
Description: Maps AMR genes to Core/Accessory genome using coordinate matching

Files:
  - acinetobacter/amr_core_accessory.csv     Per-strain AMR mapping (Acinetobacter)
  - enterobacter/amr_core_accessory.csv      Per-strain AMR mapping (Enterobacter)
  - enterococcus/amr_core_accessory.csv      Per-strain AMR mapping (Enterococcus)
  - klebsiella/amr_core_accessory.csv        Per-strain AMR mapping (Klebsiella)
  - pseudomonas/amr_core_accessory.csv       Per-strain AMR mapping (Pseudomonas)
  - staphylococcus/amr_core_accessory.csv    Per-strain AMR mapping (Staphylococcus)
  - amr_core_accessory_all.csv               Combined table (all 120 strains)
  - final_summary_core_accessory.csv         Per-species summary counts
  - core_vs_accessory_amr_genes.png          Bar plot of core vs accessory AMR hits

------------------------------------------------------------------------------
04_final_summary/
------------------------------------------------------------------------------
Script: 04_final_summary.py
Inputs: 02_mdr_analysis/mdr_from_amrfinder_class.csv + 03_core_accessory_amr/final_summary_core_accessory.csv + roary_out/*/summary_statistics.txt + checkm2_results/*/quality_report.tsv
Description: Merges MDR, pangenome, CheckM2, and core/accessory data into integrated summary

Files:
  - final_species_summary.csv           Integrated species summary table
  - final_mdr_vs_accessory_clean.png    Scatter plot (accessory genome vs MDR%)
  - checkm2_summary.csv                 Per-strain CheckM2 quality metrics
  - checkm2_species_summary.csv         Per-species CheckM2 averages

------------------------------------------------------------------------------
05_heatmaps/
------------------------------------------------------------------------------
Script: 05_heatmaps.py
Input: 01_build_amr_matrix/amr_matrix.csv
Description: Generates species-wise clustered heatmaps of AMR gene presence/absence

Files:
  - Acinetobacter_heatmap.png
  - Enterobacter_heatmap.png
  - Enterococcus_heatmap.png
  - Klebsiella_heatmap.png
  - Pseudomonas_heatmap.png
  - Staphylococcus_heatmap.png

------------------------------------------------------------------------------
06_within_species_correlation/
------------------------------------------------------------------------------
Script: 06_within_species_correlation.py
Inputs: roary_out/*/gene_presence_absence.csv + 02_mdr_analysis/mdr_from_amrfinder_class.csv
Description: Tests per-strain correlation between accessory genome size and AMR class burden

Files:
  - within_species_correlation_summary.csv   Table of Pearson/Spearman per species
  - within_species_accessory_vs_amr.png      Faceted scatter plots

------------------------------------------------------------------------------
07_fisher_core_accessory/
------------------------------------------------------------------------------
Script: 07_fisher_core_accessory.py
Inputs: 03_core_accessory_amr/*/amr_core_accessory.csv + roary_out/*/gene_presence_absence.csv
Description: Fisher's exact tests for AMR enrichment in core vs accessory genome

Files:
  - fisher_core_accessory_enrichment.csv   Full statistical results (OR, CI, p-values)
  - fisher_core_accessory_oddsratios.png   Plot of ORs with 95% CIs (log scale)

------------------------------------------------------------------------------
08_virulence_analysis/
------------------------------------------------------------------------------
Script: 08_virulence_analysis.py
Inputs: prokka_out/*/*/*.faa + databases/vfdb.dmnd
Description: DIAMOND blastp screening against VFDB for virulence factor detection

Files:
  - *_vfdb.tsv (120 files)              Raw DIAMOND output per strain
  - virulence_matrix.csv                 Binary presence/absence matrix (120 x 2160 VFs)
  - virulence_strain_summary.csv         Per-strain VF counts
  - virulence_species_summary.csv        Per-species aggregated counts
  - virulence_count_boxplot.png          Box plot of VF counts per species
  - virulence_heatmap_top40.png          Heatmap of top 40 most prevalent VFs

------------------------------------------------------------------------------
09_phylogenetic_tree/
------------------------------------------------------------------------------
Script: 09_phylogenetic_tree.py
Input: roary_out/*/core_gene_alignment.aln
Description: Core genome neighbor-joining phylogenies per species

Files:
  - acinetobacter_core_tree.nwk     Newick format tree
  - enterobacter_core_tree.nwk      Newick format tree
  - enterococcus_core_tree.nwk      Newick format tree
  - klebsiella_core_tree.nwk        Newick format tree
  - pseudomonas_core_tree.nwk       Newick format tree
  - staphylococcus_core_tree.nwk    Newick format tree
  - tree_summary.csv                 Reconstruction summary (SNP counts, etc.)

------------------------------------------------------------------------------
BACKUP
------------------------------------------------------------------------------
analysis_results_backup/ contains the original flat directory structure before
reorganization (created prior to moving files).

------------------------------------------------------------------------------
10_fragmentation_accessory_correlation/
------------------------------------------------------------------------------
Script: 10_fragmentation_accessory_correlation.py
Inputs: roary_out/*/gene_presence_absence.csv + 04_final_summary/checkm2_summary.csv
Description: Correlates per-strain accessory genome size with assembly fragmentation
             metrics (total contigs, contig N50) to rule out assembly quality bias.

Files:
  - fragmentation_accessory_correlation.csv   # Pearson/Spearman summary per species
  - per_strain_accessory_fragmentation.csv    # Per-strain accessory counts and QC metrics
  - fragmentation_accessory_correlation.png   # Supplementary Figure S9

Manuscript reference:
  - Supplementary Figure S9
  - Supplementary Table S6

Key result: No consistent evidence that fragmentation inflates accessory gene counts;
            in five species, larger accessory genomes were associated with higher
            (more contiguous) contig N50 values.
