# Zenodo Dataset Guide

This document explains how to access, verify, and cite the complete ESKAPE comparative genomics dataset archived on Zenodo.

---

## Dataset Record

| Field | Value |
|-------|-------|
| **Title** | Comparative Genomics of ESKAPE Pathogens: Integrating Pan-genome Architecture, Antimicrobial Resistance, and Virulence Factor Repertoires — Full Dataset |
| **Authors** | Pankaj Kumar, Srishti Singh |
| **Affiliation** | Indian Agricultural Research Institute (IARI), New Delhi, India |
| **DOI** | [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) |
| **URL** | [https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) |
| **License** | CC-BY 4.0 |
| **Upload type** | Dataset |
| **Size** | ~8 GB |

---

## What Is Included in the Zenodo Archive

The Zenodo deposit contains the **complete input-to-output dataset** (everything except manuscript files):

```
eskape_full_dataset/
├── analysis_scripts/              # Python pipeline scripts (01–09)
├── analysis_results/              # All result tables and figures
│   ├── 01_build_amr_matrix/
│   ├── 02_mdr_analysis/
│   ├── 03_core_accessory_amr/
│   ├── 04_final_summary/
│   ├── 05_heatmaps/
│   ├── 06_within_species_correlation/
│   ├── 07_fisher_core_accessory/
│   ├── 08_virulence_analysis/     # Includes 120 per-strain VFDB TSVs
│   └── 09_phylogenetic_tree/
├── genomes/                       # 120 genome assemblies (.fna)
├── prokka_out/                    # Prokka annotations for all 120 strains
├── roary_out/                     # Roary pan-genome outputs for all 6 species
├── checkm2_results/               # CheckM2 quality reports
├── amrfinder_results/             # Raw AMRFinderPlus outputs
├── databases/                     # VFDB reference files used for screening
├── metadata/
│   └── strain_accessions.csv      # NCBI accession numbers
├── AUTHORS.md
├── README.md
└── ZENODO_GUIDE.md
```

---

## How to Download

### Option A: Direct Download (Browser)
1. Visit the Zenodo record URL.
2. Click the **Download** button next to `eskape_full_dataset.zip`.

### Option B: Command Line (wget/curl)
```bash
# Replace XXXXXXX with the actual Zenodo record ID
wget https://zenodo.org/record/XXXXXXX/files/eskape_full_dataset.zip

# Verify checksum (provided on Zenodo page)
md5sum eskape_full_dataset.zip

# Extract
unzip eskape_full_dataset.zip
```

### Option C: Python
```python
import requests

url = "https://zenodo.org/record/XXXXXXX/files/eskape_full_dataset.zip"
r = requests.get(url)
with open("eskape_full_dataset.zip", "wb") as f:
    f.write(r.content)
```

---

## Verification

After downloading and extracting, verify the dataset integrity:

```bash
# Count genome files (should be 120)
find eskape_full_dataset/genomes/ -name "*.fna" | wc -l

# Count Prokka annotation folders (should be 120)
ls -1 eskape_full_dataset/prokka_out/*/* | wc -l

# Count Roary species folders (should be 6)
ls -1 eskape_full_dataset/roary_out/ | wc -l

# Verify CheckM2 reports exist
ls eskape_full_dataset/checkm2_results/*/quality_report.tsv

# Check total size
du -sh eskape_full_dataset/
```

Expected output:
- 120 genome files
- 120 Prokka annotation directories
- 6 Roary species directories
- Total size: ~7.9 GB

---

## How to Cite the Dataset

### APA
Kumar, P., & Singh, S. (2024). *Comparative Genomics of ESKAPE Pathogens: Integrating Pan-genome Architecture, Antimicrobial Resistance, and Virulence Factor Repertoires — Full Dataset* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

### BibTeX
```bibtex
@dataset{kumar_singh_2024_eskape,
  author       = {Kumar, Pankaj and Singh, Srishti},
  title        = {Comparative Genomics of ESKAPE Pathogens: 
                  Integrating Pan-genome Architecture, Antimicrobial 
                  Resistance, and Virulence Factor Repertoires — 
                  Full Dataset},
  month        = may,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

---

## Relationship to GitHub Repository

| Platform | Content | Size |
|----------|---------|------|
| **GitHub** | Code + summary tables + metadata | ~30 MB |
| **Zenodo** | Full input-to-output dataset (genomes, annotations, intermediates, results) | ~8 GB |

The GitHub repository ([https://github.com/pankaj357/eskape-comparative-genomics](https://github.com/pankaj357/eskape-comparative-genomics)) contains the analysis scripts and key summary tables. The Zenodo archive contains the complete dataset required to fully reproduce every figure and table from scratch.

---

## Questions?

For issues with the dataset or reproduction, please:
1. Open an issue on GitHub: [https://github.com/pankaj357/eskape-comparative-genomics/issues](https://github.com/pankaj357/eskape-comparative-genomics/issues)
2. Or contact the authors directly (see `AUTHORS.md`).
