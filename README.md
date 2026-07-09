# Comparative Genomics of ESKAPE Pathogens

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21283840.svg)](https://doi.org/10.5281/zenodo.21283840)

This repository contains the analysis pipeline and summary data for:

> **Comparative Genomics of ESKAPE Pathogens: Integrating Pan-genome Architecture, Antimicrobial Resistance, and Virulence Factor Repertoires**
 
## Authors

**Pankaj Kumar**¹* · **Srishti Singh**¹*

¹ Division of Bioinformatics, Indian Agricultural Research Institute (IARI), New Delhi, India

*These authors contributed equally to this work.

| Author | ORCID | GitHub |
|--------|-------|--------|
| Pankaj Kumar | [0009-0006-3422-8881](https://orcid.org/0009-0006-3422-8881) | [@pankaj357](https://github.com/pankaj357) |
| Srishti Singh | [0009-0000-0949-4279](https://orcid.org/0009-0000-0949-4279) | [@srishtisingh03](https://github.com/srishtisingh03) |

---

## Overview

We analyzed 120 high-quality genomes (20 strains × 6 species) from six ESKAPE pathogen groups (*Acinetobacter*, *Enterobacter*, *Enterococcus*, *Klebsiella*, *Pseudomonas*, *Staphylococcus*). The pipeline integrates:

- **Genome quality control** — CheckM2
- **Genome annotation** — Prokka v1.15.6
- **Pan-genome analysis** — Roary v1.7.8
- **AMR profiling** — AMRFinderPlus v4.2.7
- **Virulence screening** — DIAMOND v2.2.1 against VFDB
- **Statistical analysis** — Python (pandas, SciPy, seaborn)

---

## Repository Structure

```
├── analysis_scripts/          # Python pipeline scripts (01–09)
├── analysis_results/          # Key result tables and figures
│   ├── 01_build_amr_matrix/
│   ├── 02_mdr_analysis/
│   ├── 04_final_summary/
│   ├── 06_within_species_correlation/
│   ├── 07_fisher_core_accessory/
│   └── 09_phylogenetic_tree/
├── metadata/                  # Strain metadata & NCBI accessions
├── AUTHORS.md                 # Authorship and contribution statement
├── README.md                  # This file
└── ZENODO_GUIDE.md            # Instructions for accessing full dataset
```

---

## Data Availability

All data supporting the findings of this study are available as follows:

| Data Type | Location |
|-----------|----------|
| **Analysis pipeline** | This repository (`analysis_scripts/`) |
| **Main result tables** | `analysis_results/` |
| **Supplementary tables** | Archived on Zenodo (see below) |
| **Supplementary figures** | Archived on Zenodo (see below) |
| **Phylogenetic trees** | `analysis_results/09_phylogenetic_tree/` |
| **Raw genome assemblies** | NCBI Genome database (`metadata/strain_accessions.csv`) |
| **Complete input-to-output dataset** | **Zenodo** https://doi.org/10.5281/zenodo.21283840 |

### Reference Databases
- **AMRFinderPlus reference database** — NCBI ([https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/))
- **VFDB (Virulence Factor Database, Set B)** — [http://www.mgc.ac.cn/VFs/](http://www.mgc.ac.cn/VFs/)

---

## Full Dataset (Zenodo)

The complete dataset—including raw genome assemblies, Prokka annotations, Roary pan-genome outputs, CheckM2 quality reports, AMRFinderPlus raw results, and per-strain DIAMOND outputs—is archived on **Zenodo**:

> **Zenodo Record (Version 2):** https://doi.org/10.5281/zenodo.21283840
> **Size:** ~8 GB  
> **License:** CC-BY 4.0

See [`ZENODO_GUIDE.md`](ZENODO_GUIDE.md) for detailed download and verification instructions.

---

## Reproducibility

### Conda Environments
Two environments were used:
- `eskaPE_env` — Prokka, AMRFinderPlus, DIAMOND, CheckM2
- `roary_env` — Roary (requires specific Perl dependencies)

### Running the Pipeline
Scripts are numbered in execution order:

```bash
python analysis_scripts/01_build_amr_matrix.py
python analysis_scripts/02_mdr_analysis.py
python analysis_scripts/03_core_accessory_amr.py
python analysis_scripts/04_final_summary.py
python analysis_scripts/05_heatmaps.py
python analysis_scripts/06_within_species_correlation.py
python analysis_scripts/07_fisher_core_accessory.py
python analysis_scripts/08_virulence_analysis.py
python analysis_scripts/09_phylogenetic_tree.py
```

---

## Citation

If you use this code or data, please cite the manuscript (in preparation) and the Zenodo dataset:

> Kumar, P., & Singh, S. (2026). *Comparative Genomics of ESKAPE Pathogens: Integrating Pan-genome Architecture, Antimicrobial Resistance, and Virulence Factor Repertoires* (Version 2) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.21283840

---

## License

- **Code:** MIT License
- **Data:** CC-BY 4.0

---

## Contact

For questions regarding the analysis pipeline or data, please open an issue on GitHub or contact:
- Pankaj Kumar: ft.pank@gmail.com
- Srishti Singh: srishtisingh5433@gmail.com
