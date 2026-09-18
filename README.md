# Epithelial OPA1 Loss Drives IBD-Concordant Gut Microbiome Remodeling

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Docs-CC%20BY%204.0-lightgrey.svg)](LICENSE-CC-BY-4.0)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.22835154)

## Summary

OPA1-dependent mitochondrial fusion is indispensable for intestinal epithelial integrity. Using a tamoxifen-inducible, villin-creERT2-driven intestine-specific *opa1* conditional knockout mouse model, we show that epithelial OPA1 loss drives progressive microbiome dysbiosis, depletion of butyrate-producing obligate anaerobes, and enrichment of aerotolerant taxa — an IBD-concordant community shift. Dysbiosis severity correlated with host metabolic decline (weight loss and reduced adiposity). OPA1-deficient mice succumbed to DSS or *Citrobacter rodentium* challenge, positioning epithelial mitochondrial dynamics as a critical upstream determinant of gut microbial community structure and disease susceptibility.

## Citation

> Roxas JL, Roxas BAP, Rutins I, Rubinstein S, Kedia S, Holyoak-Aguirre AA, Lucarevskiy L, Cocchi K, Lindsey J, Anwar F, Sullivan A, Ghosh B, Obergh V, Scranton CE, Cooper KK, Wilson J, Vedantam G, and Viswanathan VK. **Intestine-specific OPA1 Loss Drives Progressive Intestinal Epithelial Damage, Metabolic Decline and IBD-Concordant Gut Microbiome Remodeling and Susceptibility to Colitogenic Insults.** *Journal TBD*, 2026.

See [CITATION.cff](CITATION.cff) for machine-readable citation. Zenodo DOI: `10.5281/zenodo.22835154` .

## Data Availability

| Resource | Location | Status |
|----------|----------|--------|
| Raw FASTQ (all 35 samples) | NCBI BioProject [PRJNA1523974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1523974) | Will be public upon publication |
| Kraken2 reports (18 S2 samples) | This repo: `data/kraken_reports_core_nt/` | Available now |
| Analysis code | This repo: `code/00_generate_all_figures.R` | Available now |
| Sample metadata | This repo: `config/sample_metadata/opa1_drp_meta_s2.csv` | Available now |
| Guild definitions (Table 1) | This repo: `data/reference/guild_definitions.csv` | Available now |

## Quick Start

Reproduce all 8 manuscript figure panels (Figures 4A-4C, 5A-5E) from the processed Kraken2 reports included in this repository:

```bash
git clone https://github.com/baproxas/opa1-intestinal-microbiome.git
cd opa1-intestinal-microbiome

