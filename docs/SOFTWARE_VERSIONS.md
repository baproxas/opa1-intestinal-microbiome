# Software Versions

## Sequencing & Basecalling

- Oxford Nanopore Technologies MinION Mk1B
- Flow cell: FLO-MIN114 (R10.4.1)
- Library prep: Native Barcoding Kit 24 V14 (SQK-NBD114-24)
- MinKNOW: v25.03.9
- Basecalling: High-accuracy mode, 72 hours

## Taxonomic Classification

- Kraken2 (latest via conda metagenomics environment)
- Host database: GRCm39 (mouse) + GRCh38 (human)
- Classification database: NCBI core_nt (release: October 15, 2025; ~316 GB index)
- Confidence threshold: 0.01

## R Packages (key packages — full list in renv.lock)

- R: ≥ 4.3
- phyloseq
- microbiome
- vegan
- mixOmics (for sPLS-DA)
- fmsb (for radar charts)
- pheatmap (for CLR heatmaps)
- ggplot2, patchwork, ggrepel, ggpubr, scales
- dplyr, tidyr, tibble, readr, stringr
- rstatix (for ANOVA/Tukey tests)
- glmnet (for LASSO)
- ranger (for random forest)
- compositions (for CLR transform helpers)
- paletteer (for color palettes)

## HPC Environment

- University of Arizona HPC (Puma cluster)
- High-memory nodes: 3008 GB RAM (constraint=hi_mem)
- SLURM workload manager

## Random Seed

- `set.seed(2026)` used throughout all stochastic operations
