# Methods Mapping

This document maps sentences from the Methods section of Roxas et al. (2026) to specific sections of `code/00_generate_all_figures.R`.

| Manuscript Statement | Script Section | Notes |
|---------------------|----------------|-------|
| "Taxa were prevalence-filtered prior to downstream analysis by retaining taxa present at >0.01% relative abundance in ≥10% of samples" | SECTION 0: `prev_filter()` function | Parameters in config.yml: min_prev=0.10, min_abund=0.0001 |
| "CLR transformation (pseudocount = 0.5) was applied prior to regression and feature selection" | SECTION 0: `clr_transform()` function | Pseudocount in config.yml: clr_pseudocount=0.5 |
| "Alpha diversity was calculated from rarefied counts and included observed richness, Chao1, ACE, Shannon entropy, Pielou evenness, and relative dominance" | SECTION 2: Fig 4A | Uses `rarefy_even_depth()` + `microbiome::alpha()` + `vegan::estimateR()` |
| "A dysbiosis score was calculated for each sample as the mean Bray-Curtis distance to all day 0 reference samples" | SECTION 5: Fig 5A | `dysbiosis_score()` function with ref = D0 group |
| "Dysbiosis scores were positively associated with metabolic disease severity" (ρ = −0.61, padj = 0.020; ρ = −0.47, padj = 0.070) | SECTION 6: Fig 5B | Spearman with BH correction across 3 outcomes |
| "sPLS-DA was performed using CLR-transformed taxonomic features, keepX=10, LOMO-CV" | SECTION 4: Fig 4C | `mixOmics::splsda()` |
| "Guild scores were computed as the sum of relative abundances of all detected guild member taxa per sample" | SECTION 7/8: Fig 5C/5D | Guild definitions from `guild_definitions.csv` |
| "Nine functional axes displayed in a radar plot, each index normalized to a 0–100% range" | SECTION 7: Fig 5C | `fmsb::radarchart()` |
| "Partial Spearman correlation analysis was performed... residualized on group identity" | SECTION 9: Fig 5E | Residual-on-residual approach, BH + BY + maxT correction |
| "Of 140 genera tested, 50 passed BH adjusted significance threshold (padj < 0.2)" | SECTION 9: Fig 5E | 22 promoting + 28 protecting |
| "set.seed(2026)" | SECTION 0 | Consistent across all stochastic operations |
