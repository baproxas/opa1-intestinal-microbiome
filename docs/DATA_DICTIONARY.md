# Data Dictionary

## config/sample_metadata/opa1_drp_meta_s2.csv

S2 subset metadata — 18 samples from 9 mice (OPA1 KO longitudinal: D0 / D8 / D21)

| Column | Type | Description |
|--------|------|-------------|
| id | string | Sample identifier (OPA01–OPA18); matches k2report filename prefix |
| fastq_filename | string | Original FASTQ filename (e.g., OPA01.fastq) |
| mouse_id | integer | Mouse identifier (4548–4560 for D0/D8; 3693–3793 for D21) |
| Day | integer | Experimental day (0, 8, or 21) relative to first tamoxifen injection |
| genotype | string | Mouse genotype (opa1floxed_villincre for all S2 samples) |
| trt | string | Treatment (none = no tamoxifen; tam = tamoxifen-treated) |
| phenotype | string | Phenotype classification (wt = wild-type OPA1+; opa1_ko = OPA1 knockout) |
| sex | string | Biological sex (m = male; f = female) |
| deletion | string | Target gene deletion (opa1 for all S2 samples) |
| wt_chg_pct | numeric | Body weight change from baseline (%, measured by scale) |
| fat_pct | numeric | Fat mass percentage (measured by EchoMRI) |
| lean_pct | numeric | Lean mass percentage (measured by EchoMRI) |
| group_id | string | Experimental group identifier combining genotype + treatment + day |
| output_label | string | Human-readable label for figure axis display |

## data/kraken_reports_core_nt/*.k2report

Kraken2 standard report format (one per sample, 18 files total). Tab-delimited with columns:

| Column | Description |
|--------|-------------|
| 1 | Percentage of fragments covered by the clade rooted at this taxon |
| 2 | Number of fragments covered by the clade rooted at this taxon |
| 3 | Number of fragments assigned directly to this taxon |
| 4 | Rank code (U, R, D, K, P, C, O, F, G, S) |
| 5 | NCBI taxonomy ID |
| 6 | Scientific name (indented by depth) |

## data/reference/guild_definitions.csv

| Column | Type | Description |
|--------|------|-------------|
| guild_name | string | Name of the functional guild (14 guilds from Table 1) |
| taxa_member | string | Species name belonging to this guild |
| references | string | Reference numbers from the manuscript |

## results/tables/*.csv

Stats output CSVs — one or more per figure; see [FIGURE_INDEX.md](FIGURE_INDEX.md) for the figure-to-table mapping. Column definitions below reflect `code/00_generate_all_figures.R` as currently written.

### Pairwise statistics files (Fig4A, Fig5A, Fig5C, Fig5D)

`Fig4A_alpha_diversity_pairwise_stats.csv`, `Fig5A_dysbiosis_score_pairwise_stats.csv`, `Fig5C_guild_radar_pairwise_stats.csv`, and `Fig5D_scfa_guild_pairwise_stats.csv` share the same statistical design (via the shared `run_kw_dunn_pairwise()` helper for 5C/5D):

- `omnibus_test` / `omnibus_p`: Kruskal-Wallis test across all three groups (D0, D8, D21).
- `test` / `p` / `p.adj`: Dunn post-hoc test with BH adjustment, run only when the omnibus p < 0.05; otherwise `test` records why it was skipped.
- `p.adj.signif` cutpoints: `< 0.001` = `***`, `< 0.01` = `**`, `< 0.05` = `*`, otherwise `ns`.
- `diagnostic`: reason a comparison was not tested (e.g., "Dunn not run: omnibus p >= 0.05", "Dunn not available", "Not tested: missing groups, too few values, or zero variance"); `NA` when the comparison was tested.
- All three pairwise group comparisons (D0-D8, D0-D21, D8-D21) are always reported as rows, including skipped ones.

| Column | Description |
|--------|-------------|
| metric / axis / feature | Which measured quantity this row's test applies to (name varies by file: `metric` in Fig4A, `axis` + `guild_column` in Fig5C, `feature` in Fig5D) |
| group1, group2 | The two groups being compared |
| n_group1, n_group2 | Sample counts in each group |
| omnibus_test, omnibus_p | See above |
| test, p, p.adj, p.adj.signif | See above |
| diagnostic | See above |

`Fig5A_dysbiosis_score_pairwise_stats.csv` has the same columns without a `metric`/`axis`/`feature` column (single outcome: `dysbiosis_score`).

### Fig4A_alpha_diversity_stats.csv

Per-sample alpha diversity metrics on the rarefied phyloseq object, joined to sample metadata.

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| observed | Observed richness |
| diversity_shannon | Shannon entropy |
| diversity_simpson | Simpson diversity |
| evenness_pielou | Pielou evenness |
| dominance_relative | Relative dominance |
| chao1, ace | Chao1 and ACE richness estimators (`vegan::estimateR`) |
| group_id, mouse_id, Day, ... | Joined sample metadata columns (see `opa1_drp_meta_s2.csv` above) |

### Fig4A_alpha_diversity_pairwise_stats.csv

See "Pairwise statistics files" above. `metric` is one of: observed, chao1, ace, diversity_shannon, evenness_pielou, dominance_relative.

### Fig4B_family_composition_plotted.csv

Per-sample, per-family relative abundances exactly as plotted in the Figure 4B stacked bar (top-20 families + "Other").

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| mouse_id | Mouse identifier |
| group_id | Experimental group |
| output_label | Figure axis label |
| family_label | Family name, or "Other" if outside the top 20 |
| relative_abundance | Fractional relative abundance (0-1) |
| percent | `relative_abundance * 100` |

### Fig4B_family_group_summary.csv

Per-family, per-group summary of relative abundance across samples.

| Column | Description |
|--------|-------------|
| family_label | Family name, or "Other" |
| group_id | Experimental group |
| n_samples | Number of samples in this group with this family |
| mean_relative_abundance, sd_relative_abundance, median_relative_abundance, min_relative_abundance, max_relative_abundance | Summary statistics of relative abundance within the group |
| pooled_total_abundance | Total abundance of this family pooled across all samples |
| top20_rank | Rank (1-20) among top families by pooled abundance; `NA` for "Other" |

### Fig4C_splsda_loadings.csv

sPLS-DA (`mixOmics::splsda`, keepX = 10, 2 components) taxon loadings.

| Column | Description |
|--------|-------------|
| taxon | Taxon identifier (best available rank label) |
| loading_comp1, loading_comp2 | Signed sPLS-DA loading on component 1 / 2 |
| abs_loading_comp1, abs_loading_comp2 | Absolute value of the loadings |
| selected_comp1, selected_comp2 | Logical; `TRUE` if the loading is non-zero (selected by keepX sparsity) on that component |
| in_heatmap | Logical; `TRUE` if this taxon is among the top taxa shown in the Figure 4C heatmap |

### Fig4C_splsda_heatmap_clr_matrix.csv

CLR-transformed abundance matrix underlying the Figure 4C heatmap.

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| group_id | Experimental group |
| (one column per top taxon) | CLR-transformed relative abundance for each taxon shown in the heatmap |

### Fig4C_splsda_sample_scores.csv

Per-sample sPLS-DA component scores.

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| mouse_id | Mouse identifier |
| group_id | Experimental group |
| comp1, comp2 | sPLS-DA component 1 / 2 score |

### Fig5A_dysbiosis_score.csv

Per-sample Bray-Curtis dysbiosis score (mean distance to the D0/reference-group centroid), joined to sample metadata.

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| dysbiosis_score | Mean Bray-Curtis distance to the reference (D0) group |
| group_id, mouse_id, Day, ... | Joined sample metadata columns |

### Fig5A_dysbiosis_score_stats.csv

Raw `rstatix::dunn_test()` output (only written if the omnibus Kruskal-Wallis p < 0.05), with standard rstatix columns (`.y.`, `group1`, `group2`, `n1`, `n2`, `statistic`, `p`, `p.adj`, `p.adj.signif`) for the `dysbiosis_score ~ group_id` comparison.

### Fig5A_dysbiosis_score_pairwise_stats.csv

See "Pairwise statistics files" above (single outcome: `dysbiosis_score`).

### Fig5B_dysbiosis_bodycomp_stats.csv

Spearman correlations between dysbiosis score and body-composition outcomes.

| Column | Description |
|--------|-------------|
| outcome | `wt_chg_pct` (weight change %) or `fat_pct` (fat mass %) |
| rho | Spearman correlation coefficient |
| pval | Nominal p-value |
| padj_BH | BH-adjusted p-value across the two outcomes |

### Fig5C_guild_scores_per_sample.csv

Per-sample functional guild relative-abundance scores (species-level guild membership sums), joined to sample metadata. **Renamed from `S2_guild_scores.csv`.**

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| (one column per guild in `guild_definitions.csv`, e.g. "Butyrate Producers", "Acetate Producers", "Depleted in IBD", ...) | Summed relative abundance of species belonging to that guild |
| ROS_sensitive | Relative abundance sum of the ROS-sensitive taxon set (multi-rank; see Figure 5C radar) |
| group_id, mouse_id, Day, ... | Joined sample metadata columns |

### Fig5C_guild_radar_scores.csv

D0-referenced validation table for the 9 Figure 5C radar axes.

- Normalization is D0-referenced: `normalized_value = 100 * group_mean / D0_group_mean`.
- `D0_is_raw_maximum` reports whether D0 has the highest raw mean of the three groups on that axis; it is **reported only, not enforced** (D8/D21 raw means may legitimately exceed D0).
- Values above 100% are preserved as-is, not clipped.

| Column | Description |
|--------|-------------|
| axis | Radar axis display label (e.g. "Butyrate", "IBD-depleted") |
| raw_D0_mean, raw_D8_mean, raw_D21_mean | Raw (un-normalized) group mean for that axis |
| normalized_D0, normalized_D8, normalized_D21 | D0-referenced percentage (D0 is always 100 when valid) |
| D0_is_raw_maximum | Logical/`NA`; see above |
| denominator_status | "ok", "D0 mean missing", or "D0 mean is zero" — whether the D0 denominator was usable |

### Fig5C_guild_radar_pairwise_stats.csv

See "Pairwise statistics files" above. Tests are run on the **raw** guild score per axis (not the normalized percentage used for plotting). `axis` is the radar display label; `guild_column` is the underlying `guild_df` column name.

### Fig5D_scfa_guild_plotted.csv

Per-sample SCFA-producer guild relative abundances exactly as plotted in the Figure 5D stacked bar.

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| mouse_id | Mouse identifier |
| group_id | Experimental group |
| guild | One of "Butyrate Producers", "Propionate Producers", "Acetate Producers" |
| relative_abundance | Relative abundance of that guild in that sample |
| scfa_total | Sum of all three SCFA guild relative abundances for that sample |

### Fig5D_scfa_guild_pairwise_stats.csv

See "Pairwise statistics files" above. `feature` is one of the three SCFA guilds or `scfa_total`.

### Fig5E_partial_spearman_dysbiosis.csv

Partial Spearman correlation (CLR-taxon residuals vs. dysbiosis-score residuals, both regressed on `group_id`) underlying the Figure 5E volcano plot.

- y-axis plots `-log10(BH-adj. p)`, not the nominal p-value.
- Visual labeling threshold is `padj < 0.05`.
- The `sig` column in this CSV uses `padj < 0.2` (consistent with the manuscript's "50/140 genera" reporting threshold) — a separate, less strict flag than the plot's visual/labeling threshold.
- Guild annotation (`func_guild`) uses a multi-rank lookup (Species/Genus/Family/Order/Phylum) reproduced from Module FH4 of `09_opa1_drp_core_nt_functional_analysis.R`.
- ML top taxa are sourced from the in-memory SECTION 4 sPLS-DA selection (`top_taxa_ff`), with a fallback to reading `Fig4C_splsda_loadings.csv` (`in_heatmap == TRUE`) if unavailable.
- Guild members keep their guild assignment; "ML top taxa" is assigned only to taxa with no guild match.
- Labels are placed inside the volcano's left/right arms using two-arm `ggrepel` with `xlim` constraints (not directly reflected as a column, but explains why some non-significant-looking points near the plot centre are unlabeled).

| Column | Type | Description |
|--------|------|-------------|
| taxon | character | Taxon identifier (best available rank label) |
| rho | numeric | Partial Spearman correlation coefficient |
| pval | numeric | Nominal p-value |
| padj | numeric | BH-adjusted p-value across all taxa |
| sig | logical | `TRUE` if `padj < 0.2` (see above) |
| dir | character | "Dysbiosis-promoting" (`rho > 0`) or "Dysbiosis-protecting" (`rho <= 0`) |
| func_guild | character | Functional guild name, "ML top taxa", or `NA` |
| in_ml_top_taxa | logical | `TRUE` if taxon is in the sPLS-DA top-20 |
| colour_cat | character | Final plotting category: guild name, "ML top taxa", "Promoting (other)", or "Protecting (other)" |
| label_set | character | "guild", "ML top taxa", "top5", or `NA` — which label rule (if any) selected this taxon for on-plot text |

### Statistical caveat (Figures 4A, 5A, 5C, 5D, 5E)

D0 and D8 samples are paired within the same mouse (longitudinal, pre/post tamoxifen); D21 samples come from a separate cohort. All Kruskal-Wallis/Dunn and partial-Spearman comparisons in this pipeline treat the three groups (D0, D8, D21) as independent samples for rank-based testing — they do not model the within-mouse pairing of D0/D8.
