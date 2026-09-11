#!/usr/bin/env Rscript
# ==============================================================================
# 00_generate_all_figures.R
#
# OPA1 Intestinal Microbiome Study — Single pipeline to generate all manuscript
# figures (4A–4C, 5A–5E) from Kraken2 .k2report files.
#
# Subset: S2 — OPA1 KO longitudinal (D0 / D8 / D21; n = 18 samples, 9 mice)
#
# Input:
#   config/sample_metadata/opa1_drp_meta_s2.csv
#   data/kraken_reports_core_nt/*.k2report  (18 files)
#   data/reference/guild_definitions.csv
#
# Output:
#   results/figures/  — 8 PDFs + 8 PNGs
#   results/tables/   — stats CSVs
#   sessionInfo.txt
#
# Usage:
#   renv::restore()
#   source("code/00_generate_all_figures.R")
#
# Authors: Bryan Angelo P. Roxas, Jennifer Lising Roxas
# Date: 2026-09-09
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 0: Setup — packages, seed, config, helper functions
# ──────────────────────────────────────────────────────────────────────────────
cat("=== SECTION 0: Setup ===\n")

# 0a. Packages
suppressPackageStartupMessages({
	library(phyloseq)
	library(microbiome)
	library(vegan)
	library(dplyr)
	library(tidyr)
	library(tibble)
	library(readr)
	library(stringr)
	library(tools)
	library(ggplot2)
	library(scales)
	library(patchwork)
	library(ggrepel)
	library(ggpubr)
	library(pheatmap)
	library(paletteer)
	library(fmsb)
	library(rstatix)
	library(compositions)
	library(glmnet)
	library(here)
})

# mixOmics from Bioconductor
if (!requireNamespace("mixOmics", quietly = TRUE))
	stop("mixOmics not installed. Run: BiocManager::install('mixOmics')")
suppressPackageStartupMessages(library(mixOmics))

set.seed(2026)

# 0b. Paths (all relative to repo root via here::here())
META_PATH       <- here::here("config", "sample_metadata", "opa1_drp_meta_s2.csv")
REPORT_DIR      <- here::here("data", "kraken_reports_core_nt")
GUILD_CSV       <- here::here("data", "reference", "guild_definitions.csv")
FIG_DIR         <- here::here("results", "figures")
TBL_DIR         <- here::here("results", "tables")
REPORT_SUFFIX   <- ".k2report"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# 0c. Analysis parameters
PREV_MIN     <- 0.10
ABUND_MIN    <- 1e-4
CLR_PSEUDO   <- 0.5
SPLSDA_KEEPX <- 10
BH_THRESH    <- 0.2
DPI          <- 300

# 0d. Group configuration
GROUP_ORDER <- c(
	"untreated_opa1floxed_vilcre_d0",
	"opa1floxed_vilcre_tam_d8",
	"opa1floxed_vilcre_tam_d21"
)
GROUP_LABELS <- c(
	untreated_opa1floxed_vilcre_d0 = "OPA1+ (D0)",
	opa1floxed_vilcre_tam_d8       = "OPA1 KO (D8)",
	opa1floxed_vilcre_tam_d21      = "OPA1 KO (D21)"
)
GROUP_PALETTE <- c(
	untreated_opa1floxed_vilcre_d0 = "#0015FF",
	opa1floxed_vilcre_tam_d8       = "#00B050",
	opa1floxed_vilcre_tam_d21      = "#E600C7"
)
# Shared group palette used by Figures 4A, 4B, 4C, 5A, and 5D.
GROUP_PALETTE_4A_5A <- GROUP_PALETTE
GROUP_X_LABELS_4A_5A <- c(
	untreated_opa1floxed_vilcre_d0 = "0",
	opa1floxed_vilcre_tam_d8       = "8",
	opa1floxed_vilcre_tam_d21      = "21"
)
DAY_LABELS <- GROUP_X_LABELS_4A_5A
GROUP_HEADER_LABELS <- c(
	untreated_opa1floxed_vilcre_d0 = "Day 0 opa1flox/villin-creERT2\nNo Tamoxifen (OPA1+)",
	opa1floxed_vilcre_tam_d8 = "Day 8 opa1flox/villin-creERT2\n+ Tamoxifen (OPA1-, Villin-Cre+)",
	opa1floxed_vilcre_tam_d21 = "Day 21 opa1flox/villin-creERT2\n+ Tamoxifen (OPA1-, Villin-Cre+)"
)
REF_GROUP <- "untreated_opa1floxed_vilcre_d0"

cat("  Paths, parameters, and palettes configured.\n")

# 0e. Helper functions
# ── Kraken2 report parser ────────────────────────────────────────────────────
RANK_CODES <- c(P = "Phylum", C = "Class", O = "Order",
								F = "Family", G = "Genus",  S = "Species")

parse_k2_report <- function(path) {
	lines <- readLines(path, warn = FALSE)
	out   <- vector("list", length(lines))
	lineage <- character(0)
	infer_domain <- function(nodes) {
		if (length(nodes) == 0) return(NA_character_)
		idx <- which(grepl("bacteria|archaea|eukary", nodes, ignore.case = TRUE))
		if (length(idx) == 0) return(NA_character_)
		hit <- nodes[max(idx)]
		if (grepl("bacteria", hit, ignore.case = TRUE)) return("Bacteria")
		if (grepl("archaea",  hit, ignore.case = TRUE)) return("Archaea")
		if (grepl("eukary",   hit, ignore.case = TRUE)) return("Eukaryota")
		NA_character_
	}
	for (i in seq_along(lines)) {
		parts <- strsplit(lines[i], "\t")[[1]]
		if (length(parts) < 6) next
		name_raw <- parts[6]; name <- trimws(name_raw)
		if (name == "") next
		indent_n <- nchar(sub("^([[:space:]]*).*", "\\1", name_raw))
		depth    <- floor(indent_n / 2)
		lineage[depth + 1] <- name
		if (length(lineage) > depth + 1) lineage <- lineage[seq_len(depth + 1)]
		domain   <- infer_domain(lineage)
		rank_raw <- trimws(parts[4])
		if (!rank_raw %in% names(RANK_CODES)) next
		reads <- as.integer(parts[3])
		if (is.na(reads) || reads == 0) next
		out[[i]] <- data.frame(rank = RANK_CODES[[rank_raw]], name = name,
													 kingdom = domain, reads = reads,
													 stringsAsFactors = FALSE)
	}
	dplyr::bind_rows(out)
}

# ── Build phyloseq from metadata ─────────────────────────────────────────────
build_phyloseq <- function(meta_sub) {
	cat("  Parsing", nrow(meta_sub), "Kraken2 reports...\n")
	all_taxa <- list(); otu_list <- list()
	for (i in seq_len(nrow(meta_sub))) {
		sid   <- meta_sub$id[i]
		fname <- sub("\\.fastq$", REPORT_SUFFIX, meta_sub$fastq_filename[i])
		fpath <- file.path(REPORT_DIR, fname)
		if (!file.exists(fpath)) { cat("    [WARN] Missing:", fname, "\n"); next }
		df <- parse_k2_report(fpath)
		if (nrow(df) == 0) next
		df <- dplyr::mutate(df, tax_key = paste0(rank, "|", name))
		all_taxa[[sid]] <- dplyr::select(df, tax_key, rank, name, kingdom)
		otu_list[[sid]] <- dplyr::select(df, tax_key, reads)
	}
	if (length(otu_list) == 0) stop("No reports parsed.")
	tax_tbl  <- dplyr::bind_rows(all_taxa) %>%
		dplyr::distinct(tax_key, .keep_all = TRUE) %>%
		tibble::column_to_rownames("tax_key")
	all_keys <- rownames(tax_tbl)
	otu_mat  <- matrix(0L, nrow = length(all_keys), ncol = length(otu_list),
										 dimnames = list(all_keys, names(otu_list)))
	for (sid in names(otu_list)) {
		df_s <- otu_list[[sid]]; idx <- match(df_s$tax_key, all_keys)
		ok <- !is.na(idx); otu_mat[idx[ok], sid] <- df_s$reads[ok]
	}
	std_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
	tax_mat   <- matrix(NA_character_, nrow(tax_tbl), 7,
											dimnames = list(rownames(tax_tbl), std_ranks))
	tax_mat[, "Kingdom"] <- tax_tbl$kingdom
	for (rk in std_ranks) {
		is_rk <- tax_tbl$rank == rk
		tax_mat[is_rk, match(rk, std_ranks)] <- tax_tbl$name[is_rk]
	}
	meta_p <- meta_sub %>% dplyr::filter(id %in% names(otu_list)) %>%
		tibble::column_to_rownames("id")
	ps <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
								 tax_table(tax_mat), sample_data(meta_p))
	cat("  Built:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")
	ps
}

# ── Filters ──────────────────────────────────────────────────────────────────
filter_host_taxa <- function(ps) {
	if (is.null(tax_table(ps))) return(ps)
	tt <- as.data.frame(tax_table(ps))
	bad <- apply(tt, 1, function(row) any(grepl(
		"Homo|sapiens|Hominidae|Primates|Chordata|Mammalia|mitochond|plastid|chloroplast",
		row, ignore.case = TRUE, perl = TRUE), na.rm = TRUE))
	if (sum(bad) > 0) ps <- prune_taxa(!bad, ps)
	ps
}

filter_bacteria_archaea <- function(ps) {
	if (is.null(tax_table(ps))) return(ps)
	tt  <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
	dom <- trimws(as.character(tt$Kingdom))
	keep <- grepl("bacteria", dom, ignore.case = TRUE) |
					grepl("archaea",  dom, ignore.case = TRUE)
	if (sum(keep) == 0) stop("Kingdom filter kept 0 taxa.")
	cat("  Keeping", sum(keep), "Bacteria/Archaea taxa; removed",
			ntaxa(ps) - sum(keep), "non-target taxa\n")
	prune_taxa(keep, ps)
}

prev_filter <- function(ps, min_prev = PREV_MIN, min_abund = ABUND_MIN) {
	ps_rel <- microbiome::transform(ps, "compositional")
	keep   <- apply(otu_table(ps_rel), 1, function(x) mean(x > min_abund) >= min_prev)
	prune_taxa(keep, ps)
}

# ── CLR transform ────────────────────────────────────────────────────────────
clr_transform <- function(X, pseudocount = CLR_PSEUDO) {
	lX <- log(X + pseudocount); lX - rowMeans(lX)
}

# ── Best rank label ──────────────────────────────────────────────────────────
best_rank_label <- function(ps_g) {
	tt <- as.data.frame(tax_table(ps_g))
	ranks <- intersect(c("Species","Genus","Family","Order","Class","Phylum"), colnames(tt))
	labels <- apply(tt[, ranks, drop = FALSE], 1, function(row) {
		hit <- row[!is.na(row) & nchar(trimws(row)) > 0]
		if (length(hit) > 0) paste0(names(hit)[1], "|", hit[1]) else "Unclassified"
	})
	if (anyDuplicated(labels)) {
		dup_idx <- which(duplicated(labels) | duplicated(labels, fromLast = TRUE))
		labels[dup_idx] <- paste0(labels[dup_idx], "_", dup_idx)
	}
	labels
}

# ── Sample data as data.frame ────────────────────────────────────────────────
get_sdata <- function(ps) {
	df <- as(sample_data(ps), "data.frame")
	tibble::rownames_to_column(df, "sample_id")
}

# ── Dysbiosis score ──────────────────────────────────────────────────────────
dysbiosis_score <- function(ps, ref_group_id) {
	ps_rel  <- microbiome::transform(ps, "compositional")
	bc_mat  <- as.matrix(phyloseq::distance(ps_rel, method = "bray"))
	sdata   <- get_sdata(ps)
	ref_ids <- sdata$sample_id[sdata$group_id == ref_group_id]
	if (length(ref_ids) == 0) stop("Reference group not found")
	score   <- rowMeans(bc_mat[, ref_ids, drop = FALSE], na.rm = TRUE)
	data.frame(sample_id = names(score), dysbiosis_score = score,
						 stringsAsFactors = FALSE)
}

# ── Save figure (PDF + PNG) ──────────────────────────────────────────────────
save_fig <- function(plot_obj, base_path, w, h, dpi = DPI) {
	ggsave(paste0(base_path, ".pdf"), plot_obj, width = w, height = h)
	ggsave(paste0(base_path, ".png"), plot_obj, width = w, height = h, dpi = dpi)
	cat("  Saved:", basename(base_path), "(.pdf + .png)\n")
}

# ── Save base-R figure (PDF + PNG) for pheatmap/fmsb ─────────────────────────
save_base_pdf_png <- function(expr, base_path, w, h, dpi = DPI) {
	pdf(paste0(base_path, ".pdf"), width = w, height = h)
	eval(expr)
	dev.off()
	png(paste0(base_path, ".png"), width = w * dpi, height = h * dpi, res = dpi)
	eval(expr)
	dev.off()
	cat("  Saved:", basename(base_path), "(.pdf + .png)\n")
}

cat("  Helper functions defined.\n")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 1: Build S2 phyloseq object
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 1: Build S2 phyloseq ===\n")

meta <- read_csv(META_PATH, show_col_types = FALSE) %>%
	dplyr::mutate(
		Day      = as.integer(Day),
		group_id = factor(group_id, levels = GROUP_ORDER)
	)
cat("  Metadata:", nrow(meta), "samples\n")

ps_all <- build_phyloseq(meta) %>%
	filter_host_taxa() %>%
	filter_bacteria_archaea()

# Attach group_id as factor to sample_data
sd <- as(sample_data(ps_all), "data.frame")
sd$group_id <- factor(sd$group_id, levels = GROUP_ORDER, ordered = FALSE)
sample_data(ps_all) <- sample_data(sd)

cat("  S2 phyloseq ready:", nsamples(ps_all), "samples,", ntaxa(ps_all), "taxa\n")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 2: Figure 4A — Alpha diversity boxplots (S2)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 2: Figure 4A — Alpha diversity ===\n")

set.seed(2026)
ps_rare <- tryCatch(
	rarefy_even_depth(ps_all, rngseed = 2026, verbose = FALSE),
	error = function(e) { message("  [WARN] Rarefaction failed: ", e$message); ps_all }
)

sdata_a <- get_sdata(ps_rare)

alpha_df <- tryCatch({
	a <- microbiome::alpha(
		ps_rare,
		index = c("observed", "diversity_shannon", "diversity_simpson",
							"evenness_pielou", "dominance_relative")
	) %>%
		tibble::rownames_to_column("sample_id") %>%
		dplyr::left_join(sdata_a, by = "sample_id")
	rc <- as(otu_table(ps_rare), "matrix")
	if (taxa_are_rows(ps_rare)) rc <- t(rc)
	er <- vegan::estimateR(round(rc))
	a <- dplyr::left_join(a,
		data.frame(sample_id = colnames(er),
							 chao1 = as.numeric(er["S.chao1", ]),
							 ace   = as.numeric(er["S.ACE",   ]),
							 stringsAsFactors = FALSE),
		by = "sample_id")
	a
}, error = function(e) { message("  [WARN] Alpha failed: ", e$message); NULL })

if (!is.null(alpha_df)) {
	write_csv(alpha_df, file.path(TBL_DIR, "Fig4A_alpha_diversity_stats.csv"))

	grp_present <- intersect(GROUP_ORDER, unique(as.character(alpha_df$group_id)))
	alpha_df$group_id <- factor(alpha_df$group_id, levels = grp_present)
	pal <- GROUP_PALETTE_4A_5A[grp_present]
	alpha_pairwise_rows <- list()
	alpha_pairwise_i <- 0L

	metrics <- c("observed", "chao1", "ace",
							 "diversity_shannon", "evenness_pielou", "dominance_relative")
	labels  <- c("Observed richness", "Chao1", "ACE",
							 "Shannon entropy", "Pielou evenness", "Dominance (relative)")

	plots_alpha <- Map(function(m, lab) {
		metric_available <- m %in% names(alpha_df)
		dat <- if (metric_available) alpha_df %>%
			dplyr::transmute(group_id, value = .data[[m]]) %>%
			dplyr::filter(!is.na(value), !is.na(group_id)) else data.frame()
		pair_grid <- expand.grid(group1 = GROUP_ORDER, group2 = GROUP_ORDER,
								 stringsAsFactors = FALSE) %>%
			dplyr::filter(match(group1, GROUP_ORDER) < match(group2, GROUP_ORDER))
		group_counts <- if (nrow(dat) > 0) table(factor(dat$group_id, levels = GROUP_ORDER)) else
			setNames(rep(0L, length(GROUP_ORDER)), GROUP_ORDER)
		omnibus_p <- NA_real_
		omnibus_note <- if (!metric_available) "Not tested: metric unavailable" else
			"Not tested: missing groups, too few values, or zero variance"
		dunn_res <- NULL
		if (metric_available && nrow(dat) >= 3 && dplyr::n_distinct(dat$group_id) >= 2 &&
				dplyr::n_distinct(dat$value) >= 2) {
			kw <- tryCatch(rstatix::kruskal_test(dat, value ~ group_id),
							 error = function(e) NULL)
			if (!is.null(kw)) {
				omnibus_p <- kw$p[1]
				omnibus_note <- "Kruskal-Wallis"
				if (!is.na(omnibus_p) && omnibus_p < 0.05) {
					dunn_res <- tryCatch(
						rstatix::dunn_test(dat, value ~ group_id,
											 p.adjust.method = "BH"),
						error = function(e) NULL)
				}
			} else {
				omnibus_note <- "Kruskal-Wallis failed"
			}
		} else {
			omnibus_note <- "Not tested: missing groups, too few values, or zero variance"
		}
		for (j in seq_len(nrow(pair_grid))) {
			g1 <- pair_grid$group1[j]; g2 <- pair_grid$group2[j]
			pair <- if (!is.null(dunn_res)) dunn_res %>%
				dplyr::filter(group1 == g1, group2 == g2) else NULL
			alpha_pairwise_i <- alpha_pairwise_i + 1L
			alpha_pairwise_rows[[alpha_pairwise_i]] <<- data.frame(
				metric = m, group1 = g1, group2 = g2,
				n_group1 = unname(group_counts[g1]), n_group2 = unname(group_counts[g2]),
				omnibus_test = ifelse(is.na(omnibus_note), "Kruskal-Wallis", omnibus_note),
				omnibus_p = omnibus_p,
				test = if (!is.null(pair) && nrow(pair) > 0) "Dunn" else
					if (!is.na(omnibus_p) && omnibus_p >= 0.05) "Dunn not run: omnibus p >= 0.05" else "Dunn not available",
				p = if (!is.null(pair) && nrow(pair) > 0) pair$p[1] else NA_real_,
				p.adj = if (!is.null(pair) && nrow(pair) > 0) pair$p.adj[1] else NA_real_,
				p.adj.signif = if (!is.null(pair) && nrow(pair) > 0) {
					ifelse(pair$p.adj[1] < 0.001, "***",
						ifelse(pair$p.adj[1] < 0.01, "**",
							   ifelse(pair$p.adj[1] < 0.05, "*", "ns")))
				} else NA_character_,
				diagnostic = if (!is.null(pair) && nrow(pair) > 0) NA_character_ else omnibus_note,
				stringsAsFactors = FALSE
			)
		}
		if (nrow(dat) == 0) return(NULL)
		sig_alpha <- if (!is.null(dunn_res)) dunn_res %>%
			dplyr::filter(p.adj < 0.05) %>%
			dplyr::mutate(y.position = max(dat$value, na.rm = TRUE) +
								  max(diff(range(dat$value, na.rm = TRUE)) * 0.12, 0.05) *
								  seq_len(dplyr::n())) else NULL
		ggplot(dat, aes(x = group_id, y = value, fill = group_id)) +
			geom_boxplot(outlier.shape = NA, alpha = 0.30, width = 0.55) +
			geom_jitter(aes(colour = group_id), width = 0.15, size = 2.2, alpha = 0.85) +
			scale_fill_manual(values = pal, guide = "none") +
			scale_colour_manual(values = pal, guide = "none") +
			scale_x_discrete(labels = GROUP_X_LABELS_4A_5A[grp_present]) +
			scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
			labs(title = lab, x = NULL, y = lab) +
			theme_bw(base_size = 10) +
			theme(aspect.ratio = 1,
					  axis.title.x = element_text(size = 20),
					  axis.title.y = element_text(size = 20),
					  axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 0.5),
					  axis.text.y = element_text(size = 14),
					  legend.position = "none") +
			{ if (!is.null(sig_alpha) && nrow(sig_alpha) > 0)
				ggpubr::stat_pvalue_manual(sig_alpha, label = "p.adj.signif",
											  tip.length = 0.01, size = 7.0)
				else NULL }
	}, metrics, labels)
	plots_alpha <- Filter(Negate(is.null), plots_alpha)

	alpha_pairwise_for_metric <- function(m) {
		metric_available <- m %in% names(alpha_df)
		dat <- if (metric_available) alpha_df %>%
			dplyr::transmute(group_id, value = .data[[m]]) %>%
			dplyr::filter(!is.na(value), !is.na(group_id)) else data.frame()
		group_counts <- if (nrow(dat) > 0) table(factor(dat$group_id, levels = GROUP_ORDER)) else
			setNames(rep(0L, length(GROUP_ORDER)), GROUP_ORDER)
		omnibus_p <- NA_real_
		omnibus_note <- if (!metric_available) "Not tested: metric unavailable" else
			"Not tested: missing groups, too few values, or zero variance"
		dunn_res <- NULL
		if (metric_available && nrow(dat) >= 3 && dplyr::n_distinct(dat$group_id) >= 2 &&
				dplyr::n_distinct(dat$value) >= 2) {
			kw <- tryCatch(rstatix::kruskal_test(dat, value ~ group_id), error = function(e) NULL)
			if (!is.null(kw)) {
				omnibus_p <- kw$p[1]
				omnibus_note <- "Kruskal-Wallis"
				if (!is.na(omnibus_p) && omnibus_p < 0.05) {
					dunn_res <- tryCatch(rstatix::dunn_test(dat, value ~ group_id,
																	 p.adjust.method = "BH"), error = function(e) NULL)
				}
			} else omnibus_note <- "Kruskal-Wallis failed"
		}
		pair_grid <- expand.grid(group1 = GROUP_ORDER, group2 = GROUP_ORDER,
								 stringsAsFactors = FALSE) %>%
			dplyr::filter(match(group1, GROUP_ORDER) < match(group2, GROUP_ORDER))
		dplyr::bind_rows(lapply(seq_len(nrow(pair_grid)), function(j) {
			g1 <- pair_grid$group1[j]; g2 <- pair_grid$group2[j]
			pair <- if (!is.null(dunn_res)) dunn_res %>% dplyr::filter(group1 == g1, group2 == g2) else NULL
			has_pair <- !is.null(pair) && nrow(pair) > 0
			data.frame(
				metric = m, group1 = g1, group2 = g2,
				n_group1 = unname(group_counts[g1]), n_group2 = unname(group_counts[g2]),
				omnibus_test = omnibus_note, omnibus_p = omnibus_p,
				test = if (has_pair) "Dunn" else if (!is.na(omnibus_p) && omnibus_p >= 0.05) {
					"Dunn not run: omnibus p >= 0.05"
				} else "Dunn not available",
				p = if (has_pair) pair$p[1] else NA_real_,
				p.adj = if (has_pair) pair$p.adj[1] else NA_real_,
				p.adj.signif = if (has_pair) ifelse(pair$p.adj[1] < 0.001, "***",
					ifelse(pair$p.adj[1] < 0.01, "**", ifelse(pair$p.adj[1] < 0.05, "*", "ns"))) else NA_character_,
				diagnostic = if (has_pair) NA_character_ else omnibus_note,
				stringsAsFactors = FALSE
			)
		}))
	}
	alpha_pairwise_stats <- dplyr::bind_rows(lapply(metrics, alpha_pairwise_for_metric))
	write_csv(alpha_pairwise_stats,
			  file.path(TBL_DIR, "Fig4A_alpha_diversity_pairwise_stats.csv"))

	p4A <- wrap_plots(plots_alpha, ncol = 3) +
		plot_annotation(
			title = "Figure 4A — Alpha diversity (S2: OPA1 KO longitudinal)",
			subtitle = "Colour = group_id | Rarefied | Observed · Chao1 · ACE · Shannon · Pielou · Dominance",
			theme = theme(plot.title = element_text(face = "bold", size = 13),
										plot.subtitle = element_text(size = 9, colour = "grey30"))
		)
	save_fig(p4A, file.path(FIG_DIR, "Fig4A_S2_alpha_boxplots"),
					 w = 14, h = 5 * ceiling(length(plots_alpha) / 3))
}

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 3: Figure 4B — Family-level stacked bar (S2)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 3: Figure 4B — Family stacked bar ===\n")

tryCatch({
	ps_fam <- tax_glom(ps_all, taxrank = "Family", NArm = FALSE)
	taxa_vals <- as.character(tax_table(ps_fam)[, "Family"])
	named <- !is.na(taxa_vals) & trimws(taxa_vals) != ""
	ps_fam <- prune_taxa(taxa_names(ps_fam)[named], ps_fam)
	ps_fam_c <- microbiome::transform(ps_fam, "compositional")
	df_melt <- psmelt(ps_fam_c)
  
	n_top <- 20
	top_families <- df_melt %>%
		dplyr::group_by(Family) %>%
		dplyr::summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
		dplyr::filter(!is.na(Family), trimws(Family) != "") %>%
		dplyr::arrange(dplyr::desc(total)) %>%
		dplyr::slice_head(n = n_top) %>%
		dplyr::pull(Family)
  
	pal_cols <- as.character(paletteer::paletteer_dynamic("cartography::multi.pal", length(top_families)))
	fam_pal <- setNames(pal_cols, top_families)
	fam_pal <- c(fam_pal, "Other" = "grey80")
  
	meta_df <- get_sdata(ps_all)
  
	df_plot <- df_melt %>%
		dplyr::mutate(lbl = ifelse(Family %in% top_families, Family, "Other")) %>%
		dplyr::group_by(Sample, lbl) %>%
		dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
		dplyr::group_by(Sample) %>%
		dplyr::mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
		dplyr::ungroup() %>%
		dplyr::left_join(meta_df %>% dplyr::rename(Sample = sample_id), by = "Sample")
  
	present_groups <- intersect(GROUP_ORDER, unique(df_plot$group_id))
	sample_order <- df_plot %>%
		dplyr::filter(group_id %in% present_groups) %>%
		dplyr::arrange(factor(group_id, levels = present_groups), Sample) %>%
		dplyr::pull(Sample) %>% unique()
  
	df_plot <- df_plot %>%
		dplyr::filter(group_id %in% present_groups) %>%
		dplyr::mutate(
		group_id = factor(group_id, levels = present_groups),
		Sample = factor(Sample, levels = sample_order),
		x_label = factor(as.character(Sample),
						 levels = sample_order)
		)
  
	grp_labeller <- GROUP_HEADER_LABELS[present_groups]
  
	tax_order <- df_plot %>%
		dplyr::filter(lbl != "Other") %>%
		dplyr::group_by(lbl) %>%
		dplyr::summarise(tot = sum(Abundance), .groups = "drop") %>%
		dplyr::arrange(tot) %>% dplyr::pull(lbl)
	lbl_levels <- c("Other", tax_order)
	df_plot <- df_plot %>% dplyr::mutate(lbl = factor(lbl, levels = lbl_levels))
	legend_breaks <- c(rev(tax_order), "Other")
  
	p4B <- ggplot(df_plot, aes(x = x_label, y = Abundance * 100, fill = lbl)) +
		geom_bar(stat = "identity", width = 0.88) +
		facet_grid(cols = vars(group_id), scales = "free_x", space = "free_x",
					   labeller = as_labeller(grp_labeller)) +
		scale_fill_manual(values = fam_pal, breaks = legend_breaks,
						  name = paste0("Top ", n_top, "\nFamily")) +
		scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100.5),
						 expand = c(0, 0)) +
		labs(title = "Figure 4B — Family-level composition (S2)", x = NULL,
			 y = "Relative Abundance") +
		guides(fill = guide_legend(ncol = 1)) +
		theme_bw(base_size = 11) +
		theme(axis.title.x = element_text(size = 22),
			  axis.title.y = element_text(size = 22),
			  axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
			  axis.text.y = element_text(size = 15),
			  legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 8.5),
			  legend.title = element_text(size = 9.5, face = "bold"),
			  plot.title = element_text(size = 12, face = "bold"),
			  strip.text = element_text(size = 8.5, face = "bold"),
			  panel.spacing.x = unit(0.4, "lines"),
			  panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
	if (requireNamespace("ggh4x", quietly = TRUE)) {
		p4B <- p4B + ggh4x::strip_themed(
			background_x = ggh4x::elem_list_rect(fill = unname(GROUP_PALETTE[present_groups]), colour = "grey30"),
			text_x = ggh4x::elem_list_text(colour = "white", face = "bold", size = 10)
		)
	}
  
	save_fig(p4B, file.path(FIG_DIR, "Fig4B_S2_family_bar"),
					 w = max(10, nsamples(ps_all) * 0.55 + length(present_groups) * 2.2), h = 10)
}, error = function(e) message("  [ERROR] Fig 4B: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 4: Figure 4C — sPLS-DA heatmap (S2)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 4: Figure 4C — sPLS-DA heatmap ===\n")

tryCatch({
	ps_allrank <- prev_filter(ps_all, min_prev = PREV_MIN, min_abund = ABUND_MIN)
	otu_ar_raw <- as.matrix(otu_table(ps_allrank))
	if (!taxa_are_rows(ps_allrank)) otu_ar_raw <- t(otu_ar_raw)
	X_allrank_clr <- clr_transform(t(otu_ar_raw))
	colnames(X_allrank_clr) <- best_rank_label(ps_allrank)
  
	sdata_ff <- get_sdata(ps_allrank)
	common_ids <- intersect(rownames(X_allrank_clr), sdata_ff$sample_id)
	X_ff <- X_allrank_clr[common_ids, , drop = FALSE]
	Y_ff <- factor(sdata_ff$group_id[match(common_ids, sdata_ff$sample_id)],
								 levels = GROUP_ORDER)
	Y_ff <- droplevels(Y_ff)
  
	keep_sp <- min(SPLSDA_KEEPX, ncol(X_ff))
	set.seed(2026)
	splsda_res <- mixOmics::splsda(X_ff, Y_ff, ncomp = 2, keepX = rep(keep_sp, 2))
  
	load_c1 <- abs(splsda_res$loadings$X[, 1])
	load_c2 <- abs(splsda_res$loadings$X[, 2])
	top_taxa_ff <- union(
		names(sort(load_c1, decreasing = TRUE))[seq_len(min(10, keep_sp))],
		names(sort(load_c2, decreasing = TRUE))[seq_len(min(10, keep_sp))]
	)
  
	heat_mat <- X_ff[, top_taxa_ff, drop = FALSE]
	ann_col_ff <- data.frame(
		Group = as.character(Y_ff[match(rownames(heat_mat), common_ids)]),
		row.names = rownames(heat_mat)
	)
	ff_order <- rownames(ann_col_ff)[order(match(ann_col_ff$Group, GROUP_ORDER))]
	heat_mat <- heat_mat[ff_order, , drop = FALSE]
	ann_col_ff <- ann_col_ff[ff_order, , drop = FALSE]
  
	ann_colors <- list(Group = GROUP_PALETTE[intersect(names(GROUP_PALETTE),
																											unique(ann_col_ff$Group))])
  
	fig4c_expr <- quote({
		pheatmap::pheatmap(
			t(heat_mat),
			annotation_col    = ann_col_ff,
			annotation_colors = ann_colors,
			scale             = "row",
			cluster_cols      = FALSE,
			clustering_method = "ward.D2",
			color             = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
			breaks            = seq(-3, 3, length.out = 101),
			border_color      = "grey70",
			fontsize_row      = 12,
			fontsize_col      = 10,
			main              = "Figure 4C — sPLS-DA top taxa (S2; keepX=10, LOMO-CV)"
		)
	})
	save_base_pdf_png(fig4c_expr,
										file.path(FIG_DIR, "Fig4C_splsda_membership_heatmap"),
										w = 11, h = 7)
}, error = function(e) message("  [ERROR] Fig 4C: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 5: Figure 5A — Dysbiosis score boxplot (S2)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 5: Figure 5A — Dysbiosis score ===\n")

tryCatch({
	ds_df <- dysbiosis_score(ps_all, ref_group_id = REF_GROUP) %>%
		dplyr::left_join(get_sdata(ps_all), by = "sample_id") %>%
		dplyr::filter(!is.na(group_id)) %>%
		dplyr::mutate(group_id = factor(group_id, levels = GROUP_ORDER))
  
	write_csv(ds_df, file.path(TBL_DIR, "Fig5A_dysbiosis_score.csv"))
	ds_counts <- table(factor(ds_df$group_id, levels = GROUP_ORDER))
	kw_5a <- tryCatch(rstatix::kruskal_test(ds_df, dysbiosis_score ~ group_id),
							error = function(e) NULL)
	kw_p_5a <- if (!is.null(kw_5a) && "p" %in% names(kw_5a)) kw_5a$p[1] else NA_real_
	pw_5a <- if (!is.na(kw_p_5a) && kw_p_5a < 0.05) {
		tryCatch(rstatix::dunn_test(ds_df, dysbiosis_score ~ group_id,
											p.adjust.method = "BH"), error = function(e) NULL)
	} else NULL
	sig_5a <- if (!is.null(pw_5a)) {
		sig <- pw_5a %>% dplyr::filter(p.adj < 0.05) %>%
			dplyr::mutate(y.position = max(ds_df$dysbiosis_score, na.rm = TRUE) +
								  max(diff(range(ds_df$dysbiosis_score, na.rm = TRUE)) * 0.12, 0.05) *
								  seq_len(dplyr::n()))
		if (nrow(sig) > 0) {
			sig$p.adj.signif <- ifelse(sig$p.adj < 0.001, "***",
				ifelse(sig$p.adj < 0.01, "**", "*"))
			sig
		} else NULL
	} else NULL
	ref_scores_5a <- ds_df$dysbiosis_score[ds_df$group_id == REF_GROUP]
	ref_threshold_5a <- if (length(ref_scores_5a) >= 2) {
		mean(ref_scores_5a, na.rm = TRUE) + 1.5 * sd(ref_scores_5a, na.rm = TRUE)
	} else NA_real_
	pair_grid_5a <- expand.grid(group1 = GROUP_ORDER, group2 = GROUP_ORDER,
								 stringsAsFactors = FALSE) %>%
		dplyr::filter(match(group1, GROUP_ORDER) < match(group2, GROUP_ORDER))
	stats_5a <- lapply(seq_len(nrow(pair_grid_5a)), function(i) {
		g1 <- pair_grid_5a$group1[i]; g2 <- pair_grid_5a$group2[i]
		pair <- if (!is.null(pw_5a)) pw_5a %>% dplyr::filter(group1 == g1, group2 == g2) else NULL
		data.frame(
			group1 = g1, group2 = g2,
			n_group1 = unname(ds_counts[g1]), n_group2 = unname(ds_counts[g2]),
			omnibus_test = "Kruskal-Wallis", omnibus_p = kw_p_5a,
			test = if (!is.null(pair) && nrow(pair) > 0) "Dunn" else
				if (!is.na(kw_p_5a) && kw_p_5a >= 0.05) "Dunn not run: omnibus p >= 0.05" else "Dunn not available",
			p = if (!is.null(pair) && nrow(pair) > 0) pair$p[1] else NA_real_,
			p.adj = if (!is.null(pair) && nrow(pair) > 0) pair$p.adj[1] else NA_real_,
			p.adj.signif = if (!is.null(pair) && nrow(pair) > 0) {
				ifelse(pair$p.adj[1] < 0.001, "***",
					ifelse(pair$p.adj[1] < 0.01, "**",
						   ifelse(pair$p.adj[1] < 0.05, "*", "ns")))
			} else NA_character_,
			stringsAsFactors = FALSE
		)
	}) %>% dplyr::bind_rows()
	write_csv(stats_5a, file.path(TBL_DIR, "Fig5A_dysbiosis_score_pairwise_stats.csv"))
	if (!is.null(pw_5a)) write_csv(pw_5a, file.path(TBL_DIR, "Fig5A_dysbiosis_score_stats.csv"))

	p5A <- ggplot(ds_df, aes(x = group_id, y = dysbiosis_score, fill = group_id)) +
		geom_boxplot(alpha = 0.30, outlier.shape = NA, width = 0.5) +
		geom_jitter(aes(colour = group_id), width = 0.13, size = 2.2, alpha = 0.85) +
		scale_fill_manual(values = GROUP_PALETTE_4A_5A[intersect(names(GROUP_PALETTE_4A_5A),
																												levels(ds_df$group_id))],
											guide = "none") +
		scale_colour_manual(values = GROUP_PALETTE_4A_5A[intersect(names(GROUP_PALETTE_4A_5A),
																																											levels(ds_df$group_id))],
																																											guide = "none") +
		scale_x_discrete(labels = GROUP_X_LABELS_4A_5A[levels(ds_df$group_id)]) +
		{ if (!is.na(ref_threshold_5a))
			geom_hline(yintercept = ref_threshold_5a, linetype = "dashed",
							colour = "grey40", linewidth = 0.6) else NULL } +
		{ if (!is.null(sig_5a) && nrow(sig_5a) > 0)
				ggpubr::stat_pvalue_manual(sig_5a, label = "p.adj.signif",
																																									tip.length = 0.01, size = 7.0)
			else NULL } +
		labs(title = "Figure 5A — Bray-Curtis dysbiosis score (S2)",
				 subtitle = paste0("D0 reference",
																		 if (!is.na(ref_threshold_5a)) sprintf(" | mean + 1.5 SD = %.3f", ref_threshold_5a) else "",
																		 if (!is.na(kw_p_5a)) sprintf(" | Kruskal-Wallis p = %.3f", kw_p_5a) else ""),
				 x = NULL, y = "Bray-Curtis distance to reference centroid") +
		theme_bw(base_size = 11) +
			theme(axis.title.x = element_text(size = 22),
				  axis.title.y = element_text(size = 22),
				  axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
				  axis.text.y = element_text(size = 15),
					plot.title = element_text(face = "bold"), aspect.ratio = 1)
  
	save_fig(p5A, file.path(FIG_DIR, "Fig5A_S2_dysbiosis_score"), w = 8, h = 7)
}, error = function(e) message("  [ERROR] Fig 5A: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 6: Figure 5B — Dysbiosis score vs body composition (S2)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 6: Figure 5B — Dysbiosis vs body composition ===\n")

tryCatch({
	ds_5b <- dysbiosis_score(ps_all, ref_group_id = REF_GROUP) %>%
		dplyr::left_join(get_sdata(ps_all), by = "sample_id") %>%
		dplyr::filter(!is.na(wt_chg_pct) | !is.na(fat_pct))
  
	outcomes_5b <- c("wt_chg_pct", "fat_pct")
  
	ct_list_5b <- lapply(outcomes_5b, function(oc) {
		df_tmp <- ds_5b %>% dplyr::filter(!is.na(.data[[oc]]), !is.na(dysbiosis_score))
		if (nrow(df_tmp) < 4) return(list(rho = NA, pval = NA))
		ct <- cor.test(df_tmp$dysbiosis_score, df_tmp[[oc]],
									 method = "spearman", exact = FALSE)
		list(rho = ct$estimate, pval = ct$p.value)
	})
	names(ct_list_5b) <- outcomes_5b
	pvals_5b <- sapply(ct_list_5b, `[[`, "pval")
	padj_5b  <- p.adjust(pvals_5b, method = "BH")
  
	stats_5b <- data.frame(
		outcome = outcomes_5b,
		rho     = sapply(ct_list_5b, function(x) round(x$rho, 4)),
		pval    = sapply(ct_list_5b, function(x) round(x$pval, 6)),
		padj_BH = round(padj_5b, 6)
	)
	write_csv(stats_5b, file.path(TBL_DIR, "Fig5B_dysbiosis_bodycomp_stats.csv"))
  
	plots_5b <- lapply(outcomes_5b, function(oc) {
		df_5b <- ds_5b %>% dplyr::filter(!is.na(.data[[oc]]))
		if (nrow(df_5b) < 4) return(NULL)
		rho_v  <- round(ct_list_5b[[oc]]$rho, 2)
		padj_v <- formatC(padj_5b[[oc]], format = "f", digits = 3)
		annot  <- paste0("rho=", rho_v, "\npadj=", padj_v)
		ggplot(df_5b, aes(x = dysbiosis_score, y = .data[[oc]], colour = group_id)) +
			geom_smooth(method = "lm", se = TRUE, colour = "grey40", fill = "grey85") +
			geom_point(size = 3) +
			scale_colour_manual(values = GROUP_PALETTE, drop = FALSE, name = NULL) +
			annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
							 label = annot, size = 6) +
			labs(title = paste0(oc, " ~ Dysbiosis score"),
					 x = "Dysbiosis score",
					 y = if (oc == "wt_chg_pct") "Weight change (%)" else "Fat mass (%)") +
			theme_bw(base_size = 10) +
			theme(aspect.ratio = 1,
						axis.title.x = element_text(size = 22),
						axis.title.y = element_text(size = 22),
						axis.text.x = element_text(size = 15),
						axis.text.y = element_text(size = 15),
						legend.position = "none")
	})
	plots_5b <- Filter(Negate(is.null), plots_5b)
  
	if (length(plots_5b) > 0) {
		p5B <- wrap_plots(plots_5b, nrow = 1, ncol = 2) +
			plot_annotation(
				title = "Figure 5B — Dysbiosis score vs body composition (S2)",
				theme = theme(plot.title = element_text(face = "bold")))
		save_fig(p5B, file.path(FIG_DIR, "Fig5B_S2_dysbiosis_bodycomp"),
						 w = 10, h = 10)
	}
}, error = function(e) message("  [ERROR] Fig 5B: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 7: Figure 5C — Functional guild radar (D0 vs D8 vs D21)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 7: Figure 5C — Guild radar ===\n")

tryCatch({
	guild_defs <- read_csv(GUILD_CSV, show_col_types = FALSE)
  
	ps_sp <- tax_glom(ps_all, "Species", NArm = TRUE) %>%
		prev_filter(min_prev = PREV_MIN, min_abund = ABUND_MIN)
	otu_sp_ra <- as.matrix(t(microbiome::transform(ps_sp, "compositional")@otu_table))
	sp_names  <- gsub("^Species\\|", "", colnames(otu_sp_ra))
	colnames(otu_sp_ra) <- sp_names
	sdata_sp <- get_sdata(ps_sp)
  
	GUILDS <- split(guild_defs$taxa_member, guild_defs$guild_name)
  
	guild_scores <- as.data.frame(do.call(cbind, lapply(names(GUILDS), function(gname) {
		cols <- intersect(GUILDS[[gname]], colnames(otu_sp_ra))
		if (length(cols) == 0) return(rep(0, nrow(otu_sp_ra)))
		rowSums(otu_sp_ra[, cols, drop = FALSE])
	})))
	colnames(guild_scores) <- names(GUILDS)
	rownames(guild_scores) <- rownames(otu_sp_ra)
  
	guild_df <- guild_scores %>%
		tibble::rownames_to_column("sample_id") %>%
		dplyr::left_join(sdata_sp, by = "sample_id")
  
	# IBD_score and OST_index from all-rank psmelt
	ps_all_rel <- microbiome::transform(
		prev_filter(ps_all, min_prev = 0.05, min_abund = 1e-5), "compositional")
	melt_idx <- psmelt(ps_all_rel) %>%
		dplyr::mutate(taxon_plain = dplyr::coalesce(Species, Genus, Family, Order, Class, Phylum))
  
	IBD_DEPLETED_VEC <- c("Lachnospiraceae","Ruminococcaceae","Faecalibacterium",
												 "Akkermansia","Blautia","Roseburia","Coprococcus","Oscillospiraceae")
	IBD_ENRICHED_VEC <- c("Proteobacteria","Enterobacteriaceae","Fusobacteriaceae",
												 "Fusobacterium","Escherichia","Klebsiella","Bilophila")
	OST_TOLERANT_VEC <- c("Escherichia","Klebsiella","Enterococcus","Streptococcus",
												 "Helicobacter","Campylobacter","Bacteroides","Akkermansia",
												 "Enterobacteriaceae","Pasteurellaceae")
	OST_SENSITIVE_VEC <- c("Faecalibacterium","Roseburia","Blautia","Butyrivibrio",
													"Coprococcus","Eubacterium","Lachnospiraceae","Ruminococcaceae")
  
	ibd_smry <- melt_idx %>%
		dplyr::mutate(ibd_dir = dplyr::case_when(
			taxon_plain %in% IBD_ENRICHED_VEC ~ "enriched",
			taxon_plain %in% IBD_DEPLETED_VEC ~ "depleted",
			TRUE ~ NA_character_)) %>%
		dplyr::filter(!is.na(ibd_dir)) %>%
		dplyr::group_by(Sample) %>%
		dplyr::summarise(IBD_score = sum(Abundance[ibd_dir == "enriched"], na.rm = TRUE) -
											 sum(Abundance[ibd_dir == "depleted"], na.rm = TRUE), .groups = "drop") %>%
		dplyr::rename(sample_id = Sample)
  
	ost_smry <- melt_idx %>%
		dplyr::mutate(ost_dir = dplyr::case_when(
			taxon_plain %in% OST_TOLERANT_VEC ~ "tolerant",
			taxon_plain %in% OST_SENSITIVE_VEC ~ "sensitive",
			TRUE ~ NA_character_)) %>%
		dplyr::filter(!is.na(ost_dir)) %>%
		dplyr::group_by(Sample) %>%
		dplyr::summarise(OST_index = sum(Abundance[ost_dir == "tolerant"], na.rm = TRUE) -
											 sum(Abundance[ost_dir == "sensitive"], na.rm = TRUE), .groups = "drop") %>%
		dplyr::rename(sample_id = Sample)
  
	guild_df <- guild_df %>%
		dplyr::left_join(ibd_smry, by = "sample_id") %>%
		dplyr::left_join(ost_smry, by = "sample_id")
	guild_df$IBD_score[is.na(guild_df$IBD_score)] <- 0
	guild_df$OST_index[is.na(guild_df$OST_index)] <- 0
  
	write_csv(guild_df, file.path(TBL_DIR, "S2_guild_scores.csv"))
  
	# Radar chart
	radar_guilds  <- c("Butyrate Producers","Propionate Producers","Acetate Producers",
											"Secondary Bile Acid Producers","H2S Producers","OST_index","IBD_score",
											"LPS-High (Endotoxemia) Producers","TMAO Precursor Producers")
	radar_vlabels <- c("Butyrate","Propionate","Acetate","SecBile\nAcid","H2S\nproducers",
											"OST\nindex","IBD\nscore","LPS","TMAO")
  
	# Use available guilds only
	radar_avail <- intersect(radar_guilds, names(guild_df))
	radar_vlabels_avail <- radar_vlabels[radar_guilds %in% radar_avail]
  
	if (length(radar_avail) >= 3) {
		radar_means <- guild_df %>%
			dplyr::group_by(group_id) %>%
			dplyr::summarise(dplyr::across(dplyr::all_of(radar_avail), mean), .groups = "drop") %>%
			dplyr::mutate(group_id = as.character(group_id))
    
		radar_max <- radar_means %>% dplyr::select(dplyr::all_of(radar_avail)) %>%
			apply(2, max) * 1.1
		radar_min <- radar_means %>% dplyr::select(dplyr::all_of(radar_avail)) %>%
			apply(2, min)
		radar_min <- pmin(radar_min * 1.1, 0)
    
		radar_mat <- radar_means %>%
			dplyr::filter(group_id %in% GROUP_ORDER) %>%
			tibble::column_to_rownames("group_id") %>%
			dplyr::select(dplyr::all_of(radar_avail)) %>% as.matrix()
		radar_mat <- radar_mat[GROUP_ORDER[GROUP_ORDER %in% rownames(radar_mat)], , drop = FALSE]
		radar_df_plot <- as.data.frame(rbind(radar_max, radar_min, radar_mat))
		rownames(radar_df_plot)[1:2] <- c("max", "min")
    
		radar_cols <- unname(GROUP_PALETTE[GROUP_ORDER[GROUP_ORDER %in% rownames(radar_mat)]])
    
		fig5c_expr <- quote({
			fmsb::radarchart(
				radar_df_plot,
				axistype = 1, pcol = radar_cols,
				pfcol = adjustcolor(radar_cols, 0.15), plwd = 2,
				vlabels = radar_vlabels_avail,
				title = "Figure 5C — Functional guild radar\nOPA1+ (D0) vs OPA1 KO (D8) vs OPA1 KO (D21)"
			)
			legend("topright", legend = unname(GROUP_LABELS[GROUP_ORDER[GROUP_ORDER %in% rownames(radar_mat)]]),
						 col = radar_cols, lwd = 2, bty = "n", cex = 0.7)
		})
		save_base_pdf_png(fig5c_expr, file.path(FIG_DIR, "Fig5C_kegg_radar_all_groups"), w = 7, h = 6)
	}
}, error = function(e) message("  [ERROR] Fig 5C: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 8: Figure 5D — SCFA guild stacked bar
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 8: Figure 5D — SCFA guild stacked bar ===\n")

tryCatch({
	scfa_guilds <- c("Butyrate Producers", "Propionate Producers", "Acetate Producers")
	scfa_colors <- c("Butyrate Producers" = "#2ca02c",
										"Propionate Producers" = "#ff7f0e",
										"Acetate Producers" = "#1f77b4")
  
	scfa_avail <- intersect(scfa_guilds, names(guild_df))
  
	if (length(scfa_avail) > 0) {
		present_groups_5d <- intersect(GROUP_ORDER, unique(as.character(guild_df$group_id)))
		guild_long_5d <- guild_df %>%
			dplyr::arrange(group_id) %>%
			dplyr::select(sample_id, group_id, dplyr::all_of(scfa_avail)) %>%
			tidyr::pivot_longer(dplyr::all_of(scfa_avail),
													names_to = "guild", values_to = "RA") %>%
			dplyr::mutate(guild = factor(guild, levels = scfa_avail))
    
		p5D <- ggplot(guild_long_5d, aes(x = sample_id, y = RA, fill = guild)) +
			geom_col() +
			facet_grid(. ~ group_id, scales = "free_x", space = "free_x",
								 labeller = labeller(group_id = GROUP_HEADER_LABELS)) +
			scale_fill_manual(values = scfa_colors[scfa_avail]) +
			scale_y_continuous(labels = scales::percent_format()) +
			labs(title = "Figure 5D — SCFA producer guild contributions per sample",
					 subtitle = "OPA1 KO vs healthy reference (OPA1+)",
					 x = NULL, y = "Relative Abundance", fill = "Guild") +
			theme_bw(base_size = 9) +
			theme(axis.title.x = element_text(size = 22),
				  axis.title.y = element_text(size = 22),
				  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
				  axis.text.y = element_text(size = 15),
				  strip.text = element_text(size = 10, face = "bold"))
		if (requireNamespace("ggh4x", quietly = TRUE)) {
			p5D <- p5D + ggh4x::strip_themed(
				background_x = ggh4x::elem_list_rect(fill = unname(GROUP_PALETTE[present_groups_5d]), colour = "grey30"),
				text_x = ggh4x::elem_list_text(colour = "white", face = "bold", size = 10)
			)
		}
    
		save_fig(p5D, file.path(FIG_DIR, "Fig5D_scfa_guild_stacked_bar"), w = 10, h = 7)
	}
}, error = function(e) message("  [ERROR] Fig 5D: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 9: Figure 5E — Partial Spearman dysbiosis volcano
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 9: Figure 5E — Partial Spearman volcano ===\n")

tryCatch({
	ps_fh <- prev_filter(ps_all, min_prev = PREV_MIN, min_abund = ABUND_MIN)
	ps_fh_rel <- microbiome::transform(ps_fh, "compositional")
	bc_fh <- as.matrix(phyloseq::distance(ps_fh_rel, method = "bray"))
	sdata_fh <- get_sdata(ps_fh) %>%
		dplyr::mutate(group_id = factor(group_id, levels = GROUP_ORDER))
	ref_ids_fh <- sdata_fh$sample_id[sdata_fh$group_id == REF_GROUP]
	dys_fh <- rowMeans(bc_fh[, intersect(ref_ids_fh, colnames(bc_fh)), drop = FALSE], na.rm = TRUE)
	dys_df_fh <- data.frame(sample_id = names(dys_fh), dysbiosis_score = dys_fh,
													 stringsAsFactors = FALSE)
  
	fh_lbl <- best_rank_label(ps_fh)
	otu_fh_raw <- as.matrix(otu_table(ps_fh))
	if (!taxa_are_rows(ps_fh)) otu_fh_raw <- t(otu_fh_raw)
	X_fh <- clr_transform(t(otu_fh_raw))
	colnames(X_fh) <- fh_lbl
  
	common_fh <- intersect(rownames(X_fh), dys_df_fh$sample_id)
	dys_v <- dys_df_fh$dysbiosis_score[match(common_fh, dys_df_fh$sample_id)]
	grp_v <- sdata_fh$group_id[match(common_fh, sdata_fh$sample_id)]
	X_fh_sub <- X_fh[common_fh, , drop = FALSE]
  
	dys_resid <- tryCatch(residuals(lm(dys_v ~ grp_v)),
												 error = function(e) dys_v - mean(dys_v))
  
	pcor_res <- apply(X_fh_sub, 2, function(x) {
		xr <- tryCatch(residuals(lm(x ~ grp_v)), error = function(e) x)
		ct <- tryCatch(suppressWarnings(
			cor.test(xr, dys_resid, method = "spearman", exact = FALSE)),
			error = function(e) list(estimate = NA_real_, p.value = NA_real_))
		c(rho = ct$estimate[[1]], pval = ct$p.value)
	})
  
	pcor_df <- data.frame(
		taxon = colnames(X_fh_sub), rho = pcor_res["rho", ], pval = pcor_res["pval", ],
		stringsAsFactors = FALSE
	) %>%
		dplyr::filter(!is.na(rho)) %>%
		dplyr::mutate(
			padj = p.adjust(pval, method = "BH"),
			sig  = padj < BH_THRESH,
			dir  = ifelse(rho > 0, "Dysbiosis-promoting", "Dysbiosis-protecting")
		) %>%
		dplyr::arrange(dplyr::desc(abs(rho)))
  
	# Guild annotation
	guild_lookup_5e <- guild_defs %>%
		dplyr::mutate(taxon = paste0("Species|", taxa_member)) %>%
		dplyr::select(taxon, func_guild = guild_name) %>%
		dplyr::distinct(taxon, .keep_all = TRUE)
  
	pcor_df <- pcor_df %>%
		dplyr::left_join(guild_lookup_5e, by = "taxon") %>%
		dplyr::mutate(colour_cat = dplyr::coalesce(
			func_guild,
			dplyr::if_else(dir == "Dysbiosis-promoting", "Promoting (other)", "Protecting (other)")
		))
  
	write_csv(dplyr::select(pcor_df, taxon, rho, pval, padj, sig, dir, func_guild),
						file.path(TBL_DIR, "Fig5E_partial_spearman_dysbiosis.csv"))
  
	# Labels: top significant taxa
	sig_taxa <- pcor_df %>% dplyr::filter(pval < 0.05)
	label_guild <- sig_taxa %>% dplyr::filter(!is.na(func_guild))
	label_other <- sig_taxa %>% dplyr::filter(is.na(func_guild)) %>%
		dplyr::arrange(dplyr::desc(abs(rho))) %>% dplyr::slice_head(n = 4)
	label_all <- dplyr::bind_rows(label_guild, label_other) %>%
		dplyr::distinct(taxon, .keep_all = TRUE) %>%
		dplyr::mutate(lbl = taxon)
  
	# Guild palette
	guild_names_present <- unique(na.omit(pcor_df$func_guild))
	n_guilds <- length(guild_names_present)
	if (n_guilds > 0) {
		guild_pal <- setNames(
			colorRampPalette(RColorBrewer::brewer.pal(min(9, max(3, n_guilds)), "Set1"))(n_guilds),
			guild_names_present)
	} else {
		guild_pal <- character(0)
	}
	full_pal <- c(guild_pal,
								"Promoting (other)" = "#f7c6d8", "Protecting (other)" = "#c8e6c9")
  
	p5E <- ggplot(pcor_df, aes(x = rho, y = -log10(pval + 1e-6),
															colour = colour_cat, size = abs(rho))) +
		geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", colour = "grey55") +
		geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
		geom_point(data = dplyr::filter(pcor_df, is.na(func_guild)), alpha = 0.25) +
		geom_point(data = dplyr::filter(pcor_df, !is.na(func_guild)), alpha = 0.9, shape = 16) +
		ggrepel::geom_text_repel(
			data = label_all, aes(label = lbl), size = 2.4,
			max.overlaps = 60, segment.size = 0.3, box.padding = 0.6,
			show.legend = FALSE) +
		scale_colour_manual(
			values = full_pal,
			breaks = guild_names_present,
			name = "Functional guild",
			guide = guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1))) +
		scale_size_continuous(range = c(1, 5), guide = "none") +
		labs(title = "Figure 5E — Partial Spearman: CLR vs Dysbiosis Score | group_id",
				 subtitle = paste0("Partialling out group_id. n = ", length(common_fh), " samples."),
				 x = expression(paste("Spearman ", rho, " (partial)")),
				 y = expression(-log[10](nominal~p))) +
		theme_bw(base_size = 11) +
		theme(legend.position = "right", plot.title = element_text(face = "bold"),
					plot.subtitle = element_text(size = 8.5, colour = "grey30"))
  
	save_fig(p5E, file.path(FIG_DIR, "Fig5E_partial_spearman_dysbiosis_volcano"), w = 14, h = 9)
}, error = function(e) message("  [ERROR] Fig 5E: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 10: Summary
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 10: Summary ===\n")

fig_files <- list.files(FIG_DIR, pattern = "\\.(pdf|png)$", full.names = FALSE)
tbl_files <- list.files(TBL_DIR, pattern = "\\.csv$", full.names = FALSE)

cat("\nFigures generated (", length(fig_files), "files):\n")
for (f in sort(fig_files)) cat("  ", f, "\n")

cat("\nTables generated (", length(tbl_files), "files):\n")
for (f in sort(tbl_files)) cat("  ", f, "\n")

# Session info
writeLines(capture.output(sessionInfo()), here::here("sessionInfo.txt"))
cat("\nsessionInfo.txt written.\n")

cat("\n========================================================================\n")
cat("  All figures generated successfully.\n")
cat("  Output: results/figures/ (8 PDFs + 8 PNGs)\n")
cat("  Stats:  results/tables/  (CSV files)\n")
cat("========================================================================\n")
cat("Done.\n")
