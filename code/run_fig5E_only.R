#!/usr/bin/env Rscript
# ==============================================================================
# run_fig5E_only.R
#
# Targeted runner for Figure 5E (partial Spearman dysbiosis volcano) ONLY.
# Self-contained extract of code/00_generate_all_figures.R — does not source
# or re-run the main script, and does not touch the other 7 figure sections.
#
# Produces:
#   results/figures/Fig5E_partial_spearman_dysbiosis_volcano.pdf
#   results/figures/Fig5E_partial_spearman_dysbiosis_volcano.png
#   results/tables/Fig5E_partial_spearman_dysbiosis.csv
#
# NOTE ON SCOPE: as currently written, SECTION 9 of the main script builds its
# own independent functional-guild lookup (fh_mk()/fh_guild_lookup) and does
# NOT reference guild_defs/guild_df from SECTION 7 (verified by inspection of
# 00_generate_all_figures.R). SECTION 7's guild-score setup (tax_glom to
# Species, compositional transform, all-rank psmelt for ROS_sensitive) is
# therefore intentionally OMITTED here — it is the most expensive step in the
# full pipeline and is not needed to reproduce Figure 5E. Only SECTION 0
# (setup/helpers), SECTION 1 (build phyloseq), and the SECTION 4 sPLS-DA fit
# (for top_taxa_ff) are run as prerequisites.
#
# Usage:
#   setwd("<repo_root>/opa1-intestinal-microbiome")
#   source("code/run_fig5E_only.R")
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 0 (extracted): Setup — packages, seed, config, helper functions
# ──────────────────────────────────────────────────────────────────────────────
cat("=== SECTION 0: Setup ===\n")

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

if (!requireNamespace("mixOmics", quietly = TRUE))
	stop("mixOmics not installed. Run: BiocManager::install('mixOmics')")
suppressPackageStartupMessages(library(mixOmics))

set.seed(2026)

META_PATH       <- here::here("config", "sample_metadata", "opa1_drp_meta_s2.csv")
REPORT_DIR      <- here::here("data", "kraken_reports_core_nt")
GUILD_CSV       <- here::here("data", "reference", "guild_definitions.csv")
FIG_DIR         <- here::here("results", "figures")
TBL_DIR         <- here::here("results", "tables")
REPORT_SUFFIX   <- ".k2report"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

PREV_MIN     <- 0.10
ABUND_MIN    <- 1e-4
CLR_PSEUDO   <- 0.5
SPLSDA_KEEPX <- 10
BH_THRESH    <- 0.2
DPI          <- 300

GROUP_ORDER <- c(
	"untreated_opa1floxed_vilcre_d0",
	"opa1floxed_vilcre_tam_d8",
	"opa1floxed_vilcre_tam_d21"
)
REF_GROUP <- "untreated_opa1floxed_vilcre_d0"

cat("  Paths, parameters, and palettes configured.\n")

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

clr_transform <- function(X, pseudocount = CLR_PSEUDO) {
	lX <- log(X + pseudocount); lX - rowMeans(lX)
}

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

get_sdata <- function(ps) {
	df <- as(sample_data(ps), "data.frame")
	tibble::rownames_to_column(df, "sample_id")
}

save_fig <- function(plot_obj, base_path, w, h, dpi = DPI) {
	ggsave(paste0(base_path, ".pdf"), plot_obj, width = w, height = h)
	ggsave(paste0(base_path, ".png"), plot_obj, width = w, height = h, dpi = dpi)
	cat("  Saved:", basename(base_path), "(.pdf + .png)\n")
}

cat("  Helper functions defined.\n")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 1 (extracted): Build S2 phyloseq object
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

sd <- as(sample_data(ps_all), "data.frame")
sd$group_id <- factor(sd$group_id, levels = GROUP_ORDER, ordered = FALSE)
sample_data(ps_all) <- sample_data(sd)

cat("  S2 phyloseq ready:", nsamples(ps_all), "samples,", ntaxa(ps_all), "taxa\n")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 4 (extracted, partial): sPLS-DA fit -> top_taxa_ff only
# (pheatmap and tabular exports skipped — not needed for Figure 5E)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== SECTION 4 (partial): sPLS-DA fit for ML top taxa ===\n")

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
	cat("  top_taxa_ff computed:", length(top_taxa_ff), "taxa\n")
}, error = function(e) message("  [WARN] SECTION 4 (partial) sPLS-DA: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 9 (verbatim): Figure 5E — Partial Spearman dysbiosis volcano
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

	# Guild annotation: multi-rank lookup reproduced from Module FH4 in
	# 09_opa1_drp_core_nt_functional_analysis.R (Species/Genus/Family/Phylum;
	# distinct(taxon) keeps the first, highest-priority guild match).
	fh_mk <- function(rank, nms, guild)
		data.frame(taxon = paste0(rank, "|", nms), func_guild = guild,
							 stringsAsFactors = FALSE)

	fh_guild_lookup <- dplyr::bind_rows(
		# Acetate guild
		fh_mk("Species", c(
			"Bifidobacterium longum", "Bifidobacterium bifidum",
			"Bifidobacterium adolescentis", "Bifidobacterium breve",
			"Lactobacillus reuteri", "Lactobacillus acidophilus",
			"Lactobacillus rhamnosus", "Lactobacillus johnsonii",
			"Blautia obeum", "Blautia wexlerae"
		), "Acetate guild"),
		fh_mk("Genus", c("Bifidobacterium", "Lactobacillus"),
					"Acetate guild"),
		fh_mk("Family", "Bifidobacteriaceae", "Acetate guild"),

		# Propionate guild
		fh_mk("Species", c(
			"Bacteroides thetaiotaomicron", "Bacteroides uniformis",
			"Bacteroides ovatus", "Bacteroides vulgatus", "Bacteroides dorei",
			"Phascolarctobacterium succinatutens", "Dialister succinatiphilus",
			"Veillonella parvula", "Veillonella atypica",
			"Akkermansia muciniphila", "Alistipes putredinis", "Prevotella copri"
		), "Propionate guild"),
		fh_mk("Genus", c("Bacteroides", "Phascolarctobacterium", "Dialister",
											"Veillonella", "Akkermansia", "Alistipes", "Prevotella"),
					"Propionate guild"),
		fh_mk("Family", "Veillonellaceae", "Propionate guild"),

		# Butyrate guild
		fh_mk("Species", c(
			"Faecalibacterium prausnitzii",
			"Roseburia intestinalis", "Roseburia hominis", "Roseburia inulinivorans",
			"Butyrivibrio fibrisolvens", "Butyrivibrio proteoclasticus",
			"Coprococcus comes", "Coprococcus catus", "Coprococcus eutactus",
			"Eubacterium hallii", "Eubacterium rectale",
			"Anaerostipes caccae", "Anaerostipes hadrus",
			"Clostridium butyricum", "Clostridium tyrobutyricum",
			"Subdoligranulum variabile", "Butyricicoccus pullicaecorum"
		), "Butyrate guild"),
		fh_mk("Genus", c("Faecalibacterium", "Roseburia", "Butyrivibrio",
											"Coprococcus", "Anaerostipes", "Subdoligranulum"),
					"Butyrate guild"),

		# Lachnospiraceae / Butyrate
		fh_mk("Species", c(
			"Blautia obeum", "Blautia producta", "Blautia hansenii", "Blautia wexlerae",
			"Lachnoclostridium phytofermentans",
			"Dorea longicatena", "Dorea formicigenerans"
		), "Lachnospiraceae / Butyrate"),
		fh_mk("Genus", c("Blautia", "Dorea", "Lachnoclostridium"),
					"Lachnospiraceae / Butyrate"),
		fh_mk("Family", "Lachnospiraceae", "Lachnospiraceae / Butyrate"),

		# SecBileAcid (7-alpha-dehydroxylation; secondary bile acid producers)
		fh_mk("Species", c(
			"Clostridium scindens", "Clostridium hylemonae",
			"Clostridium hiranonis", "Clostridium absonum",
			"Ruminococcus gnavus",
			"Eggerthella lenta"
		), "SecBileAcid"),
		fh_mk("Genus", c("Eggerthella"), "SecBileAcid"),
		fh_mk("Family", "Eggerthellaceae", "SecBileAcid"),

		# IBD-depleted
		fh_mk("Species", c(
			"Ruminococcus bromii", "Ruminococcus champanellensis",
			"Ruminococcus albus", "Ruminococcus callidus",
			"Fusicatenibacter saccharivorans",
			"Anaerotruncus colihominis", "Roseburia faecis"
		), "IBD-depleted"),
		fh_mk("Genus", c("Ruminococcus", "Fusicatenibacter", "Anaerotruncus"),
					"IBD-depleted"),
		fh_mk("Family", c("Ruminococcaceae", "Oscillospiraceae"), "IBD-depleted"),

		# IBD-enriched
		fh_mk("Species", c(
			"Escherichia coli", "Klebsiella pneumoniae", "Klebsiella oxytoca",
			"Fusobacterium nucleatum", "Fusobacterium varium",
			"Peptostreptococcus anaerobius", "Bilophila wadsworthensis",
			"Enterobacter cloacae"
		), "IBD-enriched"),
		fh_mk("Genus", c("Escherichia", "Klebsiella", "Fusobacterium",
											"Bilophila", "Peptostreptococcus", "Enterobacter",
											"Desulfovibrio"),
					"IBD-enriched"),
		fh_mk("Family", c("Fusobacteriaceae", "Enterobacteriaceae"), "IBD-enriched"),
		fh_mk("Phylum", c("Fusobacteria"), "IBD-enriched"),

		# ROS sensitive: strict anaerobes from Oscillospiraceae not assigned above
		fh_mk("Species", c(
			"Oscillibacter valericigenes", "Oscillibacter ruminantium",
			"Eubacterium limosum", "Eubacterium yurii",
			"Pseudobutyrivibrio ruminis", "Pseudobutyrivibrio xylanivorans",
			"Intestinimonas butyriciproducens"
		), "ROS sensitive"),
		fh_mk("Genus", c("Oscillibacter", "Pseudobutyrivibrio", "Intestinimonas",
											"Eubacterium"), "ROS sensitive"),
		fh_mk("Family", "Oscillospiraceae", "ROS sensitive"),

		# ROS tolerant: aerotolerant taxa not in IBD-enriched
		fh_mk("Species", c(
			"Enterococcus faecalis", "Enterococcus faecium",
			"Helicobacter hepaticus", "Helicobacter pylori",
			"Campylobacter jejuni", "Campylobacter coli",
			"Staphylococcus aureus", "Staphylococcus epidermidis",
			"Streptococcus thermophilus"
		), "ROS tolerant"),
		fh_mk("Genus", c("Enterococcus", "Helicobacter", "Campylobacter",
											"Staphylococcus", "Streptococcus"), "ROS tolerant"),
		fh_mk("Family", c("Enterococcaceae", "Helicobacteraceae",
											"Campylobacteraceae"), "ROS tolerant"),

		# DysbiosisBloom: pathobionts not already in IBD-enriched or ROS tolerant
		fh_mk("Species", c(
			"Proteus mirabilis", "Morganella morganii",
			"Pseudomonas aeruginosa",
			"Citrobacter rodentium", "Citrobacter freundii",
			"Enterobacter aerogenes", "Enterobacter hormaechei",
			"Klebsiella variicola"
		), "DysbiosisBloom"),
		fh_mk("Genus", c("Proteus", "Morganella", "Pseudomonas", "Citrobacter",
											"Serratia", "Acinetobacter"), "DysbiosisBloom"),
		fh_mk("Family", "Pseudomonadaceae", "DysbiosisBloom"),
		fh_mk("Order",  "Pseudomonadales",  "DysbiosisBloom"),

		# H2S_producers: sulfate/sulfite reducers
		fh_mk("Species", c(
			"Desulfovibrio piger", "Desulfovibrio desulfuricans",
			"Desulfobulbus propionicus",
			"Desulfotomaculum ruminis",
			"Bilophila wadsworthensis",
			"Fusobacterium nucleatum"
		), "H2S_producers"),
		fh_mk("Genus", c("Desulfovibrio", "Desulfobulbus", "Desulfotomaculum"), "H2S_producers"),
		fh_mk("Family", "Desulfovibrionaceae", "H2S_producers")
	) %>%
		dplyr::distinct(taxon, .keep_all = TRUE)

	# Legend order/colour: "ML top taxa" first, then the 11 FH4 guilds
	fh_guild_order <- c(
		"ML top taxa",
		"Acetate guild", "Propionate guild", "Butyrate guild",
		"Lachnospiraceae / Butyrate", "SecBileAcid", "H2S_producers",
		"IBD-depleted", "IBD-enriched",
		"ROS sensitive", "ROS tolerant", "DysbiosisBloom"
	)
	fh_guild_pal <- c(
		"ML top taxa" = "#3f3f3f",
		setNames(
			colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(fh_guild_order) - 1),
			fh_guild_order[-1]
		)
	)
	fh_guild_pal["IBD-depleted"] <- "#E67300"

	# ML top taxa: prefer the in-memory SECTION 4 sPLS-DA selection; fall back
	# to the persisted Fig4C loadings table if SECTION 4 did not run/complete.
	if (exists("top_taxa_ff", inherits = TRUE) && length(top_taxa_ff) > 0) {
		ml_top_taxa <- top_taxa_ff
		ml_taxa_source <- "in-memory top_taxa_ff (SECTION 4)"
	} else {
		fig4c_loadings_path <- file.path(TBL_DIR, "Fig4C_splsda_loadings.csv")
		if (file.exists(fig4c_loadings_path)) {
			fig4c_loadings <- read_csv(fig4c_loadings_path, show_col_types = FALSE)
			ml_top_taxa <- fig4c_loadings$taxon[fig4c_loadings$in_heatmap == TRUE]
			ml_taxa_source <- "Fig4C_splsda_loadings.csv (in_heatmap == TRUE) fallback"
		} else {
			ml_top_taxa <- character(0)
			ml_taxa_source <- "unavailable: no in-memory top_taxa_ff and no Fig4C_splsda_loadings.csv"
		}
	}
	cat("  [Fig 5E] ML top taxa source:", ml_taxa_source,
			"(", length(ml_top_taxa), "taxa)\n")

	# Join guilds; "ML top taxa" is assigned only where no guild already matched.
	pcor_df <- pcor_df %>%
		dplyr::left_join(fh_guild_lookup, by = "taxon") %>%
		dplyr::mutate(
			func_guild = dplyr::if_else(is.na(func_guild) & taxon %in% ml_top_taxa,
																	 "ML top taxa", func_guild),
			in_ml_top_taxa = taxon %in% ml_top_taxa,
			colour_cat = dplyr::coalesce(
				func_guild,
				dplyr::if_else(dir == "Dysbiosis-promoting", "Promoting (other)", "Protecting (other)")
			)
		)

	# Labeling uses BH-adjusted p (padj); BH_THRESH stays 0.2 and is only used
	# for the `sig` column above, not for which points get a text label here.
	SIG_THRESH_LABEL <- 0.05

	# func_guild is only "ML top taxa" for points with no prior guild match
	# (see the mutate above), so label_ml uses that value directly rather than
	# is.na(func_guild), which would now always be FALSE for those points.
	label_guild <- pcor_df %>%
		dplyr::filter(padj < SIG_THRESH_LABEL, !is.na(func_guild), func_guild != "ML top taxa") %>%
		dplyr::mutate(lbl = taxon)
	label_ml <- pcor_df %>%
		dplyr::filter(padj < SIG_THRESH_LABEL, func_guild == "ML top taxa") %>%
		dplyr::mutate(lbl = taxon)
	label_top5_promoting <- pcor_df %>%
		dplyr::filter(padj < SIG_THRESH_LABEL, is.na(func_guild), !in_ml_top_taxa, rho > 0) %>%
		dplyr::arrange(dplyr::desc(rho)) %>% dplyr::slice_head(n = 5)
	label_top5_protecting <- pcor_df %>%
		dplyr::filter(padj < SIG_THRESH_LABEL, is.na(func_guild), !in_ml_top_taxa, rho < 0) %>%
		dplyr::arrange(rho) %>% dplyr::slice_head(n = 5)
	label_top5 <- dplyr::bind_rows(label_top5_promoting, label_top5_protecting) %>%
		dplyr::mutate(lbl = taxon)

	taxa_to_exclude <- c("Species|Enterocloster clostridioformis", "Species|Vescimonas coprocola")
	label_guild <- label_guild %>% dplyr::filter(!taxon %in% taxa_to_exclude)
	label_ml <- label_ml %>% dplyr::filter(!taxon %in% taxa_to_exclude)
	label_top5 <- label_top5 %>% dplyr::filter(!taxon %in% taxa_to_exclude)

	pcor_df <- pcor_df %>%
		dplyr::mutate(label_set = dplyr::case_when(
			taxon %in% label_guild$taxon ~ "guild",
			taxon %in% label_ml$taxon    ~ "ML top taxa",
			taxon %in% label_top5$taxon  ~ "top5",
			TRUE ~ NA_character_
		))

	write_csv(dplyr::select(pcor_df, taxon, rho, pval, padj, sig, dir, func_guild,
											 in_ml_top_taxa, colour_cat, label_set),
						file.path(TBL_DIR, "Fig5E_partial_spearman_dysbiosis.csv"))

	# Palette: fh_guild_pal (11 guilds + "ML top taxa") plus subdued background colours
	full_pal <- c(fh_guild_pal,
								"Promoting (other)" = "#f7c6d8", "Protecting (other)" = "#c8e6c9")

	# Two-arm interior labels: split each label set by sign of rho so labels
	# repel outward from the plot centre instead of overlapping near x = 0.
	label_guild_L <- dplyr::filter(label_guild, rho < 0)
	label_guild_R <- dplyr::filter(label_guild, rho > 0)
	label_ml_L    <- dplyr::filter(label_ml, rho < 0)
	label_ml_R    <- dplyr::filter(label_ml, rho > 0)
	label_top5_L  <- dplyr::filter(label_top5, rho < 0)
	label_top5_R  <- dplyr::filter(label_top5, rho > 0)

	p5E <- ggplot(pcor_df, aes(x = rho, y = -log10(padj + 1e-6),
															colour = colour_cat, size = abs(rho))) +
		geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", colour = "grey55") +
		geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
		geom_point(data = dplyr::filter(pcor_df, is.na(func_guild), padj >= SIG_THRESH_LABEL),
							 alpha = 0.15) +
		geom_point(data = dplyr::filter(pcor_df, is.na(func_guild), padj < SIG_THRESH_LABEL),
							 alpha = 0.45) +
		geom_point(data = dplyr::filter(pcor_df, func_guild == "ML top taxa",
																		 padj < SIG_THRESH_LABEL),
							 alpha = 0.92, shape = 16) +
		geom_point(data = dplyr::filter(pcor_df, !is.na(func_guild), func_guild != "ML top taxa",
																		 padj < SIG_THRESH_LABEL),
							 alpha = 0.92, shape = 16) +
		ggrepel::geom_text_repel(
			data = label_guild_L, aes(label = lbl),
			nudge_x = 0.10, hjust = 0, direction = "both", xlim = c(-0.45, -0.10),
			force = 4, box.padding = unit(0.8, "lines"),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		ggrepel::geom_text_repel(
			data = label_guild_R, aes(label = lbl),
			nudge_x = -0.10, hjust = 1, direction = "y", xlim = c(0.10, 0.45),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		ggrepel::geom_text_repel(
			data = label_ml_L, aes(label = lbl), colour = "#3f3f3f",
			nudge_x = 0.10, hjust = 0, direction = "both", xlim = c(-0.45, -0.10),
			force = 4, box.padding = unit(0.8, "lines"),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		ggrepel::geom_text_repel(
			data = label_ml_R, aes(label = lbl), colour = "#3f3f3f",
			nudge_x = -0.10, hjust = 1, direction = "y", xlim = c(0.10, 0.45),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		ggrepel::geom_text_repel(
			data = label_top5_L, aes(label = lbl), colour = "black",
			nudge_x = 0.10, hjust = 0, direction = "both", xlim = c(-0.45, -0.10),
			force = 4, box.padding = unit(0.8, "lines"),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		ggrepel::geom_text_repel(
			data = label_top5_R, aes(label = lbl), colour = "black",
			nudge_x = -0.10, hjust = 1, direction = "y", xlim = c(0.10, 0.45),
			force_pull = 0, min.segment.length = 0, max.overlaps = Inf, seed = 2026,
			size = 2.4, show.legend = FALSE) +
		scale_colour_manual(
			values = full_pal,
			breaks = fh_guild_order,
			name = "Functional guild",
			guide = guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1))) +
		scale_size_continuous(range = c(1, 5), guide = "none") +
		labs(title = "Figure 5E — Partial Spearman: CLR vs Dysbiosis Score | group_id",
				 subtitle = paste0(
					 "CLR-taxon and dysbiosis-score residuals after regressing out group_id (partial Spearman). ",
					 "y = -log10(BH-adj. p); labels = guild members + ML top-20 taxa + top 5 promoting/protecting ",
					 "(visual threshold padj < 0.05; CSV `sig` column uses padj < 0.2). n = ", length(common_fh),
					 " samples."),
				 x = expression(paste("Spearman ", rho, " (partial)")),
				 y = expression(-log[10](BH~adj.~p))) +
		theme_bw(base_size = 11) +
		theme(legend.position = "right",
					legend.key.size = unit(0.65, "cm"),
					legend.text = element_text(size = 12),
					legend.title = element_text(size = 14, face = "bold"),
					axis.title.x = element_text(size = 22),
					axis.title.y = element_text(size = 22),
					axis.text.x = element_text(size = 16),
					axis.text.y = element_text(size = 16),
					plot.title = element_text(face = "bold"),
					plot.subtitle = element_text(size = 8.5, colour = "grey30"))

	save_fig(p5E, file.path(FIG_DIR, "Fig5E_partial_spearman_dysbiosis_volcano"), w = 18, h = 9)
}, error = function(e) message("  [ERROR] Fig 5E: ", e$message))

# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")

fig5e_files <- c(
	file.path(FIG_DIR, "Fig5E_partial_spearman_dysbiosis_volcano.pdf"),
	file.path(FIG_DIR, "Fig5E_partial_spearman_dysbiosis_volcano.png"),
	file.path(TBL_DIR, "Fig5E_partial_spearman_dysbiosis.csv")
)
for (f in fig5e_files) {
	status <- if (file.exists(f) && file.info(f)$size > 0) "OK" else "MISSING/EMPTY"
	cat("  [", status, "]", f, "\n")
}

writeLines(capture.output(sessionInfo()), here::here("sessionInfo.txt"))
cat("\nsessionInfo.txt written.\nDone.\n")
