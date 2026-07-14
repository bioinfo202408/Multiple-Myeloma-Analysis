#!/usr/bin/env Rscript
# ============================================================================
# State 2 / State 3 two-stage analysis in GSE117847
#
# Biological hypotheses (pre-specified from scRNA-seq programs):
#   H1: State 2 (adhesion/niche-interaction program) is elevated at SMM
#       diagnosis in patients who subsequently progress to MM (P-SMM).
#   H2: State 3 (proliferation/replication-stress program) increases within
#       progressor patients during the transition from baseline SMM to paired MM.
#
# IMPORTANT:
# - Do not merge baseline P-SMM and paired MM as independent observations.
# - Baseline comparison uses NP-SMM vs P-SMM only.
# - Longitudinal comparison uses matched P-SMM -> Paired MM pairs only.
# - This is a discovery analysis in a small GEO cohort, not a final clinical model.
# ============================================================================
write.csv(meta117847,"/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE117847_sample_annotation.csv",row.names = F)
saveRDS(
   expr_gene117847,
   file = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE117847_expression_gene_symbol.rds"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(pheatmap)
  library(scales)
  library(grid)
})

# -----------------------------------------------------------------------------
# 0. User settings
# -----------------------------------------------------------------------------
# Option A (recommended): run this script AFTER your existing GSE117847 script,
# which should leave the following objects in the R environment:
#   expr117847 : numeric matrix, gene symbol x GSM/sample
#   sample_info: data.frame with gsm_id, group, patient_id
#
# Option B: save the objects and set the two file paths below.
expr_rds_path <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE117847_expression_gene_symbol.rds"
meta_csv_path <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE117847_sample_annotation.csv"

outdir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260702)

# Expected metadata format:
# gsm_id,group,patient_id
# GSM123...,NP-SMM,NP_01
# GSM123...,P-SMM,P_01
# GSM123...,Paired MM,P_01
#
# patient_id may be arbitrary for NP-SMM, but P-SMM and Paired MM must share
# exactly the same patient_id for genuine longitudinal pairs.

# -----------------------------------------------------------------------------
# 1. Program definitions fixed BEFORE viewing GSE117847 results
# -----------------------------------------------------------------------------
State2_genes <- c(
  "MCC", "CADPS2", "SFMBT2", "COBLL1", "SAMD12", "KIAA1217",
  "COL6A3", "ADGRB3", "COL4A5", "EPHA6", "NEU3", "NBEA",
  "RASSF6", "STARD9", "CNTN5", "NCAM1", "SYT1", "FMN1",
  "ATP10B", "PCDH9", "RELN", "MEF2C", "KHDRBS2", "DDX31",
  "RAPGEF5", "TMTC2", "HOMER1", "NDNF", "AJAP1", "NEB",
  "ESRRG", "MAP2", "TTC28", "PRKG1", "DCC", "SP4"
)

State3_genes <- c(
   "MAD2L1", "CLSPN", "ASPM", "SGO1", "CDCA5", "WDR76", "TEDC2", "ZNF367",
   "DTL", "FAM111B", "NUP155", "FANCD2", "CKAP2L", "ECT2", "WDR62", "RAD54L",
   "CGAS", "LIN9", "INTS7", "BDH1", "STIL", "E2F7", "FAM133A", "CHAF1B",
   "TSPAN5", "SLC35F2", "TRIM59", "GINS3", "ORC1", "CGREF1", "DMRT2",
   "FOXRED2", "DUSP14", "MSRA"
)

State2_genes <- unique(toupper(State2_genes))
State3_genes <- unique(toupper(State3_genes))

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------
assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install missing package(s): ", paste(missing, collapse = ", "))
  }
}

# Collapse duplicated gene symbols (e.g., multiple microarray probes) by mean.
# Input must be a numeric matrix with genes in rows and samples in columns.
collapse_duplicate_genes <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  rn <- toupper(trimws(rownames(mat)))
  keep <- !is.na(rn) & rn != "" & rn != "---"
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  collapsed <- rowsum(mat, group = rn, reorder = FALSE, na.rm = TRUE)
  n_per_gene <- as.numeric(table(rn)[rownames(collapsed)])
  collapsed <- sweep(collapsed, 1, n_per_gene, "/")
  collapsed
}

# Program score = mean within-sample percentile rank across signature genes.
# This rank-based score is robust to different microarray/RNA-seq dynamic ranges.
rank_program_score <- function(expr, signature) {
  signature <- intersect(toupper(signature), rownames(expr))
  if (length(signature) < 5) {
    stop("Fewer than 5 signature genes are detected; check gene-symbol mapping.")
  }
  rank_mat <- apply(expr, 2, function(z) {
    rank(z, ties.method = "average", na.last = "keep") / (sum(!is.na(z)) + 1)
  })
  if (is.null(dim(rank_mat))) rank_mat <- matrix(rank_mat, ncol = 1)
  rownames(rank_mat) <- rownames(expr)
  colnames(rank_mat) <- colnames(expr)
  score <- colMeans(rank_mat[signature, , drop = FALSE], na.rm = TRUE)
  attr(score, "detected_genes") <- signature
  score
}

cliffs_delta <- function(x, y) {
  # Positive values indicate x > y.
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  mean(sign(outer(x, y, "-")))
}

bootstrap_cliffs_delta <- function(x, y, B = 5000, seed = 20260702) {
  set.seed(seed)
  boot <- replicate(B, {
    cliffs_delta(sample(x, length(x), replace = TRUE), sample(y, length(y), replace = TRUE))
  })
  unname(stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = 3))))
}

make_label <- function(raw, n) paste0(raw, "\n(n = ", n, ")")

pub_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.55, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(colour = "grey30", size = base_size - 1),
      plot.margin = margin(8, 10, 8, 8),
      legend.title = element_text(face = "bold"),
      legend.background = element_blank()
    )
}

save_plot_dual <- function(plot, filename, width, height) {
   # 改用原生pdf设备，移除cairo_pdf
   ggsave(file.path(outdir, paste0(filename, ".pdf")), plot, width = width, height = height,
          units = "in", device = "pdf", bg = "white")
   ggsave(file.path(outdir, paste0(filename, ".png")), plot, width = width, height = height,
          units = "in", dpi = 600, bg = "white")
}

# -----------------------------------------------------------------------------
# 3. Load expression matrix and metadata
# -----------------------------------------------------------------------------
if (!exists("expr117847", inherits = TRUE)) {
  if (!file.exists(expr_rds_path)) {
    stop(
      "Cannot find expr117847 in the R environment or ", expr_rds_path, ".\n",
      "Run your previous GSE117847 extraction script first, then assign:\n",
      "  expr117847 <- <gene_symbol_by_sample_matrix>\n",
      "  sample_info <- <metadata_with_gsm_id_group_patient_id>"
    )
  }
  expr117847 <- readRDS(expr_rds_path)
}

if (!exists("sample_info", inherits = TRUE)) {
  if (!file.exists(meta_csv_path)) {
    stop(
      "Cannot find sample_info in the R environment or ", meta_csv_path, ".\n",
      "Metadata must contain: gsm_id, group, patient_id"
    )
  }
  sample_info <- read.csv(meta_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
}

expr <- as.matrix(expr117847)
if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("expr117847 needs row names (gene symbols) and column names (GSM/sample identifiers).")
}

required_meta <- c("gsm_id", "group", "patient_id")
if (!all(required_meta %in% colnames(sample_info))) {
  stop("sample_info is missing required column(s): ",
       paste(setdiff(required_meta, colnames(sample_info)), collapse = ", "))
}

sample_info <- sample_info %>%
  mutate(
    gsm_id = as.character(gsm_id),
    group = as.character(group),
    patient_id = as.character(patient_id)
  )

valid_groups <- c("NP-SMM", "P-SMM", "Paired MM")
if (!all(sample_info$group %in% valid_groups)) {
  bad <- unique(sample_info$group[!sample_info$group %in% valid_groups])
  stop("group must be exactly one of: ", paste(valid_groups, collapse = ", "),
       ". Invalid value(s): ", paste(bad, collapse = ", "))
}

# Orient matrix if needed, then harmonize sample order.
if (!all(sample_info$gsm_id %in% colnames(expr)) && all(sample_info$gsm_id %in% rownames(expr))) {
  expr <- t(expr)
}
if (!all(sample_info$gsm_id %in% colnames(expr))) {
  stop("Some gsm_id values are not found among expression-matrix columns.\nMissing: ",
       paste(setdiff(sample_info$gsm_id, colnames(expr)), collapse = ", "))
}

expr <- expr[, sample_info$gsm_id, drop = FALSE]
expr <- collapse_duplicate_genes(expr)

# Metadata QC
if (anyDuplicated(sample_info$gsm_id)) stop("gsm_id must be unique.")
if (anyNA(expr)) warning("Expression matrix contains NA values; ranking ignores NAs sample-wise.")

# Pair QC. NP-SMM IDs are irrelevant. Each P-SMM and Paired MM must have one sample.
# p_counts <- sample_info %>% filter(group %in% c("P-SMM", "Paired MM")) %>% count(patient_id, group) %>%
#   tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)

p_counts <- sample_info %>%
   dplyr::filter(.data$group %in% c("P-SMM", "Paired MM")) %>%
   dplyr::count(.data$patient_id, .data$group, name = "n") %>%
   tidyr::pivot_wider(
      names_from = .data$group,
      values_from = .data$n,
      values_fill = list(n = 0)
   )

if (!all(c("P-SMM", "Paired MM") %in% colnames(p_counts))) {
  stop("Both P-SMM and Paired MM groups are required.")
}
if (any(p_counts$`P-SMM` != 1 | p_counts$`Paired MM` != 1)) {
  print(p_counts)
  stop("Every progressor patient must have exactly one baseline P-SMM and one Paired MM sample.")
}

# Save normalized inputs for reproducibility.
saveRDS(expr, file.path(outdir, "GSE117847_gene_symbol_expression_collapsed.rds"))
write.csv(sample_info, file.path(outdir, "GSE117847_metadata_used.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Score State 2 and State 3
# -----------------------------------------------------------------------------
state2_score <- rank_program_score(expr, State2_genes)
state3_score <- rank_program_score(expr, State3_genes)

coverage <- bind_rows(
  tibble(
    program = "State 2",
    genes_requested = length(State2_genes),
    genes_detected = length(attr(state2_score, "detected_genes")),
    genes_missing = paste(setdiff(State2_genes, attr(state2_score, "detected_genes")), collapse = "; ")
  ),
  tibble(
    program = "State 3",
    genes_requested = length(State3_genes),
    genes_detected = length(attr(state3_score, "detected_genes")),
    genes_missing = paste(setdiff(State3_genes, attr(state3_score, "detected_genes")), collapse = "; ")
  )
)
write.csv(coverage, file.path(outdir, "Table_gene_coverage.csv"), row.names = FALSE)
print(coverage)

score_df <- sample_info %>%
  mutate(
    State2_rank = unname(state2_score[gsm_id]),
    State3_rank = unname(state3_score[gsm_id]),
    group = factor(group, levels = valid_groups)
  )
write.csv(score_df, file.path(outdir, "Table_patient_level_program_scores.csv"), row.names = FALSE)

# Exact pair table for longitudinal inference.
pair_df <- score_df %>%
  filter(group == "P-SMM") %>%
  select(patient_id, gsm_id, State2_rank, State3_rank) %>%
  rename(P_SMM_gsm = gsm_id, State2_P = State2_rank, State3_P = State3_rank) %>%
  inner_join(
    score_df %>% filter(group == "Paired MM") %>%
      select(patient_id, gsm_id, State2_rank, State3_rank) %>%
      rename(MM_gsm = gsm_id, State2_MM = State2_rank, State3_MM = State3_rank),
    by = "patient_id"
  ) %>%
  mutate(
    Delta_State2 = State2_MM - State2_P,
    Delta_State3 = State3_MM - State3_P
  )
write.csv(pair_df, file.path(outdir, "Table_matched_P_SMM_to_MM_deltas.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Pre-specified statistical tests
# -----------------------------------------------------------------------------
baseline <- score_df %>% filter(group %in% c("NP-SMM", "P-SMM"))

run_baseline_test <- function(score_name) {
  x_np <- baseline %>% filter(group == "NP-SMM") %>% pull(.data[[score_name]])
  x_p  <- baseline %>% filter(group == "P-SMM") %>% pull(.data[[score_name]])
  wt <- wilcox.test(x_p, x_np, alternative = "two.sided", exact = FALSE)
  cd <- cliffs_delta(x_p, x_np)
  cd_ci <- bootstrap_cliffs_delta(x_p, x_np)
  tibble(
    analysis = "Baseline P-SMM versus NP-SMM",
    program = sub("_rank", "", score_name),
    n_NP_SMM = length(x_np),
    n_P_SMM = length(x_p),
    statistic_W = unname(wt$statistic),
    p_value = wt$p.value,
    effect = "Cliff's delta (P-SMM minus NP-SMM)",
    effect_estimate = cd,
    ci_low = cd_ci[1],
    ci_high = cd_ci[2]
  )
}

run_paired_test <- function(program) {
  x_p <- pair_df[[paste0(program, "_P")]]
  x_m <- pair_df[[paste0(program, "_MM")]]
  wt <- wilcox.test(x_m, x_p, paired = TRUE, alternative = "two.sided", exact = FALSE, conf.int = TRUE)
  tibble(
    analysis = "Paired P-SMM to MM",
    program = program,
    n_pairs = length(x_p),
    statistic_V = unname(wt$statistic),
    p_value = wt$p.value,
    effect = "Hodges-Lehmann paired shift (MM minus P-SMM)",
    effect_estimate = unname(wt$estimate),
    ci_low = unname(wt$conf.int[1]),
    ci_high = unname(wt$conf.int[2])
  )
}

stat_tbl <- bind_rows(
  run_baseline_test("State2_rank"),
  run_baseline_test("State3_rank"),
  run_paired_test("State2"),
  run_paired_test("State3")
) %>%
  mutate(p_BH_within_analysis = p.adjust(p_value, method = "BH"))
write.csv(stat_tbl, file.path(outdir, "Table_pre_specified_statistics.csv"), row.names = FALSE)
print(stat_tbl)

# ROC analyses use BASELINE SMM samples only. These are exploratory because N is small.
roc_state2 <- pROC::roc(
  response = baseline$group,
  predictor = baseline$State2_rank,
  levels = c("NP-SMM", "P-SMM"),
  direction = "<",
  ci = TRUE, quiet = TRUE
)
roc_state3 <- pROC::roc(
  response = baseline$group,
  predictor = baseline$State3_rank,
  levels = c("NP-SMM", "P-SMM"),
  direction = "<",
  ci = TRUE, quiet = TRUE
)

roc_tbl <- bind_rows(
  tibble(program = "State 2", AUC = as.numeric(auc(roc_state2)),
         CI_low = as.numeric(ci.auc(roc_state2)[1]), CI_high = as.numeric(ci.auc(roc_state2)[3])),
  tibble(program = "State 3", AUC = as.numeric(auc(roc_state3)),
         CI_low = as.numeric(ci.auc(roc_state3)[1]), CI_high = as.numeric(ci.auc(roc_state3)[3]))
)
write.csv(roc_tbl, file.path(outdir, "Table_exploratory_ROC.csv"), row.names = FALSE)
print(roc_tbl)

# -----------------------------------------------------------------------------
# 6. Publication-quality plots
# -----------------------------------------------------------------------------
group_cols <- c("NP-SMM" = "#4C78A8", "P-SMM" = "#E45756", "Paired MM" = "#54A24B")
state2_base_stat <- stat_tbl %>% filter(analysis == "Baseline P-SMM versus NP-SMM", program == "State2")
state3_base_stat <- stat_tbl %>% filter(analysis == "Baseline P-SMM versus NP-SMM", program == "State3")
state2_pair_stat <- stat_tbl %>% filter(analysis == "Paired P-SMM to MM", program == "State2")
state3_pair_stat <- stat_tbl %>% filter(analysis == "Paired P-SMM to MM", program == "State3")

plot_baseline <- function(score_var, title, stat_row) {
  plot_dat <- baseline %>%
    mutate(group_plot = factor(group, levels = c("NP-SMM", "P-SMM")))
  ns <- plot_dat %>% dplyr::count(group_plot)
  x_labels <- setNames(make_label(as.character(ns$group_plot), ns$n), as.character(ns$group_plot))
  y_max <- max(plot_dat[[score_var]], na.rm = TRUE)
  y_min <- min(plot_dat[[score_var]], na.rm = TRUE)
  y_pad <- max((y_max - y_min) * 0.19, 0.025)
  ann <- paste0(
    fmt_p(stat_row$p_value), "\n",
    "Cliff's delta = ", sprintf("%.2f", stat_row$effect_estimate)
  )

  ggplot(plot_dat, aes(x = group_plot, y = .data[[score_var]], fill = group_plot, colour = group_plot)) +
    geom_violin(trim = FALSE, alpha = 0.22, linewidth = 0.75, width = 0.88) +
    geom_boxplot(width = 0.27, alpha = 0.42, outlier.shape = NA, linewidth = 0.65) +
    geom_jitter(width = 0.075, size = 2.55, alpha = 0.92, show.legend = FALSE) +
    annotate("text", x = 1.5, y = y_max + y_pad, label = ann, fontface = "bold", size = 4.0) +
    scale_fill_manual(values = group_cols[c("NP-SMM", "P-SMM")]) +
    scale_colour_manual(values = group_cols[c("NP-SMM", "P-SMM")]) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(y_min - y_pad * 0.15, y_max + y_pad * 1.55), clip = "off") +
    labs(title = title, x = NULL, y = "Rank-based program score") +
    pub_theme()
}

plot_paired <- function(program, title, stat_row) {
  pcol <- paste0(program, "_P")
  mcol <- paste0(program, "_MM")
  long <- pair_df %>%
    select(patient_id, all_of(c(pcol, mcol))) %>%
    pivot_longer(-patient_id, names_to = "time", values_to = "score") %>%
    mutate(
      time = recode(time, !!pcol := "P-SMM", !!mcol := "Paired MM"),
      time = factor(time, levels = c("P-SMM", "Paired MM"))
    )
  ns <- long %>% dplyr::count(.data$time)
  x_labels <- setNames(make_label(as.character(ns$time), ns$n), as.character(ns$time))
  y_max <- max(long$score, na.rm = TRUE)
  y_min <- min(long$score, na.rm = TRUE)
  y_pad <- max((y_max - y_min) * 0.18, 0.025)
  delta_name <- paste0("Delta_", program)
  ann <- paste0(
    "Paired ", fmt_p(stat_row$p_value), "\n",
    "Median delta = ", sprintf("%.3f", median(pair_df[[delta_name]], na.rm = TRUE))
  )

  ggplot(long, aes(x = time, y = score, group = patient_id)) +
    geom_line(colour = "grey65", linewidth = 0.72, alpha = 0.88) +
    geom_point(aes(colour = time), size = 3.0, alpha = 0.96) +
    stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black",
                 linetype = "dashed", linewidth = 0.88) +
    annotate("text", x = 1.5, y = y_max + y_pad, label = ann, fontface = "bold", size = 4.0) +
    scale_colour_manual(values = group_cols[c("P-SMM", "Paired MM")]) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(y_min - y_pad * 0.12, y_max + y_pad * 1.55), clip = "off") +
    labs(title = title, x = NULL, y = "Rank-based program score") +
    pub_theme() +
    theme(legend.position = "none")
}

p_A <- plot_baseline("State2_rank", "A  State 2 is elevated in progressor SMM at diagnosis", state2_base_stat)
p_B <- plot_baseline("State3_rank", "B  State 3 does not discriminate SMM progression at diagnosis", state3_base_stat)
p_C <- plot_paired("State2", "C  State 2 is retained during SMM-to-MM transition", state2_pair_stat)
p_D <- plot_paired("State3", "D  State 3 intensifies at overt-MM conversion", state3_pair_stat)

main_figure <- (p_A | p_B) / (p_C | p_D) +
  plot_annotation(
    title = "A two-stage transcriptional model of SMM progression",
    subtitle = "State 2 marks a high-risk precursor phenotype, whereas State 3 marks late proliferative escalation during clinical conversion"
  )
save_plot_dual(main_figure, "Figure_State2_State3_two_stage_main", width = 13, height = 8)

# Exploratory ROC figure: both AUCs are displayed so State 2's diagnostic value
# is interpreted relative to State 3 rather than shown in isolation.
plot_roc <- function(roc_obj, table_row, title, color) {
  roc_dat <- tibble(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  label <- paste0("AUC = ", sprintf("%.3f", table_row$AUC),
                  "\n95% CI ", sprintf("%.3f", table_row$CI_low), "–", sprintf("%.3f", table_row$CI_high))
  ggplot(roc_dat, aes(x = fpr, y = tpr)) +
    geom_abline(linetype = "dashed", colour = "grey70", linewidth = 0.7) +
    geom_path(colour = color, linewidth = 1.15) +
    annotate("text", x = 0.62, y = 0.15, label = label, hjust = 0, fontface = "bold", size = 4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title, x = "False-positive rate", y = "True-positive rate") +
    pub_theme()
}

p_roc2 <- plot_roc(roc_state2, roc_tbl %>% filter(program == "State 2"),
                   "State 2: baseline P-SMM versus NP-SMM", "#7A5195")
p_roc3 <- plot_roc(roc_state3, roc_tbl %>% filter(program == "State 3"),
                   "State 3: baseline P-SMM versus NP-SMM", "#7A5195")
roc_figure <- p_roc2 | p_roc3
save_plot_dual(roc_figure, "Figure_State2_State3_baseline_ROC", width = 11, height = 5)

# Joint trajectory plot: gives an intuitive state-space view of the two-stage model.
score_z <- score_df %>%
  mutate(
    State2_z = as.numeric(scale(State2_rank)),
    State3_z = as.numeric(scale(State3_rank))
  )
traj_segments <- score_z %>% filter(group == "P-SMM") %>%
  select(patient_id, State2_z, State3_z) %>%
  rename(x = State2_z, y = State3_z) %>%
  inner_join(
    score_z %>% filter(group == "Paired MM") %>%
      select(patient_id, State2_z, State3_z) %>% rename(xend = State2_z, yend = State3_z),
    by = "patient_id"
  )

p_space <- ggplot(score_z, aes(x = State2_z, y = State3_z)) +
  geom_segment(data = traj_segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey60", linewidth = 0.65,
               arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_point(aes(colour = group, shape = group), size = 3.1, alpha = 0.94) +
  scale_colour_manual(values = group_cols) +
  scale_shape_manual(values = c("NP-SMM" = 16, "P-SMM" = 16, "Paired MM" = 17)) +
  labs(
    title = "State-space representation of disease evolution",
    subtitle = "Arrows link matched P-SMM and MM samples; axes are z-scored within program",
    x = "State 2 score (z)", y = "State 3 score (z)", colour = NULL, shape = NULL
  ) +
  pub_theme() +
  theme(legend.position = "right")
save_plot_dual(p_space, "Figure_State2_State3_state_space", width = 7, height = 6)

# Heatmaps: program scores and gene-level expression patterns.
all_order <- score_df %>%
  arrange(group, patient_id) %>%
  pull(gsm_id)
ann_col <- score_df %>%
  select(gsm_id, group) %>%
  arrange(match(gsm_id, all_order)) %>%
  as.data.frame()
rownames(ann_col) <- ann_col$gsm_id
ann_col$gsm_id <- NULL

program_mat <- rbind(
  `State 2` = score_df$State2_rank[match(all_order, score_df$gsm_id)],
  `State 3` = score_df$State3_rank[match(all_order, score_df$gsm_id)]
)
program_z <- t(scale(t(program_mat)))
colnames(program_z) <- all_order

pdf(file.path(outdir, "Figure_State2_State3_program_score_heatmap.pdf"), width = 12, height = 3.2)
pheatmap(
  program_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = list(group = group_cols),
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 7,
  main = "Relative State 2 and State 3 activity across the GSE117847 trajectory"
)
dev.off()

# Marker heatmaps are supplementary evidence only; scoring remains defined by all detected genes.
# make_marker_heatmap <- function(signature, samples, filename, title, width, height) {
#   genes <- intersect(signature, rownames(expr))
#   mat <- expr[genes, samples, drop = FALSE]
#   mat_z <- t(scale(t(mat)))
#   # Drop invariant genes, which otherwise create NA rows in heatmaps.
#   mat_z <- mat_z[apply(mat_z, 1, function(z) all(is.finite(z)) && sd(z) > 0), , drop = FALSE]
#   if (nrow(mat_z) < 2) return(invisible(NULL))
#   ann <- score_df %>% filter(gsm_id %in% samples) %>%
#     arrange(match(gsm_id, samples)) %>% select(gsm_id, group) %>% as.data.frame()
#   rownames(ann) <- ann$gsm_id
#   ann$gsm_id <- NULL
#   pdf(file.path(outdir, filename), width = width, height = height)
#   pheatmap(
#     mat_z,
#     cluster_rows = TRUE,
#     cluster_cols = FALSE,
#     annotation_col = ann,
#     annotation_colors = list(group = group_cols),
#     color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
#     border_color = NA,
#     fontsize_row = 6.5,
#     fontsize_col = 7,
#     main = title
#   )
#   dev.off()
# }
# 
# baseline_order <- baseline %>% arrange(group, patient_id) %>% pull(gsm_id)
# pair_order <- pair_df %>% arrange(patient_id) %>%
#   select(P_SMM_gsm, MM_gsm) %>% unlist(use.names = FALSE)
# make_marker_heatmap(State2_genes, baseline_order,
#                     "FigureS_State2_baseline_marker_heatmap.pdf",
#                     "State 2 marker expression at SMM diagnosis", 10, 8)
# make_marker_heatmap(State3_genes, pair_order,
#                     "FigureS_State3_paired_marker_heatmap.pdf",
#                     "State 3 marker expression across paired SMM-to-MM transition", 12, 10)
# 
# # Record session information for reproducibility.
# writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))


make_marker_heatmap <- function(signature, samples, filename, title, width, height) {
   
   samples <- intersect(samples, colnames(expr))
   genes <- intersect(toupper(signature), rownames(expr))
   
   message("Heatmap: ", filename)
   message("Detected genes: ", length(genes))
   message("Detected samples: ", length(samples))
   
   if (length(genes) < 2 || length(samples) < 2) {
      warning("Too few genes or samples for heatmap: ", filename)
      return(invisible(NULL))
   }
   
   mat <- expr[genes, samples, drop = FALSE]
   
   keep_gene <- apply(mat, 1, function(x) {
      all(is.finite(x)) && sd(x, na.rm = TRUE) > 0
   })
   mat <- mat[keep_gene, , drop = FALSE]
   
   message("Genes after variance filter: ", nrow(mat))
   
   if (nrow(mat) < 2) {
      warning("Too few variable genes after filtering: ", filename)
      return(invisible(NULL))
   }
   
   mat_z <- t(scale(t(mat)))
   mat_z[mat_z > 2] <- 2
   mat_z[mat_z < -2] <- -2
   
   ann <- score_df %>%
      dplyr::filter(.data$gsm_id %in% samples) %>%
      dplyr::arrange(match(.data$gsm_id, samples)) %>%
      dplyr::select(gsm_id, group) %>%
      as.data.frame()
   
   rownames(ann) <- ann$gsm_id
   ann$gsm_id <- NULL
   
   ann <- ann[samples, , drop = FALSE]
   
   out_file <- file.path(outdir, filename)
   
   p <- pheatmap::pheatmap(
      mat_z,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = ann,
      annotation_colors = list(group = group_cols),
      color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
      border_color = NA,
      fontsize_row = 6.5,
      fontsize_col = 7,
      main = title,
      silent = TRUE
   )
   
   grDevices::pdf(out_file, width = width, height = height, useDingbats = FALSE)
   grid::grid.newpage()
   grid::grid.draw(p$gtable)
   grDevices::dev.off()
   
   message("Saved: ", out_file)
}

baseline_order <- baseline %>%
   dplyr::arrange(group, patient_id) %>%
   dplyr::pull(gsm_id)

pair_order <- pair_df %>%
   dplyr::arrange(patient_id) %>%
   dplyr::select(P_SMM_gsm, MM_gsm) %>%
   unlist(use.names = FALSE)

make_marker_heatmap(
   State2_genes,
   baseline_order,
   "FigureS_State2_baseline_marker_heatmap.pdf",
   "State 2 marker expression at SMM diagnosis",
   10,
   8
)

make_marker_heatmap(
   State3_genes,
   pair_order,
   "FigureS_State3_paired_marker_heatmap.pdf",
   "State 3 marker expression across paired SMM-to-MM transition",
   12,
   10
)
