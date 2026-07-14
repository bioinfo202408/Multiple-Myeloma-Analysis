# ============================================================
# Marker-based relative deconvolution of MM plasma-cell states
# and overall-survival analysis in MMRF-CoMMpass
#
# Purpose:
#   1. Estimate relative State1/State2/State3 abundance from bulk RNA-seq
#      using state-specific marker genes.
#   2. Aggregate repeated samples from the same patient.
#   3. Divide patients into State3-high and State3-low groups by the
#      median estimated State3 abundance.
#   4. Perform Kaplan-Meier and Cox regression analyses for overall survival.
#
# Important methodological note:
#   Because only marker-gene lists are supplied, this script performs a
#   marker-based relative abundance estimation rather than reference-matrix
#   deconvolution such as BayesPrism or CIBERSORTx. The resulting values should
#   be described as "estimated relative state abundance" or
#   "State3-like transcriptional contribution", not an absolute cell fraction.
# ============================================================

# ============================================================
# 0. Required packages
# ============================================================
required_packages <- c("survival", "survminer", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(survival)
library(survminer)
library(ggplot2)

set.seed(20260708)

# ============================================================
# 1. Input and output paths
# ============================================================
input_file <- paste0(
  "/home/yjliu/mmProj/clinical/MMRF/UCSC/",
  "MMRF-COMMPASS.star_tpm_行859样本_列基因_SymbolID_包含OS等其他临床信息.csv"
)

output_dir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/State_marker_deconvolution_OS"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. State-specific marker genes
# ============================================================
MM1_specific <- c(
  "AIM2", "ALKBH2", "BYSL", "CADM1", "DKK1", "ETV2", "FHL1", "FZD3",
  "GNB1L", "HEY2", "KBTBD3", "LAMP5", "MTO1", "MYEOV", "NOMO2",
  "PDCD2L", "PFDN4", "PLEKHF2", "PRXL2A", "RSU1", "SELENOP",
  "SYNC", "TFAP4", "TNFRSF13B", "ZNF296", "ZNF749"
)

MM2_specific <- c(
  "ADGRB3", "AJAP1", "ATP10B", "CADPS2", "CNTN5", "COBLL1", "COL4A5",
  "COL6A3", "DCC", "DDX31", "EPHA6", "ESRRG", "FMN1", "HOMER1",
  "KHDRBS2", "KIAA1217", "MAP2", "MCC", "MEF2C", "NBEA", "NCAM1",
  "NDNF", "NEB", "NEU3", "PCDH9", "PRKG1", "RAPGEF5", "RASSF6",
  "RELN", "SAMD12", "SFMBT2", "SP4", "STARD9", "SYT1", "TMTC2", "TTC28"
)

MM3_specific <- c(
  "ASPM", "CDCA5", "CGAS", "CGREF1", "CHAF1B", "CKAP2L", "CLSPN",
  "DMRT2", "DTL", "DUSP14", "E2F7", "ECT2", "FAM111B", "FANCD2",
  "FOXRED2", "GINS3", "INTS7", "LIN9", "MAD2L1", "MSRA", "ORC1",
  "RAD54L", "SGO1", "SLC35F2", "STIL", "TEDC2", "TRIM59", "TSPAN5",
  "WDR62", "WDR76", "ZNF367"
)

state_markers <- list(
  State1 = unique(MM1_specific),
  State2 = unique(MM2_specific),
  State3 = unique(MM3_specific)
)

# Confirm that the marker sets are mutually exclusive.
marker_overlap <- combn(names(state_markers), 2, simplify = FALSE)
overlap_table <- do.call(
  rbind,
  lapply(marker_overlap, function(x) {
    overlap_genes <- intersect(state_markers[[x[1]]], state_markers[[x[2]]])
    data.frame(
      state_pair = paste(x, collapse = "_vs_"),
      n_overlap = length(overlap_genes),
      overlap_genes = paste(overlap_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  overlap_table,
  file.path(output_dir, "state_marker_overlap_check.csv"),
  row.names = FALSE
)

if (any(overlap_table$n_overlap > 0)) {
  stop("The state marker sets are not mutually exclusive. Please check the marker lists.")
}

# ============================================================
# 3. Read CoMMpass data
# ============================================================
data <- read.table(
  input_file,
  header = TRUE,
  sep = ",",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

required_survival_columns <- c("OS.time", "OS")
missing_survival_columns <- setdiff(required_survival_columns, colnames(data))

if (length(missing_survival_columns) > 0) {
  stop(
    "The following required survival columns were not found: ",
    paste(missing_survival_columns, collapse = ", ")
  )
}

# ============================================================
# 4. Detect sample identifiers
# ============================================================
detect_sample_ids <- function(df) {
  candidate_names <- c(
    "sample_id", "SampleID", "sample", "Sample", "sample_name",
    "...1", "X", ""
  )

  matched_names <- candidate_names[candidate_names %in% colnames(df)]

  if (length(matched_names) > 0) {
    idx <- match(matched_names[1], colnames(df))
    candidate_values <- as.character(df[[idx]])
    if (length(unique(candidate_values)) == nrow(df)) {
      return(candidate_values)
    }
  }

  # Use the first column if it is not a known clinical variable and is unique.
  known_clinical <- c(
    "OS.time", "OS", "Gender", "Age", "iss_stage", "treatment_type"
  )

  first_name <- colnames(df)[1]
  first_values <- as.character(df[[1]])

  if (
    !(first_name %in% known_clinical) &&
    length(unique(first_values)) == nrow(df)
  ) {
    return(first_values)
  }

  # Otherwise, use row names if informative.
  rn <- rownames(df)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) {
    return(rn)
  }

  warning("No explicit sample-ID column was detected; artificial sample IDs were created.")
  paste0("sample_", seq_len(nrow(df)))
}

data$sample_id_internal <- detect_sample_ids(data)

# Extract patient IDs from identifiers such as:
# MMRF_1024_1_BM_CD138pos -> MMRF_1024
extract_patient_id <- function(sample_id) {
  sample_id <- as.character(sample_id)
  is_mmrf <- grepl("^MMRF_[^_]+", sample_id)

  patient_id <- sample_id
  patient_id[is_mmrf] <- sub(
    "^(MMRF_[^_]+).*$",
    "\\1",
    sample_id[is_mmrf],
    perl = TRUE
  )

  # Generic fallback for non-MMRF identifiers.
  patient_id[!is_mmrf] <- sub(
    "_[0-9]+_BM.*$",
    "",
    sample_id[!is_mmrf],
    perl = TRUE
  )

  patient_id
}

data$patient_id_internal <- extract_patient_id(data$sample_id_internal)

# ============================================================
# 5. Check marker availability
# ============================================================
marker_check <- do.call(
  rbind,
  lapply(names(state_markers), function(state_name) {
    genes <- state_markers[[state_name]]
    present <- genes[genes %in% colnames(data)]
    missing <- setdiff(genes, present)

    data.frame(
      state = state_name,
      total_markers = length(genes),
      present_markers = length(present),
      missing_markers = length(missing),
      present_gene_names = paste(present, collapse = ";"),
      missing_gene_names = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_check,
  file.path(output_dir, "state_marker_availability.csv"),
  row.names = FALSE
)

print(marker_check[, c("state", "total_markers", "present_markers", "missing_markers")])

if (any(marker_check$present_markers < 5)) {
  stop(
    "Fewer than five marker genes were detected for at least one state. ",
    "The marker-based abundance estimation would be unreliable."
  )
}

present_markers_by_state <- lapply(
  state_markers,
  function(x) intersect(x, colnames(data))
)
all_present_markers <- unique(unlist(present_markers_by_state, use.names = FALSE))

# ============================================================
# 6. Prepare marker-expression matrix
# ============================================================
expression_matrix <- data[, all_present_markers, drop = FALSE]
expression_matrix[] <- lapply(
  expression_matrix,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)
expression_matrix <- as.matrix(expression_matrix)
rownames(expression_matrix) <- data$sample_id_internal

# Remove genes that are entirely missing.
all_na_genes <- colnames(expression_matrix)[
  apply(expression_matrix, 2, function(x) all(is.na(x)))
]

if (length(all_na_genes) > 0) {
  expression_matrix <- expression_matrix[, !colnames(expression_matrix) %in% all_na_genes, drop = FALSE]
}

# Median-impute occasional missing expression values gene by gene.
for (j in seq_len(ncol(expression_matrix))) {
  missing_idx <- is.na(expression_matrix[, j])
  if (any(missing_idx)) {
    gene_median <- median(expression_matrix[, j], na.rm = TRUE)
    expression_matrix[missing_idx, j] <- gene_median
  }
}

# Automatically determine whether log2 transformation is needed.
# Values such as TPM/counts usually have a high upper quantile, whereas
# log2(TPM + 1) values are usually below approximately 20-30.
expression_q99 <- as.numeric(
  quantile(expression_matrix, probs = 0.99, na.rm = TRUE)
)

if (is.finite(expression_q99) && expression_q99 > 50) {
  message("The expression matrix appears to be on a linear scale; applying log2(x + 1).")
  expression_matrix <- log2(pmax(expression_matrix, 0) + 1)
  expression_scale_used <- "log2(x + 1) applied by script"
} else {
  message("The expression matrix appears to be already log-transformed; no additional log transformation was applied.")
  expression_scale_used <- "input scale retained"
}

# Remove zero-variance genes because they provide no discriminatory information.
gene_sd <- apply(expression_matrix, 2, sd, na.rm = TRUE)
valid_genes <- names(gene_sd)[is.finite(gene_sd) & gene_sd > 0]
removed_constant_genes <- setdiff(colnames(expression_matrix), valid_genes)
expression_matrix <- expression_matrix[, valid_genes, drop = FALSE]

writeLines(
  c(
    paste0("Expression scale: ", expression_scale_used),
    paste0("99th percentile before optional transformation: ", signif(expression_q99, 6)),
    paste0("Number of marker genes retained: ", ncol(expression_matrix)),
    paste0(
      "Zero-variance/all-missing genes removed: ",
      ifelse(
        length(c(all_na_genes, removed_constant_genes)) == 0,
        "None",
        paste(unique(c(all_na_genes, removed_constant_genes)), collapse = ";")
      )
    )
  ),
  file.path(output_dir, "deconvolution_preprocessing_log.txt")
)

# Update state marker sets after filtering.
present_markers_by_state <- lapply(
  present_markers_by_state,
  function(x) intersect(x, colnames(expression_matrix))
)

if (any(vapply(present_markers_by_state, length, integer(1)) < 5)) {
  stop("After removing non-informative genes, fewer than five markers remained for at least one state.")
}

# ============================================================
# 7. Marker-based relative deconvolution
# ============================================================
# Standardize each marker gene across samples so that genes with large absolute
# expression values do not dominate the state score.
z_expression <- scale(expression_matrix)
z_expression <- as.matrix(z_expression)

# Winsorize extreme z-scores to reduce sensitivity to isolated outliers.
z_expression[z_expression > 3] <- 3
z_expression[z_expression < -3] <- -3

# Mean standardized marker expression for each state.
state_score_matrix <- sapply(
  names(present_markers_by_state),
  function(state_name) {
    genes <- present_markers_by_state[[state_name]]
    rowMeans(z_expression[, genes, drop = FALSE], na.rm = TRUE)
  }
)

state_score_matrix <- as.matrix(state_score_matrix)
colnames(state_score_matrix) <- paste0(names(present_markers_by_state), "_score")
rownames(state_score_matrix) <- data$sample_id_internal

# Convert the three state scores into positive relative abundances summing to 1
# in each sample using a softmax transformation.
softmax_rows <- function(x) {
  t(
    apply(x, 1, function(v) {
      v <- v - max(v, na.rm = TRUE)
      ev <- exp(v)
      ev / sum(ev)
    })
  )
}

state_fraction_matrix <- softmax_rows(state_score_matrix)
colnames(state_fraction_matrix) <- c(
  "State1_fraction", "State2_fraction", "State3_fraction"
)
rownames(state_fraction_matrix) <- data$sample_id_internal

# Confirm that estimated fractions sum to 1.
fraction_sum <- rowSums(state_fraction_matrix)
if (any(abs(fraction_sum - 1) > 1e-8)) {
  stop("Estimated state fractions do not sum to 1. Please check the calculation.")
}

sample_results <- data.frame(
  sample_id = data$sample_id_internal,
  patient_id = data$patient_id_internal,
  state_score_matrix,
  state_fraction_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Add available clinical variables.
clinical_columns <- intersect(
  c("OS.time", "OS", "Gender", "Age", "iss_stage", "treatment_type"),
  colnames(data)
)

for (clinical_name in clinical_columns) {
  sample_results[[clinical_name]] <- data[[clinical_name]]
}

write.csv(
  sample_results,
  file.path(output_dir, "sample_level_state_abundance.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Aggregate repeated samples to the patient level
# ============================================================
# The input data may contain multiple CD138-positive samples for the same
# patient. Survival analysis must use one row per patient to avoid artificial
# sample-size inflation. State scores/fractions are averaged across samples;
# clinical and survival variables are taken from the first non-missing value.
first_non_missing <- function(x) {
  valid <- !is.na(x)
  if (is.character(x)) {
    valid <- valid & nzchar(x)
  }

  idx <- which(valid)
  if (length(idx) == 0) {
    return(NA)
  }
  x[idx[1]]
}

patient_split <- split(seq_len(nrow(sample_results)), sample_results$patient_id)

patient_results <- do.call(
  rbind,
  lapply(names(patient_split), function(pid) {
    idx <- patient_split[[pid]]
    x <- sample_results[idx, , drop = FALSE]

    output <- data.frame(
      patient_id = pid,
      sample_ids = paste(unique(x$sample_id), collapse = ";"),
      n_samples = nrow(x),
      State1_score = mean(x$State1_score, na.rm = TRUE),
      State2_score = mean(x$State2_score, na.rm = TRUE),
      State3_score = mean(x$State3_score, na.rm = TRUE),
      State1_fraction = mean(x$State1_fraction, na.rm = TRUE),
      State2_fraction = mean(x$State2_fraction, na.rm = TRUE),
      State3_fraction = mean(x$State3_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    for (clinical_name in clinical_columns) {
      output[[clinical_name]] <- first_non_missing(x[[clinical_name]])
    }

    output
  })
)

rownames(patient_results) <- NULL

# Re-normalize fractions after aggregation to guard against numerical drift.
patient_fraction_sum <- rowSums(
  patient_results[, c("State1_fraction", "State2_fraction", "State3_fraction")]
)
patient_results$State1_fraction <- patient_results$State1_fraction / patient_fraction_sum
patient_results$State2_fraction <- patient_results$State2_fraction / patient_fraction_sum
patient_results$State3_fraction <- patient_results$State3_fraction / patient_fraction_sum

write.csv(
  patient_results,
  file.path(output_dir, "patient_level_state_abundance_before_survival_filter.csv"),
  row.names = FALSE
)

# ============================================================
# 9. Prepare patient-level survival data
# ============================================================
patient_results$OS.time <- suppressWarnings(
  as.numeric(as.character(patient_results$OS.time))
)
patient_results$OS <- suppressWarnings(
  as.numeric(as.character(patient_results$OS))
)
patient_results$State3_fraction <- suppressWarnings(
  as.numeric(as.character(patient_results$State3_fraction))
)

surv_data <- patient_results[
  complete.cases(
    patient_results[, c("OS.time", "OS", "State3_fraction")]
  ),
  ,
  drop = FALSE
]

surv_data <- surv_data[
  is.finite(surv_data$OS.time) &
    surv_data$OS.time > 0 &
    surv_data$OS %in% c(0, 1) &
    is.finite(surv_data$State3_fraction),
  ,
  drop = FALSE
]

if (nrow(surv_data) < 20) {
  stop("Fewer than 20 patients remained for survival analysis.")
}

cat("Number of unique patients used for survival analysis:", nrow(surv_data), "\n")
cat("Number of deaths:", sum(surv_data$OS == 1), "\n")

# ============================================================
# 10. Divide patients into State3-high and State3-low groups
# ============================================================
state3_cutoff <- median(surv_data$State3_fraction, na.rm = TRUE)

surv_data$State3_group <- ifelse(
  surv_data$State3_fraction >= state3_cutoff,
  "State3-high",
  "State3-low"
)

surv_data$State3_group <- factor(
  surv_data$State3_group,
  levels = c("State3-low", "State3-high")
)

# Continuous variable scaled to a 10-percentage-point increase for Cox models.
surv_data$State3_per10pct <- surv_data$State3_fraction / 0.10

write.csv(
  surv_data,
  file.path(output_dir, "patient_level_state_abundance_with_State3_group.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("State3 median cutoff = ", signif(state3_cutoff, 8)),
    paste0("State3-low n = ", sum(surv_data$State3_group == "State3-low")),
    paste0("State3-high n = ", sum(surv_data$State3_group == "State3-high"))
  ),
  file.path(output_dir, "State3_group_cutoff.txt")
)

# ============================================================
# 11. Kaplan-Meier overall-survival analysis
# ============================================================
fit_km <- survfit(
  Surv(OS.time, OS) ~ State3_group,
  data = surv_data
)

km_plot <- ggsurvplot(
  fit_km,
  data = surv_data,
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377EB8", "#E41A1C"),
  legend.title = "Estimated State 3 abundance",
  legend.labs = c("State3-low", "State3-high"),
  xlab = "Overall survival time (days)",
  ylab = "Overall survival probability",
  title = "Overall survival according to estimated State 3 abundance",
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 12),
  tables.theme = theme_cleantable(),
  censor = TRUE
)

pdf(
  file.path(output_dir, "State3_OS_KM_curve.pdf"),
  width = 8,
  height = 7.5,
  onefile = FALSE
)
print(km_plot)
dev.off()

png(
  file.path(output_dir, "State3_OS_KM_curve.png"),
  width = 2400,
  height = 2250,
  res = 300
)
print(km_plot)
dev.off()

# Median OS and group summary.
median_os <- survminer::surv_median(fit_km)
median_os$group <- sub("^State3_group=", "", median_os$strata)
median_os$median_OS_months <- median_os$median / 30.44
median_os$lower_95CI_months <- median_os$lower / 30.44
median_os$upper_95CI_months <- median_os$upper / 30.44

state3_group_summary <- do.call(
  rbind,
  lapply(levels(surv_data$State3_group), function(group_name) {
    x <- surv_data[surv_data$State3_group == group_name, , drop = FALSE]
    data.frame(
      group = group_name,
      n_patients = nrow(x),
      n_deaths = sum(x$OS == 1),
      median_State3_fraction = median(x$State3_fraction, na.rm = TRUE),
      mean_State3_fraction = mean(x$State3_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

state3_group_summary <- merge(
  state3_group_summary,
  median_os[, c(
    "group", "median", "lower", "upper",
    "median_OS_months", "lower_95CI_months", "upper_95CI_months"
  )],
  by = "group",
  all.x = TRUE,
  sort = FALSE
)

colnames(state3_group_summary)[
  colnames(state3_group_summary) == "median"
] <- "median_OS_days"
colnames(state3_group_summary)[
  colnames(state3_group_summary) == "lower"
] <- "lower_95CI_days"
colnames(state3_group_summary)[
  colnames(state3_group_summary) == "upper"
] <- "upper_95CI_days"

write.csv(
  state3_group_summary,
  file.path(output_dir, "State3_group_OS_summary.csv"),
  row.names = FALSE
)

# Log-rank test result.
logrank_test <- survdiff(
  Surv(OS.time, OS) ~ State3_group,
  data = surv_data
)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)

logrank_result <- data.frame(
  comparison = "State3-high vs State3-low",
  chisq = unname(logrank_test$chisq),
  degrees_of_freedom = length(logrank_test$n) - 1,
  pvalue = logrank_p,
  stringsAsFactors = FALSE
)

write.csv(
  logrank_result,
  file.path(output_dir, "State3_logrank_test.csv"),
  row.names = FALSE
)

# ============================================================
# 12. Cox regression analyses
# ============================================================
extract_cox_result <- function(model) {
  model_summary <- summary(model)

  data.frame(
    variable = rownames(model_summary$coefficients),
    HR = model_summary$conf.int[, "exp(coef)"],
    lower95CI = model_summary$conf.int[, "lower .95"],
    upper95CI = model_summary$conf.int[, "upper .95"],
    pvalue = model_summary$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

# 12.1 Univariate Cox: State3 abundance as a continuous variable.
cox_continuous <- coxph(
  Surv(OS.time, OS) ~ State3_per10pct,
  data = surv_data,
  ties = "efron"
)

cox_continuous_result <- extract_cox_result(cox_continuous)
cox_continuous_result$interpretation <- "HR per 10-percentage-point increase in estimated State3 abundance"

write.csv(
  cox_continuous_result,
  file.path(output_dir, "State3_continuous_univariate_cox.csv"),
  row.names = FALSE
)

# 12.2 Univariate Cox: State3-high versus State3-low.
cox_group <- coxph(
  Surv(OS.time, OS) ~ State3_group,
  data = surv_data,
  ties = "efron"
)

cox_group_result <- extract_cox_result(cox_group)

write.csv(
  cox_group_result,
  file.path(output_dir, "State3_group_univariate_cox.csv"),
  row.names = FALSE
)

# 12.3 Multivariable Cox adjusted for available clinical variables.
adjustment_candidates <- c("Gender", "Age", "iss_stage", "treatment_type")
adjustment_variables <- intersect(adjustment_candidates, colnames(surv_data))

# Remove variables with fewer than two observed levels.
adjustment_variables <- adjustment_variables[
  vapply(
    adjustment_variables,
    function(v) length(unique(na.omit(surv_data[[v]]))) >= 2,
    logical(1)
  )
]

for (v in adjustment_variables) {
  surv_data[[v]] <- factor(surv_data[[v]])
}

multi_variables <- c(
  "OS.time", "OS", "State3_per10pct", adjustment_variables
)
cox_multi_data <- surv_data[
  complete.cases(surv_data[, multi_variables, drop = FALSE]),
  multi_variables,
  drop = FALSE
]

if (length(adjustment_variables) > 0 && nrow(cox_multi_data) >= 20) {
  multi_formula <- as.formula(
    paste(
      "Surv(OS.time, OS) ~ State3_per10pct +",
      paste(adjustment_variables, collapse = " + ")
    )
  )

  cox_multi <- coxph(
    multi_formula,
    data = cox_multi_data,
    ties = "efron",
    x = TRUE
  )

  cox_multi_result <- extract_cox_result(cox_multi)

  write.csv(
    cox_multi_result,
    file.path(output_dir, "State3_multivariable_cox.csv"),
    row.names = FALSE
  )

  # Proportional-hazards assumption test.
  cox_ph_test <- cox.zph(cox_multi)
  cox_ph_table <- data.frame(
    variable = rownames(cox_ph_test$table),
    chisq = cox_ph_test$table[, "chisq"],
    df = cox_ph_test$table[, "df"],
    pvalue = cox_ph_test$table[, "p"],
    row.names = NULL,
    check.names = FALSE
  )

  write.csv(
    cox_ph_table,
    file.path(output_dir, "State3_multivariable_cox_PH_test.csv"),
    row.names = FALSE
  )
} else {
  warning(
    "The multivariable Cox model was not fitted because no usable adjustment ",
    "variables were available or too few complete patients remained."
  )
}

# ============================================================
# 13. Save the state-fraction distribution plot
# ============================================================
state_fraction_long <- rbind(
  data.frame(
    patient_id = surv_data$patient_id,
    state = "State1",
    fraction = surv_data$State1_fraction,
    stringsAsFactors = FALSE
  ),
  data.frame(
    patient_id = surv_data$patient_id,
    state = "State2",
    fraction = surv_data$State2_fraction,
    stringsAsFactors = FALSE
  ),
  data.frame(
    patient_id = surv_data$patient_id,
    state = "State3",
    fraction = surv_data$State3_fraction,
    stringsAsFactors = FALSE
  )
)

state_fraction_long$state <- factor(
  state_fraction_long$state,
  levels = c("State1", "State2", "State3")
)

fraction_plot <- ggplot(
  state_fraction_long,
  aes(x = state, y = fraction, fill = state)
) +
  geom_violin(trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(
    x = NULL,
    y = "Estimated relative abundance",
    title = "Estimated plasma-cell state abundance in CoMMpass"
  ) +
  scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#E41A1C")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(output_dir, "State_fraction_distribution.pdf"),
  plot = fraction_plot,
  width = 6.5,
  height = 5.5
)

# ============================================================
# 14. Session information and completion message
# ============================================================
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

cat("\nAnalysis completed successfully.\n")
cat("Output directory:\n", output_dir, "\n")
cat("Main files:\n")
cat("  - patient_level_state_abundance_with_State3_group.csv\n")
cat("  - State3_OS_KM_curve.pdf\n")
cat("  - State3_group_OS_summary.csv\n")
cat("  - State3_continuous_univariate_cox.csv\n")
cat("  - State3_multivariable_cox.csv (when adjustment variables are usable)\n")
