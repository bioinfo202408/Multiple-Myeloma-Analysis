# =============================================================================
# Figure 2 add-on analysis
# Linking diagnostic eRNAs to immune-inflammatory pathway activity
#
# Two analyses are implemented:
#   A) Top SHAP eRNA expression vs immune/inflammatory ssGSEA pathway activity
#   B) Sample-level eRNA SHAP burden vs immune/inflammatory ssGSEA pathway activity
#
# This script is designed to reuse outputs from:
#   Fig2_MMFinder_GSVA_immune_correlation_v3_glmnet_fixed.R
#
# Main outputs:
#   - CSV correlation tables
#   - PDF heatmap: Top SHAP eRNAs x immune-inflammatory pathways
#   - PDF dotplot: eRNA SHAP burden x immune-inflammatory pathways
#   - PDF representative scatterplots
# =============================================================================

# ------------------------------
# 0. User configuration
# ------------------------------

set.seed(2024)

# Core project paths, following your previous scripts.
SHAP_ANALYSIS_DIR <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/SHAP_analysis"
TRAIN_MATRIX_PATH <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data.csv"
FEATURE_ANNOTATION_PATH <- file.path(SHAP_ANALYSIS_DIR, "特征类型.csv")

# Reuse ssGSEA outputs from the previous MMFinder score-GSVA script.
GSVA_OUTPUT_DIR <- file.path(SHAP_ANALYSIS_DIR, "GSVA_immune_MMFinder_correlation")
IMMUNE_SCORE_LONG_PATH <- file.path(GSVA_OUTPUT_DIR, "Immune_inflammatory_ssGSEA_scores_long.csv")
MMFINDER_COR_PATH <- file.path(GSVA_OUTPUT_DIR, "MMFinder_score_immune_pathway_correlation_all.csv")

# Output directory for the new eRNA-focused analyses.
OUTPUT_DIR <- file.path(SHAP_ANALYSIS_DIR, "eRNA_immune_GSVA_correlation")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Optional manual paths. Leave NA to auto-detect from SHAP_ANALYSIS_DIR.
# The summary file should contain one row per feature and a mean/consensus SHAP importance column.
CONSENSUS_SHAP_SUMMARY_PATH <- NA_character_

# The sample-level SHAP matrix should contain samples x features, or long format with sample/feature/shap columns.
# This is required for Analysis B. If auto-detection fails, set this path manually.
SAMPLE_SHAP_MATRIX_PATH <- NA_character_

# Analysis parameters.
TOP_ERNA_N <- 100                    # Top SHAP eRNAs used for correlations
TOP_ERNA_HEATMAP_N <- 50            # Number of eRNAs displayed in heatmap
PATHWAY_N <- 14                     # Number of immune pathways displayed
COR_METHOD <- "spearman"
FDR_CUTOFF <- 0.05
TOP_SCATTER_N <- 6

# If TRUE, the script attempts to install missing packages. For HPC/server use, FALSE is safer.
INSTALL_MISSING_PACKAGES <- FALSE

# ------------------------------
# 1. Package loading
# ------------------------------

cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "ggplot2", "stringr", "readr",
  "forcats", "scales", "patchwork", "matrixStats", "purrr"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!INSTALL_MISSING_PACKAGES) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      "\nInstall them first, or set INSTALL_MISSING_PACKAGES <- TRUE."
    )
  }
  install.packages(missing, repos = "https://cloud.r-project.org")
}

install_if_missing(cran_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(forcats)
  library(scales)
  library(patchwork)
  library(matrixStats)
  library(purrr)
})

message("Packages loaded.")

# ------------------------------
# 2. Helper functions
# ------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

clean_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

clean_pathway_label <- function(x) {
  x %>%
    stringr::str_replace("^HALLMARK_", "") %>%
    stringr::str_replace("^REACTOME_", "") %>%
    stringr::str_replace("^GOBP_", "") %>%
    stringr::str_replace("^GO_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_sentence()
}

short_label <- function(x, max_len = 58) {
  x <- clean_pathway_label(x)
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 1), "…"), x)
}

standardize_group <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    stringr::str_detect(tolower(x), "tumor|mm|case|patient|disease") ~ "MM",
    stringr::str_detect(tolower(x), "health|normal|control|ctrl") ~ "Healthy",
    TRUE ~ x
  )
}

standardize_rna_class <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    stringr::str_detect(tolower(x), "^mrna$|messenger") ~ "mRNA",
    stringr::str_detect(tolower(x), "erna|enhancer") ~ "eRNA",
    stringr::str_detect(tolower(x), "lncrna|long") ~ "lncRNA",
    stringr::str_detect(tolower(x), "mirna|mir-") ~ "miRNA",
    TRUE ~ x
  )
}

infer_rna_class <- function(transcript_ids) {
  ids <- as.character(transcript_ids)
  classes <- rep("mRNA", length(ids))
  names(classes) <- ids
  classes[grepl("^Ens[0-9]|enhancer|eRNA|FANTOM", ids, ignore.case = TRUE)] <- "eRNA"
  classes[grepl("^hsa-miR|^hsa-let", ids, ignore.case = TRUE)] <- "miRNA"
  classes[grepl("^NON|^LINC|^MALAT|^HOTAIR|^XIST|^NEAT|^GAS5|^H19|^MEG3|^PVT1|lnc|antisense", ids, ignore.case = TRUE)] <- "lncRNA"
  classes
}

load_feature_classes <- function(feature_ids, anno_path) {
  feature_ids <- as.character(feature_ids)

  if (!is.null(anno_path) && !is.na(anno_path) && file.exists(anno_path)) {
    anno <- read.csv(anno_path, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(anno) <- trimws(colnames(anno))
    id_col <- intersect(c("transcript_id", "feature", "gene", "id", "Feature", "ID", "feature_id"), colnames(anno))[1]
    class_col <- intersect(c("rna_class", "RNA_class", "type", "class", "Type", "Class", "RNA_type", "feature_type"), colnames(anno))[1]

    if (!is.na(id_col) && !is.na(class_col)) {
      anno_ids <- as.character(anno[[id_col]])
      fc <- stats::setNames(standardize_rna_class(anno[[class_col]]), anno_ids)

      # Exact matching first.
      out <- fc[feature_ids]

      # Clean-name matching for any missing values.
      miss <- is.na(out)
      if (any(miss)) {
        clean_map <- stats::setNames(fc, clean_name(names(fc)))
        out[miss] <- clean_map[clean_name(feature_ids[miss])]
      }

      # Fall back to rule-based inference for remaining missing values.
      miss <- is.na(out)
      if (any(miss)) out[miss] <- infer_rna_class(feature_ids[miss])
      names(out) <- feature_ids
      return(out)
    }
  }

  infer_rna_class(feature_ids)
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey12"),
      axis.line = ggplot2::element_line(linewidth = 0.42, color = "grey25"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, color = "grey25"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", color = "grey15")
    )
}

save_pdf <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    device = "pdf"
  )
}

safe_cor_test <- function(x, y, method = COR_METHOD) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(c(rho = NA_real_, pvalue = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(rho = unname(ct$estimate), pvalue = ct$p.value, n = sum(ok))
}

read_csv_flexible <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# ------------------------------
# 3. Load immune-inflammatory ssGSEA scores
# ------------------------------

if (!file.exists(IMMUNE_SCORE_LONG_PATH)) {
  stop(
    "Cannot find immune ssGSEA score file: ", IMMUNE_SCORE_LONG_PATH,
    "\nRun the previous MMFinder score-GSVA script first."
  )
}

message("Loading immune/inflammatory ssGSEA scores...")
score_long <- read_csv_flexible(IMMUNE_SCORE_LONG_PATH)
colnames(score_long) <- trimws(colnames(score_long))

required_score_cols <- c("sample_id", "pathway", "pathway_score")
missing_score_cols <- setdiff(required_score_cols, colnames(score_long))
if (length(missing_score_cols) > 0) {
  stop("ssGSEA long file is missing required columns: ", paste(missing_score_cols, collapse = ", "))
}

if ("group" %in% colnames(score_long)) score_long$group <- standardize_group(score_long$group)

score_wide <- score_long %>%
  dplyr::select(.data$sample_id, .data$pathway, .data$pathway_score) %>%
  tidyr::pivot_wider(names_from = .data$pathway, values_from = .data$pathway_score) %>%
  as.data.frame()

rownames(score_wide) <- score_wide$sample_id
score_mat <- as.matrix(score_wide[, setdiff(colnames(score_wide), "sample_id"), drop = FALSE])
storage.mode(score_mat) <- "double"

sample_metadata <- score_long %>%
  dplyr::select(dplyr::any_of(c("sample_id", "group", "MMFinder_score"))) %>%
  dplyr::distinct(.data$sample_id, .keep_all = TRUE)

message("Immune pathway score matrix: ", nrow(score_mat), " samples x ", ncol(score_mat), " pathways.")

# Select representative immune pathways for display. Prefer pathways positively associated with MMFinder score.
if (file.exists(MMFINDER_COR_PATH)) {
  mmfinder_cor <- read_csv_flexible(MMFINDER_COR_PATH)
  colnames(mmfinder_cor) <- trimws(colnames(mmfinder_cor))
  if ("FDR" %in% colnames(mmfinder_cor)) {
    mmfinder_cor$sort_p <- ifelse(is.finite(mmfinder_cor$FDR), mmfinder_cor$FDR, 1)
  } else if ("pvalue" %in% colnames(mmfinder_cor)) {
    mmfinder_cor$sort_p <- ifelse(is.finite(mmfinder_cor$pvalue), mmfinder_cor$pvalue, 1)
  } else {
    mmfinder_cor$sort_p <- 1
  }

  selected_pathways <- mmfinder_cor %>%
    dplyr::filter(.data$pathway %in% colnames(score_mat)) %>%
    dplyr::filter(is.finite(.data$rho), .data$rho > 0) %>%
    dplyr::arrange(.data$sort_p, dplyr::desc(.data$rho)) %>%
    dplyr::slice_head(n = PATHWAY_N) %>%
    dplyr::pull(.data$pathway)

  if (length(selected_pathways) < 3) {
    selected_pathways <- mmfinder_cor %>%
      dplyr::filter(.data$pathway %in% colnames(score_mat)) %>%
      dplyr::arrange(.data$sort_p, dplyr::desc(abs(.data$rho))) %>%
      dplyr::slice_head(n = PATHWAY_N) %>%
      dplyr::pull(.data$pathway)
  }
} else {
  selected_pathways <- colnames(score_mat)[seq_len(min(PATHWAY_N, ncol(score_mat)))]
}

selected_pathways <- unique(selected_pathways[selected_pathways %in% colnames(score_mat)])
if (length(selected_pathways) == 0) stop("No immune pathways selected for display.")

readr::write_csv(
  tibble::tibble(pathway = selected_pathways, pathway_label = clean_pathway_label(selected_pathways)),
  file.path(OUTPUT_DIR, "Selected_immune_inflammatory_pathways_for_eRNA_analysis.csv")
)

message("Selected pathways: ", paste(selected_pathways, collapse = "; "))

# ------------------------------
# 4. Load expression matrix and identify eRNA features
# ------------------------------

if (!file.exists(TRAIN_MATRIX_PATH)) stop("Training expression matrix not found: ", TRAIN_MATRIX_PATH)

message("Loading training expression matrix...")
expr_raw <- read.csv(TRAIN_MATRIX_PATH, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

numeric_cols <- vapply(expr_raw, is.numeric, logical(1))
if (sum(!numeric_cols) > 0) {
  message("Dropping non-numeric columns from training matrix: ", paste(colnames(expr_raw)[!numeric_cols], collapse = ", "))
}
expr <- as.matrix(expr_raw[, numeric_cols, drop = FALSE])
storage.mode(expr) <- "double"

# Align expression rows with ssGSEA samples.
common_samples <- intersect(rownames(expr), rownames(score_mat))
if (length(common_samples) >= 5) {
  expr <- expr[common_samples, , drop = FALSE]
  score_mat_aligned <- score_mat[common_samples, , drop = FALSE]
  sample_metadata_aligned <- sample_metadata %>% dplyr::filter(.data$sample_id %in% common_samples)
  sample_metadata_aligned <- sample_metadata_aligned[match(common_samples, sample_metadata_aligned$sample_id), , drop = FALSE]
  message("Expression and ssGSEA scores aligned by sample_id: ", length(common_samples), " samples.")
} else if (nrow(expr) == nrow(score_mat)) {
  rownames(expr) <- rownames(score_mat)
  score_mat_aligned <- score_mat
  sample_metadata_aligned <- sample_metadata[match(rownames(score_mat), sample_metadata$sample_id), , drop = FALSE]
  message("Expression and ssGSEA scores aligned by row order. Please verify if sample IDs are unavailable.")
} else {
  stop(
    "Cannot align expression matrix with ssGSEA score matrix. ",
    "Expression rows = ", nrow(expr), "; score rows = ", nrow(score_mat), "; common IDs = ", length(common_samples), "."
  )
}

feature_ids <- colnames(expr)
feature_classes <- load_feature_classes(feature_ids, FEATURE_ANNOTATION_PATH)
feature_classes <- standardize_rna_class(feature_classes)
names(feature_classes) <- feature_ids

erna_features <- feature_ids[feature_classes == "eRNA"]
if (length(erna_features) == 0) {
  stop(
    "No eRNA features were detected. Check FEATURE_ANNOTATION_PATH or feature naming rules.\n",
    "Feature annotation path used: ", FEATURE_ANNOTATION_PATH
  )
}
message("Detected eRNA features in expression matrix: ", length(erna_features))

readr::write_csv(
  tibble::tibble(feature = feature_ids, rna_class = feature_classes),
  file.path(OUTPUT_DIR, "Feature_RNA_class_used_for_eRNA_analysis.csv")
)

# ------------------------------
# 5. Load or infer SHAP feature importance; rank Top SHAP eRNAs
# ------------------------------

find_candidate_file <- function(root_dir, include_patterns, exclude_patterns = character(), extensions = c("csv", "tsv", "txt", "rds")) {
  if (!dir.exists(root_dir)) return(NA_character_)
  files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files)]
  files <- files[grepl(paste0("\\.(", paste(extensions, collapse = "|"), ")$"), files, ignore.case = TRUE)]
  if (length(files) == 0) return(NA_character_)

  base <- basename(files)
  keep <- rep(TRUE, length(files))
  for (pat in include_patterns) keep <- keep & grepl(pat, base, ignore.case = TRUE)
  for (pat in exclude_patterns) keep <- keep & !grepl(pat, base, ignore.case = TRUE)
  cand <- files[keep]
  if (length(cand) == 0) return(NA_character_)
  cand[order(nchar(basename(cand)))][1]
}

normalize_shap_summary <- function(df) {
  colnames(df) <- trimws(colnames(df))
  feature_col <- intersect(
    c("feature", "Feature", "feature_id", "transcript_id", "gene", "Gene", "id", "ID", "variable", "Variable"),
    colnames(df)
  )[1]

  shap_col <- intersect(
    c(
      "consensus_shap", "Consensus_SHAP", "mean_abs_shap", "mean_abs_SHAP", "avg_abs_shap",
      "average_abs_shap", "abs_mean_shap", "mean_SHAP", "mean_shap", "importance",
      "Importance", "shap_importance", "SHAP_importance", "mean_abs", "MeanDecreaseGini"
    ),
    colnames(df)
  )[1]

  if (is.na(feature_col)) return(NULL)

  if (is.na(shap_col)) {
    numeric_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, feature_col)
    if (length(numeric_cols) == 0) return(NULL)
    # Prefer columns whose names contain shap/importance.
    pref <- numeric_cols[grepl("shap|importance|mean|abs", numeric_cols, ignore.case = TRUE)]
    shap_col <- if (length(pref) > 0) pref[1] else numeric_cols[1]
  }

  out <- df %>%
    dplyr::transmute(
      feature = as.character(.data[[feature_col]]),
      shap_importance = as.numeric(.data[[shap_col]])
    ) %>%
    dplyr::filter(!is.na(.data$feature), .data$feature != "", is.finite(.data$shap_importance)) %>%
    dplyr::group_by(.data$feature) %>%
    dplyr::summarise(shap_importance = max(abs(.data$shap_importance), na.rm = TRUE), .groups = "drop")

  if (nrow(out) == 0) return(NULL)
  out
}

load_shap_summary_auto <- function(path, root_dir) {
  tried <- character()

  if (!is.na(path) && file.exists(path)) {
    tried <- c(tried, path)
    obj <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else read_csv_flexible(path)
    if (is.matrix(obj)) obj <- as.data.frame(obj)
    if (is.data.frame(obj)) {
      out <- normalize_shap_summary(obj)
      if (!is.null(out)) return(list(summary = out, source = path, tried = tried))
    }
  }

  # Try common summary/ranking file names first.
  candidates <- c(
    find_candidate_file(root_dir, c("shap", "summary"), c("sample|matrix|long|burden|correlation")),
    find_candidate_file(root_dir, c("shap", "rank"), c("sample|matrix|long|burden|correlation")),
    find_candidate_file(root_dir, c("consensus", "shap"), c("sample|matrix|long|burden|correlation")),
    find_candidate_file(root_dir, c("feature", "importance"), c("sample|matrix|long|burden|correlation")),
    find_candidate_file(root_dir, c("mean", "shap"), c("sample|matrix|long|burden|correlation"))
  )
  candidates <- unique(candidates[!is.na(candidates)])

  for (f in candidates) {
    tried <- c(tried, f)
    obj <- tryCatch({
      if (grepl("\\.rds$", f, ignore.case = TRUE)) readRDS(f) else read_csv_flexible(f)
    }, error = function(e) NULL)
    if (is.null(obj)) next
    if (is.matrix(obj)) obj <- as.data.frame(obj)
    if (!is.data.frame(obj)) next
    out <- normalize_shap_summary(obj)
    if (!is.null(out)) return(list(summary = out, source = f, tried = tried))
  }

  list(summary = NULL, source = NA_character_, tried = tried)
}

# Forward declaration: sample-level SHAP matrix loading also allows fallback ranking.
load_sample_shap_matrix <- function(path, root_dir, expected_samples, known_features) {
  candidate_paths <- character()

  if (!is.na(path) && file.exists(path)) candidate_paths <- c(candidate_paths, path)

  auto_candidates <- c(
    find_candidate_file(root_dir, c("shap", "matrix"), c("summary|rank|importance|burden|correlation")),
    find_candidate_file(root_dir, c("shap", "values"), c("summary|rank|importance|burden|correlation")),
    find_candidate_file(root_dir, c("shap", "sample"), c("summary|rank|importance|burden|correlation")),
    find_candidate_file(root_dir, c("shap", "train"), c("summary|rank|importance|burden|correlation")),
    find_candidate_file(root_dir, c("sample", "shap"), c("summary|rank|importance|burden|correlation"))
  )
  candidate_paths <- unique(c(candidate_paths, auto_candidates[!is.na(auto_candidates)]))

  tried <- character()

  for (f in candidate_paths) {
    tried <- c(tried, f)
    message("Trying sample-level SHAP file: ", f)

    obj <- tryCatch({
      if (grepl("\\.rds$", f, ignore.case = TRUE)) readRDS(f) else read_csv_flexible(f)
    }, error = function(e) {
      message("  Failed to read: ", e$message)
      NULL
    })
    if (is.null(obj)) next

    if (is.matrix(obj)) obj <- as.data.frame(obj)

    # Long format: sample_id / feature / shap_value.
    if (is.data.frame(obj)) {
      colnames(obj) <- trimws(colnames(obj))
      sample_col <- intersect(c("sample_id", "sample", "Sample", "SampleID", "ID", "X"), colnames(obj))[1]
      feature_col <- intersect(c("feature", "Feature", "feature_id", "transcript_id", "gene", "Gene", "variable"), colnames(obj))[1]
      shap_col <- intersect(c("shap", "SHAP", "shap_value", "SHAP_value", "value", "phi", "contribution"), colnames(obj))[1]

      if (!is.na(sample_col) && !is.na(feature_col) && !is.na(shap_col)) {
        long_df <- obj %>%
          dplyr::transmute(
            sample_id = as.character(.data[[sample_col]]),
            feature = as.character(.data[[feature_col]]),
            shap_value = as.numeric(.data[[shap_col]])
          ) %>%
          dplyr::filter(.data$sample_id %in% expected_samples, .data$feature %in% known_features, is.finite(.data$shap_value))

        if (nrow(long_df) > 0) {
          wide <- long_df %>%
            tidyr::pivot_wider(names_from = .data$feature, values_from = .data$shap_value, values_fill = 0) %>%
            as.data.frame()
          rownames(wide) <- wide$sample_id
          mat <- as.matrix(wide[, setdiff(colnames(wide), "sample_id"), drop = FALSE])
          storage.mode(mat) <- "double"
          return(list(matrix = mat[intersect(expected_samples, rownames(mat)), , drop = FALSE], source = f, tried = tried))
        }
      }

      # Wide format: samples x features, maybe first column is sample ID.
      df <- obj
      if (ncol(df) >= 2 && !is.numeric(df[[1]]) && !is.logical(df[[1]])) {
        first_col <- colnames(df)[1]
        if (any(as.character(df[[first_col]]) %in% expected_samples)) {
          rownames(df) <- as.character(df[[first_col]])
          df[[first_col]] <- NULL
        }
      }

      # Keep only numeric columns when features are in columns.
      numeric_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
      row_ids <- rownames(df)

      # Case 1: rows are samples and columns are features.
      feature_overlap <- intersect(numeric_cols, known_features)
      sample_overlap <- intersect(row_ids, expected_samples)
      if (length(sample_overlap) >= 5 && length(feature_overlap) >= 2) {
        mat <- as.matrix(df[sample_overlap, feature_overlap, drop = FALSE])
        storage.mode(mat) <- "double"
        return(list(matrix = mat, source = f, tried = tried))
      }

      # Case 2: row order equals expected samples and columns are known features.
      if (nrow(df) == length(expected_samples) && length(feature_overlap) >= 2) {
        mat <- as.matrix(df[, feature_overlap, drop = FALSE])
        rownames(mat) <- expected_samples
        storage.mode(mat) <- "double"
        return(list(matrix = mat, source = f, tried = tried))
      }

      # Case 3: rows are features and columns are samples; transpose.
      col_sample_overlap <- intersect(colnames(df), expected_samples)
      row_feature_overlap <- intersect(row_ids, known_features)
      if (length(col_sample_overlap) >= 5 && length(row_feature_overlap) >= 2) {
        mat <- t(as.matrix(df[row_feature_overlap, col_sample_overlap, drop = FALSE]))
        storage.mode(mat) <- "double"
        return(list(matrix = mat, source = f, tried = tried))
      }
    }
  }

  list(matrix = NULL, source = NA_character_, tried = tried)
}

message("Loading SHAP feature-importance summary...")
shap_summary_obj <- load_shap_summary_auto(CONSENSUS_SHAP_SUMMARY_PATH, SHAP_ANALYSIS_DIR)
shap_summary <- shap_summary_obj$summary
sample_shap_obj <- NULL

if (is.null(shap_summary)) {
  message("No valid SHAP summary file auto-detected. Trying sample-level SHAP matrix to derive mean absolute SHAP ranking...")
  sample_shap_obj <- load_sample_shap_matrix(
    SAMPLE_SHAP_MATRIX_PATH,
    SHAP_ANALYSIS_DIR,
    expected_samples = rownames(expr),
    known_features = feature_ids
  )
  if (!is.null(sample_shap_obj$matrix)) {
    shap_summary <- tibble::tibble(
      feature = colnames(sample_shap_obj$matrix),
      shap_importance = matrixStats::colMeans2(abs(sample_shap_obj$matrix), na.rm = TRUE)
    )
    shap_summary_obj$source <- paste0(sample_shap_obj$source, " [mean(abs(sample-level SHAP))]")
  }
}

if (is.null(shap_summary)) {
  stop(
    "Could not detect a SHAP summary or sample-level SHAP matrix.\n",
    "Set CONSENSUS_SHAP_SUMMARY_PATH or SAMPLE_SHAP_MATRIX_PATH manually near the top of this script.\n",
    "Files tried for SHAP summary: ", paste(shap_summary_obj$tried, collapse = "; ")
  )
}

# Match summary features to expression features.
feature_clean_to_actual <- stats::setNames(feature_ids, clean_name(feature_ids))
shap_summary <- shap_summary %>%
  dplyr::mutate(
    feature_raw = .data$feature,
    feature = dplyr::if_else(.data$feature %in% feature_ids, .data$feature, feature_clean_to_actual[clean_name(.data$feature)]),
    feature = as.character(.data$feature)
  ) %>%
  dplyr::filter(!is.na(.data$feature), .data$feature %in% feature_ids) %>%
  dplyr::mutate(rna_class = feature_classes[.data$feature]) %>%
  dplyr::arrange(dplyr::desc(.data$shap_importance)) %>%
  dplyr::distinct(.data$feature, .keep_all = TRUE)

erna_shap_rank <- shap_summary %>%
  dplyr::filter(.data$rna_class == "eRNA", .data$feature %in% erna_features) %>%
  dplyr::arrange(dplyr::desc(.data$shap_importance)) %>%
  dplyr::mutate(shap_rank_erNA = dplyr::row_number())

if (nrow(erna_shap_rank) == 0) {
  stop("SHAP ranking was loaded, but no ranked eRNA features overlapped with expression matrix.")
}

top_erna <- erna_shap_rank %>%
  dplyr::slice_head(n = TOP_ERNA_N) %>%
  dplyr::pull(.data$feature)

top_erna_heatmap_n <- min(TOP_ERNA_HEATMAP_N, nrow(erna_shap_rank))

top_erna_heatmap <- erna_shap_rank %>%
   dplyr::slice_head(n = top_erna_heatmap_n) %>%
   dplyr::pull(.data$feature)

readr::write_csv(shap_summary, file.path(OUTPUT_DIR, "SHAP_feature_importance_used_all_features.csv"))
readr::write_csv(erna_shap_rank, file.path(OUTPUT_DIR, "Top_SHAP_eRNA_ranking_used.csv"))
writeLines(
  c(
    paste0("SHAP summary source: ", shap_summary_obj$source),
    paste0("Number of eRNAs ranked by SHAP: ", nrow(erna_shap_rank)),
    paste0("Top eRNAs used for correlation: ", length(top_erna)),
    paste0("Top eRNAs shown in heatmap: ", length(top_erna_heatmap))
  ),
  file.path(OUTPUT_DIR, "SHAP_source_and_eRNA_ranking_summary.txt")
)

message("SHAP source: ", shap_summary_obj$source)
message("Top SHAP eRNAs selected: ", length(top_erna))

# ------------------------------
# 6. Analysis A: Top SHAP eRNA expression vs immune pathway activity
# ------------------------------

message("Running Analysis A: Top SHAP eRNA expression vs immune pathway activity...")

erna_expr <- expr[, top_erna, drop = FALSE]
pathway_scores_selected <- score_mat_aligned[, selected_pathways, drop = FALSE]

cor_grid <- tidyr::expand_grid(
  eRNA = colnames(erna_expr),
  pathway = colnames(pathway_scores_selected)
)

cor_results_erna_expr <- purrr::pmap_dfr(cor_grid, function(eRNA, pathway) {
  res <- safe_cor_test(erna_expr[, eRNA], pathway_scores_selected[, pathway], method = COR_METHOD)
  tibble::tibble(
    eRNA = eRNA,
    pathway = pathway,
    rho = as.numeric(res["rho"]),
    pvalue = as.numeric(res["pvalue"]),
    n_samples = as.numeric(res["n"])
  )
}) %>%
  dplyr::mutate(
    FDR = p.adjust(.data$pvalue, method = "BH"),
    direction = dplyr::case_when(.data$rho > 0 ~ "Positive", .data$rho < 0 ~ "Negative", TRUE ~ "Zero"),
    pathway_label = short_label(.data$pathway, max_len = 54),
    eRNA_rank = match(.data$eRNA, erna_shap_rank$feature),
    eRNA_label = paste0(.data$eRNA, "  (#", .data$eRNA_rank, ")"),
    sig_label = dplyr::case_when(
      .data$FDR < 0.001 ~ "***",
      .data$FDR < 0.01 ~ "**",
      .data$FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    minus_log10_p = -log10(pmax(.data$pvalue, .Machine$double.xmin)),
    minus_log10_FDR = -log10(pmax(.data$FDR, .Machine$double.xmin))
  ) %>%
  dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$rho)))

readr::write_csv(cor_results_erna_expr, file.path(OUTPUT_DIR, "Top_SHAP_eRNA_expression_vs_immune_pathway_correlation_all.csv"))
readr::write_csv(
  cor_results_erna_expr %>% dplyr::filter(.data$FDR < FDR_CUTOFF),
  file.path(OUTPUT_DIR, "Top_SHAP_eRNA_expression_vs_immune_pathway_correlation_FDR_lt_0.05.csv")
)

# eRNA-level summary for interpretation.
erna_expr_summary <- cor_results_erna_expr %>%
  dplyr::group_by(.data$eRNA) %>%
  dplyr::summarise(
    eRNA_rank = min(.data$eRNA_rank, na.rm = TRUE),
    shap_importance = erna_shap_rank$shap_importance[match(dplyr::first(.data$eRNA), erna_shap_rank$feature)],
    n_positive_FDR_lt_0.05 = sum(.data$rho > 0 & .data$FDR < FDR_CUTOFF, na.rm = TRUE),
    n_negative_FDR_lt_0.05 = sum(.data$rho < 0 & .data$FDR < FDR_CUTOFF, na.rm = TRUE),
    max_positive_rho = suppressWarnings(max(.data$rho, na.rm = TRUE)),
    min_negative_rho = suppressWarnings(min(.data$rho, na.rm = TRUE)),
    strongest_pathway = .data$pathway[which.max(abs(.data$rho))],
    strongest_rho = .data$rho[which.max(abs(.data$rho))],
    strongest_FDR = .data$FDR[which.max(abs(.data$rho))],
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$eRNA_rank)

readr::write_csv(erna_expr_summary, file.path(OUTPUT_DIR, "Top_SHAP_eRNA_immune_association_summary_by_eRNA.csv"))

# Plot A1: dot heatmap of Top SHAP eRNAs x immune pathways.
heatmap_df <- cor_results_erna_expr %>%
  dplyr::filter(.data$eRNA %in% top_erna_heatmap, .data$pathway %in% selected_pathways)

# Order pathways by MMFinder correlation if available; otherwise by mean eRNA correlation.
if (exists("mmfinder_cor") && all(c("pathway", "rho") %in% colnames(mmfinder_cor))) {
  pathway_order <- mmfinder_cor %>%
    dplyr::filter(.data$pathway %in% selected_pathways) %>%
    dplyr::arrange(.data$rho) %>%
    dplyr::pull(.data$pathway)
} else {
  pathway_order <- heatmap_df %>%
    dplyr::group_by(.data$pathway) %>%
    dplyr::summarise(mean_rho = mean(.data$rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(.data$mean_rho) %>%
    dplyr::pull(.data$pathway)
}
pathway_order <- unique(pathway_order[pathway_order %in% selected_pathways])

# Order eRNAs by hierarchical clustering of rho profile when possible.
rho_mat <- heatmap_df %>%
  dplyr::select(.data$eRNA, .data$pathway, .data$rho) %>%
  tidyr::pivot_wider(names_from = .data$pathway, values_from = .data$rho) %>%
  as.data.frame()
rownames(rho_mat) <- rho_mat$eRNA
rho_mat$eRNA <- NULL
rho_mat <- as.matrix(rho_mat[, intersect(pathway_order, colnames(rho_mat)), drop = FALSE])
rho_mat[!is.finite(rho_mat)] <- 0

if (nrow(rho_mat) >= 2 && ncol(rho_mat) >= 2) {
  erna_order <- rownames(rho_mat)[stats::hclust(stats::dist(rho_mat), method = "ward.D2")$order]
} else {
  erna_order <- top_erna_heatmap
}

heatmap_df <- heatmap_df %>%
  dplyr::mutate(
    eRNA_label = paste0(.data$eRNA, "  #", .data$eRNA_rank),
    eRNA_label = factor(.data$eRNA_label, levels = paste0(rev(erna_order), "  #", match(rev(erna_order), erna_shap_rank$feature))),
    pathway_label = factor(short_label(.data$pathway, 54), levels = short_label(pathway_order, 54))
  )

p_heat <- ggplot(heatmap_df, aes(x = pathway_label, y = eRNA_label)) +
  geom_tile(aes(fill = rho), color = "white", linewidth = 0.28) +
  geom_point(aes(size = pmin(minus_log10_FDR, 10)), shape = 21, fill = "white", color = "grey20", stroke = 0.25, alpha = 0.85) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-max(abs(heatmap_df$rho), na.rm = TRUE), max(abs(heatmap_df$rho), na.rm = TRUE)),
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(0.2, 3.1), breaks = c(1.3, 2, 3, 5), limits = c(0, 10)) +
  labs(
    title = "High-contribution eRNAs are linked to immune-inflammatory pathway activity",
    subtitle = paste0("Top SHAP-ranked eRNA expression vs ssGSEA pathway scores; ", COR_METHOD, " correlation"),
    x = NULL,
    y = "Top SHAP eRNAs"
  ) +
  theme_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.1),
    axis.text.y = element_text(size = 7.6),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.35),
    legend.position = "right"
  )

save_pdf(
  p_heat,
  file.path(OUTPUT_DIR, "Fig2M_Top_SHAP_eRNA_expression_immune_pathway_correlation_heatmap.pdf"),
  width = max(10, 0.46 * length(pathway_order) + 5.6),
  height = max(7, 0.24 * length(top_erna_heatmap) + 2.3)
)

# Plot A2: rank eRNAs by number of significant positive immune-pathway associations.
p_erna_summary_n <- min(20, nrow(erna_expr_summary))

p_erna_summary_df <- erna_expr_summary %>%
   dplyr::slice_head(n = p_erna_summary_n) %>%
   dplyr::mutate(
      eRNA_label = paste0(.data$eRNA, "  #", .data$eRNA_rank),
      eRNA_label = forcats::fct_reorder(
         .data$eRNA_label,
         .data$n_positive_FDR_lt_0.05 + 0.1 * .data$max_positive_rho
      )
   )

p_erna_summary <- ggplot(p_erna_summary_df, aes(x = n_positive_FDR_lt_0.05, y = eRNA_label)) +
  geom_segment(aes(x = 0, xend = n_positive_FDR_lt_0.05, yend = eRNA_label), linewidth = 0.45, color = "grey70") +
  geom_point(aes(size = shap_importance, color = max_positive_rho), alpha = 0.92) +
  scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "Max positive\nrho") +
  scale_size_continuous(name = "SHAP\nimportance", range = c(2, 6)) +
  labs(
    title = "Immune-inflammatory coupling of top diagnostic eRNAs",
    subtitle = paste0("Number of positively correlated pathways at FDR < ", FDR_CUTOFF),
    x = "Number of significant positive immune-pathway associations",
    y = NULL
  ) +
  theme_pub(base_size = 10) +
  theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25))

save_pdf(
  p_erna_summary,
  file.path(OUTPUT_DIR, "Fig2N_Top_SHAP_eRNA_immune_association_lollipop.pdf"),
  width = 8.5,
  height = max(5.2, 0.24 * nrow(p_erna_summary_df) + 1.7)
)

# Plot A3: representative eRNA-pathway scatter plots.
top_pairs_expr <- cor_results_erna_expr %>%
  dplyr::filter(.data$rho > 0, .data$pvalue < 0.05) %>%
  dplyr::arrange(.data$FDR, dplyr::desc(.data$rho)) %>%
  dplyr::slice_head(n = TOP_SCATTER_N)

if (nrow(top_pairs_expr) > 0) {
  scatter_expr_df <- purrr::map_dfr(seq_len(nrow(top_pairs_expr)), function(i) {
    e <- top_pairs_expr$eRNA[i]
    p <- top_pairs_expr$pathway[i]
    tibble::tibble(
      sample_id = rownames(expr),
      eRNA_expression = expr[, e],
      pathway_score = score_mat_aligned[, p],
      group = sample_metadata_aligned$group %||% NA_character_,
      pair_label = paste0(
        e, "\n", short_label(p, 46),
        "\nrho=", sprintf("%.2f", top_pairs_expr$rho[i]),
        ", FDR=", scales::scientific(top_pairs_expr$FDR[i], digits = 2)
      )
    )
  })

  p_scatter_expr <- ggplot(scatter_expr_df, aes(x = eRNA_expression, y = pathway_score)) +
    geom_point(aes(fill = group), shape = 21, size = 2.2, color = "white", stroke = 0.25, alpha = 0.82) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "grey20", fill = "grey82") +
    facet_wrap(~ pair_label, scales = "free", ncol = 2) +
    labs(
      title = "Representative eRNA–immune pathway associations",
      subtitle = "Top positive correlations between high-SHAP eRNA expression and immune-inflammatory pathway activity",
      x = "eRNA expression",
      y = "ssGSEA pathway score",
      fill = "Group"
    ) +
    theme_pub(base_size = 9) +
    theme(legend.position = "bottom")

  save_pdf(
    p_scatter_expr,
    file.path(OUTPUT_DIR, "Fig2O_representative_Top_SHAP_eRNA_expression_immune_pathway_scatterplots.pdf"),
    width = 9.5,
    height = 3.15 * ceiling(nrow(top_pairs_expr) / 2)
  )
}

# ------------------------------
# 7. Analysis B: eRNA SHAP burden vs immune pathway activity
# ------------------------------

message("Running Analysis B: eRNA SHAP burden vs immune pathway activity...")

if (is.null(sample_shap_obj)) {
  sample_shap_obj <- load_sample_shap_matrix(
    SAMPLE_SHAP_MATRIX_PATH,
    SHAP_ANALYSIS_DIR,
    expected_samples = rownames(expr),
    known_features = feature_ids
  )
}

if (is.null(sample_shap_obj$matrix)) {
  warning(
    "Sample-level SHAP matrix could not be loaded. Analysis B was skipped.\n",
    "Set SAMPLE_SHAP_MATRIX_PATH manually near the top of this script.\n",
    "Files tried: ", paste(sample_shap_obj$tried, collapse = "; ")
  )

  writeLines(
    c(
      "Analysis B skipped: sample-level SHAP matrix could not be detected.",
      "To enable eRNA SHAP burden analysis, set SAMPLE_SHAP_MATRIX_PATH manually.",
      "Expected formats:",
      "  1) wide samples x features SHAP matrix, with sample IDs as rownames or first column; or",
      "  2) long table with columns sample_id, feature, shap_value.",
      paste0("Files tried: ", paste(sample_shap_obj$tried, collapse = "; "))
    ),
    file.path(OUTPUT_DIR, "Analysis_B_eRNA_SHAP_burden_SKIPPED_readme.txt")
  )
} else {
  shap_mat <- sample_shap_obj$matrix
  # Align sample-level SHAP with expression/score samples.
  common_shap_samples <- intersect(rownames(expr), rownames(shap_mat))
  if (length(common_shap_samples) >= 5) {
    shap_mat <- shap_mat[common_shap_samples, , drop = FALSE]
    score_for_burden <- score_mat_aligned[common_shap_samples, , drop = FALSE]
    meta_for_burden <- sample_metadata_aligned[match(common_shap_samples, sample_metadata_aligned$sample_id), , drop = FALSE]
  } else if (nrow(shap_mat) == nrow(expr)) {
    rownames(shap_mat) <- rownames(expr)
    score_for_burden <- score_mat_aligned
    meta_for_burden <- sample_metadata_aligned
  } else {
    stop(
      "Sample-level SHAP matrix loaded but cannot be aligned with expression/ssGSEA samples. ",
      "SHAP rows = ", nrow(shap_mat), "; expression rows = ", nrow(expr), "."
    )
  }

  # Clean/match SHAP feature names to expression feature names.
  shap_features_raw <- colnames(shap_mat)
  matched_features <- ifelse(shap_features_raw %in% feature_ids, shap_features_raw, feature_clean_to_actual[clean_name(shap_features_raw)])
  keep_features <- !is.na(matched_features) & matched_features %in% feature_ids
  shap_mat <- shap_mat[, keep_features, drop = FALSE]
  colnames(shap_mat) <- matched_features[keep_features]

  # Collapse duplicated features if any.
  if (anyDuplicated(colnames(shap_mat)) > 0) {
    split_idx <- split(seq_len(ncol(shap_mat)), colnames(shap_mat))
    collapsed <- do.call(cbind, lapply(split_idx, function(ii) {
      if (length(ii) == 1) shap_mat[, ii] else rowSums(shap_mat[, ii, drop = FALSE], na.rm = TRUE)
    }))
    rownames(collapsed) <- rownames(shap_mat)
    shap_mat <- as.matrix(collapsed)
  }
  storage.mode(shap_mat) <- "double"

  erna_shap_features <- intersect(erna_features, colnames(shap_mat))
  top_erna_shap_features <- intersect(top_erna, colnames(shap_mat))

  if (length(erna_shap_features) < 2) {
    stop("Sample-level SHAP matrix loaded, but fewer than 2 eRNA SHAP columns overlap with feature annotation.")
  }

  total_abs_shap <- rowSums(abs(shap_mat), na.rm = TRUE)
  erna_abs_burden <- rowSums(abs(shap_mat[, erna_shap_features, drop = FALSE]), na.rm = TRUE)
  erna_signed_burden <- rowSums(shap_mat[, erna_shap_features, drop = FALSE], na.rm = TRUE)

  if (length(top_erna_shap_features) >= 2) {
    top_erna_abs_burden <- rowSums(abs(shap_mat[, top_erna_shap_features, drop = FALSE]), na.rm = TRUE)
    top_erna_signed_burden <- rowSums(shap_mat[, top_erna_shap_features, drop = FALSE], na.rm = TRUE)
  } else {
    top_erna_abs_burden <- rep(NA_real_, nrow(shap_mat))
    top_erna_signed_burden <- rep(NA_real_, nrow(shap_mat))
    warning("Fewer than 2 Top SHAP eRNAs overlap with sample-level SHAP matrix; Top eRNA burden metrics set to NA.")
  }

  burden_df <- tibble::tibble(
    sample_id = rownames(shap_mat),
    group = meta_for_burden$group %||% NA_character_,
    MMFinder_score = meta_for_burden$MMFinder_score %||% NA_real_,
    eRNA_abs_SHAP_burden = erna_abs_burden,
    eRNA_signed_SHAP_burden = erna_signed_burden,
    Top_eRNA_abs_SHAP_burden = top_erna_abs_burden,
    Top_eRNA_signed_SHAP_burden = top_erna_signed_burden,
    eRNA_SHAP_fraction = erna_abs_burden / pmax(total_abs_shap, .Machine$double.eps)
  )

  readr::write_csv(burden_df, file.path(OUTPUT_DIR, "Sample_level_eRNA_SHAP_burden_scores.csv"))

  burden_metrics <- c(
    "eRNA_abs_SHAP_burden",
    "Top_eRNA_abs_SHAP_burden",
    "eRNA_signed_SHAP_burden",
    "eRNA_SHAP_fraction"
  )
  burden_metrics <- burden_metrics[vapply(burden_metrics, function(m) any(is.finite(burden_df[[m]])), logical(1))]

  burden_grid <- tidyr::expand_grid(
    metric = burden_metrics,
    pathway = selected_pathways
  )

  cor_results_burden <- purrr::pmap_dfr(burden_grid, function(metric, pathway) {
    res <- safe_cor_test(burden_df[[metric]], score_for_burden[, pathway], method = COR_METHOD)
    tibble::tibble(
      metric = metric,
      pathway = pathway,
      rho = as.numeric(res["rho"]),
      pvalue = as.numeric(res["pvalue"]),
      n_samples = as.numeric(res["n"])
    )
  }) %>%
    dplyr::group_by(.data$metric) %>%
    dplyr::mutate(FDR = p.adjust(.data$pvalue, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pathway_label = short_label(.data$pathway, max_len = 54),
      metric_label = dplyr::recode(
        .data$metric,
        eRNA_abs_SHAP_burden = "All eRNA absolute SHAP burden",
        Top_eRNA_abs_SHAP_burden = "Top eRNA absolute SHAP burden",
        eRNA_signed_SHAP_burden = "All eRNA signed SHAP burden",
        Top_eRNA_signed_SHAP_burden = "Top eRNA signed SHAP burden",
        eRNA_SHAP_fraction = "eRNA fraction of total SHAP burden"
      ),
      minus_log10_p = -log10(pmax(.data$pvalue, .Machine$double.xmin)),
      minus_log10_FDR = -log10(pmax(.data$FDR, .Machine$double.xmin))
    ) %>%
    dplyr::arrange(.data$metric, .data$FDR, dplyr::desc(abs(.data$rho)))

  readr::write_csv(cor_results_burden, file.path(OUTPUT_DIR, "eRNA_SHAP_burden_vs_immune_pathway_correlation_all.csv"))
  readr::write_csv(
    cor_results_burden %>% dplyr::filter(.data$FDR < FDR_CUTOFF),
    file.path(OUTPUT_DIR, "eRNA_SHAP_burden_vs_immune_pathway_correlation_FDR_lt_0.05.csv")
  )

  # Plot B1: dot plot of burden metrics vs pathways.
  p_burden_df <- cor_results_burden %>%
    dplyr::filter(.data$pathway %in% selected_pathways) %>%
    dplyr::mutate(
      pathway_label = factor(.data$pathway_label, levels = rev(short_label(pathway_order, 54))),
      metric_label = factor(
        .data$metric_label,
        levels = c(
          "All eRNA absolute SHAP burden",
          "Top eRNA absolute SHAP burden",
          "All eRNA signed SHAP burden",
          "eRNA fraction of total SHAP burden"
        )
      )
    )

  p_burden <- ggplot(p_burden_df, aes(x = rho, y = pathway_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey58") +
    geom_point(aes(size = pmin(minus_log10_FDR, 12), color = rho), alpha = 0.94) +
    facet_wrap(~ metric_label, ncol = 2) +
    scale_color_gradient2(low = "#2166AC", mid = "grey88", high = "#B2182B", midpoint = 0, name = "Spearman\nrho") +
    scale_size_continuous(name = expression(-log[10](FDR)), range = c(1.8, 6.2), limits = c(0, 12)) +
    labs(
      title = "eRNA-derived SHAP burden is coupled to immune-inflammatory pathway activity",
      subtitle = paste0("Sample-level eRNA SHAP burden vs ssGSEA pathway scores; ", COR_METHOD, " correlation"),
      x = "Spearman correlation with eRNA SHAP burden",
      y = NULL
    ) +
    theme_pub(base_size = 9.5) +
    theme(
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      strip.text = element_text(size = 8.5, face = "bold")
    )

  save_pdf(
    p_burden,
    file.path(OUTPUT_DIR, "Fig2P_eRNA_SHAP_burden_immune_pathway_correlation_dotplot.pdf"),
    width = 11,
    height = max(7.2, 0.24 * length(selected_pathways) + 4.8)
  )

  # Plot B2: representative scatterplots for Top eRNA absolute SHAP burden if available.
  metric_for_scatter <- if ("Top_eRNA_abs_SHAP_burden" %in% burden_metrics && any(is.finite(burden_df$Top_eRNA_abs_SHAP_burden))) {
    "Top_eRNA_abs_SHAP_burden"
  } else {
    "eRNA_abs_SHAP_burden"
  }

  top_burden_paths <- cor_results_burden %>%
    dplyr::filter(.data$metric == metric_for_scatter, .data$rho > 0, .data$pvalue < 0.05) %>%
    dplyr::arrange(.data$FDR, dplyr::desc(.data$rho)) %>%
    dplyr::slice_head(n = TOP_SCATTER_N)

  if (nrow(top_burden_paths) == 0) {
    top_burden_paths <- cor_results_burden %>%
      dplyr::filter(.data$metric == metric_for_scatter, is.finite(.data$rho)) %>%
      dplyr::arrange(dplyr::desc(.data$rho)) %>%
      dplyr::slice_head(n = TOP_SCATTER_N)
  }

  if (nrow(top_burden_paths) > 0) {
    scatter_burden_df <- purrr::map_dfr(seq_len(nrow(top_burden_paths)), function(i) {
      p <- top_burden_paths$pathway[i]
      tibble::tibble(
        sample_id = burden_df$sample_id,
        burden = burden_df[[metric_for_scatter]],
        pathway_score = score_for_burden[, p],
        group = burden_df$group,
        facet_label = paste0(
          short_label(p, 48),
          "\nrho=", sprintf("%.2f", top_burden_paths$rho[i]),
          ", FDR=", scales::scientific(top_burden_paths$FDR[i], digits = 2)
        )
      )
    })

    p_scatter_burden <- ggplot(scatter_burden_df, aes(x = burden, y = pathway_score)) +
      geom_point(aes(fill = group), shape = 21, size = 2.25, color = "white", stroke = 0.25, alpha = 0.82) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "grey20", fill = "grey82") +
      facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
      labs(
        title = "Representative immune pathways associated with eRNA SHAP burden",
        subtitle = dplyr::recode(
          metric_for_scatter,
          Top_eRNA_abs_SHAP_burden = "Top eRNA absolute SHAP burden",
          eRNA_abs_SHAP_burden = "All eRNA absolute SHAP burden"
        ),
        x = dplyr::recode(
          metric_for_scatter,
          Top_eRNA_abs_SHAP_burden = "Top eRNA absolute SHAP burden",
          eRNA_abs_SHAP_burden = "All eRNA absolute SHAP burden"
        ),
        y = "ssGSEA pathway score",
        fill = "Group"
      ) +
      theme_pub(base_size = 9) +
      theme(legend.position = "bottom")

    save_pdf(
      p_scatter_burden,
      file.path(OUTPUT_DIR, "Fig2Q_representative_eRNA_SHAP_burden_immune_pathway_scatterplots.pdf"),
      width = 9.5,
      height = 3.15 * ceiling(nrow(top_burden_paths) / 2)
    )
  }
}

# ------------------------------
# 8. Summary
# ------------------------------

summary_lines <- c(
  "============================================================",
  "eRNA-immune pathway linkage analysis summary",
  "============================================================",
  paste0("Analysis time: ", Sys.time()),
  paste0("Samples used for expression-pathway correlations: ", nrow(expr)),
  paste0("eRNA features detected: ", length(erna_features)),
  paste0("Top SHAP eRNAs tested: ", length(top_erna)),
  paste0("Immune-inflammatory pathways displayed/tested: ", length(selected_pathways)),
  paste0("SHAP summary source: ", shap_summary_obj$source),
  "",
  "Analysis A outputs:",
  paste0("  ", file.path(OUTPUT_DIR, "Top_SHAP_eRNA_expression_vs_immune_pathway_correlation_all.csv")),
  paste0("  ", file.path(OUTPUT_DIR, "Top_SHAP_eRNA_immune_association_summary_by_eRNA.csv")),
  paste0("  ", file.path(OUTPUT_DIR, "Fig2M_Top_SHAP_eRNA_expression_immune_pathway_correlation_heatmap.pdf")),
  paste0("  ", file.path(OUTPUT_DIR, "Fig2N_Top_SHAP_eRNA_immune_association_lollipop.pdf")),
  paste0("  ", file.path(OUTPUT_DIR, "Fig2O_representative_Top_SHAP_eRNA_expression_immune_pathway_scatterplots.pdf")),
  "",
  "Analysis B outputs, if sample-level SHAP matrix was found:",
  paste0("  ", file.path(OUTPUT_DIR, "Sample_level_eRNA_SHAP_burden_scores.csv")),
  paste0("  ", file.path(OUTPUT_DIR, "eRNA_SHAP_burden_vs_immune_pathway_correlation_all.csv")),
  paste0("  ", file.path(OUTPUT_DIR, "Fig2P_eRNA_SHAP_burden_immune_pathway_correlation_dotplot.pdf")),
  paste0("  ", file.path(OUTPUT_DIR, "Fig2Q_representative_eRNA_SHAP_burden_immune_pathway_scatterplots.pdf")),
  "============================================================"
)

writeLines(summary_lines, file.path(OUTPUT_DIR, "eRNA_immune_pathway_linkage_analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

message("Done.")
