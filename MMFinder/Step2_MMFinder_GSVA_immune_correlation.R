# =============================================================================
# Figure 2 add-on analysis
# MMFinder score–GSVA/ssGSEA immune-inflammatory pathway correlation
#
# Purpose:
#   Test whether MMFinder diagnostic probability captures sample-level
#   immune/inflammatory transcriptional remodeling.
#
# Inputs reused from the MMFinder SHAP analysis workflow:
#   1) Ridge meta-learner model
#   2) Base learner training predictions
#   3) Training expression matrix
#   4) Feature annotation file, if available
#
# Main outputs:
#   - CSV: pathway scores and correlation statistics
#   - PDF: publication-style dot plot of pathway correlations
#   - PDF: selected scatter plots for the strongest correlations
#   - PDF: MM vs healthy pathway-score comparison for significant pathways
# =============================================================================

# ------------------------------
# 0. User configuration
# ------------------------------

set.seed(2024)

# Paths copied from your SHAP aggregation script. Modify only if your project paths changed.
RIDGE_MODEL_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/ridge/model_package/ridge_regression_model.rds"
BASE_PRED_TRAIN_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/10model_train_predictions.csv"
TRAIN_MATRIX_PATH <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data.csv"
FEATURE_ANNOTATION_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/SHAP_analysis/特征类型.csv"

OUTPUT_DIR <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/SHAP_analysis/GSVA_immune_MMFinder_correlation"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ssGSEA/GSVA settings
MIN_GENESET_SIZE <- 5
MAX_GENESET_SIZE <- 500
COR_METHOD <- "spearman"
PLOT_FDR_CUTOFF <- 0.10       # for highlighting; all terms are saved regardless of FDR
TOP_SCATTER_N <- 6            # number of pathways shown as individual scatter plots

# Optional: If TRUE, try to install missing packages automatically.
# For a server/HPC environment, FALSE is safer; install packages manually if needed.
INSTALL_MISSING_PACKAGES <- FALSE

# ------------------------------
# 1. Package loading
# ------------------------------

cran_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "stringr", "readr", "patchwork", "forcats", "scales", "glmnet")
bioc_pkgs <- c("GSVA", "GSEABase", "msigdbr", "org.Hs.eg.db", "AnnotationDbi")

install_if_missing <- function(pkgs, bioc = FALSE) {
   missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
   if (length(missing) == 0) return(invisible(TRUE))
   
   if (!INSTALL_MISSING_PACKAGES) {
      stop(
         "Missing required package(s): ", paste(missing, collapse = ", "),
         "\nInstall them first, or set INSTALL_MISSING_PACKAGES <- TRUE."
      )
   }
   
   if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(missing, ask = FALSE, update = FALSE)
   } else {
      install.packages(missing, repos = "https://cloud.r-project.org")
   }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
   library(dplyr)
   library(tidyr)
   library(tibble)
   library(ggplot2)
   library(stringr)
   library(readr)
   library(patchwork)
   library(forcats)
   library(scales)
   library(glmnet)
   library(GSVA)
   library(msigdbr)
   library(org.Hs.eg.db)
   library(AnnotationDbi)
})

message("Packages loaded.")

# ------------------------------
# 2. Helper functions
# ------------------------------

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
   classes <- rep("mRNA", length(transcript_ids))
   names(classes) <- transcript_ids
   classes[grepl("^Ens[0-9]|enhancer|eRNA|FANTOM", transcript_ids, ignore.case = TRUE)] <- "eRNA"
   classes[grepl("^hsa-miR|^hsa-let", transcript_ids, ignore.case = TRUE)] <- "miRNA"
   classes[grepl("^NON|^LINC|^MALAT|^HOTAIR|^XIST|^NEAT|^GAS5|^H19|^MEG3|^PVT1|lnc|antisense", transcript_ids, ignore.case = TRUE)] <- "lncRNA"
   classes
}

# Robustly load feature annotation if available.
load_feature_classes <- function(feature_ids, anno_path) {
   if (!is.null(anno_path) && file.exists(anno_path)) {
      anno <- read.csv(anno_path, check.names = FALSE, stringsAsFactors = FALSE)
      id_col <- intersect(c("transcript_id", "feature", "gene", "id", "Feature", "ID"), colnames(anno))[1]
      class_col <- intersect(c("rna_class", "RNA_class", "type", "class", "Type", "Class"), colnames(anno))[1]
      
      if (!is.na(id_col) && !is.na(class_col)) {
         fc <- stats::setNames(standardize_rna_class(anno[[class_col]]), as.character(anno[[id_col]]))
         missing_ids <- setdiff(feature_ids, names(fc))
         if (length(missing_ids) > 0) {
            inferred <- infer_rna_class(missing_ids)
            fc <- c(fc, inferred)
         }
         return(fc[feature_ids])
      }
   }
   
   infer_rna_class(feature_ids)
}

# Convert feature IDs to gene symbols. Handles common cases:
#   - already gene symbols
#   - Ensembl IDs, with or without version suffix
#   - ENTREZ IDs
# Unmapped IDs are returned as original IDs, then filtered by overlap with pathway genes.
map_features_to_symbols <- function(feature_ids) {
   ids <- as.character(feature_ids)
   ids_no_version <- sub("\\.\\d+$", "", ids)
   
   # Start with direct symbols.
   out <- ids
   
   # Ensembl mapping.
   ensembl_like <- grepl("^ENSG", ids_no_version, ignore.case = TRUE)
   if (any(ensembl_like)) {
      ens_map <- AnnotationDbi::select(
         org.Hs.eg.db,
         keys = unique(ids_no_version[ensembl_like]),
         keytype = "ENSEMBL",
         columns = c("SYMBOL")
      ) %>%
         dplyr::filter(!is.na(.data$SYMBOL)) %>%
         dplyr::distinct(.data$ENSEMBL, .keep_all = TRUE)
      
      idx <- match(ids_no_version, ens_map$ENSEMBL)
      out[!is.na(idx)] <- ens_map$SYMBOL[idx[!is.na(idx)]]
   }
   
   # Entrez mapping for purely numeric IDs.
   numeric_like <- grepl("^[0-9]+$", ids)
   if (any(numeric_like)) {
      entrez_map <- AnnotationDbi::select(
         org.Hs.eg.db,
         keys = unique(ids[numeric_like]),
         keytype = "ENTREZID",
         columns = c("SYMBOL")
      ) %>%
         dplyr::filter(!is.na(.data$SYMBOL)) %>%
         dplyr::distinct(.data$ENTREZID, .keep_all = TRUE)
      
      idx <- match(ids, entrez_map$ENTREZID)
      out[!is.na(idx)] <- entrez_map$SYMBOL[idx[!is.na(idx)]]
   }
   
   toupper(out)
}

# Collapse duplicated gene symbols by mean expression.
collapse_expression_by_symbol <- function(expr_samples_by_features, symbols) {
   stopifnot(ncol(expr_samples_by_features) == length(symbols))
   
   keep <- !is.na(symbols) & symbols != "" & is.finite(colSums(expr_samples_by_features, na.rm = TRUE))
   expr <- expr_samples_by_features[, keep, drop = FALSE]
   symbols <- symbols[keep]
   
   expr_t <- t(expr) # genes/features x samples
   rownames(expr_t) <- symbols
   
   # Collapse duplicated symbols using mean.
   if (anyDuplicated(rownames(expr_t)) > 0) {
      split_idx <- split(seq_len(nrow(expr_t)), rownames(expr_t))
      collapsed <- do.call(rbind, lapply(split_idx, function(ii) {
         if (length(ii) == 1) expr_t[ii, , drop = FALSE] else matrix(colMeans(expr_t[ii, , drop = FALSE], na.rm = TRUE), nrow = 1)
      }))
      rownames(collapsed) <- names(split_idx)
      colnames(collapsed) <- colnames(expr_t)
      expr_t <- collapsed
   }
   
   storage.mode(expr_t) <- "double"
   expr_t
}

# msigdbr API changed in recent versions. This helper supports both old and new argument names.
msigdbr_safe <- function(species = "Homo sapiens", category = NULL, subcategory = NULL) {
   f <- msigdbr::msigdbr
   args <- names(formals(f))
   
   call_args <- list(species = species)
   if ("category" %in% args) {
      call_args$category <- category
      if (!is.null(subcategory) && "subcategory" %in% args) call_args$subcategory <- subcategory
   } else {
      call_args$collection <- category
      if (!is.null(subcategory) && "subcollection" %in% args) call_args$subcollection <- subcategory
   }
   
   do.call(f, call_args)
}

# Run ssGSEA in a way that works across old and new GSVA versions.
run_ssgsea <- function(expr_gene_by_sample, gene_sets) {
   gene_sets <- gene_sets[lengths(gene_sets) >= MIN_GENESET_SIZE]
   gene_sets <- gene_sets[lengths(gene_sets) <= MAX_GENESET_SIZE]
   gene_sets <- lapply(gene_sets, function(g) intersect(unique(toupper(g)), rownames(expr_gene_by_sample)))
   gene_sets <- gene_sets[lengths(gene_sets) >= MIN_GENESET_SIZE]
   
   if (length(gene_sets) == 0) stop("No gene set has enough genes after intersecting with expression matrix.")
   
   # New GSVA versions use parameter objects.
   if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      param <- GSVA::ssgseaParam(
         exprData = expr_gene_by_sample,
         geneSets = gene_sets,
         minSize = MIN_GENESET_SIZE,
         maxSize = MAX_GENESET_SIZE,
         normalize = TRUE
      )
      scores <- GSVA::gsva(param, verbose = FALSE)
   } else {
      scores <- GSVA::gsva(
         expr = expr_gene_by_sample,
         gset.idx.list = gene_sets,
         method = "ssgsea",
         kcdf = "Gaussian",
         min.sz = MIN_GENESET_SIZE,
         max.sz = MAX_GENESET_SIZE,
         abs.ranking = FALSE,
         verbose = FALSE
      )
   }
   
   as.matrix(scores)
}

clean_pathway_label <- function(x) {
   x %>%
      stringr::str_replace("^HALLMARK_", "") %>%
      stringr::str_replace("^REACTOME_", "") %>%
      stringr::str_replace("^GOBP_", "") %>%
      stringr::str_replace("^GO_", "") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence()
}

# ------------------------------
# 3. Load model predictions and compute final MMFinder score
# ------------------------------

message("Loading ridge meta-learner and base-learner predictions...")

if (!file.exists(RIDGE_MODEL_PATH)) stop("Ridge model not found: ", RIDGE_MODEL_PATH)
if (!file.exists(BASE_PRED_TRAIN_PATH)) stop("Base prediction file not found: ", BASE_PRED_TRAIN_PATH)
if (!file.exists(TRAIN_MATRIX_PATH)) stop("Training expression matrix not found: ", TRAIN_MATRIX_PATH)

ridge_model <- readRDS(RIDGE_MODEL_PATH)
train_preds <- read.csv(BASE_PRED_TRAIN_PATH, check.names = FALSE, stringsAsFactors = FALSE)

# Clean column names immediately. Some exported CSV files contain leading/trailing spaces.
colnames(train_preds) <- trimws(colnames(train_preds))

sample_id_col <- intersect(c("X", "sample", "sample_id", "Sample", "SampleID", "ID"), colnames(train_preds))[1]
group_col <- intersect(c("group", "Group", "condition", "Condition", "label", "Label", "class", "Class"), colnames(train_preds))[1]

# ---- Robustly identify base-learner prediction columns ----
# Do NOT subset train_preds by an unchecked character vector. In some R sessions,
# NA/empty/old column names can cause: undefined columns selected.
exclude_cols <- unique(stats::na.omit(c(sample_id_col, group_col)))
candidate_cols <- setdiff(colnames(train_preds), exclude_cols)
candidate_cols <- candidate_cols[!is.na(candidate_cols) & candidate_cols != "" & candidate_cols %in% colnames(train_preds)]

# Convert numeric-like prediction columns stored as character/factor back to numeric.
is_numeric_like <- function(x) {
   if (is.numeric(x)) return(TRUE)
   x_chr <- as.character(x)
   x_chr <- x_chr[!is.na(x_chr) & x_chr != ""]
   if (length(x_chr) == 0) return(FALSE)
   all(grepl("^-?[0-9.]+([eE][-+]?[0-9]+)?$", x_chr))
}

numeric_like_cols <- candidate_cols[vapply(candidate_cols, function(cc) is_numeric_like(train_preds[[cc]]), logical(1))]
if (length(numeric_like_cols) == 0) {
   stop(
      "No numeric base-learner prediction columns found in BASE_PRED_TRAIN_PATH.
",
      "Detected columns were: ", paste(colnames(train_preds), collapse = ", "), "
",
"Excluded columns were: ", paste(exclude_cols, collapse = ", ")
   )
}

for (cc in numeric_like_cols) {
   train_preds[[cc]] <- as.numeric(as.character(train_preds[[cc]]))
}

# If the ridge/glmnet model exposes coefficient names, use them to order meta-features.
get_model_feature_names <- function(model) {
   rn <- NULL
   cf <- tryCatch(stats::coef(model), error = function(e) NULL)
   if (!is.null(cf)) rn <- rownames(as.matrix(cf))
   rn <- setdiff(rn, "(Intercept)")
   rn[!is.na(rn) & rn != ""]
}

model_feature_names <- get_model_feature_names(ridge_model)
clean_name <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", x))

if (length(model_feature_names) > 0) {
   # First try exact matching; then try cleaned-name matching.
   exact_cols <- intersect(model_feature_names, numeric_like_cols)
   if (length(exact_cols) == length(model_feature_names)) {
      meta_feature_cols <- model_feature_names
   } else {
      actual_clean <- setNames(numeric_like_cols, clean_name(numeric_like_cols))
      matched_actual <- unname(actual_clean[clean_name(model_feature_names)])
      keep <- !is.na(matched_actual) & matched_actual %in% numeric_like_cols
      if (sum(keep) >= 2) {
         meta_feature_cols <- matched_actual[keep]
         # Rename columns in the matrix later to the model's expected names.
         model_feature_names <- model_feature_names[keep]
      } else {
         warning(
            "Could not match ridge coefficient names to prediction-table columns. ",
            "Using all numeric-like prediction columns in file order."
         )
         meta_feature_cols <- numeric_like_cols
         model_feature_names <- meta_feature_cols
      }
   }
} else {
   meta_feature_cols <- numeric_like_cols
   model_feature_names <- meta_feature_cols
}

# Final safety check before subsetting.
meta_feature_cols <- meta_feature_cols[!is.na(meta_feature_cols) & meta_feature_cols %in% colnames(train_preds)]
if (length(meta_feature_cols) == 0) {
   stop("After matching/cleaning, no valid base-learner prediction columns remain.")
}

message("Base-learner prediction columns used for MMFinder score: ", paste(meta_feature_cols, collapse = ", "))

X_meta <- as.matrix(train_preds[, meta_feature_cols, drop = FALSE])
storage.mode(X_meta) <- "double"
colnames(X_meta) <- model_feature_names[seq_len(ncol(X_meta))]

# Predict final MMFinder probability.
# Important: the saved ridge meta-learner is often a glmnet/lognet object.
# Generic predict() may fail if glmnet is not attached or if the object is a raw glmnet object.
# This wrapper explicitly dispatches to glmnet-compatible prediction first.
predict_mmfinder_score <- function(model, X_meta) {
   if (any(c("cv.glmnet", "glmnet", "lognet", "elnet") %in% class(model))) {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
         stop("Package 'glmnet' is required to predict from the saved ridge/lognet model.")
      }
      
      # For cv.glmnet, use lambda.min by default. For raw glmnet/lognet objects,
      # do not pass s unless the object provides lambda; predict.glmnet will return
      # one column per lambda if s is omitted, so we use the last lambda column by default.
      pred <- NULL
      
      if ("cv.glmnet" %in% class(model)) {
         pred <- glmnet::predict.cv.glmnet(model, newx = X_meta, s = "lambda.min", type = "response")
      } else {
         # Prefer a scalar lambda if available. In most glmnet fits, the last lambda
         # is the least penalized model and usually closest to the exported model behavior.
         if (!is.null(model$lambda) && length(model$lambda) > 0) {
            lambda_use <- model$lambda[length(model$lambda)]
            pred <- glmnet::predict.glmnet(model, newx = X_meta, s = lambda_use, type = "response")
         } else {
            pred <- glmnet::predict.glmnet(model, newx = X_meta, type = "response")
         }
      }
      
      pred <- as.matrix(pred)
      
      # Binary lognet may return n x 1, or n x lambda_count. Use the last column
      # if multiple lambdas remain.
      if (ncol(pred) > 1) pred <- pred[, ncol(pred), drop = FALSE]
      return(as.numeric(pred[, 1]))
   }
   
   # Non-glmnet fallback.
   tryCatch({
      as.numeric(stats::predict(model, newx = X_meta, type = "response"))
   }, error = function(e1) {
      tryCatch({
         as.numeric(stats::predict(model, newdata = as.data.frame(X_meta), type = "response"))
      }, error = function(e2) {
         stop("Could not predict final MMFinder score from ridge_model. First error: ", e1$message,
              "\nSecond error: ", e2$message)
      })
   })
}

final_score <- predict_mmfinder_score(ridge_model, X_meta)

metadata <- tibble::tibble(
   row_index = seq_len(nrow(train_preds)),
   sample_id = if (!is.na(sample_id_col)) as.character(train_preds[[sample_id_col]]) else paste0("Sample_", seq_len(nrow(train_preds))),
   group = if (!is.na(group_col)) standardize_group(train_preds[[group_col]]) else NA_character_,
   MMFinder_score = final_score
)

readr::write_csv(metadata, file.path(OUTPUT_DIR, "MMFinder_training_scores.csv"))
message("MMFinder final scores computed for ", nrow(metadata), " samples.")

# ------------------------------
# 4. Load training expression matrix and align samples
# ------------------------------

message("Loading training expression matrix...")
train_matrix_raw <- read.csv(TRAIN_MATRIX_PATH, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Drop non-numeric columns such as group/label if present.
numeric_cols <- vapply(train_matrix_raw, is.numeric, logical(1))
if (sum(!numeric_cols) > 0) {
   message("Dropping non-numeric columns from training matrix: ", paste(colnames(train_matrix_raw)[!numeric_cols], collapse = ", "))
}
train_matrix <- as.matrix(train_matrix_raw[, numeric_cols, drop = FALSE])
storage.mode(train_matrix) <- "double"

# Align train_matrix rows to metadata.
if (!is.na(sample_id_col) && all(metadata$sample_id %in% rownames(train_matrix))) {
   train_matrix <- train_matrix[metadata$sample_id, , drop = FALSE]
   message("Samples aligned by sample_id column: ", sample_id_col)
} else if (nrow(train_matrix) == nrow(metadata)) {
   # Use row order if sample IDs do not match.
   rownames(train_matrix) <- metadata$sample_id
   message("Samples aligned by row order. Please verify if sample IDs are unavailable or inconsistent.")
} else {
   stop("Cannot align training expression matrix with base prediction table. Matrix rows = ", nrow(train_matrix),
        "; prediction rows = ", nrow(metadata), ".")
}

message("Training matrix dimensions after cleaning: ", nrow(train_matrix), " samples x ", ncol(train_matrix), " features.")

# ------------------------------
# 5. Build mRNA gene-symbol expression matrix
# ------------------------------

feature_ids <- colnames(train_matrix)
feature_classes <- load_feature_classes(feature_ids, FEATURE_ANNOTATION_PATH)
feature_classes <- standardize_rna_class(feature_classes)

# Use only mRNAs for pathway activity scoring.
mrna_features <- feature_ids[feature_classes == "mRNA"]
if (length(mrna_features) < 50) {
   warning("Only ", length(mrna_features), " mRNA features detected. Check FEATURE_ANNOTATION_PATH and feature names.")
}

mrna_matrix <- train_matrix[, mrna_features, drop = FALSE]
mrna_symbols <- map_features_to_symbols(mrna_features)
expr_gene_by_sample <- collapse_expression_by_symbol(mrna_matrix, mrna_symbols)

message("mRNA expression matrix for ssGSEA: ", nrow(expr_gene_by_sample), " unique gene symbols x ", ncol(expr_gene_by_sample), " samples.")
readr::write_csv(
   tibble::tibble(feature = mrna_features, mapped_symbol = mrna_symbols),
   file.path(OUTPUT_DIR, "mRNA_feature_to_symbol_mapping.csv")
)

# ------------------------------
# 6. Build predefined immune/inflammatory gene sets
# ------------------------------

message("Loading immune/inflammatory gene sets from MSigDB via msigdbr...")

hallmark <- msigdbr_safe(species = "Homo sapiens", category = "H")
reactome <- msigdbr_safe(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
go_bp <- msigdbr_safe(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

# Normalize column names across msigdbr versions.
get_gs_col <- function(df) {
   intersect(c("gs_name", "geneset", "gene_set_name"), colnames(df))[1]
}
get_symbol_col <- function(df) {
   intersect(c("gene_symbol", "human_gene_symbol", "symbol"), colnames(df))[1]
}

h_gs <- get_gs_col(hallmark); h_sym <- get_symbol_col(hallmark)
r_gs <- get_gs_col(reactome); r_sym <- get_symbol_col(reactome)
g_gs <- get_gs_col(go_bp); g_sym <- get_symbol_col(go_bp)

# Hallmark immune/inflammatory programs. Curated and intentionally small.
hallmark_keep <- c(
   "HALLMARK_INFLAMMATORY_RESPONSE",
   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
   "HALLMARK_COMPLEMENT",
   "HALLMARK_ALLOGRAFT_REJECTION",
   "HALLMARK_IL2_STAT5_SIGNALING"
)

hallmark_sub <- hallmark %>%
   dplyr::filter(.data[[h_gs]] %in% hallmark_keep) %>%
   dplyr::transmute(pathway = .data[[h_gs]], gene = toupper(.data[[h_sym]]))

# Reactome immune/inflammatory pathways selected by keyword.
reactome_pattern <- paste(
   c(
      "CYTOKINE", "INTERLEUKIN", "INTERFERON", "TNF", "NFKB", "NF_KB",
      "ANTIGEN", "COMPLEMENT", "TOLL", "TLR", "INNATE_IMMUNE", "ADAPTIVE_IMMUNE",
      "CHEMOKINE", "JAK_STAT", "IL_6", "NEUTROPHIL", "MACROPHAGE", "LYMPHOCYTE"
   ),
   collapse = "|"
)

reactome_sub <- reactome %>%
   dplyr::filter(stringr::str_detect(.data[[r_gs]], reactome_pattern)) %>%
   dplyr::transmute(pathway = .data[[r_gs]], gene = toupper(.data[[r_sym]]))

# GO BP immune/inflammatory pathways selected by keyword, then restrict to reasonably interpretable terms.
go_pattern <- paste(
   c(
      "IMMUNE_RESPONSE", "INFLAMMATORY_RESPONSE", "LEUKOCYTE_ACTIVATION", "LYMPHOCYTE_ACTIVATION",
      "T_CELL_ACTIVATION", "B_CELL_ACTIVATION", "MYELOID_LEUKOCYTE", "CYTOKINE",
      "CHEMOKINE", "INTERFERON", "ANTIGEN_PROCESSING", "ANTIGEN_PRESENTATION",
      "COMPLEMENT", "NF_KAPPAB", "NFKB", "JAK_STAT", "INTERLEUKIN"
   ),
   collapse = "|"
)

go_sub <- go_bp %>%
   dplyr::filter(stringr::str_detect(.data[[g_gs]], go_pattern)) %>%
   dplyr::transmute(pathway = .data[[g_gs]], gene = toupper(.data[[g_sym]]))

# Combine, remove overly broad/very small sets after intersection.
immune_gene_sets_df <- dplyr::bind_rows(hallmark_sub, reactome_sub, go_sub) %>%
   dplyr::distinct(.data$pathway, .data$gene) %>%
   dplyr::filter(!is.na(.data$gene), .data$gene != "")

immune_gene_sets <- split(immune_gene_sets_df$gene, immune_gene_sets_df$pathway)
immune_gene_sets <- lapply(immune_gene_sets, unique)

# Intersect gene sets with expression genes and filter by size.
gene_set_size_df <- tibble::tibble(
   pathway = names(immune_gene_sets),
   original_size = lengths(immune_gene_sets),
   detected_size = lengths(lapply(immune_gene_sets, intersect, rownames(expr_gene_by_sample)))
) %>%
   dplyr::filter(.data$detected_size >= MIN_GENESET_SIZE, .data$detected_size <= MAX_GENESET_SIZE) %>%
   dplyr::arrange(dplyr::desc(.data$detected_size))

immune_gene_sets <- immune_gene_sets[gene_set_size_df$pathway]
immune_gene_sets <- lapply(immune_gene_sets, function(g) intersect(unique(toupper(g)), rownames(expr_gene_by_sample)))

readr::write_csv(gene_set_size_df, file.path(OUTPUT_DIR, "Immune_inflammatory_gene_set_detected_sizes.csv"))
saveRDS(immune_gene_sets, file.path(OUTPUT_DIR, "Immune_inflammatory_gene_sets_used.rds"))

message("Immune/inflammatory gene sets retained for ssGSEA: ", length(immune_gene_sets))

# ------------------------------
# 7. Compute ssGSEA pathway activity scores
# ------------------------------

message("Running ssGSEA/GSVA pathway scoring...")
ssgsea_scores <- run_ssgsea(expr_gene_by_sample, immune_gene_sets)

# rows = pathways, columns = samples
readr::write_csv(
   as.data.frame(ssgsea_scores) %>% tibble::rownames_to_column("pathway"),
   file.path(OUTPUT_DIR, "Immune_inflammatory_ssGSEA_scores_pathway_by_sample.csv")
)

score_long <- as.data.frame(t(ssgsea_scores)) %>%
   tibble::rownames_to_column("sample_id") %>%
   tidyr::pivot_longer(-.data$sample_id, names_to = "pathway", values_to = "pathway_score") %>%
   dplyr::left_join(metadata, by = "sample_id")

readr::write_csv(score_long, file.path(OUTPUT_DIR, "Immune_inflammatory_ssGSEA_scores_long.csv"))

# ------------------------------
# 8. Correlate MMFinder probability with pathway activity
# ------------------------------

message("Computing MMFinder score–pathway activity correlations...")

cor_results <- score_long %>%
   dplyr::group_by(.data$pathway) %>%
   dplyr::summarise(
      n_samples = sum(is.finite(.data$MMFinder_score) & is.finite(.data$pathway_score)),
      rho = suppressWarnings(stats::cor(.data$MMFinder_score, .data$pathway_score, method = COR_METHOD, use = "complete.obs")),
      pvalue = suppressWarnings(stats::cor.test(.data$MMFinder_score, .data$pathway_score, method = COR_METHOD, exact = FALSE)$p.value),
      .groups = "drop"
   ) %>%
   dplyr::mutate(
      FDR = p.adjust(.data$pvalue, method = "BH"),
      abs_rho = abs(.data$rho),
      direction = ifelse(.data$rho >= 0, "Positive", "Negative"),
      pathway_label = clean_pathway_label(.data$pathway),
      minus_log10_FDR = -log10(pmax(.data$FDR, .Machine$double.xmin)),
      minus_log10_pvalue = -log10(pmax(.data$pvalue, .Machine$double.xmin))
   ) %>%
   dplyr::arrange(.data$FDR, dplyr::desc(.data$abs_rho))

readr::write_csv(cor_results, file.path(OUTPUT_DIR, "MMFinder_score_immune_pathway_correlation_all.csv"))
readr::write_csv(
   cor_results %>% dplyr::filter(.data$pvalue < 0.05),
   file.path(OUTPUT_DIR, "MMFinder_score_immune_pathway_correlation_pvalue_lt_0.05.csv")
)
readr::write_csv(
   cor_results %>% dplyr::filter(.data$FDR < 0.10),
   file.path(OUTPUT_DIR, "MMFinder_score_immune_pathway_correlation_FDR_lt_0.10.csv")
)

message("Top correlations:")
print(utils::head(cor_results, 20))

# ------------------------------
# 9. Optional MM vs Healthy pathway-score comparison
# ------------------------------

if (!all(is.na(metadata$group)) && length(unique(na.omit(metadata$group))) >= 2) {
   message("Computing MM vs Healthy pathway score comparisons...")
   
   group_results <- score_long %>%
      dplyr::filter(.data$group %in% c("MM", "Healthy")) %>%
      dplyr::group_by(.data$pathway) %>%
      dplyr::summarise(
         n_MM = sum(.data$group == "MM"),
         n_Healthy = sum(.data$group == "Healthy"),
         median_MM = median(.data$pathway_score[.data$group == "MM"], na.rm = TRUE),
         median_Healthy = median(.data$pathway_score[.data$group == "Healthy"], na.rm = TRUE),
         delta_median = median_MM - median_Healthy,
         pvalue = suppressWarnings(stats::wilcox.test(pathway_score ~ group, data = dplyr::cur_data(), exact = FALSE)$p.value),
         .groups = "drop"
      ) %>%
      dplyr::mutate(
         FDR = p.adjust(.data$pvalue, method = "BH"),
         pathway_label = clean_pathway_label(.data$pathway)
      ) %>%
      dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$delta_median)))
   
   readr::write_csv(group_results, file.path(OUTPUT_DIR, "Immune_pathway_scores_MM_vs_Healthy_all.csv"))
} else {
   group_results <- NULL
   message("Group labels unavailable or insufficient; skipping MM vs Healthy comparison.")
}

# ------------------------------
# 10. Publication-style plots
# ------------------------------

message("Generating plots...")

theme_pub <- function(base_size = 11) {
   ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
         plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
         plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30", hjust = 0),
         axis.title = ggplot2::element_text(face = "bold"),
         axis.text = ggplot2::element_text(color = "grey15"),
         axis.line = ggplot2::element_line(linewidth = 0.45, color = "grey20"),
         axis.ticks = ggplot2::element_line(linewidth = 0.35, color = "grey20"),
         legend.title = ggplot2::element_text(face = "bold"),
         legend.position = "right",
         strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
         strip.text = ggplot2::element_text(face = "bold")
      )
}

# 10.1 Main dot plot: correlation between MMFinder score and pathway activity.
plot_cor_df <- cor_results %>%
   dplyr::filter(is.finite(.data$rho), is.finite(.data$pvalue)) %>%
   dplyr::arrange(.data$rho) %>%
   dplyr::mutate(
      pathway_label = factor(.data$pathway_label, levels = unique(.data$pathway_label)),
      sig_shape = ifelse(.data$FDR < PLOT_FDR_CUTOFF, "FDR < cutoff", "Nominal/NS")
   )

# If there are too many pathways, show top 30 by FDR/absolute rho in the main figure.
MAX_PATHWAYS_IN_DOTPLOT <- 30
plot_cor_df_main <- plot_cor_df %>%
   dplyr::arrange(.data$FDR, dplyr::desc(.data$abs_rho)) %>%
   dplyr::slice_head(n = MAX_PATHWAYS_IN_DOTPLOT) %>%
   dplyr::arrange(.data$rho) %>%
   dplyr::mutate(pathway_label = factor(.data$pathway_label, levels = unique(.data$pathway_label)))

p_cor <- ggplot(plot_cor_df_main, aes(x = rho, y = pathway_label)) +
   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey55") +
   geom_point(aes(size = minus_log10_pvalue, color = rho), alpha = 0.92) +
   scale_color_gradient2(low = "#2166AC", mid = "grey88", high = "#B2182B", midpoint = 0, name = "Spearman rho") +
   scale_size_continuous(name = expression(-log[10](italic(P))), range = c(2.2, 7.0)) +
   labs(
      title = "MMFinder score is associated with immune-inflammatory pathway activity",
      subtitle = paste0("ssGSEA pathway scores from mRNA expression; ", COR_METHOD, " correlation; top pathways shown"),
      x = "Spearman correlation with MMFinder diagnostic probability",
      y = NULL
   ) +
   theme_pub(base_size = 10) +
   theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25))

ggsave(
   filename = file.path(OUTPUT_DIR, "Fig2J_MMFinder_score_immune_pathway_correlation_dotplot.pdf"),
   plot = p_cor,
   width = 15,
   height = max(6, 0.26 * nrow(plot_cor_df_main) + 1.8),
   device = "pdf"
)


# 10.2 Scatter plots for top positively correlated pathways.
scatter_pathways <- cor_results %>%
   dplyr::filter(rho > 0, pvalue < 0.05) %>%
   dplyr::arrange(.data$FDR, dplyr::desc(.data$rho)) %>%
   dplyr::slice_head(n = TOP_SCATTER_N) %>%
   dplyr::pull(.data$pathway)

# If no nominal positive pathways, use top positive rho terms.
if (length(scatter_pathways) == 0) {
   scatter_pathways <- cor_results %>%
      dplyr::filter(rho > 0) %>%
      dplyr::arrange(dplyr::desc(.data$rho)) %>%
      dplyr::slice_head(n = TOP_SCATTER_N) %>%
      dplyr::pull(.data$pathway)
}

scatter_df <- score_long %>%
   dplyr::filter(.data$pathway %in% scatter_pathways) %>%
   dplyr::left_join(cor_results %>% dplyr::select(.data$pathway, .data$pathway_label, .data$rho, .data$pvalue, .data$FDR), by = "pathway") %>%
   dplyr::mutate(
      facet_label = paste0(pathway_label, "\n", "rho=", sprintf("%.2f", rho), ", P=", scales::scientific(pvalue, digits = 2))
   )

if (nrow(scatter_df) > 0) {
   p_scatter <- ggplot(scatter_df, aes(x = MMFinder_score, y = pathway_score)) +
      geom_point(aes(fill = group), shape = 21, size = 2.3, alpha = 0.78, color = "white", stroke = 0.25) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "grey20", fill = "grey80") +
      facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
      labs(
         title = "Representative immune-inflammatory pathways associated with MMFinder score",
         x = "MMFinder diagnostic probability",
         y = "ssGSEA pathway score",
         fill = "Group"
      ) +
      theme_pub(base_size = 10) +
      theme(legend.position = "bottom")
   
   ggsave(
      filename = file.path(OUTPUT_DIR, "Fig2K_representative_MMFinder_score_pathway_scatterplots.pdf"),
      plot = p_scatter,
      width = 9.5,
      height = 3.1 * ceiling(length(scatter_pathways) / 2)
   )
}

# 10.3 MM vs Healthy pathway-score boxplot for significant positively correlated pathways.
if (!is.null(group_results)) {
   box_pathways <- cor_results %>%
      dplyr::filter(rho > 0, pvalue < 0.05) %>%
      dplyr::arrange(.data$FDR, dplyr::desc(.data$rho)) %>%
      dplyr::slice_head(n = 12) %>%
      dplyr::pull(.data$pathway)
   
   if (length(box_pathways) > 0) {
      box_df <- score_long %>%
         dplyr::filter(.data$pathway %in% box_pathways, .data$group %in% c("Healthy", "MM")) %>%
         dplyr::left_join(cor_results %>% dplyr::select(.data$pathway, .data$pathway_label, .data$rho, .data$pvalue), by = "pathway") %>%
         dplyr::mutate(
            pathway_label = forcats::fct_reorder(.data$pathway_label, .data$rho),
            group = factor(.data$group, levels = c("Healthy", "MM"))
         )
      
      p_box <- ggplot(box_df, aes(x = group, y = pathway_score, fill = group)) +
         geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
         geom_jitter(width = 0.15, size = 0.75, alpha = 0.35, color = "grey25") +
         facet_wrap(~ pathway_label, scales = "free_y", ncol = 3) +
         labs(
            title = "Immune-inflammatory pathway activity in MM and healthy samples",
            subtitle = "Pathways selected from nominally significant positive MMFinder-score correlations",
            x = NULL,
            y = "ssGSEA pathway score"
         ) +
         theme_pub(base_size = 9) +
         theme(legend.position = "none")
      
      ggsave(
         filename = file.path(OUTPUT_DIR, "Fig2L_immune_pathway_activity_MM_vs_Healthy_boxplots.pdf"),
         plot = p_box,
         width = 10,
         height = 7.2
      )
   }
}

# ------------------------------
# 11. Diagnostics and manuscript-ready summary
# ------------------------------

summary_lines <- c(
   "============================================================",
   "MMFinder score–immune/inflammatory pathway correlation summary",
   "============================================================",
   paste0("Analysis time: ", Sys.time()),
   paste0("Samples: ", nrow(metadata)),
   paste0("Unique mRNA gene symbols used for ssGSEA: ", nrow(expr_gene_by_sample)),
   paste0("Immune/inflammatory gene sets scored: ", length(immune_gene_sets)),
   paste0("Correlation method: ", COR_METHOD),
   "",
   "Top 10 pathways by FDR:",
   paste(utils::capture.output(
      cor_results %>%
         dplyr::select(.data$pathway, .data$rho, .data$pvalue, .data$FDR) %>%
         dplyr::slice_head(n = 10)
   ), collapse = "\n"),
   "",
   paste0("Nominal P < 0.05 pathways: ", sum(cor_results$pvalue < 0.05, na.rm = TRUE)),
   paste0("FDR < 0.10 pathways: ", sum(cor_results$FDR < 0.10, na.rm = TRUE)),
   "",
   "Output files:",
   paste0("  ", file.path(OUTPUT_DIR, "MMFinder_score_immune_pathway_correlation_all.csv")),
   paste0("  ", file.path(OUTPUT_DIR, "Immune_inflammatory_ssGSEA_scores_pathway_by_sample.csv")),
   paste0("  ", file.path(OUTPUT_DIR, "Fig2J_MMFinder_score_immune_pathway_correlation_dotplot.pdf")),
   paste0("  ", file.path(OUTPUT_DIR, "Fig2K_representative_MMFinder_score_pathway_scatterplots.pdf")),
   paste0("  ", file.path(OUTPUT_DIR, "Fig2L_immune_pathway_activity_MM_vs_Healthy_boxplots.pdf")),
   "============================================================"
)

writeLines(summary_lines, file.path(OUTPUT_DIR, "MMFinder_GSVA_immune_correlation_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

message("Done.")
