# ========== 0. Environment Configuration ==========
library(shapviz)       # Tree SHAP for XGBoost/LightGBM/CatBoost/RF
library(glmnet)       # Ridge coefficient extraction
library(xgboost)      # XGBoost model handling
library(lightgbm)     # LightGBM model handling
library(catboost)     # CatBoost model handling (if needed)
library(ranger)       # Random Forest model handling
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)    # Plot assembly
library(tibble)
library(forcats)

set.seed(2024)

# ========== 0.1 User Configuration Area ==========

# Meta-learner (Ridge) paths
RIDGE_MODEL_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/ridge/model_package/ridge_regression_model.rds"
RIDGE_COEF_PATH  <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/ridge/ridge_regression_feature_importance_full.csv"

# Base learner model paths (comment unused models)
BASE_MODEL_PATHS <- list(
  lightgbm   = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/LightGBM/model_package/lightgbm_model.txt",
  catboost   = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/model_package/catboost_model.cbm",
  rf         = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/RandomForest/model_package/randomforest_model.rds",
  lasso      = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/lasso_regression/model_package/lasso_regression_model.rds",
  Ridge_Regression= "/home/yjliu/mmProj/data_process/Human/Machine_Learning/ridge_regression/model_package/ridge_regression_model.rds",
  Elastic_Net       = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/elastic_net/model_package/elastic_net_model.rds",
  adaboost   = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/model_package/adaboost_model.rds",
  svm        = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/SVM_polynomial/model_package/svm_model.rds",
  Neural_Network       = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/model_package/neuralnetwork_model.rds"
)

# Training expression matrix for SHAP calculation
TRAIN_MATRIX_PATH <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data.csv"

# Base learner prediction files for meta-level analysis
BASE_PRED_TRAIN_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/10model_train_predictions.csv"
BASE_PRED_TEST_PATH  <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/predictions/10model_test_predictions.csv"

# Feature to RNA type annotation file
FEATURE_ANNOTATION_PATH <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/SHAP_analysis/特征类型.csv"

# Output directory
OUTPUT_DIR <- "/home/yjliu/mmProj/data_process/Human/Ensemble_Model/Stacking/SHAP_analysis/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========== 0.2 Utility Functions ==========

# Infer RNA biotype from transcript ID naming rules
infer_rna_class <- function(transcript_ids) {
  classes <- rep("mRNA", length(transcript_ids))
  names(classes) <- transcript_ids
  
  # eRNA: EnsXXXXX / ENSG with enhancer tag / FANTOM5 eRNA
  eRNA_pattern <- "^Ens[0-9]|enhancer|eRNA|FANTOM"
  classes[grepl(eRNA_pattern, transcript_ids, ignore.case = TRUE)] <- "eRNA"
  
  # miRNA: hsa-miR / hsa-let
  miRNA_pattern <- "^hsa-miR|^hsa-let"
  classes[grepl(miRNA_pattern, transcript_ids, ignore.case = TRUE)] <- "miRNA"
  
  # lncRNA: Known lncRNA prefixes / NONCODE IDs
  lncRNA_pattern <- "^NON|^LINC|^MALAT|^HOTAIR|^XIST|^NEAT|^GAS5|^H19|^MEG3|^PVT1|lnc|antisense"
  classes[grepl(lncRNA_pattern, transcript_ids, ignore.case = TRUE)] <- "lncRNA"
  
  return(classes)
}

cat("✅ Environment configuration completed\n")
# =============================================================================
# SECTION 1: Load Ridge Meta-learner & Basic Data
# =============================================================================

cat("\n========== 1. Load Ridge Meta-learner ==========\n")

# Load Ridge model
ridge_model <- readRDS(RIDGE_MODEL_PATH)
cat("Ridge model class:", class(ridge_model), "\n")

# Load Ridge coefficients (base learner weights)
ridge_coef_raw <- read.csv(RIDGE_COEF_PATH)
cat("Ridge coefficients:\n")
print(head(ridge_coef_raw, 12))

# Extract base learner names and coefficients
ridge_coef_df <- ridge_coef_raw %>%
  filter(Feature != "(Intercept)") %>%
  rename(
    base_learner = Feature,
    coef = Coefficient,
    abs_coef = Abs_Coefficient,
    odds_ratio = Odds_Ratio
  ) %>%
  mutate(
    weight = abs_coef / sum(abs_coef),
    weight_pct = round(weight * 100, 1)
  )

cat("\nBase learner weights from Ridge meta-learner:\n")
ridge_coef_df %>%
  select(base_learner, coef, weight_pct) %>%
  arrange(desc(weight_pct)) %>%
  print()

# =============================================================================
# SECTION 2: Meta-level Analysis — Base Learner Contribution
# =============================================================================

cat("\n========== 2. Meta-level Analysis ==========\n")

# Load base learner predictions
train_preds <- read.csv(BASE_PRED_TRAIN_PATH)

base_learner_names <- setdiff(names(train_preds), c("X","group"))
cat("Detected base learners (", length(base_learner_names), "):\n")
cat(paste(base_learner_names, collapse = ", "), "\n")

stopifnot(all(base_learner_names %in% ridge_coef_df$base_learner))

# Calculate meta-level SHAP for linear Ridge: SHAP = coef * (pred - mean_pred)
meta_shap_list <- list()
for (bl in base_learner_names) {
  coef_bl <- ridge_coef_df$coef[ridge_coef_df$base_learner == bl]
  pred_bl <- train_preds[[bl]]
  mean_pred_bl <- mean(pred_bl)
  meta_shap_list[[bl]] <- coef_bl * (pred_bl - mean_pred_bl)
}
meta_shap_matrix <- as.data.frame(meta_shap_list)

# Meta contribution statistics
meta_contribution <- data.frame(
  base_learner = base_learner_names,
  mean_abs_shap = sapply(meta_shap_list, function(x) mean(abs(x))),
  sd_abs_shap   = sapply(meta_shap_list, function(x) sd(abs(x)))
) %>%
  left_join(ridge_coef_df[, c("base_learner", "weight_pct")], by = "base_learner") %>%
  mutate(
    meta_contribution = mean_abs_shap * weight_pct / 100
  ) %>%
  arrange(desc(meta_contribution))

cat("\nMeta-level base learner contribution ranking:\n")
print(meta_contribution)

# Save meta-level results
write.csv(meta_contribution,
          file.path(OUTPUT_DIR, "meta_level_base_learner_contribution.csv"),
          row.names = FALSE)

# Meta contribution bar plot
p_meta <- ggplot(meta_contribution, 
                 aes(x = reorder(base_learner, meta_contribution), 
                     y = meta_contribution,
                     fill = meta_contribution)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#6baed6", high = "#08519c") +
  labs(x = "Base Learner", y = "Meta-level Contribution",
       title = "Base Learner Contribution in MMFinder Meta-learner (Ridge)",
       subtitle = "SHAP-based contribution weighted by Ridge coefficients") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "Fig_meta_base_learner_contribution.pdf"),
       p_meta, width = 8, height = 5)
cat("✅ Meta-level analysis completed\n")
# =============================================================================
cat("\n========== 3. Load Data & Base Learner Models ==========\n")

# 3.1 Load training expression matrix
if (grepl("\\.csv$|.csv", TRAIN_MATRIX_PATH)) {
   train_matrix_raw <- read.csv(TRAIN_MATRIX_PATH, row.names = 1)
} else if (grepl("\\.rds$|.rds", TRAIN_MATRIX_PATH)) {
   train_matrix_raw <- readRDS(TRAIN_MATRIX_PATH)
} else {
   stop("Unsupported file format: ", TRAIN_MATRIX_PATH)
}
train_matrix_raw <- train_matrix_raw[, -ncol(train_matrix_raw)]
train_matrix <- as.matrix(train_matrix_raw)
cat("Training set dimension:", dim(train_matrix), "(samples × features)\n")

all_features <- colnames(train_matrix)
n_features_total <- length(all_features)

# 3.2 Feature to RNA biotype mapping
if (!is.null(FEATURE_ANNOTATION_PATH) && file.exists(FEATURE_ANNOTATION_PATH)) {
   feature_anno <- read.csv(FEATURE_ANNOTATION_PATH)
   id_col <- intersect(c("transcript_id", "gene", "feature", "id"), names(feature_anno))[1]
   class_col <- intersect(c("rna_class", "type", "class", "RNA_class"), names(feature_anno))[1]
   
   if (!is.na(id_col) && !is.na(class_col)) {
      feature_classes <- setNames(as.character(feature_anno[[class_col]]),
                                  as.character(feature_anno[[id_col]]))
      cat("✅ RNA biotype mapping loaded from annotation file\n")
   } else {
      cat("⚠️ Annotation file column names mismatch, fallback to name-based inference\n")
      feature_classes <- infer_rna_class(all_features)
   }
} else {
   cat("⚠️ Feature annotation file not found\n")
   feature_classes <- infer_rna_class(all_features)
}
cat("RNA biotype distribution:\n")
print(table(feature_classes))

# =========================================================================
# 3.3 General Model Loader
# =========================================================================

load_base_learner <- function(path, model_name) {
   
   if (!file.exists(path)) {
      cat("  ❌ File not found:", path, "\n")
      return(NULL)
   }
   
   ext <- tolower(tools::file_ext(path))
   cat("  Load", model_name, "(", ext, "):", path, "\n")
   
   model <- NULL
   model_type <- "unknown"
   predict_fn <- NULL
   
   # XGBoost native format
   if (ext %in% c("model") || (!grepl("\\.", basename(path)) && model_name == "xgboost")) {
      model <- xgboost::xgb.load(path)
      model_type <- "tree"
      predict_fn <- function(m, X) predict(m, newdata = as.matrix(X))
   }
   
   # LightGBM text format
   else if (ext %in% c("txt") || model_name == "lightgbm") {
      model <- lightgbm::lgb.load(path)
      model_type <- "tree"
      predict_fn <- function(m, X) predict(m, data = as.matrix(X))
   }
   
   # CatBoost format
   else if (ext %in% c("cbm") || model_name == "catboost") {
      model <- catboost::catboost.load_model(path)
      model_type <- "tree"
      predict_fn <- function(m, X) {
         catboost::catboost.predict(m, catboost::catboost.load_pool(as.data.frame(X)),
                                    prediction_type = "Probability")
      }
   }
   
   # RDS file objects
   else {
      model <- readRDS(path)
      obj_class <- class(model)[1]
      cat("    RDS object type:", obj_class, "\n")
      
      # Extract finalModel from caret train object
      if (any(c("train", "train.formula") %in% class(model))) {
         cat("    → Detected caret train object, extract $finalModel\n")
         model <- model$finalModel
         obj_class <- class(model)[1]
         cat("    Underlying model:", obj_class, "\n")
      }
      
      # Dispatch by model class
      if ("xgb.Booster" %in% class(model)) {
         model_type <- "tree"
         predict_fn <- function(m, X) predict(m, newdata = as.matrix(X))
      } else if ("lgb.Booster" %in% class(model)) {
         model_type <- "tree"
         predict_fn <- function(m, X) predict(m, data = as.matrix(X))
      } else if ("ranger" %in% class(model)) {
         model_type <- "tree"
         predict_fn <- function(m, X) predict(m, data = as.data.frame(X))$predictions
      } else if ("randomForest" %in% class(model)) {
         model_type <- "tree"
         predict_fn <- function(m, X) predict(m, newdata = as.data.frame(X), type = "prob")[, 2]
      } else if (any(c("cv.glmnet", "glmnet", "lognet", "elnet") %in% class(model))) {
         # LASSO / Ridge / Elastic Net
         model_type <- "linear"
         predict_fn <- function(m, X) as.numeric(predict(m, newx = as.matrix(X), type = "response"))
         
      } else if ("svm" %in% class(model) || "svm.formula" %in% class(model)) {
         model_type <- "svm"
         predict_fn <- function(m, X) {
            pred <- predict(m, newdata = as.data.frame(X), probability = TRUE)
            prob_mat <- attr(pred, "probabilities")
            if (!is.null(prob_mat)) {
               tumor_col <- grep("tumor", colnames(prob_mat), ignore.case = TRUE, value = TRUE)
               if (length(tumor_col) > 0) prob_mat[, tumor_col[1]] else prob_mat[, 2]
            } else as.numeric(pred == "tumor")
         }
         
      } else if ("ksvm" %in% class(model)) {
         model_type <- "svm"
         predict_fn <- function(m, X) {
            as.numeric(kernlab::predict(m, as.data.frame(X), type = "probabilities")[, "tumor"])
         }
         
      } else if ("nnet" %in% class(model) || "nnet.formula" %in% class(model)) {
         model_type <- "nnet"
         predict_fn <- function(m, X) as.numeric(predict(m, newdata = as.data.frame(X), type = "raw"))
         
      } else if ("nn" %in% class(model)) {
         model_type <- "nnet"
         predict_fn <- function(m, X) as.numeric(neuralnet::compute(m, as.data.frame(X))$net.result)
         
      } else if (any(c("boosting", "bagging") %in% class(model))) {
         model_type <- "adaboost"
         predict_fn <- function(m, X) {
            pred <- predict(m, newdata = as.data.frame(X))
            if (is.list(pred)) pred$prob[, "tumor"] else pred
         }
         
      } else if (any(c("gbm", "mboost", "glmboost") %in% class(model))) {
         model_type <- "tree"
         predict_fn <- function(m, X) as.numeric(predict(m, newdata = as.data.frame(X), type = "response"))
         
      } else {
         cat("    ⚠️ Unrecognized model type:", obj_class, "→ Use generic predict()\n")
         model_type <- "unknown"
         predict_fn <- function(m, X) {
            p <- tryCatch(
               predict(m, newdata = as.data.frame(X), type = "prob"),
               error = function(e) predict(m, newdata = as.data.frame(X))
            )
            if (is.data.frame(p) || is.matrix(p)) p[, 2] else as.numeric(p)
         }
      }
   }
   
   return(list(model = model, model_type = model_type, predict_fn = predict_fn))
}

# =========================================================================
# 3.4 Batch Model Loading
# =========================================================================

available_models  <- list()
model_types       <- list()
model_predictors  <- list()

for (bl_name in names(BASE_MODEL_PATHS)) {
   path <- BASE_MODEL_PATHS[[bl_name]]
   if (is.null(path) || path == "") {
      cat("⚠️", bl_name, "path is empty, skip\n")
      next
   }
   result <- load_base_learner(path, bl_name)
   if (!is.null(result)) {
      available_models[[bl_name]]  <- result$model
      model_types[[bl_name]]       <- result$model_type
      model_predictors[[bl_name]]  <- result$predict_fn
      cat("  ✅ Load successfully, type:", result$model_type, "\n")
   }
}

n_available <- length(available_models)
cat("\nAvailable base learners:", n_available, "/ 10\n")
cat("Tree-based models:", paste(names(model_types)[model_types == "tree"], collapse = ", "), "\n")

# =========================================================================
# 3.5 Name Matching
# =========================================================================
ridge_bl_names <- ridge_coef_df$base_learner

match_bl_name <- function(model_name, ridge_names) {
   # Exact match
   if (model_name %in% ridge_names) return(model_name)
   
   # Manual name mapping
   manual_map <- list(
      "nnet"       = "Neural_Network",
      "ridge_base" = "Ridge_Regression",
      "enet"       = "Elastic_Net",
      "lasso"      = "Lasso_Regression"
   )
   if (model_name %in% names(manual_map)) {
      candidate <- manual_map[[model_name]]
      if (candidate %in% ridge_names) return(candidate)
   }
   
   # Fuzzy match (remove underscore & case insensitive)
   model_clean <- tolower(gsub("[_-]", "", model_name))
   for (rn in ridge_names) {
      rn_clean <- tolower(gsub("[_-]", "", rn))
      if (model_clean == rn_clean ||
          grepl(model_clean, rn_clean, fixed = TRUE) ||
          grepl(rn_clean, model_clean, fixed = TRUE)) {
         return(rn)
      }
   }
   return(NA)
}

name_map <- list()
for (mname in names(available_models)) {
   matched <- match_bl_name(mname, ridge_bl_names)
   name_map[[mname]] <- matched
   if (!is.na(matched)) {
      cat("  ", mname, "→", matched, "\n")
   } else {
      cat("  ⚠️", mname, "→ No match found in Ridge coefficients\n")
   }
}

cat("✅ Section 3 completed\n")

# =============================================================================
# SECTION 4: Calculate Feature-level SHAP for Each Base Learner
# =============================================================================

cat("\n========== 4. Calculate Base Learner SHAP ==========\n")

# Store SHAP matrices and feature importances
all_shap_matrices <- list()
all_shap_importances <- list()

# Background samples for Kernel SHAP
n_bg <- min(150, nrow(train_matrix))
bg_indices <- sample(1:nrow(train_matrix), n_bg)
bg_matrix <- train_matrix[bg_indices, ]

library(kernelshap)

for (mname in names(available_models)) {
  model <- available_models[[mname]]
  mtype <- model_types[[mname]]
  cat("\n---", mname, "(", mtype, ") ---\n")
  
  shp_obj <- NULL
  
  # 4a. Tree SHAP for XGBoost / LightGBM / CatBoost / RF
  if (mtype == "tree") {
    tryCatch({
      shp_obj <- shapviz(model, X_pred = train_matrix, X = train_matrix)
      cat("  ✅ Tree SHAP calculated successfully\n")
    }, error = function(e) {
      cat("  ❌ Tree SHAP failed:", e$message, "\n")
      cat("  → Fallback to Kernel SHAP approximation (", n_bg, " background samples)\n")
      tryCatch({
        predict_fn <- function(X) predict(model, as.matrix(X))
        shp_obj <- kernelshap(model, X = train_matrix[1:min(500, nrow(train_matrix)), ],
                              bg_X = bg_matrix)
      }, error = function(e2) {
        cat("  ❌ Kernel SHAP also failed:", e2$message, "→ Skip this model\n")
      })
    })
  }
  
  # 4b. Linear model SHAP for LASSO / Ridge / Elastic Net
  else if (mtype == "linear") {
    tryCatch({
      coefs <- as.matrix(coef(model))
      coefs <- coefs[-1, ]
      
      common_feat <- intersect(names(coefs), all_features)
      
      X_centered <- scale(train_matrix[, common_feat], center = TRUE, scale = FALSE)
      shap_matrix <- sweep(X_centered, 2, coefs[common_feat], "*")
      
      shp_obj <- list(
        S = shap_matrix,
        X = train_matrix[, common_feat],
        baseline = attr(X_centered, "scaled:center") %*% coefs[common_feat]
      )
      class(shp_obj) <- "shapviz"
      cat("  ✅ Linear SHAP calculated successfully (", length(common_feat), "features)\n")
    }, error = function(e) {
      cat("  ❌ Linear SHAP failed:", e$message, "\n")
    })
  }
  
  # 4c. Kernel SHAP for SVM / NNet / AdaBoost
  else if (mtype %in% c("svm", "nnet", "adaboost")) {
    tryCatch({
      predict_fn <- function(X) {
        if (mtype == "svm") {
          as.numeric(attr(predict(model, newdata = as.data.frame(X), probability = TRUE), 
                          "probabilities")[, "tumor"])
        } else {
          as.numeric(predict(model, newdata = as.data.frame(X), type = "prob")[, "tumor"])
        }
      }
      shp_obj <- kernelshap(predict_fn, 
                            X = train_matrix[1:min(500, nrow(train_matrix)), ],
                            bg_X = bg_matrix)
      cat("  ✅ Kernel SHAP calculated successfully\n")
    }, error = function(e) {
      cat("  ❌ Kernel SHAP failed:", e$message, "→ Skip this model\n")
    })
  }
  
  # 4d. Extract SHAP importance
  if (!is.null(shp_obj)) {
    all_shap_matrices[[mname]] <- shp_obj$S
    
    importance <- colMeans(abs(shp_obj$S))
    all_shap_importances[[mname]] <- data.frame(
      feature = names(importance),
      mean_abs_shap = as.numeric(importance),
      base_learner = mname,
      rna_class = feature_classes[names(importance)]
    ) %>% arrange(desc(mean_abs_shap))
    
    cat("  Top-5 features (", mname, "):\n")
    print(head(all_shap_importances[[mname]], 5))
  }
}

cat("\nSuccessfully calculated SHAP for:", length(all_shap_matrices), "/", n_available, "base learners\n")
# =============================================================================
# SECTION 5: Weighted SHAP Aggregation → Consensus SHAP
# =============================================================================

cat("\n========== 5. Weighted SHAP Aggregation ==========\n")

# 5.1 Get weights from Ridge meta-learner
available_bl_names <- names(all_shap_matrices)

ridge_weights <- sapply(available_bl_names, function(mname) {
  ridge_name <- name_map[[mname]]
  if (is.na(ridge_name)) return(NA)
  ridge_coef_df$weight[ridge_coef_df$base_learner == ridge_name]
})

# Re-normalize weights
ridge_weights <- ridge_weights / sum(ridge_weights, na.rm = TRUE)
cat("SHAP aggregation weights:\n")
for (i in seq_along(ridge_weights)) {
  cat(sprintf("  %-15s → %-20s  weight = %.3f\n", 
              names(ridge_weights)[i], name_map[[names(ridge_weights)[i]]], 
              ridge_weights[i]))
}

# 5.2 Unify feature space
all_shap_features <- unique(unlist(lapply(all_shap_matrices, colnames)))
cat("\nTotal features covered across all models:", length(all_shap_features), "\n")

# 5.3 Weighted aggregation
consensus_shap_list <- list()

for (feat in all_shap_features) {
  shap_values <- c()
  weights <- c()
  
  for (mname in available_bl_names) {
    if (feat %in% colnames(all_shap_matrices[[mname]])) {
      shap_values <- c(shap_values, 
                       all_shap_importances[[mname]]$mean_abs_shap[
                         all_shap_importances[[mname]]$feature == feat])
      weights <- c(weights, ridge_weights[mname])
    }
  }
  
  if (length(shap_values) > 0) {
    consensus_shap_list[[feat]] <- list(
      feature = feat,
      consensus_shap = weighted.mean(shap_values, weights),
      n_models = length(shap_values),
      shap_values = shap_values,
      rna_class = feature_classes[feat]
    )
  }
}

# Convert to data frame
consensus_shap_df <- do.call(rbind, lapply(consensus_shap_list, function(x) {
   len <- length(x$feature)
   if (length(x$consensus_shap) != len || length(x$n_models) != len) {
      stop("Inconsistent vector lengths in list elements")
   }
   
   data.frame(
      feature         = x$feature,
      consensus_shap  = x$consensus_shap,
      n_models        = x$n_models,
      rna_class       = ifelse(is.na(x$rna_class), "mRNA", x$rna_class),
      row.names       = NULL,
      stringsAsFactors = FALSE
   )
})) %>% arrange(desc(consensus_shap))

rownames(consensus_shap_df) <- NULL

cat("\n====== Consensus SHAP Top-20 Features ======\n")
print(head(consensus_shap_df, 20))

# 5.4 Save full ranking
write.csv(consensus_shap_df,
          file.path(OUTPUT_DIR, "consensus_shap_feature_importance.csv"),
          row.names = FALSE)

# 5.5 Model coverage analysis
coverage_dist <- table(consensus_shap_df$n_models)
cat("\nFeature coverage distribution across models:\n")
print(coverage_dist)

# Robust features present in >= 60% of models
min_models <- max(2, ceiling(length(available_bl_names) * 0.6))
robust_features <- consensus_shap_df %>% filter(n_models >= min_models)
cat("\nRobust features covered by ≥", min_models, "models:", nrow(robust_features), "\n")

cat("✅ SHAP aggregation completed\n")

# =============================================================================
# SECTION 6: Per RNA Biotype Aggregation & Statistical Test
# =============================================================================

cat("\n========== 6. Per-class SHAP Aggregation & Statistical Test ==========\n")

class_summary <- consensus_shap_df %>%
   group_by(rna_class) %>%
   summarise(
      n_features        = n(),
      total_shap        = sum(consensus_shap),
      mean_shap         = mean(consensus_shap),
      median_shap       = median(consensus_shap),
      sd_shap           = sd(consensus_shap),
      se_shap           = sd_shap / sqrt(n()),
      max_shap          = max(consensus_shap),
      top5_features     = paste(head(feature, 5), collapse = "; ")
   ) %>%
   mutate(
      shap_share_pct    = total_shap / sum(total_shap) * 100,
      feature_share_pct = n_features / sum(n_features) * 100,
      efficiency_ratio  = shap_share_pct / feature_share_pct
   ) %>%
   arrange(desc(mean_shap))

cat("\n====== RNA Biotype SHAP Contribution Summary ======\n")
print(class_summary)
write.csv(class_summary, file.path(OUTPUT_DIR, "per_class_shap_summary.csv"), row.names = FALSE)

# Wilcoxon rank-sum test: eRNA SHAP > other biotypes
eRNA_shap  <- consensus_shap_df$consensus_shap[consensus_shap_df$rna_class == "eRNA"]
mRNA_shap  <- consensus_shap_df$consensus_shap[consensus_shap_df$rna_class == "mRNA"]
lncRNA_shap <- consensus_shap_df$consensus_shap[consensus_shap_df$rna_class == "lncRNA"]
miRNA_shap <- consensus_shap_df$consensus_shap[consensus_shap_df$rna_class == "miRNA"]

mw_e_m  <- wilcox.test(eRNA_shap, mRNA_shap,  alternative = "greater")
mw_e_l  <- wilcox.test(eRNA_shap, lncRNA_shap, alternative = "greater")
mw_e_mi <- wilcox.test(eRNA_shap, miRNA_shap, alternative = "greater")

test_results <- data.frame(
   comparison   = c("eRNA > mRNA", "eRNA > lncRNA", "eRNA > miRNA"),
   W_statistic  = c(mw_e_m$statistic, mw_e_l$statistic, mw_e_mi$statistic),
   P_value      = c(mw_e_m$p.value,  mw_e_l$p.value,  mw_e_mi$p.value)
) %>%
   mutate(
      FDR = p.adjust(P_value, method = "BH"),
      significance = case_when(
         FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "ns"
      )
   )

cat("\nStatistical test results:\n")
print(test_results)
write.csv(test_results, file.path(OUTPUT_DIR, "statistical_tests.csv"), row.names = FALSE)

# =============================================================================
# SECTION 7: Top-N Enrichment Analysis
# =============================================================================

cat("\n========== 7. Top-N Enrichment Analysis ==========\n")

N_values <- c(5, 10, 20, 30, 50, 100, 200, 500)
enrichment_results <- list()
total_n <- nrow(consensus_shap_df)

for (N in N_values) {
   topN <- consensus_shap_df[1:min(N, total_n), ]
   for (cls in names(table(consensus_shap_df$rna_class))) {
      n_topN      <- sum(topN$rna_class == cls)
      n_total_cls <- sum(consensus_shap_df$rna_class == cls)
      expected    <- n_total_cls / total_n * N
      
      mat <- matrix(c(n_topN, N - n_topN,
                      n_total_cls - n_topN, total_n - N - (n_total_cls - n_topN)),
                    nrow = 2)
      ft <- fisher.test(mat)
      
      enrichment_results[[length(enrichment_results) + 1]] <- data.frame(
         N = N, rna_class = cls,
         n_in_topN = n_topN, expected = round(expected, 2),
         enrichment_ratio = (n_topN / N) / (n_total_cls / total_n),
         odds_ratio = ft$estimate, p_value = ft$p.value
      )
   }
}

enrichment_df <- bind_rows(enrichment_results) %>%
   mutate(
      FDR = p.adjust(p_value, method = "BH"),
      sig_label = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "")
   )

cat("\n====== Top-N Enrichment (eRNA focus) ======\n")
print(enrichment_df %>% filter(rna_class == "eRNA"))
write.csv(enrichment_df, file.path(OUTPUT_DIR, "topN_enrichment.csv"), row.names = FALSE)
# =============================================================================
# SECTION 8: MMER1 Specific Analysis
# =============================================================================

cat("\n========== 8. MMER1 Specific Analysis ==========\n")

mmer1_match <- grep("MMER1|Ens223489|223489", consensus_shap_df$feature,
                    ignore.case = TRUE, value = FALSE)

if (length(mmer1_match) > 0) {
   mmer1_rank <- mmer1_match[1]
   mmer1_row  <- consensus_shap_df[mmer1_rank, ]
   mmer1_feature_name <- mmer1_row$feature
   
   cat(sprintf("MMER1 global rank: %d/%d (top %.1f%%)\n",
               mmer1_rank, nrow(consensus_shap_df),
               mmer1_rank / nrow(consensus_shap_df) * 100))
   cat(sprintf("  Consensus SHAP = %.6f | Covered models = %d | Biotype = %s\n",
               mmer1_row$consensus_shap, mmer1_row$n_models, mmer1_row$rna_class))
   
   # Rank within eRNA
   eRNA_only <- consensus_shap_df %>% filter(rna_class == "eRNA")
   mmer1_eRNA_rank <- which(eRNA_only$feature == mmer1_feature_name)
   cat(sprintf("  Rank among %d eRNAs: %d/%d\n",
               nrow(eRNA_only), mmer1_eRNA_rank, nrow(eRNA_only)))
   
   # Co-expression analysis between MMER1 and other eRNAs
   if (mmer1_feature_name %in% all_features) {
      mmer1_expr <- train_matrix[, mmer1_feature_name]
      eRNA_features <- intersect(names(feature_classes)[feature_classes == "eRNA"], all_features)
      
      eRNA_cor <- data.frame(
         feature = eRNA_features,
         pearson_r = sapply(eRNA_features, function(f) cor(mmer1_expr, train_matrix[, f]))
      ) %>% mutate(abs_r = abs(pearson_r)) %>% arrange(desc(abs_r))
      
      coexpressed <- eRNA_cor %>% filter(abs_r > 0.5)
      cat(sprintf("eRNAs with |Pearson r| > 0.5 with MMER1: %d\n", nrow(coexpressed)))
      if (nrow(coexpressed) > 0) print(head(coexpressed, 10))
      write.csv(coexpressed, file.path(OUTPUT_DIR, "MMER1_coexpressed_eRNAs.csv"), row.names = FALSE)
   }
} else {
   cat("⚠️ MMER1 not found in consensus SHAP results. Please check feature names in training matrix\n")
   cat("   Preview of top 20 feature names:\n")
   print(head(consensus_shap_df$feature, 20))
   mmer1_feature_name <- NULL
}

# Top-15 eRNAs ranked by SHAP
top_eRNAs <- consensus_shap_df %>% filter(rna_class == "eRNA") %>% head(15)
cat("\n====== Top-15 eRNAs by SHAP ======\n")
print(top_eRNAs)
write.csv(top_eRNAs, file.path(OUTPUT_DIR, "top_eRNAs_shap_importance.csv"), row.names = FALSE)

cat("✅ Section 8 completed\n")
# =============================================================================
# SECTION 9: Visualization — Fig D1 / D2 / D5
# =============================================================================

cat("\n========== 9. Generate Visualizations ==========\n")

# --- Fig D1: RNA Biotype Contribution (3-panel) ---

p_share <- ggplot(class_summary, aes(x = "", y = shap_share_pct, fill = rna_class)) +
   geom_bar(stat = "identity", width = 1, color = "white") +
   coord_polar("y", start = 0) +
   geom_text(aes(label = sprintf("%s\n%.1f%%", rna_class, shap_share_pct)),
             position = position_stack(vjust = 0.5), size = 3.2) +
   scale_fill_manual(values = c("eRNA" = "#d95f02", "mRNA" = "#1b9e77",
                                "lncRNA" = "#7570b3", "miRNA" = "#e7298a")) +
   labs(title = "Total SHAP Contribution") +
   theme_void() + theme(legend.position = "none")

p_mean <- ggplot(class_summary, aes(x = reorder(rna_class, mean_shap),
                                    y = mean_shap, fill = rna_class)) +
   geom_bar(stat = "identity", width = 0.6) +
   geom_errorbar(aes(ymin = mean_shap - se_shap, ymax = mean_shap + se_shap), width = 0.15) +
   geom_text(aes(label = sprintf("%.4f", mean_shap)), vjust = -0.5, size = 3) +
   scale_fill_manual(values = c("eRNA" = "#d95f02", "mRNA" = "#1b9e77",
                                "lncRNA" = "#7570b3", "miRNA" = "#e7298a")) +
   labs(x = "", y = "Mean |SHAP| per Transcript",
        subtitle = sprintf("eRNA vs mRNA: P = %.2e", mw_e_m$p.value)) +
   theme_minimal() + theme(legend.position = "none")

p_eff <- ggplot(class_summary, aes(x = reorder(rna_class, efficiency_ratio),
                                   y = efficiency_ratio, fill = rna_class)) +
   geom_bar(stat = "identity", width = 0.6) +
   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
   geom_text(aes(label = sprintf("%.2f", efficiency_ratio)), vjust = -0.5, size = 3.5) +
   scale_fill_manual(values = c("eRNA" = "#d95f02", "mRNA" = "#1b9e77",
                                "lncRNA" = "#7570b3", "miRNA" = "#e7298a")) +
   labs(x = "", y = "Efficiency (SHAP% / Feature%)",
        subtitle = "> 1 = contributes more than expected") +
   theme_minimal() + theme(legend.position = "none")

p_d1 <- (p_share | p_mean | p_eff) +
   plot_annotation(title = "Fig D1: RNA Class SHAP Contribution to MMFinder",
                   subtitle = "eRNAs show disproportionately high per-transcript importance",
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUTPUT_DIR, "Fig_D1_RNA_class_SHAP_contribution.pdf"),
       p_d1, width = 15, height = 5.5)

# --- Fig D2: Top-N Enrichment Bar Plot ---

library(scales)

p_d2_data <- enrichment_df %>%
   mutate(N_label = factor(paste0("Top-", N), levels = paste0("Top-", rev(N_values))))

p_d2_labels <- p_d2_data %>%
   filter(rna_class == "eRNA", FDR < 0.05) %>%
   mutate(y_pos = 1.02)

p_d2 <- ggplot(p_d2_data, 
               aes(x = N_label, y = n_in_topN / N, fill = rna_class)) +
   geom_bar(stat = "identity", position = "fill", width = 0.7) +
   geom_text(data = p_d2_labels,
             aes(x = N_label, y = y_pos, label = sig_label),
             size = 5, color = "#d95f02", inherit.aes = FALSE) +
   scale_fill_manual(values = c("eRNA" = "#d95f02", "mRNA" = "#1b9e77",
                                "lncRNA" = "#7570b3", "miRNA" = "#e7298a")) +
   scale_y_continuous(labels = percent_format()) +
   coord_flip() +
   labs(x = "", y = "Proportion", fill = "RNA Class",
        title = "Fig D2: RNA Class Composition in Top-N SHAP Features",
        subtitle = "*** FDR < 0.001 (Fisher) for eRNA enrichment") +
   theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "Fig_D2_topN_enrichment.pdf"), p_d2, width = 9, height = 6)

# --- Fig D5: Cross Base Learner SHAP Consistency Heatmap ---

if (length(all_shap_importances) >= 2) {
   top50_by_model <- lapply(all_shap_importances, function(df) head(df$feature, 50))
   mnames <- names(top50_by_model)
   jaccard_mat <- matrix(NA, nrow = length(mnames), ncol = length(mnames),
                         dimnames = list(mnames, mnames))
   
   for (i in seq_along(mnames))
      for (j in seq_along(mnames))
         jaccard_mat[i, j] <- length(intersect(top50_by_model[[i]], top50_by_model[[j]])) /
      length(union(top50_by_model[[i]], top50_by_model[[j]]))
   
   jaccard_long <- as.data.frame(as.table(jaccard_mat))
   names(jaccard_long) <- c("Model1", "Model2", "Jaccard")
   
   p_d5 <- ggplot(jaccard_long, aes(Model1, Model2, fill = Jaccard)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3) +
      scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1)) +
      labs(title = "Fig D5: Cross-Model SHAP Consistency (Top-50 Jaccard)") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   ggsave(file.path(OUTPUT_DIR, "Fig_D5_cross_model_consistency.pdf"),
          p_d5, width = 8, height = 6)
}

cat("✅ Visualizations saved\n")
# =============================================================================
# SECTION 10: Waterfall Plot & MMER1 SHAP Dependence Plot
# =============================================================================

cat("\n========== 10. Waterfall Plot & Dependence Plot ==========\n")

# Select representative samples
train_preds <- read.csv(BASE_PRED_TRAIN_PATH)
X_train_meta <- train_preds[, setdiff(names(train_preds), c("X", "group")), drop = FALSE]
final_scores_train <- as.numeric(predict(ridge_model, newx = as.matrix(X_train_meta), type = "response"))

has_train_group <- "group" %in% names(train_preds)
if (has_train_group) {
   train_group <- train_preds$group
   mm_idx_train      <- which(train_group == "tumor"  & final_scores_train > 0.5)
   healthy_idx_train <- which(train_group == "health" & final_scores_train < 0.5)
   sample_mm      <- if (length(mm_idx_train) > 0) 
      mm_idx_train[which.max(final_scores_train[mm_idx_train])] else which.max(final_scores_train)
   sample_healthy <- if (length(healthy_idx_train) > 0) 
      healthy_idx_train[which.min(final_scores_train[healthy_idx_train])] else which.min(final_scores_train)
} else {
   sample_mm      <- which.max(final_scores_train)
   sample_healthy <- which.min(final_scores_train)
}

cat(sprintf("Representative MM sample: Training row %d (score=%.4f)\n", sample_mm, final_scores_train[sample_mm]))
cat(sprintf("Representative Healthy sample: Training row %d (score=%.4f)\n", sample_healthy, final_scores_train[sample_healthy]))

# =============================================================================
# 10.2 Waterfall Plot Function
# =============================================================================

draw_waterfall <- function(shap_values, feature_names, sample_id, 
                           label, max_display = 15, color_pos = "#d95f02", 
                           color_neg = "#1b9e77", outpath = NULL) {
   
   s <- sort(shap_values, decreasing = TRUE)
   n <- min(max_display, length(s))
   s <- s[1:n]
   feat_names <- names(s)
   
   pdf(outpath, width = 10, height = 0.5 + 0.35 * n)
   
   cols <- ifelse(s > 0, color_pos, color_neg)
   
   par(mar = c(4, 12, 3, 2))
   bp <- barplot(rev(s), horiz = TRUE, las = 1, col = rev(cols),
                 border = NA, xlab = "SHAP value",
                 main = paste0("SHAP Waterfall - ", label, " (row ", sample_id, ")"))
   
   text(x = ifelse(rev(s) > 0, rev(s) + max(abs(s)) * 0.02, rev(s) - max(abs(s)) * 0.1),
        y = bp, labels = rev(feat_names), pos = ifelse(rev(s) > 0, 4, 2),
        cex = 0.7, col = "gray30")
   
   text(x = rev(s), y = bp, 
        labels = sprintf("%.4f", rev(s)),
        pos = ifelse(rev(s) > 0, 2, 4), cex = 0.6)
   
   abline(v = 0, lty = 2, col = "gray50")
   
   dev.off()
   cat("✅", label, "waterfall plot saved to:", outpath, "\n")
}

# =============================================================================
# 10.3 Dependence Plot Function
# =============================================================================

draw_dependence <- function(feature_values, shap_values, feat_name,
                            color_var_values = NULL, color_var_name = NULL,
                            outpath = NULL) {
   
   pdf(outpath, width = 7, height = 5)
   par(mar = c(4, 4, 3, 1))
   
   if (!is.null(color_var_values)) {
      n_colors <- 100
      pal <- colorRampPalette(c("#2171b5", "#f7f7f7", "#d95f02"))(n_colors)
      color_idx <- as.numeric(cut(color_var_values, breaks = n_colors))
      color_idx[is.na(color_idx)] <- 1
      
      plot(feature_values, shap_values,
           pch = 16, cex = 0.5, col = adjustcolor(pal[color_idx], alpha.f = 0.6),
           xlab = paste0(feat_name, " expression"),
           ylab = "SHAP value",
           main = paste0("SHAP Dependence: ", feat_name))
      
      legend("topright", legend = c("high", "mid", "low"),
             fill = pal[c(n_colors, 50, 1)], 
             title = color_var_name, cex = 0.7, border = NA)
   } else {
      plot(feature_values, shap_values,
           pch = 16, cex = 0.5, col = adjustcolor("#d95f02", alpha.f = 0.3),
           xlab = paste0(feat_name, " expression"),
           ylab = "SHAP value",
           main = paste0("SHAP Dependence: ", feat_name))
   }
   
   lo <- lowess(feature_values[is.finite(feature_values)], 
                shap_values[is.finite(feature_values)])
   lines(lo, col = "red", lwd = 2)
   abline(h = 0, lty = 2, col = "gray50")
   
   dev.off()
   cat("✅", feat_name, "dependence plot saved to:", outpath, "\n")
}

# =============================================================================
# 10.4 Execute Plotting
# =============================================================================

tree_available <- names(all_shap_matrices)[model_types[names(all_shap_matrices)] == "tree"]
pref_order <- c(intersect(tree_available, names(all_shap_matrices)), 
                setdiff(names(all_shap_matrices), tree_available))

if (length(pref_order) > 0) {
   best_model <- pref_order[1]
   S_use <- all_shap_matrices[[best_model]]
   cat("Use SHAP data from:", best_model, "(", ncol(S_use), "features )\n")
   
   if (sample_mm <= nrow(S_use)) {
      draw_waterfall(
         shap_values  = S_use[sample_mm, ],
         feature_names = colnames(S_use),
         sample_id    = sample_mm,
         label        = "MM",
         outpath      = file.path(OUTPUT_DIR, "Fig_D4_waterfall_MM.pdf")
      )
   } else {
      cat("⚠️ MM sample index out of SHAP matrix range, skip\n")
   }
   
   if (sample_healthy <= nrow(S_use)) {
      draw_waterfall(
         shap_values  = S_use[sample_healthy, ],
         feature_names = colnames(S_use),
         sample_id    = sample_healthy,
         label        = "Healthy",
         outpath      = file.path(OUTPUT_DIR, "Fig_D4_waterfall_healthy.pdf")
      )
   } else {
      cat("⚠️ Healthy sample index out of SHAP matrix range, skip\n")
   }
   
   # MMER1 dependence plot
   if (!is.null(mmer1_feature_name) && mmer1_feature_name %in% colnames(S_use)) {
      draw_dependence(
         feature_values = train_matrix[1:nrow(S_use), mmer1_feature_name],
         shap_values    = S_use[, mmer1_feature_name],
         feat_name      = mmer1_feature_name,
         outpath        = file.path(OUTPUT_DIR, "Fig_D3_MMER1_dependence.pdf")
      )
      
      ikzf1_match <- grep("IKZF1|IKAROS", all_features,