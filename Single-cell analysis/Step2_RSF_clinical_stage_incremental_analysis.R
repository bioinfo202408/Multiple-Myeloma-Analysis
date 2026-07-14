############################################################
## RSF prognostic model vs ISS clinical staging
## Author: generated for RSF clinical incremental analyses
## Purpose:
##   1) Head-to-head C-index comparison: RSF vs ISS and clinical models
##   2) Multivariable Cox regression: whether RSF is independent of ISS/age/etc.
##   3) Incremental analysis: whether RSF refines high/low risk within ISS stage
##
## Input:
##   /home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis/RSF_full_workspace.RData
##
## Expected objects in workspace:
##   risk_train, risk_valid, risk_test
##   Each risk_* table should contain: OS.time, OS, riskscore, risk_group,
##   iss_stage, Age, Gender, treatment_type, if available.
##
## Output:
##   All tables and publication-ready PDF figures are saved under:
##   Clinical_Incremental_Analysis/
##
## Notes:
##   - This script does NOT retrain the RSF model.
##   - Validation and test cohorts are merged for the primary incremental analysis.
##   - All figures are saved using the traditional pdf() device.
############################################################

##############################
## 0. Global settings
##############################
set.seed(123)
options(stringsAsFactors = FALSE)

WORKSPACE_PATH <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis/RSF_full_workspace.RData"
DEFAULT_BASE_DIR <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis"

## Bootstrap iterations for C-index confidence intervals and paired delta tests.
## For final manuscript output, 1000 is recommended. For quick testing, set to 200.
BOOTSTRAP_B <- 1000

## Time points for DCA and optional fixed-time analyses
DCA_TIMES <- c(365, 730, 1095)
DCA_THRESHOLDS <- seq(0.01, 0.80, by = 0.01)

##############################
## 1. Packages
##############################
required_pkgs <- c(
   "survival", "survminer", "ggplot2", "dplyr", "scales"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
   stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running this script."
   )
}

suppressPackageStartupMessages({
   library(survival)
   library(survminer)
   library(ggplot2)
   library(dplyr)
   library(scales)
})

##############################
## 2. Load workspace and create output folders
##############################
if (!file.exists(WORKSPACE_PATH)) {
   stop("Workspace file does not exist: ", WORKSPACE_PATH)
}

load(WORKSPACE_PATH)

if (!exists("base_dir")) base_dir <- DEFAULT_BASE_DIR
out_dir <- file.path(base_dir, "Clinical_Incremental_Analysis")
fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
log_dir <- file.path(out_dir, "Logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Cindex"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Multivariable_Cox"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "ISS_Incremental"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "DCA"), recursive = TRUE, showWarnings = FALSE)

required_objects <- c("risk_train", "risk_valid", "risk_test")
missing_objects <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
   stop(
      "The following required objects were not found in the workspace: ",
      paste(missing_objects, collapse = ", ")
   )
}

##############################
## 3. Helper functions
##############################
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_filename <- function(x) {
   x <- as.character(x)
   x <- gsub("[[:space:]]+", "_", x)
   x <- gsub("[/\\\\|:*?\"<>]", "_", x)
   x <- gsub("[^A-Za-z0-9_.-]", "_", x)
   x
}

format_p <- function(p) {
   ifelse(
      is.na(p), "NA",
      ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
   )
}

smart_numeric <- function(x) {
   x0 <- as.character(x)
   x0 <- trimws(x0)
   x0[x0 %in% c("", "NA", "N/A", "NaN", "NULL", "null", "Unknown", "unknown")] <- NA
   out <- suppressWarnings(as.numeric(x0))
   out
}

is_numeric_like <- function(x, min_prop = 0.70) {
   x0 <- as.character(x)
   x0 <- trimws(x0)
   x0[x0 %in% c("", "NA", "N/A", "NaN", "NULL", "null", "Unknown", "unknown")] <- NA
   nonmiss <- !is.na(x0)
   if (sum(nonmiss) == 0) return(FALSE)
   val <- suppressWarnings(as.numeric(x0[nonmiss]))
   mean(!is.na(val)) >= min_prop
}

extract_iss_ord <- function(x) {
   z <- toupper(trimws(as.character(x)))
   z <- gsub("[_-]", " ", z)
   z[z %in% c("", "NA", "N/A", "UNKNOWN", "NULL")] <- NA

   out <- rep(NA_real_, length(z))
   out[!is.na(z) & grepl("III|\\b3\\b|STAGE[[:space:]]*3|ISS[[:space:]]*3", z)] <- 3
   out[!is.na(z) & is.na(out) & grepl("II|\\b2\\b|STAGE[[:space:]]*2|ISS[[:space:]]*2", z)] <- 2
   out[!is.na(z) & is.na(out) & grepl("\\bI\\b|\\b1\\b|STAGE[[:space:]]*1|ISS[[:space:]]*1", z)] <- 1
   out
}

make_iss_factor <- function(x) {
   ord <- extract_iss_ord(x)
   factor(
      ifelse(is.na(ord), NA, paste0("ISS ", c("I", "II", "III")[ord])),
      levels = c("ISS I", "ISS II", "ISS III")
   )
}

cohort_title <- function(x) {
   x <- gsub("_", " + ", x)
   tools::toTitleCase(x)
}

save_pdf_plot <- function(plot_obj, filename, width = 7, height = 6) {
   pdf(filename, width = width, height = height, onefile = FALSE, family = "Helvetica")
   print(plot_obj)
   dev.off()
}

## Publication themes
pub_theme <- theme_classic(base_size = 14) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      axis.title = element_text(face = "bold", size = 13, color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "gray95", color = "gray50"),
      strip.text = element_text(face = "bold", size = 11),
      plot.margin = margin(8, 12, 8, 8)
   )

risk_table_theme <- theme_classic(base_size = 10) +
   theme(
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
   )

##############################
## 4. Clean and harmonize risk data
##############################
prepare_risk_df <- function(df, set_name) {
   df <- as.data.frame(df)

   if (!"Sample" %in% colnames(df)) {
      df$Sample <- rownames(df)
   }
   if (!"set" %in% colnames(df)) {
      df$set <- set_name
   }

   needed <- c("OS.time", "OS", "riskscore")
   missing_cols <- setdiff(needed, colnames(df))
   if (length(missing_cols) > 0) {
      stop(set_name, " is missing columns: ", paste(missing_cols, collapse = ", "))
   }

   df$OS.time <- smart_numeric(df$OS.time)
   df$OS <- smart_numeric(df$OS)
   df$OS <- ifelse(is.na(df$OS), NA, ifelse(df$OS > 0, 1, 0))
   df$riskscore <- smart_numeric(df$riskscore)

   if (!"risk_group" %in% colnames(df)) {
      if (exists("cutoff")) {
         use_cutoff <- cutoff
      } else {
         use_cutoff <- median(risk_train$riskscore, na.rm = TRUE)
      }
      df$risk_group <- ifelse(df$riskscore > use_cutoff, "High risk", "Low risk")
   }
   df$risk_group <- factor(as.character(df$risk_group), levels = c("Low risk", "High risk"))

   if ("iss_stage" %in% colnames(df)) {
      df$iss_ord <- extract_iss_ord(df$iss_stage)
      df$iss_stage_f <- make_iss_factor(df$iss_stage)
   } else {
      df$iss_ord <- NA_real_
      df$iss_stage_f <- factor(NA, levels = c("ISS I", "ISS II", "ISS III"))
   }

   if ("Age" %in% colnames(df)) {
      if (is_numeric_like(df$Age)) {
         df$Age_num <- smart_numeric(df$Age)
         df$Age_f <- NA
      } else {
         df$Age_num <- NA_real_
         df$Age_f <- factor(as.character(df$Age))
      }
   } else {
      df$Age_num <- NA_real_
      df$Age_f <- NA
   }

   if ("Gender" %in% colnames(df)) {
      df$Gender_f <- factor(as.character(df$Gender))
   } else {
      df$Gender_f <- NA
   }

   if ("treatment_type" %in% colnames(df)) {
      df$treatment_type_f <- factor(as.character(df$treatment_type))
   } else {
      df$treatment_type_f <- NA
   }

   df$set <- set_name
   df
}

risk_train <- prepare_risk_df(risk_train, "training")
risk_valid <- prepare_risk_df(risk_valid, "validation")
risk_test  <- prepare_risk_df(risk_test,  "test")

## Risk-score standardization uses training-set mean and SD.
risk_mean <- mean(risk_train$riskscore, na.rm = TRUE)
risk_sd <- sd(risk_train$riskscore, na.rm = TRUE)
if (is.na(risk_sd) || risk_sd == 0) risk_sd <- 1

add_risk_z <- function(df) {
   df$risk_z <- (df$riskscore - risk_mean) / risk_sd
   df
}

risk_train <- add_risk_z(risk_train)
risk_valid <- add_risk_z(risk_valid)
risk_test  <- add_risk_z(risk_test)

## Harmonize factor levels between training and validation/test for Cox prediction.
harmonize_levels <- function(train_df, new_df, factor_vars) {
   for (v in factor_vars) {
      if (v %in% colnames(train_df) && v %in% colnames(new_df)) {
         train_df[[v]] <- factor(as.character(train_df[[v]]))
         new_df[[v]] <- factor(as.character(new_df[[v]]), levels = levels(train_df[[v]]))
      }
   }
   list(train = train_df, new = new_df)
}

risk_validtest <- bind_rows(risk_valid, risk_test)
risk_validtest$set <- "validation_test"

##############################
## 5. Build Cox models on training set for model-based C-index/DCA
##############################
choose_age_var <- function(df) {
   if ("Age_num" %in% colnames(df) && sum(!is.na(df$Age_num)) >= 20 && length(unique(na.omit(df$Age_num))) >= 5) {
      return("Age_z")
   }
   if ("Age_f" %in% colnames(df) && length(unique(na.omit(df$Age_f))) >= 2) {
      return("Age_f")
   }
   NULL
}

## Scale Age by training mean/SD if available.
age_mean <- mean(risk_train$Age_num, na.rm = TRUE)
age_sd <- sd(risk_train$Age_num, na.rm = TRUE)
if (is.na(age_mean)) age_mean <- 0
if (is.na(age_sd) || age_sd == 0) age_sd <- 1

add_age_z <- function(df) {
   df$Age_z <- ifelse(is.na(df$Age_num), NA_real_, (df$Age_num - age_mean) / age_sd)
   df
}

risk_train <- add_age_z(risk_train)
risk_valid <- add_age_z(risk_valid)
risk_test <- add_age_z(risk_test)
risk_validtest <- add_age_z(risk_validtest)

age_var <- choose_age_var(risk_train)
clinical_covars <- c("iss_stage_f", age_var, "Gender_f", "treatment_type_f")
clinical_covars <- clinical_covars[!is.na(clinical_covars) & !is.null(clinical_covars)]

valid_covariate <- function(df, v) {
   if (!v %in% colnames(df)) return(FALSE)
   x <- df[[v]]
   x <- x[!is.na(x)]
   if (length(x) < 10) return(FALSE)
   if (is.factor(x) || is.character(x)) return(length(unique(as.character(x))) >= 2)
   length(unique(x)) >= 2
}

clinical_covars <- clinical_covars[sapply(clinical_covars, function(v) valid_covariate(risk_train, v))]

fit_cox_safely <- function(df, covars, model_name = "cox") {
   covars <- covars[sapply(covars, function(v) valid_covariate(df, v))]
   if (length(covars) == 0) return(NULL)

   dat <- df[, c("OS.time", "OS", covars), drop = FALSE]
   dat <- dat[complete.cases(dat), , drop = FALSE]

   for (v in covars) {
      if (is.factor(dat[[v]]) || is.character(dat[[v]])) {
         dat[[v]] <- droplevels(factor(as.character(dat[[v]])))
      }
   }

   covars <- covars[sapply(covars, function(v) valid_covariate(dat, v))]
   if (length(covars) == 0) return(NULL)
   if (nrow(dat) < 30 || sum(dat$OS == 1, na.rm = TRUE) < 10) {
      message("Skip ", model_name, ": insufficient samples or events.")
      return(NULL)
   }

   fml <- as.formula(paste0("Surv(OS.time, OS) ~ ", paste(covars, collapse = " + ")))
   fit <- tryCatch(
      coxph(fml, data = dat, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
      error = function(e) {
         message("Failed to fit ", model_name, ": ", e$message)
         NULL
      }
   )
   fit
}

## Models used for head-to-head comparison and DCA.
fit_iss_train <- fit_cox_safely(risk_train, c("iss_stage_f"), "ISS-only model")
fit_clinical_train <- fit_cox_safely(risk_train, clinical_covars, "clinical model")
fit_iss_rsf_train <- fit_cox_safely(risk_train, c("iss_stage_f", "risk_z"), "ISS + RSF model")
fit_clinical_rsf_train <- fit_cox_safely(risk_train, unique(c(clinical_covars, "risk_z")), "clinical + RSF model")

predict_lp_safely <- function(fit, newdata) {
   if (is.null(fit)) return(rep(NA_real_, nrow(newdata)))
   nd <- as.data.frame(newdata)
   if (!is.null(fit$xlevels)) {
      for (v in names(fit$xlevels)) {
         if (v %in% colnames(nd)) {
            nd[[v]] <- factor(as.character(nd[[v]]), levels = fit$xlevels[[v]])
         }
      }
   }
   lp <- tryCatch(
      as.numeric(predict(fit, newdata = nd, type = "lp")),
      error = function(e) {
         message("Prediction failed: ", e$message)
         rep(NA_real_, nrow(nd))
      }
   )
   lp
}

add_model_markers <- function(df) {
   df$marker_RSF <- df$risk_z
   df$marker_ISS <- df$iss_ord
   df$marker_ISS_Cox <- predict_lp_safely(fit_iss_train, df)
   df$marker_Clinical <- predict_lp_safely(fit_clinical_train, df)
   df$marker_ISS_RSF <- predict_lp_safely(fit_iss_rsf_train, df)
   df$marker_Clinical_RSF <- predict_lp_safely(fit_clinical_rsf_train, df)
   df
}

risk_train_m <- add_model_markers(risk_train)
risk_valid_m <- add_model_markers(risk_valid)
risk_test_m <- add_model_markers(risk_test)
risk_validtest_m <- add_model_markers(risk_validtest)

##############################
## 6. Analysis 1: C-index head-to-head comparison
##############################
calc_cindex <- function(df, marker_col) {
   dat <- df[, c("OS.time", "OS", marker_col), drop = FALSE]
   colnames(dat) <- c("time", "status", "marker")
   dat <- dat[complete.cases(dat), , drop = FALSE]

   if (nrow(dat) < 10 || sum(dat$status == 1, na.rm = TRUE) < 3) return(NA_real_)
   if (length(unique(dat$marker)) < 2) return(NA_real_)

   fit <- tryCatch(
      survival::concordance(Surv(time, status) ~ marker, data = dat, reverse = TRUE),
      error = function(e) NULL
   )
   if (is.null(fit)) return(NA_real_)
   as.numeric(fit$concordance)
}

bootstrap_cindex <- function(df, marker_cols, reference_col = "marker_ISS", B = 1000) {
   dat0 <- df[, unique(c("OS.time", "OS", marker_cols)), drop = FALSE]
   dat0 <- dat0[complete.cases(dat0[, c("OS.time", "OS")]), , drop = FALSE]

   point <- sapply(marker_cols, function(m) calc_cindex(dat0, m))

   boot_mat <- matrix(NA_real_, nrow = B, ncol = length(marker_cols))
   colnames(boot_mat) <- marker_cols

   n <- nrow(dat0)
   if (n < 20) {
      return(list(point = point, boot = boot_mat))
   }

   for (b in seq_len(B)) {
      idx <- sample(seq_len(n), size = n, replace = TRUE)
      boot_df <- dat0[idx, , drop = FALSE]
      boot_mat[b, ] <- sapply(marker_cols, function(m) calc_cindex(boot_df, m))
   }
   list(point = point, boot = boot_mat)
}

make_cindex_table <- function(df, cohort_name, B = 1000) {
   model_map <- c(
      "RSF risk score" = "marker_RSF",
      "ISS stage" = "marker_ISS",
      "ISS Cox" = "marker_ISS_Cox",
      "Clinical model" = "marker_Clinical",
      "ISS + RSF" = "marker_ISS_RSF",
      "Clinical + RSF" = "marker_Clinical_RSF"
   )

   ## Keep only markers with at least two non-missing unique values.
   keep <- sapply(model_map, function(m) {
      m %in% colnames(df) && length(unique(na.omit(df[[m]]))) >= 2
   })
   model_map <- model_map[keep]

   boot_res <- bootstrap_cindex(df, unname(model_map), reference_col = "marker_ISS", B = B)

   out <- data.frame(
      cohort = cohort_name,
      model = names(model_map),
      marker = unname(model_map),
      cindex = as.numeric(boot_res$point[unname(model_map)]),
      lower95 = NA_real_,
      upper95 = NA_real_,
      delta_vs_ISS = NA_real_,
      delta_lower95 = NA_real_,
      delta_upper95 = NA_real_,
      p_vs_ISS = NA_real_,
      stringsAsFactors = FALSE
   )

   for (i in seq_len(nrow(out))) {
      m <- out$marker[i]
      bvals <- boot_res$boot[, m]
      bvals <- bvals[is.finite(bvals)]
      if (length(bvals) > 20) {
         out$lower95[i] <- as.numeric(quantile(bvals, 0.025, na.rm = TRUE))
         out$upper95[i] <- as.numeric(quantile(bvals, 0.975, na.rm = TRUE))
      }

      if ("marker_ISS" %in% colnames(boot_res$boot) && m != "marker_ISS") {
         d <- boot_res$boot[, m] - boot_res$boot[, "marker_ISS"]
         d <- d[is.finite(d)]
         if (length(d) > 20) {
            out$delta_vs_ISS[i] <- out$cindex[i] - out$cindex[out$marker == "marker_ISS"][1]
            out$delta_lower95[i] <- as.numeric(quantile(d, 0.025, na.rm = TRUE))
            out$delta_upper95[i] <- as.numeric(quantile(d, 0.975, na.rm = TRUE))
            out$p_vs_ISS[i] <- 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))
            out$p_vs_ISS[i] <- min(out$p_vs_ISS[i], 1)
         }
      }
   }
   out
}

message("Running C-index bootstrap comparisons. B = ", BOOTSTRAP_B)

cindex_train <- make_cindex_table(risk_train_m, "training", B = BOOTSTRAP_B)
cindex_valid <- make_cindex_table(risk_valid_m, "validation", B = BOOTSTRAP_B)
cindex_test <- make_cindex_table(risk_test_m, "test", B = BOOTSTRAP_B)
cindex_validtest <- make_cindex_table(risk_validtest_m, "validation_test", B = BOOTSTRAP_B)

cindex_all <- bind_rows(cindex_train, cindex_valid, cindex_test, cindex_validtest)
write.csv(cindex_all, file.path(tab_dir, "Cindex_head_to_head_all_cohorts.csv"), row.names = FALSE)

## C-index figure
cindex_plot_df <- cindex_all %>%
   mutate(
      cohort = factor(cohort, levels = c("training", "validation", "test", "validation_test")),
      model = factor(
         model,
         levels = c("ISS stage", "ISS Cox", "RSF risk score", "ISS + RSF", "Clinical model", "Clinical + RSF")
      ),
      label = sprintf("%.3f", cindex)
   )

p_cindex <- ggplot(cindex_plot_df, aes(x = model, y = cindex, color = model)) +
   geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.7, color = "gray55") +
   geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.16, linewidth = 0.8, na.rm = TRUE) +
   geom_point(size = 3.2, shape = 18, na.rm = TRUE) +
   geom_text(aes(label = label), vjust = -1.0, size = 3.2, color = "black", na.rm = TRUE) +
   facet_wrap(~ cohort, nrow = 1, labeller = as_labeller(function(x) cohort_title(x))) +
   scale_y_continuous(limits = c(0.45, 1.00), breaks = seq(0.5, 1.0, 0.1)) +
   scale_color_manual(
      values = c(
         "ISS stage" = "#4D4D4D",
         "ISS Cox" = "#7F7F7F",
         "RSF risk score" = "#C00000",
         "ISS + RSF" = "#2F5597",
         "Clinical model" = "#70AD47",
         "Clinical + RSF" = "#7030A0"
      )
   ) +
   labs(
      title = "Head-to-head prognostic discrimination",
      subtitle = "C-index with bootstrap 95% confidence intervals",
      x = NULL,
      y = "Harrell's C-index"
   ) +
   pub_theme +
   theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
   )

save_pdf_plot(
   p_cindex,
   file.path(fig_dir, "Cindex", "Cindex_head_to_head_all_cohorts.pdf"),
   width = 13,
   height = 5.8
)

## Delta C-index figure in merged validation+test cohort
cindex_delta_df <- cindex_validtest %>%
   filter(model != "ISS stage") %>%
   mutate(
      model = factor(
         model,
         levels = c("ISS Cox", "RSF risk score", "ISS + RSF", "Clinical model", "Clinical + RSF")
      ),
      p_label = paste0("P = ", format_p(p_vs_ISS)),
      delta_label = sprintf("%+.3f", delta_vs_ISS)
   )

p_delta <- ggplot(cindex_delta_df, aes(x = model, y = delta_vs_ISS, fill = model)) +
   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "gray50") +
   geom_col(width = 0.65, color = "black", linewidth = 0.25, alpha = 0.95) +
   geom_errorbar(aes(ymin = delta_lower95, ymax = delta_upper95), width = 0.18, linewidth = 0.75, na.rm = TRUE) +
   geom_text(aes(label = paste0(delta_label, "\n", p_label)), vjust = -0.35, size = 3.3, na.rm = TRUE) +
   scale_fill_manual(
      values = c(
         "ISS Cox" = "#7F7F7F",
         "RSF risk score" = "#C00000",
         "ISS + RSF" = "#2F5597",
         "Clinical model" = "#70AD47",
         "Clinical + RSF" = "#7030A0"
      )
   ) +
   labs(
      title = "Incremental discrimination beyond ISS",
      subtitle = "Merged validation + test cohort; delta C-index relative to ISS stage",
      x = NULL,
      y = expression(Delta*" C-index vs ISS")
   ) +
   pub_theme +
   theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
   )

save_pdf_plot(
   p_delta,
   file.path(fig_dir, "Cindex", "Delta_Cindex_vs_ISS_validation_test.pdf"),
   width = 7.8,
   height = 6.0
)

##############################
## 7. Analysis 2: Multivariable Cox regression
##############################
make_multivariable_cox <- function(df, cohort_name, risk_var = c("risk_z", "risk_group")) {
   risk_var <- match.arg(risk_var)
   covars <- unique(c(risk_var, clinical_covars))
   if (risk_var == "risk_group") {
      covars <- unique(c("risk_group", clinical_covars))
   }
   fit <- fit_cox_safely(df, covars, paste0(cohort_name, " multivariable ", risk_var))
   if (is.null(fit)) return(NULL)

   s <- summary(fit)
   coef_df <- as.data.frame(s$coefficients)
   ci_df <- as.data.frame(s$conf.int)

   out <- data.frame(
      cohort = cohort_name,
      model = risk_var,
      term = rownames(coef_df),
      HR = ci_df[, "exp(coef)"],
      lower95 = ci_df[, "lower .95"],
      upper95 = ci_df[, "upper .95"],
      pvalue = coef_df[, "Pr(>|z|)"],
      stringsAsFactors = FALSE,
      row.names = NULL
   )

   ph <- tryCatch(cox.zph(fit), error = function(e) NULL)
   list(fit = fit, table = out, ph = ph)
}

## Primary independence analysis in merged validation + test cohort.
cox_validtest_score <- make_multivariable_cox(risk_validtest, "validation_test", risk_var = "risk_z")
cox_validtest_group <- make_multivariable_cox(risk_validtest, "validation_test", risk_var = "risk_group")

## Sensitivity analyses in training, validation, and test cohorts separately.
cox_train_score <- make_multivariable_cox(risk_train, "training", risk_var = "risk_z")
cox_valid_score <- make_multivariable_cox(risk_valid, "validation", risk_var = "risk_z")
cox_test_score <- make_multivariable_cox(risk_test, "test", risk_var = "risk_z")

cox_tables <- bind_rows(
   if (!is.null(cox_validtest_score)) cox_validtest_score$table,
   if (!is.null(cox_validtest_group)) cox_validtest_group$table,
   if (!is.null(cox_train_score)) cox_train_score$table,
   if (!is.null(cox_valid_score)) cox_valid_score$table,
   if (!is.null(cox_test_score)) cox_test_score$table
)

write.csv(cox_tables, file.path(tab_dir, "Multivariable_Cox_all_results.csv"), row.names = FALSE)

if (!is.null(cox_validtest_score$ph)) {
   capture.output(cox_validtest_score$ph, file = file.path(tab_dir, "PH_test_multivariable_cox_validtest_risk_score.txt"))
}
if (!is.null(cox_validtest_group$ph)) {
   capture.output(cox_validtest_group$ph, file = file.path(tab_dir, "PH_test_multivariable_cox_validtest_risk_group.txt"))
}

clean_term_label <- function(term) {
   out <- term
   out <- gsub("risk_z", "RSF risk score (per 1-SD increase)", out)
   out <- gsub("risk_groupHigh risk", "RSF high risk vs low risk", out)
   out <- gsub("iss_stage_f", "", out)
   out <- gsub("Age_z", "Age (per 1-SD increase)", out)
   out <- gsub("Age_f", "Age: ", out)
   out <- gsub("Gender_f", "Gender: ", out)
   out <- gsub("treatment_type_f", "Treatment: ", out)
   out <- gsub("`", "", out)
   out
}

plot_cox_forest <- function(cox_df, title_text, out_pdf) {
   if (is.null(cox_df) || nrow(cox_df) == 0) return(NULL)

   plot_df <- cox_df %>%
      mutate(
         label = clean_term_label(term),
         label = factor(label, levels = rev(clean_term_label(term))),
         p_label = paste0("P = ", format_p(pvalue)),
         hr_label = sprintf("%.2f (%.2f-%.2f)", HR, lower95, upper95)
      ) %>%
      filter(is.finite(HR), is.finite(lower95), is.finite(upper95))

   if (nrow(plot_df) == 0) return(NULL)

   xmax <- max(plot_df$upper95, na.rm = TRUE)
   xmin <- min(plot_df$lower95, na.rm = TRUE)
   text_x <- xmax * 1.25
   p_x <- xmax * 2.05

   p <- ggplot(plot_df, aes(x = HR, y = label)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, color = "gray45") +
      geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.18, linewidth = 0.8, color = "gray25") +
      geom_point(size = 3.2, shape = 18, color = "#C00000") +
      geom_text(aes(x = text_x, label = hr_label), hjust = 0, size = 3.2) +
      geom_text(aes(x = p_x, label = p_label), hjust = 0, size = 3.2) +
      annotate("text", x = text_x, y = nrow(plot_df) + 0.75, label = "HR (95% CI)", hjust = 0, fontface = "bold", size = 3.4) +
      annotate("text", x = p_x, y = nrow(plot_df) + 0.75, label = "P value", hjust = 0, fontface = "bold", size = 3.4) +
      scale_x_log10() +
      coord_cartesian(xlim = c(xmin * 0.8, xmax * 2.7), clip = "off") +
      labs(
         title = title_text,
         x = "Hazard ratio (log scale)",
         y = NULL
      ) +
      pub_theme +
      theme(
         plot.margin = margin(10, 140, 10, 10),
         legend.position = "none"
      )

   save_pdf_plot(p, out_pdf, width = 9.5, height = max(5.2, 0.50 * nrow(plot_df) + 2.2))
}

if (!is.null(cox_validtest_score)) {
   plot_cox_forest(
      cox_validtest_score$table,
      "Multivariable Cox analysis: RSF risk score",
      file.path(fig_dir, "Multivariable_Cox", "Multivariable_Cox_validtest_RSF_score.pdf")
   )
}

if (!is.null(cox_validtest_group)) {
   plot_cox_forest(
      cox_validtest_group$table,
      "Multivariable Cox analysis: RSF risk group",
      file.path(fig_dir, "Multivariable_Cox", "Multivariable_Cox_validtest_RSF_group.pdf")
   )
}

##############################
## 8. Analysis 3: ISS incremental stratification by RSF high/low risk
##############################
validtest_iss <- risk_validtest %>%
   filter(
      !is.na(OS.time), !is.na(OS), !is.na(risk_group),
      !is.na(iss_stage_f)
   ) %>%
   mutate(
      risk_group = factor(risk_group, levels = c("Low risk", "High risk")),
      iss_stage_f = droplevels(iss_stage_f)
   )

write.csv(validtest_iss, file.path(tab_dir, "Validation_test_data_for_ISS_incremental_analysis.csv"), row.names = FALSE)

## 8.1 KM curves within each ISS stage
plot_km_within_iss <- function(df, iss_level, out_pdf) {
   dat <- df %>% filter(iss_stage_f == iss_level)
   dat$risk_group <- droplevels(dat$risk_group)

   if (nrow(dat) < 15 || length(unique(dat$risk_group)) < 2 || sum(dat$OS == 1, na.rm = TRUE) < 3) {
      message("Skip KM for ", iss_level, ": insufficient data.")
      return(NULL)
   }

   km_fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = dat)
   cox_fit <- tryCatch(coxph(Surv(OS.time, OS) ~ risk_group, data = dat), error = function(e) NULL)
   survdiff_fit <- survdiff(Surv(OS.time, OS) ~ risk_group, data = dat)
   p_logrank <- 1 - pchisq(survdiff_fit$chisq, df = length(survdiff_fit$n) - 1)

   hr_text <- ""
   if (!is.null(cox_fit)) {
      cox_sum <- summary(cox_fit)
      hr <- cox_sum$conf.int[1, "exp(coef)"]
      lo <- cox_sum$conf.int[1, "lower .95"]
      hi <- cox_sum$conf.int[1, "upper .95"]
      hr_text <- sprintf("HR = %.2f (95%% CI: %.2f-%.2f)", hr, lo, hi)
   }

   max_time <- max(dat$OS.time, na.rm = TRUE)
   break_by <- signif(max_time / 5, 1)
   if (is.na(break_by) || break_by <= 0) break_by <- 500

   n_low <- sum(dat$risk_group == "Low risk", na.rm = TRUE)
   n_high <- sum(dat$risk_group == "High risk", na.rm = TRUE)

   g <- ggsurvplot(
      km_fit,
      data = dat,
      title = paste0("RSF risk stratification within ", iss_level),
      xlab = "Time (days)",
      ylab = "Overall survival probability",
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.y.text = FALSE,
      risk.table.height = 0.22,
      risk.table.theme = risk_table_theme,
      legend.title = "RSF group",
      legend.labs = c(paste0("Low risk (n = ", n_low, ")"), paste0("High risk (n = ", n_high, ")")),
      palette = c("#2F5597", "#C00000"),
      censor = TRUE,
      censor.shape = 124,
      censor.size = 2.8,
      conf.int = FALSE,
      pval = FALSE,
      break.time.by = break_by,
      xlim = c(0, max_time * 1.02),
      size = 1.15,
      ggtheme = pub_theme
   )

   g$plot <- g$plot +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
      labs(subtitle = "Merged validation + test cohort") +
      annotate(
         "text",
         x = max_time * 0.52,
         y = 0.18,
         hjust = 0,
         size = 4,
         label = paste(hr_text, paste0("Log-rank P = ", format_p(p_logrank)), sep = "\n")
      ) +
      theme(legend.position = c(0.78, 0.84))

   pdf(out_pdf, width = 7.3, height = 6.3, onefile = FALSE, family = "Helvetica")
   print(g)
   dev.off()
}

for (lv in levels(droplevels(validtest_iss$iss_stage_f))) {
   plot_km_within_iss(
      validtest_iss,
      lv,
      file.path(fig_dir, "ISS_Incremental", paste0("KM_within_", safe_filename(lv), "_by_RSF.pdf"))
   )
}

## 8.2 Six-group joint stratification: ISS x RSF
validtest_iss <- validtest_iss %>%
   mutate(
      ISS_RSF_group = interaction(iss_stage_f, risk_group, sep = " / ", drop = TRUE)
   )

joint_levels <- as.vector(outer(
   c("ISS I", "ISS II", "ISS III"),
   c("Low risk", "High risk"),
   paste,
   sep = " / "
))
validtest_iss$ISS_RSF_group <- factor(as.character(validtest_iss$ISS_RSF_group), levels = joint_levels)
validtest_iss$ISS_RSF_group <- droplevels(validtest_iss$ISS_RSF_group)

if (length(unique(validtest_iss$ISS_RSF_group)) >= 2 && sum(validtest_iss$OS == 1, na.rm = TRUE) >= 5) {
   joint_fit <- survfit(Surv(OS.time, OS) ~ ISS_RSF_group, data = validtest_iss)
   joint_palette <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6", "#F79646")
   joint_palette <- joint_palette[seq_len(length(levels(validtest_iss$ISS_RSF_group)))]
   max_time <- max(validtest_iss$OS.time, na.rm = TRUE)
   break_by <- signif(max_time / 5, 1)
   if (is.na(break_by) || break_by <= 0) break_by <- 500

   g_joint <- ggsurvplot(
      joint_fit,
      data = validtest_iss,
      title = "Joint risk stratification by ISS and RSF",
      xlab = "Time (days)",
      ylab = "Overall survival probability",
      risk.table = TRUE,
      risk.table.height = 0.30,
      risk.table.y.text = FALSE,
      risk.table.theme = risk_table_theme,
      legend.title = "ISS / RSF group",
      legend.labs = levels(validtest_iss$ISS_RSF_group),
      palette = joint_palette,
      censor = TRUE,
      conf.int = FALSE,
      pval = TRUE,
      break.time.by = break_by,
      xlim = c(0, max_time * 1.02),
      size = 1.0,
      ggtheme = pub_theme
   )
   g_joint$plot <- g_joint$plot +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
      labs(subtitle = "Merged validation + test cohort") +
      theme(legend.position = "right")

   pdf(
      file.path(fig_dir, "ISS_Incremental", "KM_joint_ISS_RSF_six_groups.pdf"),
      width = 9.5,
      height = 7.2,
      onefile = FALSE,
      family = "Helvetica"
   )
   print(g_joint)
   dev.off()
}

## 8.3 Cox HR for high vs low RSF risk within each ISS stage
cox_within_iss <- function(df, iss_level) {
   dat <- df %>% filter(iss_stage_f == iss_level)
   dat$risk_group <- droplevels(dat$risk_group)
   if (nrow(dat) < 15 || length(unique(dat$risk_group)) < 2 || sum(dat$OS == 1, na.rm = TRUE) < 3) return(NULL)

   fit <- tryCatch(coxph(Surv(OS.time, OS) ~ risk_group, data = dat), error = function(e) NULL)
   if (is.null(fit)) return(NULL)
   s <- summary(fit)

   data.frame(
      ISS_stage = as.character(iss_level),
      n = nrow(dat),
      events = sum(dat$OS == 1, na.rm = TRUE),
      n_low = sum(dat$risk_group == "Low risk", na.rm = TRUE),
      n_high = sum(dat$risk_group == "High risk", na.rm = TRUE),
      HR_high_vs_low = s$conf.int[1, "exp(coef)"],
      lower95 = s$conf.int[1, "lower .95"],
      upper95 = s$conf.int[1, "upper .95"],
      pvalue = s$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
   )
}

iss_hr_table <- bind_rows(lapply(levels(droplevels(validtest_iss$iss_stage_f)), function(lv) cox_within_iss(validtest_iss, lv)))
write.csv(iss_hr_table, file.path(tab_dir, "ISS_stratified_RSF_high_vs_low_Cox.csv"), row.names = FALSE)

if (nrow(iss_hr_table) > 0) {
   iss_hr_plot <- iss_hr_table %>%
      mutate(
         ISS_stage = factor(ISS_stage, levels = rev(c("ISS I", "ISS II", "ISS III"))),
         label = sprintf("%.2f (%.2f-%.2f)", HR_high_vs_low, lower95, upper95),
         p_label = paste0("P = ", format_p(pvalue))
      )

   xmax <- max(iss_hr_plot$upper95, na.rm = TRUE)
   xmin <- min(iss_hr_plot$lower95, na.rm = TRUE)
   text_x <- xmax * 1.20
   p_x <- xmax * 2.00

   p_iss_forest <- ggplot(iss_hr_plot, aes(x = HR_high_vs_low, y = ISS_stage)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, color = "gray45") +
      geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.18, linewidth = 0.9, color = "gray25") +
      geom_point(size = 3.4, shape = 18, color = "#C00000") +
      geom_text(aes(x = text_x, label = label), hjust = 0, size = 3.4) +
      geom_text(aes(x = p_x, label = p_label), hjust = 0, size = 3.4) +
      annotate("text", x = text_x, y = nrow(iss_hr_plot) + 0.55, label = "HR (95% CI)", hjust = 0, fontface = "bold", size = 3.5) +
      annotate("text", x = p_x, y = nrow(iss_hr_plot) + 0.55, label = "P value", hjust = 0, fontface = "bold", size = 3.5) +
      scale_x_log10() +
      coord_cartesian(xlim = c(xmin * 0.8, xmax * 2.7), clip = "off") +
      labs(
         title = "RSF high-risk effect within ISS strata",
         subtitle = "Merged validation + test cohort",
         x = "Hazard ratio for high vs low RSF risk (log scale)",
         y = NULL
      ) +
      pub_theme +
      theme(plot.margin = margin(10, 145, 10, 10))

   save_pdf_plot(
      p_iss_forest,
      file.path(fig_dir, "ISS_Incremental", "Forest_RSF_high_vs_low_within_ISS.pdf"),
      width = 9.0,
      height = 5.2
   )
}

## 8.4 Composition plot: percentage of RSF high/low risk in each ISS stage

validtest_iss <- dplyr::bind_rows(risk_valid, risk_test)
cutoff <- median(risk_train$riskscore, na.rm = TRUE)

validtest_iss$risk_group <- ifelse(
   validtest_iss$riskscore > cutoff,
   "High risk",
   "Low risk"
)

validtest_iss$risk_group <- factor(validtest_iss$risk_group,
                                   levels = c("Low risk", "High risk"))
validtest_iss$iss_stage_f <- as.factor(validtest_iss$iss_stage)

composition_df <- validtest_iss %>%
   dplyr::count(iss_stage_f, risk_group) %>%
   dplyr::group_by(iss_stage_f) %>%
   dplyr::mutate(prop = n / sum(n)) %>%
   dplyr::ungroup()

write.csv(composition_df, file.path(tab_dir, "ISS_RSF_group_composition.csv"), row.names = FALSE)

p_comp <- ggplot(composition_df, aes(x = iss_stage_f, y = prop, fill = risk_group)) +
   geom_col(width = 0.68, color = "black", linewidth = 0.25) +
   geom_text(aes(label = paste0(n, "\n", percent(prop, accuracy = 1))),
             position = position_stack(vjust = 0.5), size = 3.4, color = "white") +
   scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
   scale_fill_manual(values = c("Low risk" = "#2F5597", "High risk" = "#C00000")) +
   labs(
      title = "Distribution of RSF risk groups within ISS stages",
      subtitle = "Merged validation + test cohort",
      x = "ISS stage",
      y = "Proportion of patients",
      fill = "RSF group"
   ) +
   pub_theme +
   theme(legend.position = "top")

save_pdf_plot(
   p_comp,
   file.path(fig_dir, "ISS_Incremental", "ISS_stage_RSF_risk_group_composition.pdf"),
   width = 6.8,
   height = 5.6
)

## 8.5 Nested model comparison in merged validation + test cohort
fit_nested_validtest <- function(df) {
   dat <- df[, c("OS.time", "OS", "iss_stage_f", "risk_z"), drop = FALSE]
   dat <- dat[complete.cases(dat), , drop = FALSE]
   dat$iss_stage_f <- droplevels(dat$iss_stage_f)
   if (nrow(dat) < 30 || sum(dat$OS == 1, na.rm = TRUE) < 10 || length(unique(dat$iss_stage_f)) < 2) return(NULL)

   fit_base <- coxph(Surv(OS.time, OS) ~ iss_stage_f, data = dat, ties = "efron", x = TRUE, y = TRUE)
   fit_plus <- coxph(Surv(OS.time, OS) ~ iss_stage_f + risk_z, data = dat, ties = "efron", x = TRUE, y = TRUE)
   lrt <- anova(fit_base, fit_plus, test = "LRT")

   data.frame(
      cohort = "validation_test",
      n = nrow(dat),
      events = sum(dat$OS == 1, na.rm = TRUE),
      model_base = "ISS",
      model_enhanced = "ISS + RSF risk score",
      AIC_base = AIC(fit_base),
      AIC_enhanced = AIC(fit_plus),
      delta_AIC = AIC(fit_plus) - AIC(fit_base),
      logLik_base = as.numeric(logLik(fit_base)),
      logLik_enhanced = as.numeric(logLik(fit_plus)),
      LRT_Chisq = lrt$Chisq[2],
      LRT_df = lrt$Df[2],
      LRT_pvalue = lrt$`Pr(>|Chi|)`[2],
      stringsAsFactors = FALSE
   )
}

nested_table <- fit_nested_validtest(validtest_iss)
if (!is.null(nested_table)) {
   write.csv(nested_table, file.path(tab_dir, "Nested_model_comparison_ISS_vs_ISS_plus_RSF_validtest.csv"), row.names = FALSE)
}

##############################
## 9. Decision curve analysis: ISS vs ISS + RSF
##############################
## DCA uses Cox models fitted on the training set and predicted event probabilities in validation + test.
basehaz_at_time <- function(fit, t) {
   bh <- basehaz(fit, centered = FALSE)
   if (nrow(bh) == 0) return(NA_real_)
   approx(bh$time, bh$hazard, xout = t, method = "constant", f = 0, rule = 2)$y
}

predict_event_prob_cox <- function(fit, newdata, t) {
   if (is.null(fit)) return(rep(NA_real_, nrow(newdata)))
   nd <- as.data.frame(newdata)
   if (!is.null(fit$xlevels)) {
      for (v in names(fit$xlevels)) {
         if (v %in% colnames(nd)) {
            nd[[v]] <- factor(as.character(nd[[v]]), levels = fit$xlevels[[v]])
         }
      }
   }
   H0 <- basehaz_at_time(fit, t)
   if (is.na(H0)) return(rep(NA_real_, nrow(nd)))
   lp <- predict_lp_safely(fit, nd)
   pred <- 1 - exp(-H0 * exp(lp))
   as.numeric(pred)
}

net_benefit_curve <- function(event, pred, thresholds) {
   ok <- !is.na(event) & !is.na(pred)
   event <- event[ok]
   pred <- pred[ok]
   n <- length(event)
   if (n < 20 || length(unique(event)) < 2) return(NULL)

   data.frame(
      threshold = thresholds,
      net_benefit = sapply(thresholds, function(pt) {
         positive <- pred >= pt
         TP <- sum(positive & event == 1)
         FP <- sum(positive & event == 0)
         TP / n - FP / n * pt / (1 - pt)
      }),
      stringsAsFactors = FALSE
   )
}

make_dca_compare <- function(df, times, thresholds) {
   out_list <- list()
   for (tt in times) {
      dat <- df[, c("OS.time", "OS", "iss_stage_f", "risk_z"), drop = FALSE]
      dat <- dat[complete.cases(dat[, c("OS.time", "OS")]), , drop = FALSE]
      dat$event_t <- NA_integer_
      dat$event_t[dat$OS == 1 & dat$OS.time <= tt] <- 1
      dat$event_t[dat$OS.time > tt] <- 0
      dat <- dat[!is.na(dat$event_t), , drop = FALSE]

      if (nrow(dat) < 30 || length(unique(dat$event_t)) < 2) {
         message("Skip DCA at ", tt, " days: insufficient eligible samples/outcomes.")
         next
      }

      pred_iss <- predict_event_prob_cox(fit_iss_train, dat, tt)
      pred_iss_rsf <- predict_event_prob_cox(fit_iss_rsf_train, dat, tt)

      nb_iss <- net_benefit_curve(dat$event_t, pred_iss, thresholds)
      nb_iss_rsf <- net_benefit_curve(dat$event_t, pred_iss_rsf, thresholds)

      prevalence <- mean(dat$event_t == 1, na.rm = TRUE)
      nb_all <- data.frame(
         threshold = thresholds,
         net_benefit = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
         model = "Treat all",
         time_day = tt,
         time_label = paste0(round(tt / 365), "-year"),
         stringsAsFactors = FALSE
      )
      nb_none <- data.frame(
         threshold = thresholds,
         net_benefit = 0,
         model = "Treat none",
         time_day = tt,
         time_label = paste0(round(tt / 365), "-year"),
         stringsAsFactors = FALSE
      )

      if (!is.null(nb_iss)) {
         nb_iss$model <- "ISS"
         nb_iss$time_day <- tt
         nb_iss$time_label <- paste0(round(tt / 365), "-year")
      }
      if (!is.null(nb_iss_rsf)) {
         nb_iss_rsf$model <- "ISS + RSF"
         nb_iss_rsf$time_day <- tt
         nb_iss_rsf$time_label <- paste0(round(tt / 365), "-year")
      }

      out_list[[as.character(tt)]] <- bind_rows(nb_iss, nb_iss_rsf, nb_all, nb_none)
   }
   bind_rows(out_list)
}

dca_df <- make_dca_compare(risk_validtest, DCA_TIMES, DCA_THRESHOLDS)

if (nrow(dca_df) > 0) {
   dca_df$time_label <- factor(dca_df$time_label, levels = paste0(c(1, 2, 3), "-year"))
   dca_df$model <- factor(dca_df$model, levels = c("ISS + RSF", "ISS", "Treat all", "Treat none"))
   write.csv(dca_df, file.path(tab_dir, "DCA_ISS_vs_ISS_plus_RSF_validtest.csv"), row.names = FALSE)

   p_dca <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = model, linetype = model)) +
      geom_line(linewidth = 1.0, na.rm = TRUE) +
      facet_wrap(~ time_label, nrow = 1) +
      scale_color_manual(
         values = c(
            "ISS + RSF" = "#C00000",
            "ISS" = "#2F5597",
            "Treat all" = "gray45",
            "Treat none" = "black"
         )
      ) +
      scale_linetype_manual(
         values = c(
            "ISS + RSF" = "solid",
            "ISS" = "solid",
            "Treat all" = "dashed",
            "Treat none" = "dotted"
         )
      ) +
      coord_cartesian(ylim = c(min(dca_df$net_benefit, na.rm = TRUE), max(dca_df$net_benefit, na.rm = TRUE))) +
      labs(
         title = "Decision curve analysis: ISS vs ISS + RSF",
         subtitle = "Predictions from training-set Cox models; evaluation in merged validation + test cohort",
         x = "Threshold probability",
         y = "Net benefit",
         color = "Model",
         linetype = "Model"
      ) +
      pub_theme +
      theme(legend.position = "top")

   save_pdf_plot(
      p_dca,
      file.path(fig_dir, "DCA", "DCA_ISS_vs_ISS_plus_RSF_validtest.pdf"),
      width = 11,
      height = 5.6
   )
}

##############################
## 10. Save analysis objects and session info
##############################
save(
   cindex_all,
   cox_tables,
   iss_hr_table,
   composition_df,
   nested_table,
   dca_df,
   file = file.path(out_dir, "Clinical_Incremental_Analysis_results.RData")
)
save.image(
   file = file.path(out_dir,  "Clinical_Incremental_Analysis_results_full_workspace.RData")
)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("\nClinical incremental analyses finished.\n")
cat("Output directory:\n", out_dir, "\n")
cat("Key tables saved in:\n", tab_dir, "\n")
cat("Key figures saved in:\n", fig_dir, "\n")
