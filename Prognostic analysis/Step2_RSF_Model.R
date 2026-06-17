##############################
## 0. Packages
##############################
library(survival)
library(randomForestSRC)
library(survminer)
library(ggplot2)
library(dplyr)

set.seed(123)
options(stringsAsFactors = FALSE)

##############################
## 1. Read data
##############################
training <- read.table(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/training_sets_symbol_clinic.csv",
   header = TRUE, sep = ",", quote = "", check.names = FALSE, row.names = 1
)

valid <- read.table(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/valid_sets_symbol.csv",
   header = TRUE, sep = ",", quote = "", check.names = FALSE, row.names = 1
)

test <- read.table(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/test_sets_symbol.csv",
   header = TRUE, sep = ",", quote = "", check.names = FALSE, row.names = 1
)

mm <- list(training = training, valid = valid, test = test)

##############################
## 2. Define columns
##############################

gene_cols <- colnames(training)[-c(1:7)]

##############################
## 3. Standardize gene features using training statistics
##############################
mu  <- sapply(training[, gene_cols, drop = FALSE], mean, na.rm = TRUE)
sdv <- sapply(training[, gene_cols, drop = FALSE], sd,   na.rm = TRUE)
sdv[sdv == 0 | is.na(sdv)] <- 1

scale_with_train <- function(df, gene_cols, mu, sdv) {
   x <- sweep(df[, gene_cols, drop = FALSE], 2, mu, "-")
   x <- sweep(x, 2, sdv, "/")
   df[, gene_cols] <- x
   df[, gene_cols][is.na(df[, gene_cols])] <- 0
   df
}

mm <- lapply(mm, scale_with_train, gene_cols = gene_cols, mu = mu, sdv = sdv)
training <- mm$training
valid    <- mm$valid
test     <- mm$test

##############################
## 4. Fix clinical variable types
##############################
fix_clin <- function(df) {
   df$Gender         <- as.factor(df$Gender)
   df$Age            <- as.factor(df$Age)
   df$iss_stage      <- as.factor(df$iss_stage)
   df$treatment_type <- as.factor(df$treatment_type)
   df
}

mm <- lapply(mm, fix_clin)
training <- mm$training
valid    <- mm$valid
test     <- mm$test

##############################
## 5. Train RSF using gene features only
##############################
rsf_formula <- as.formula(
   paste0("Surv(OS.time, OS) ~ ", paste(gene_cols, collapse = " + "))
)

rsf_fit <- rfsrc(
   formula    = rsf_formula,
   data       = training,
   ntree      = 1000,
   nodesize   = 10,
   splitrule  = "logrank",
   importance = TRUE,
   forest     = TRUE,
   na.action  = "na.impute"
)

##############################
## 6. Predict risk score
##############################
get_risk_df <- function(df, set_name, fit) {
   pred <- predict(fit, newdata = df, na.action = "na.impute")
   score <- as.numeric(pred$predicted)
   
   out <- data.frame(
      Sample         = rownames(df),
      set            = set_name,
      OS.time        = df$OS.time,
      OS             = df$OS,
      riskscore      = score,
      Gender         = df$Gender,
      Age            = df$Age,
      iss_stage      = df$iss_stage,
      treatment_type = df$treatment_type,
      stringsAsFactors = FALSE
   )
   out
}

risk_train <- get_risk_df(training, "training", rsf_fit)
risk_valid <- get_risk_df(valid,    "validation", rsf_fit)
risk_test  <- get_risk_df(test,     "test", rsf_fit)

##############################
## 7. Use training median as cutoff
##############################
cutoff <- median(risk_train$riskscore, na.rm = TRUE)

assign_group <- function(risk_df, cutoff) {
   risk_df$risk_group <- factor(
      ifelse(risk_df$riskscore > cutoff, "High risk", "Low risk"),
      levels = c("Low risk", "High risk")
   )
   risk_df
}

risk_train <- assign_group(risk_train, cutoff)
risk_valid <- assign_group(risk_valid, cutoff)
risk_test  <- assign_group(risk_test, cutoff)

##############################
## 8. Output directories
##############################
base_dir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "Overall_KM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "Subgroup_KM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "ForestPlot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

##############################
## 9. Plot overall KM (publication-ready)
##############################
library(survival)
library(survminer)
library(ggplot2)
library(scales)

## -------- 1. Publication-style themes --------
pub_theme <- theme_classic(base_size = 14) +
   theme(
      plot.title      = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.title      = element_text(size = 13, face = "bold", color = "black"),
      axis.text       = element_text(size = 11, color = "black"),
      axis.line       = element_line(linewidth = 0.8, color = "black"),
      axis.ticks      = element_line(linewidth = 0.7, color = "black"),
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10, color = "black"),
      legend.position = c(0.80, 0.85),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin     = margin(8, 12, 5, 8)
   )

risk_table_theme <- theme_classic(base_size = 11) +
   theme(
      axis.title.y    = element_blank(),
      axis.title.x    = element_text(size = 11, face = "bold", color = "black"),
      axis.text.x     = element_text(size = 10, color = "black"),
      axis.text.y     = element_text(size = 10, color = "black"),
      axis.line       = element_blank(),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      plot.margin     = margin(0, 12, 0, 8)
   )

## -------- 2. Main plotting function --------
plot_km_sci <- function(risk_df, title_text, out_pdf) {
   
   ## keep only complete samples for survival analysis
   risk_df <- risk_df[!is.na(risk_df$OS.time) & !is.na(risk_df$OS) & !is.na(risk_df$risk_group), ]
   risk_df$risk_group <- factor(risk_df$risk_group, levels = c("Low risk", "High risk"))
   
   ## sample size in each group
   n_low  <- sum(risk_df$risk_group == "Low risk")
   n_high <- sum(risk_df$risk_group == "High risk")
   
   ## KM and Cox
   km_fit  <- survfit(Surv(OS.time, OS) ~ risk_group, data = risk_df)
   cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = risk_df)
   cox_sum <- summary(cox_fit)
   
   ## HR and 95% CI
   hr       <- cox_sum$coefficients[1, "exp(coef)"]
   hr_low   <- cox_sum$conf.int[1, "lower .95"]
   hr_high  <- cox_sum$conf.int[1, "upper .95"]
   p_cox    <- cox_sum$coefficients[1, "Pr(>|z|)"]
   
   ## log-rank p value
   survdiff_fit <- survdiff(Surv(OS.time, OS) ~ risk_group, data = risk_df)
   p_logrank <- 1 - pchisq(survdiff_fit$chisq, df = length(survdiff_fit$n) - 1)
   
   ## pretty p text
   p_text <- ifelse(p_logrank < 0.001, "Log-rank P < 0.001",
                    paste0("Log-rank P = ", formatC(p_logrank, format = "f", digits = 3)))
   
   hr_text <- sprintf("HR = %.2f (95%% CI: %.2f-%.2f)", hr, hr_low, hr_high)
   
   ## dynamic x-axis break
   max_time <- max(risk_df$OS.time, na.rm = TRUE)
   break_by <- signif(max_time / 5, 1)
   if (break_by <= 0 || is.na(break_by)) break_by <- 500
   
   ## legend labels with n
   legend_labels <- c(
      paste0("Low risk (n = ", n_low, ")"),
      paste0("High risk (n = ", n_high, ")")
   )
   
   ## colors suitable for publication
   km_cols <- c("#2F5597", "#C00000")
   
   ## draw plot
   g <- ggsurvplot(
      fit               = km_fit,
      data              = risk_df,
      title             = title_text,
      xlab              = "Time (days)",
      ylab              = "Overall survival probability",
      palette           = km_cols,
      legend.title      = "Risk group",
      legend.labs       = legend_labels,
      risk.table        = TRUE,
      risk.table.col    = "strata",
      risk.table.y.text = FALSE,
      risk.table.height = 0.22,
      risk.table.fontsize = 3.5,
      risk.table.theme  = risk_table_theme,
      break.time.by     = break_by,
      xlim              = c(0, max_time * 1.02),
      conf.int          = FALSE,
      censor            = TRUE,
      censor.shape      = 124,
      censor.size       = 3,
      pval              = FALSE,
      surv.median.line  = "none",
      linetype          = "solid",
      size              = 1.2,
      ggtheme           = pub_theme
   )
   
   ## refine main KM panel
   g$plot <- g$plot +
      scale_y_continuous(
         limits = c(0, 1),
         labels = percent_format(accuracy = 1),
         expand = expansion(mult = c(0, 0.02))
      ) +
      annotate(
         "text",
         x = max_time * 0.55,
         y = 0.20,
         label = paste(hr_text, p_text, sep = "\n"),
         hjust = 0,
         size = 4.2,
         fontface = "plain"
      )
   
   ## refine risk table
   g$table <- g$table +
      theme(
         legend.position = "none"
      )
   
   ## save as vector PDF
   pdf(out_pdf, width = 7.2, height = 6.2, onefile = FALSE)
   print(g)
   dev.off()
}

## -------- 3. Output figures --------
plot_km_sci(
   risk_train,
   "Training cohort",
   file.path(base_dir, "Overall_KM", "KM_train_SCI.pdf")
)

plot_km_sci(
   risk_valid,
   "Validation cohort",
   file.path(base_dir, "Overall_KM", "KM_valid_SCI.pdf")
)

plot_km_sci(
   risk_test,
   "Test cohort",
   file.path(base_dir, "Overall_KM", "KM_test_SCI.pdf")
)


##############################
## 10. Save risk tables
##############################
write.csv(risk_train, file.path(base_dir, "Tables", "risk_train.csv"), row.names = FALSE)
write.csv(risk_valid, file.path(base_dir, "Tables", "risk_valid.csv"), row.names = FALSE)
write.csv(risk_test,  file.path(base_dir, "Tables", "risk_test.csv"), row.names = FALSE)

##############################
## 10A. Time-dependent ROC (publication-ready)
##############################
library(timeROC)
library(ggplot2)
library(dplyr)

dir.create(file.path(base_dir, "ROC"), recursive = TRUE, showWarnings = FALSE)

## 1/2/3年，对应天数
roc_times <- c(365, 730, 1095)

## -------- 1. publication theme --------
roc_theme <- theme_classic(base_size = 14) +
   theme(
      plot.title      = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.title      = element_text(size = 13, face = "bold", color = "black"),
      axis.text       = element_text(size = 11, color = "black"),
      axis.line       = element_line(linewidth = 0.8, color = "black"),
      axis.ticks      = element_line(linewidth = 0.7, color = "black"),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10.5, color = "black"),
      legend.position = c(0.73, 0.20),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin     = margin(8, 10, 8, 8)
   )

## -------- 2. main plotting function --------
plot_timeROC_sci <- function(risk_df, set_name, out_pdf, roc_times = c(365, 730, 1095)) {
   
   ## keep complete cases
   df <- risk_df[complete.cases(risk_df[, c("OS.time", "OS", "riskscore")]), , drop = FALSE]
   
   ## basic QC
   if (nrow(df) < 10 || sum(df$OS == 1, na.rm = TRUE) < 5) {
      message("Skip ROC for ", set_name, ": insufficient samples or events.")
      return(NULL)
   }
   
   ## keep only time points within observable follow-up range
   max_time <- max(df$OS.time, na.rm = TRUE)
   valid_times <- roc_times[roc_times < max_time]
   
   if (length(valid_times) == 0) {
      message("Skip ROC for ", set_name, ": follow-up time shorter than requested ROC times.")
      return(NULL)
   }
   
   ## fit timeROC
   roc_obj <- tryCatch(
      timeROC(
         T      = df$OS.time,
         delta  = df$OS,
         marker = df$riskscore,
         cause  = 1,
         times  = valid_times,
         iid    = TRUE
      ),
      error = function(e) {
         message("Skip ROC for ", set_name, ": ", e$message)
         return(NULL)
      }
   )
   
   if (is.null(roc_obj)) return(NULL)
   
   ## colors for 1/2/3-year ROC
   roc_cols <- c("#C00000", "#ED7D31", "#2F5597")[seq_along(valid_times)]
   
   ## labels
   time_labels <- paste0(round(valid_times / 365), "-year")
   auc_labels <- paste0(
      time_labels, " AUC = ",
      formatC(roc_obj$AUC, format = "f", digits = 3)
   )
   
   ## build ROC data.frame for ggplot
   roc_df_list <- lapply(seq_along(valid_times), function(i) {
      data.frame(
         FPR  = roc_obj$FP[, i],
         TPR  = roc_obj$TP[, i],
         Time = factor(auc_labels[i], levels = auc_labels)
      )
   })
   roc_plot_df <- bind_rows(roc_df_list)
   
   ## title style
   title_text <- paste0(tools::toTitleCase(set_name), " cohort")
   
   ## plot
   p <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = Time)) +
      geom_line(linewidth = 1.2) +
      geom_abline(
         slope = 1, intercept = 0,
         linetype = "dashed", linewidth = 0.7, color = "gray55"
      ) +
      scale_color_manual(values = roc_cols) +
      scale_x_continuous(
         limits = c(0, 1),
         breaks = seq(0, 1, 0.2),
         expand = c(0, 0)
      ) +
      scale_y_continuous(
         limits = c(0, 1),
         breaks = seq(0, 1, 0.2),
         expand = c(0, 0)
      ) +
      coord_fixed(ratio = 1) +
      labs(
         title = title_text,
         x = "1 - Specificity",
         y = "Sensitivity"
      ) +
      roc_theme +
      theme(
         legend.position = c(0.72, 0.22)
      )
   
   ## save publication-quality PDF
   ggsave(
      filename = out_pdf,
      plot = p,
      width = 6.5,
      height = 6.0,
      units = "in"
   )
   
   ## save AUC table
   auc_df <- data.frame(
      dataset = set_name,
      time = time_labels,
      days = valid_times,
      AUC = roc_obj$AUC,
      stringsAsFactors = FALSE
   )
   
   return(list(
      roc = roc_obj,
      auc = auc_df,
      plot = p
   ))
}

## -------- 3. batch output --------
roc_train <- plot_timeROC_sci(
   risk_train,
   "training",
   file.path(base_dir, "ROC", "ROC_train_1y_2y_3y_SCI.pdf"),
   roc_times = roc_times
)

roc_valid <- plot_timeROC_sci(
   risk_valid,
   "validation",
   file.path(base_dir, "ROC", "ROC_validation_1y_2y_3y_SCI.pdf"),
   roc_times = roc_times
)

roc_test <- plot_timeROC_sci(
   risk_test,
   "test",
   file.path(base_dir, "ROC", "ROC_test_1y_2y_3y_SCI.pdf"),
   roc_times = roc_times
)

## -------- 4. merge AUC summary --------
auc_all <- dplyr::bind_rows(
   if (!is.null(roc_train)) roc_train$auc,
   if (!is.null(roc_valid)) roc_valid$auc,
   if (!is.null(roc_test))  roc_test$auc
)

write.csv(
   auc_all,
   file.path(base_dir, "Tables", "timeROC_AUC_summary_SCI.csv"),
   row.names = FALSE
)

##############################
## 11. Subgroup KM plots (publication-ready)
##############################

subgroup_vars <- c("Gender", "Age", "iss_stage", "treatment_type")

## -------- 1. helper: safe filename --------
safe_filename <- function(x) {
   x <- as.character(x)
   x <- gsub("[[:space:]]+", "_", x)
   x <- gsub("[/\\|:*?\"<>]", "_", x)
   x <- gsub("[^A-Za-z0-9_.-]", "_", x)   # ← 修复这里
   x
}

## -------- 2. helper: subgroup label prettify --------
pretty_subgroup_label <- function(var, lv) {
   paste0(var, ": ", lv)
}

## -------- 3. main subgroup KM function --------
plot_km_subgroup_sci <- function(risk_df, dataset_name, subgroup_var, out_dir,
                                 min_n = 10, min_events = 3, min_group_n = 3) {
   
   ## basic check
   if (!subgroup_var %in% colnames(risk_df)) return(NULL)
   
   ## remove missing values in key columns
   sub_dat <- risk_df[
      !is.na(risk_df[[subgroup_var]]) &
         !is.na(risk_df$OS.time) &
         !is.na(risk_df$OS) &
         !is.na(risk_df$risk_group),
      , drop = FALSE
   ]
   
   if (nrow(sub_dat) == 0) return(NULL)
   
   sub_dat[[subgroup_var]] <- as.factor(sub_dat[[subgroup_var]])
   levs <- levels(sub_dat[[subgroup_var]])
   
   for (lv in levs) {
      df_lv <- sub_dat[sub_dat[[subgroup_var]] == lv, , drop = FALSE]
      df_lv$risk_group <- factor(df_lv$risk_group, levels = c("Low risk", "High risk"))
      
      ## -------- QC --------
      if (nrow(df_lv) < min_n) next
      if (length(unique(df_lv$risk_group)) < 2) next
      if (sum(df_lv$OS == 1, na.rm = TRUE) < min_events) next
      
      tab_group <- table(df_lv$risk_group)
      if (length(tab_group) < 2) next
      if (any(tab_group < min_group_n)) next
      
      ## -------- fit model safely --------
      fit_ok <- TRUE
      
      km_fit <- tryCatch(
         survfit(Surv(OS.time, OS) ~ risk_group, data = df_lv),
         error = function(e) { fit_ok <<- FALSE; NULL }
      )
      
      cox_fit <- tryCatch(
         coxph(Surv(OS.time, OS) ~ risk_group, data = df_lv),
         error = function(e) { fit_ok <<- FALSE; NULL }
      )
      
      survdiff_fit <- tryCatch(
         survdiff(Surv(OS.time, OS) ~ risk_group, data = df_lv),
         error = function(e) { fit_ok <<- FALSE; NULL }
      )
      
      if (!fit_ok || is.null(km_fit) || is.null(cox_fit) || is.null(survdiff_fit)) next
      
      ## -------- extract stats --------
      cox_sum <- tryCatch(summary(cox_fit), error = function(e) NULL)
      if (is.null(cox_sum)) next
      
      hr <- tryCatch(cox_sum$coefficients[1, "exp(coef)"], error = function(e) NA)
      hr_low <- tryCatch(cox_sum$conf.int[1, "lower .95"], error = function(e) NA)
      hr_high <- tryCatch(cox_sum$conf.int[1, "upper .95"], error = function(e) NA)
      
      p_logrank <- tryCatch(
         1 - pchisq(survdiff_fit$chisq, df = length(survdiff_fit$n) - 1),
         error = function(e) NA
      )
      
      if (is.na(hr) || is.na(hr_low) || is.na(hr_high)) next
      
      p_text <- if (is.na(p_logrank)) {
         "Log-rank P = NA"
      } else if (p_logrank < 0.001) {
         "Log-rank P < 0.001"
      } else {
         paste0("Log-rank P = ", formatC(p_logrank, format = "f", digits = 3))
      }
      
      hr_text <- sprintf("HR = %.2f (95%% CI: %.2f-%.2f)", hr, hr_low, hr_high)
      
      ## -------- sample size --------
      n_low <- sum(df_lv$risk_group == "Low risk", na.rm = TRUE)
      n_high <- sum(df_lv$risk_group == "High risk", na.rm = TRUE)
      
      legend_labels <- c(
         paste0("Low risk (n = ", n_low, ")"),
         paste0("High risk (n = ", n_high, ")")
      )
      
      ## -------- axis control --------
      max_time <- max(df_lv$OS.time, na.rm = TRUE)
      break_by <- signif(max_time / 5, 1)
      if (is.na(break_by) || break_by <= 0) break_by <- 500
      
      ## -------- output filename --------
      out_pdf <- file.path(
         out_dir,
         paste0(
            safe_filename(dataset_name), "_",
            safe_filename(subgroup_var), "_",
            safe_filename(lv),
            "_KM_SCI.pdf"
         )
      )
      
      ## -------- title --------
      title_text <- paste0(
         tools::toTitleCase(dataset_name), " cohort | ",
         pretty_subgroup_label(subgroup_var, lv)
      )
      
      ## -------- draw figure --------
      g <- ggsurvplot(
         fit               = km_fit,
         data              = df_lv,
         title             = title_text,
         xlab              = "Time (days)",
         ylab              = "Overall survival probability",
         palette           = c("#2F5597", "#C00000"),
         legend.title      = "Risk group",
         legend.labs       = legend_labels,
         risk.table        = TRUE,
         risk.table.col    = "strata",
         risk.table.y.text = FALSE,
         risk.table.height = 0.22,
         risk.table.fontsize = 3.2,
         risk.table.theme  = risk_table_theme,
         break.time.by     = break_by,
         xlim              = c(0, max_time * 1.02),
         conf.int          = FALSE,
         censor            = TRUE,
         censor.shape      = 124,
         censor.size       = 2.8,
         pval              = FALSE,
         surv.median.line  = "none",
         linetype          = "solid",
         size              = 1.1,
         ggtheme           = pub_theme
      )
      
      ## -------- annotation --------
      g$plot <- g$plot +
         scale_y_continuous(
            limits = c(0, 1),
            labels = scales::percent_format(accuracy = 1),
            expand = expansion(mult = c(0, 0.02))
         ) +
         annotate(
            "text",
            x = max_time * 0.52,
            y = 0.18,
            label = paste(hr_text, p_text, sep = "\n"),
            hjust = 0,
            size = 3.8
         ) +
         theme(
            legend.position = c(0.80, 0.84)
         )
      
      g$table <- g$table +
         theme(
            legend.position = "none"
         )
      
      ## -------- save pdf --------
      pdf(out_pdf, width = 7.2, height = 6.2, onefile = FALSE)
      print(g)
      dev.off()
   }
}

## -------- 4. batch output --------
for (v in subgroup_vars) {
   plot_km_subgroup_sci(
      risk_train, "training", v,
      file.path(base_dir, "Subgroup_KM")
   )
   plot_km_subgroup_sci(
      risk_valid, "validation", v,
      file.path(base_dir, "Subgroup_KM")
   )
   plot_km_subgroup_sci(
      risk_test, "test", v,
      file.path(base_dir, "Subgroup_KM")
   )
}


##############################
## 12. Subgroup Cox + interaction test
##############################
safe_cox_subgroup <- function(df, subgroup_var, level_name) {
   sub_df <- df[df[[subgroup_var]] == level_name, , drop = FALSE]
   
   if (nrow(sub_df) < 10) return(NULL)
   if (length(unique(sub_df$risk_group)) < 2) return(NULL)
   if (sum(sub_df$OS == 1, na.rm = TRUE) < 3) return(NULL)
   
   tab_group <- table(sub_df$risk_group)
   if (any(tab_group < 3)) return(NULL)
   
   fit <- tryCatch(
      coxph(Surv(OS.time, OS) ~ risk_group, data = sub_df),
      error = function(e) NULL
   )
   if (is.null(fit)) return(NULL)
   
   s <- summary(fit)
   ci <- tryCatch(confint(fit), error = function(e) NULL)
   if (is.null(ci)) return(NULL)
   
   data.frame(
      subgroup = subgroup_var,
      level    = as.character(level_name),
      n        = nrow(sub_df),
      events   = sum(sub_df$OS == 1, na.rm = TRUE),
      HR       = as.numeric(exp(coef(fit))[1]),
      lower    = as.numeric(exp(ci[1, 1])),
      upper    = as.numeric(exp(ci[1, 2])),
      pvalue   = as.numeric(s$coefficients[1, 5]),
      stringsAsFactors = FALSE
   )
}

get_interaction_p <- function(df, subgroup_var) {
   df2 <- df
   
   df2$risk_group <- as.factor(df2$risk_group)
   df2[[subgroup_var]] <- as.factor(df2[[subgroup_var]])
   
   ## 去掉缺失
   df2 <- df2[complete.cases(df2[, c("OS.time", "OS", "risk_group", subgroup_var)]), , drop = FALSE]
   
   ## 基本检查
   if (nrow(df2) < 10) return(NA_real_)
   if (length(unique(df2$risk_group)) < 2) return(NA_real_)
   if (length(unique(df2[[subgroup_var]])) < 2) return(NA_real_)
   if (sum(df2$OS == 1, na.rm = TRUE) < 5) return(NA_real_)
   
   ## 主效应模型
   fit_main <- tryCatch(
      coxph(
         as.formula(paste0("Surv(OS.time, OS) ~ risk_group + ", subgroup_var)),
         data = df2,
         x = TRUE
      ),
      error = function(e) NULL
   )
   
   ## 交互模型
   fit_int <- tryCatch(
      coxph(
         as.formula(paste0("Surv(OS.time, OS) ~ risk_group * ", subgroup_var)),
         data = df2,
         x = TRUE
      ),
      error = function(e) NULL
   )
   
   if (is.null(fit_main) || is.null(fit_int)) return(NA_real_)
   
   ## 用 logLik 做 LRT
   ll_main <- tryCatch(as.numeric(logLik(fit_main)), error = function(e) NA_real_)
   ll_int  <- tryCatch(as.numeric(logLik(fit_int)),  error = function(e) NA_real_)
   
   df_main <- tryCatch(attr(logLik(fit_main), "df"), error = function(e) NA_integer_)
   df_int  <- tryCatch(attr(logLik(fit_int),  "df"), error = function(e) NA_integer_)
   
   if (is.na(ll_main) || is.na(ll_int) || is.na(df_main) || is.na(df_int)) return(NA_real_)
   if (df_int <= df_main) return(NA_real_)
   
   chisq <- 2 * (ll_int - ll_main)
   df_diff <- df_int - df_main
   
   if (!is.finite(chisq) || !is.finite(df_diff) || df_diff <= 0) return(NA_real_)
   
   pval <- pchisq(chisq, df = df_diff, lower.tail = FALSE)
   
   if (length(pval) == 0 || is.na(pval) || !is.finite(pval)) return(NA_real_)
   
   return(as.numeric(pval))
}
make_subgroup_table <- function(risk_df, dataset_name, subgroup_vars) {
   res_list <- list()
   
   for (sg in subgroup_vars) {
      risk_df[[sg]] <- as.factor(risk_df[[sg]])
      levs <- levels(risk_df[[sg]])
      
      tmp_list <- lapply(levs, function(lv) safe_cox_subgroup(risk_df, sg, lv))
      tmp <- dplyr::bind_rows(tmp_list)
      
      if (nrow(tmp) > 0) {
         int_p <- get_interaction_p(risk_df, sg)
         if (length(int_p) == 0 || is.na(int_p) || !is.finite(int_p)) int_p <- NA_real_
         
         tmp$interaction_p <- NA_real_
         tmp$interaction_p[1] <- int_p
         tmp$dataset <- dataset_name
         
         res_list[[sg]] <- tmp
      }
   }
   
   dplyr::bind_rows(res_list)
}
sub_train <- make_subgroup_table(risk_train, "training", subgroup_vars)
sub_valid <- make_subgroup_table(risk_valid, "validation", subgroup_vars)
sub_test  <- make_subgroup_table(risk_test,  "test", subgroup_vars)

write.csv(sub_train, file.path(base_dir, "Tables", "subgroup_cox_training.csv"), row.names = FALSE)
write.csv(sub_valid, file.path(base_dir, "Tables", "subgroup_cox_validation.csv"), row.names = FALSE)
write.csv(sub_test,  file.path(base_dir, "Tables", "subgroup_cox_test.csv"), row.names = FALSE)

##############################
## 13. Forest plot
##############################
plot_forest <- function(sub_df, dataset_name, out_pdf) {
   if (nrow(sub_df) == 0) return(NULL)
   
   sub_df <- sub_df %>%
      mutate(
         subgroup = factor(subgroup, levels = rev(unique(subgroup))),
         label = paste0(level, " (n=", n, ", events=", events, ")"),
         p_txt = ifelse(is.na(pvalue), "",
                        ifelse(pvalue < 0.001, "<0.001", sprintf("%.3f", pvalue))),
         int_txt = ifelse(is.na(interaction_p), "",
                          ifelse(interaction_p < 0.001, "<0.001", sprintf("%.3f", interaction_p)))
      )
   
   ## 为了让每个 subgroup 内部按原顺序显示
   sub_df$display <- seq(nrow(sub_df), 1)
   
   p <- ggplot(sub_df, aes(x = HR, y = display)) +
      geom_vline(xintercept = 1, linetype = 2, color = "gray40") +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, size = 0.7) +
      geom_point(size = 2.4, shape = 18, color = "#C0392B") +
      scale_x_log10() +
      scale_y_continuous(
         breaks = sub_df$display,
         labels = paste0(sub_df$subgroup, " : ", sub_df$label)
      ) +
      theme_bw() +
      theme(
         panel.grid = element_blank(),
         axis.title.y = element_blank(),
         axis.text.y = element_text(size = 10, color = "black"),
         axis.text.x = element_text(size = 10, color = "black"),
         plot.title = element_text(hjust = 0.5, face = "bold"),
         panel.border = element_rect(color = "black", linewidth = 0.8)
      ) +
      labs(
         title = paste0("Subgroup analysis of RSF risk model (", dataset_name, ")"),
         x = "Hazard ratio for High risk vs Low risk (log scale)"
      )
   
   ## 右侧加文字
   x_max <- max(sub_df$upper, na.rm = TRUE) * 2.5
   x_hr  <- max(sub_df$upper, na.rm = TRUE) * 1.15
   x_p   <- max(sub_df$upper, na.rm = TRUE) * 1.65
   x_int <- max(sub_df$upper, na.rm = TRUE) * 2.20
   
   p <- p +
      geom_text(aes(x = x_hr, label = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper)),
                hjust = 0, size = 3.1) +
      geom_text(aes(x = x_p, label = p_txt), hjust = 0, size = 3.1) +
      geom_text(aes(x = x_int, label = int_txt), hjust = 0, size = 3.1) +
      annotate("text", x = x_hr,  y = max(sub_df$display) + 1.0, label = "HR (95% CI)", hjust = 0, size = 3.3, fontface = 2) +
      annotate("text", x = x_p,   y = max(sub_df$display) + 1.0, label = "P value", hjust = 0, size = 3.3, fontface = 2) +
      annotate("text", x = x_int, y = max(sub_df$display) + 1.0, label = "P for interaction", hjust = 0, size = 3.3, fontface = 2) +
      coord_cartesian(xlim = c(min(sub_df$lower, na.rm = TRUE) * 0.8, x_max), clip = "off") +
      theme(plot.margin = margin(10, 220, 10, 10))
   
   pdf(out_pdf, width = 11, height = max(5, 0.55 * nrow(sub_df) + 2))
   print(p)
   dev.off()
}

plot_forest(sub_train, "training",   file.path(base_dir, "ForestPlot", "forest_training.pdf"))
plot_forest(sub_valid, "validation", file.path(base_dir, "ForestPlot", "forest_validation.pdf"))
plot_forest(sub_test,  "test",       file.path(base_dir, "ForestPlot", "forest_test.pdf"))

##############################
## 14. Optional: combined table
##############################
sub_all <- bind_rows(sub_train, sub_valid, sub_test)
write.csv(sub_all, file.path(base_dir, "Tables", "subgroup_cox_all.csv"), row.names = FALSE)

##############################
## 15. Message
##############################
cat("All analyses finished.\n")
cat("Results saved in:\n", base_dir, "\n")

##############################
## 18. Calibration curve with 95% CI
## one figure per cohort: 1/2/3-year together
##############################
library(survival)
library(ggplot2)
library(dplyr)

dir.create(file.path(base_dir, "Calibration"), recursive = TRUE, showWarnings = FALSE)

## 1/2/3年
times_day <- c(365, 730, 1095)

## -------- 1. fit mapping on training set --------
cox_score <- coxph(
   Surv(OS.time, OS) ~ riskscore,
   data = risk_train,
   x = TRUE
)

## -------- 2. predict survival probability at time t --------
pred_surv_at_t <- function(cox_fit, newdata, t_day) {
   sf <- survfit(cox_fit, newdata = newdata)
   sm <- summary(sf, times = t_day)
   
   ## 若某些个体在该时间点无法返回，补成 NA
   s <- sm$surv
   as.numeric(s)
}

## -------- 3. safely extract KM estimate at time t --------
km_at_t <- function(dfg, t_day) {
   km <- survfit(Surv(OS.time, OS) ~ 1, data = dfg)
   sm <- summary(km, times = t_day)
   
   if (length(sm$surv) == 0) {
      return(c(surv = NA, lower = NA, upper = NA))
   }
   
   c(
      surv  = as.numeric(sm$surv),
      lower = as.numeric(sm$lower),
      upper = as.numeric(sm$upper)
   )
}

## -------- 4. build calibration data for one cohort at one time --------
get_calibration_df <- function(df, set_name, t_day, n_group = 5) {
   df <- df[, c("OS.time", "OS", "riskscore"), drop = FALSE]
   df <- df[complete.cases(df), , drop = FALSE]
   
   ## basic QC
   if (nrow(df) < 20) {
      message(set_name, " ", t_day, " days: too few samples.")
      return(NULL)
   }
   if (sum(df$OS == 1, na.rm = TRUE) < 5) {
      message(set_name, " ", t_day, " days: too few events.")
      return(NULL)
   }
   
   ## if follow-up too short, skip
   max_time <- max(df$OS.time, na.rm = TRUE)
   if (t_day >= max_time) {
      message(set_name, " ", t_day, " days: follow-up too short.")
      return(NULL)
   }
   
   ## predicted survival
   df$predS <- pred_surv_at_t(cox_score, df, t_day)
   df <- df[!is.na(df$predS), , drop = FALSE]
   
   if (nrow(df) < 20) {
      message(set_name, " ", t_day, " days: too few predictable samples.")
      return(NULL)
   }
   
   ## quantile grouping
   qs <- quantile(df$predS, probs = seq(0, 1, length.out = n_group + 1), na.rm = TRUE)
   qs[1] <- qs[1] - 1e-12
   qs[length(qs)] <- qs[length(qs)] + 1e-12
   qs <- unique(as.numeric(qs))
   
   if (length(qs) < 3) {
      message(set_name, " ", t_day, " days: prediction variation too low.")
      return(NULL)
   }
   
   df$grp <- cut(df$predS, breaks = qs, include.lowest = TRUE, labels = FALSE)
   
   ## mean predicted survival in each bin
   pred_mean <- tapply(df$predS, df$grp, mean, na.rm = TRUE)
   
   ## observed KM survival in each bin
   grp_ids <- sort(unique(df$grp))
   km_mat <- sapply(grp_ids, function(g) {
      dfg <- df[df$grp == g, , drop = FALSE]
      km_at_t(dfg, t_day)
   })
   
   plot_df <- data.frame(
      pred  = as.numeric(pred_mean),
      obs   = as.numeric(km_mat["surv", ]),
      lower = as.numeric(km_mat["lower", ]),
      upper = as.numeric(km_mat["upper", ])
   )
   
   plot_df <- plot_df[complete.cases(plot_df), , drop = FALSE]
   
   if (nrow(plot_df) < 2) {
      message(set_name, " ", t_day, " days: insufficient valid groups.")
      return(NULL)
   }
   
   plot_df$dataset <- set_name
   plot_df$time_day <- t_day
   plot_df$time_lab <- factor(
      paste0(round(t_day / 365), "-year"),
      levels = c("1-year", "2-year", "3-year")
   )
   
   plot_df
}

## -------- 5. publication theme --------
calib_theme <- theme_classic(base_size = 14) +
   theme(
      plot.title      = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.title      = element_text(size = 13, face = "bold", color = "black"),
      axis.text       = element_text(size = 11, color = "black"),
      axis.line       = element_line(linewidth = 0.8, color = "black"),
      axis.ticks      = element_line(linewidth = 0.7, color = "black"),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10.5, color = "black"),
      legend.position = c(0.80, 0.20),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin     = margin(8, 10, 8, 8)
   )

## -------- 6. draw one combined calibration plot for one cohort --------
plot_calibration_combined <- function(df, set_name, out_pdf, n_group = 5) {
   cal_list <- lapply(times_day, function(tt) {
      get_calibration_df(df, set_name, tt, n_group = n_group)
   })
   
   cal_all <- bind_rows(cal_list)
   
   if (nrow(cal_all) == 0) {
      message("Skip calibration plot for ", set_name, ": no valid calibration data.")
      return(NULL)
   }
   
   ## publication colors
   calib_cols <- c(
      "1-year" = "#C00000",
      "2-year" = "#ED7D31",
      "3-year" = "#2F5597"
   )
   
   ## only keep colors for available time points
   calib_cols <- calib_cols[levels(droplevels(cal_all$time_lab))]
   
   title_text <- paste0(tools::toTitleCase(set_name), " cohort")
   
   p <- ggplot(cal_all, aes(x = pred, y = obs, color = time_lab, group = time_lab)) +
      geom_abline(
         slope = 1, intercept = 0,
         linetype = "dashed", linewidth = 0.8, color = "gray55"
      ) +
      geom_errorbar(
         aes(ymin = lower, ymax = upper),
         width = 0.015,
         linewidth = 0.6,
         alpha = 0.85
      ) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.8, stroke = 0.3) +
      scale_color_manual(values = calib_cols, drop = FALSE) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      labs(
         title = title_text,
         x = "Predicted survival probability",
         y = "Observed survival probability"
      ) +
      calib_theme
   
   ggsave(
      filename = out_pdf,
      plot = p,
      width = 6.4,
      height = 6.0,
      units = "in"
   )
   
   return(cal_all)
}

## -------- 7. output one figure per cohort --------
cal_train <- plot_calibration_combined(
   risk_train,
   "training",
   file.path(base_dir, "Calibration", "calibration_training_1y_2y_3y_SCI.pdf"),
   n_group = 5
)

cal_valid <- plot_calibration_combined(
   risk_valid,
   "validation",
   file.path(base_dir, "Calibration", "calibration_validation_1y_2y_3y_SCI.pdf"),
   n_group = 5
)

cal_test <- plot_calibration_combined(
   risk_test,
   "test",
   file.path(base_dir, "Calibration", "calibration_test_1y_2y_3y_SCI.pdf"),
   n_group = 5
)

## -------- 8. save underlying calibration data --------
if (!is.null(cal_train)) {
   write.csv(
      cal_train,
      file.path(base_dir, "Tables", "calibration_training_1y_2y_3y_data.csv"),
      row.names = FALSE
   )
}

if (!is.null(cal_valid)) {
   write.csv(
      cal_valid,
      file.path(base_dir, "Tables", "calibration_validation_1y_2y_3y_data.csv"),
      row.names = FALSE
   )
}

if (!is.null(cal_test)) {
   write.csv(
      cal_test,
      file.path(base_dir, "Tables", "calibration_test_1y_2y_3y_data.csv"),
      row.names = FALSE
   )
}
##############################
## 19. Decision Curve Analysis (DCA)
##############################
library(rmda)

dir.create(file.path(base_dir, "DCA"), recursive = TRUE, showWarnings = FALSE)

## -------------------------
# Function to prepare dataset for DCA analysis
make_dca_data <- function(df, fit, u) {
   dat <- df[, c("OS.time", "OS", "riskscore"), drop = FALSE]
   dat <- dat[complete.cases(dat), , drop = FALSE]

   ## Predict probability of event at specified time point
   dat$pred <- get_pred_event_prob(fit, dat, u)

   ## Define fixed time endpoint label
   dat$event_u <- NA_integer_
   dat$event_u[dat$OS == 1 & dat$OS.time <= u] <- 1
   dat$event_u[dat$OS.time > u] <- 0

   ## Remove patients with indeterminate status: censored before time u
   dat <- dat[!is.na(dat$event_u), , drop = FALSE]

   dat
}

## -------------------------
## 19.2 DCA plotting wrapper function
## -------------------------
plot_dca <- function(df, fit, u, set_name, out_pdf) {
   dat <- make_dca_data(df, fit, u)

   if (nrow(dat) < 30) {
      message("Skip DCA for ", set_name, " at ", u, ": too few eligible samples.")
      return(NULL)
   }

   if (length(unique(dat$event_u)) < 2) {
      message("Skip DCA for ", set_name, " at ", u, ": only one outcome class.")
      return(NULL)
   }

   dca_fit <- decision_curve(
      event_u ~ pred,
      data = dat,
      family = binomial(link = "logit"),
      thresholds = seq(0.01, 0.80, by = 0.01),
      confidence.intervals = 0.95,
      study.design = "cohort",
      bootstraps = 200
   )

   pdf(out_pdf, width = 7, height = 8)
   plot_decision_curve(
      dca_fit,
      curve.names = paste0(set_name, " model"),
      cost.benefit.axis = FALSE,
      standardize = FALSE,
      confidence.intervals = FALSE,
      col = "#C0392B",
      lwd = 2
   )
   title(main = paste0(set_name, ": DCA at ", round(u / 365, 1), "-year"))
   dev.off()

   nb_df <- data.frame(
      threshold = dca_fit$derived.data$thresholds,
      net_benefit = dca_fit$derived.data$NB,
      model = set_name,
      time = u
   )

   return(nb_df)
}

## Training cohort DCA
dca_train_1 <- plot_dca(
   risk_train, cox_cal_fit, 365, "Training set",
   file.path(base_dir, "DCA", "DCA_train_1year.pdf")
)
dca_train_2 <- plot_dca(
   risk_train, cox_cal_fit, 730, "Training set",
   file.path(base_dir, "DCA", "DCA_train_2year.pdf")
)
dca_train_3 <- plot_dca(
   risk_train, cox_cal_fit, 1095, "Training set",
   file.path(base_dir, "DCA", "DCA_train_3year.pdf")
)

## Validation cohort DCA
dca_valid_1 <- plot_dca(
   risk_valid, cox_cal_fit, 365, "Validation set",
   file.path(base_dir, "DCA", "DCA_valid_1year.pdf")
)
dca_valid_2 <- plot_dca(
   risk_valid, cox_cal_fit, 730, "Validation set",
   file.path(base_dir, "DCA", "DCA_valid_2year.pdf")
)
dca_valid_3 <- plot_dca(
   risk_valid, cox_cal_fit, 1095, "Validation set",
   file.path(base_dir, "DCA", "DCA_valid_3year.pdf")
)

## Test cohort DCA
dca_test_1 <- plot_dca(
   risk_test, cox_cal_fit, 365, "Test set",
   file.path(base_dir, "DCA", "DCA_test_1year.pdf")
)
dca_test_2 <- plot_dca(
   risk_test, cox_cal_fit, 730, "Test set",
   file.path(base_dir, "DCA", "DCA_test_2year.pdf")
)
dca_test_3 <- plot_dca(
   risk_test, cox_cal_fit, 1095, "Test set",
   file.path(base_dir, "DCA", "DCA_test_3year.pdf")
)

## Export DCA net benefit data tables
save_dca_df <- function(x, name) {
   if (!is.null(x)) write.csv(x, file.path(base_dir, "Tables", paste0(name, ".csv")), row.names = FALSE)
}

save_dca_df(dca_train_1, "dca_train_1year")
save_dca_df(dca_train_2, "dca_train_2year")
save_dca_df(dca_train_3, "dca_train_3year")
save_dca_df(dca_valid_1, "dca_valid_1year")
save_dca_df(dca_valid_2, "dca_valid_2year")
save_dca_df(dca_valid_3, "dca_valid_3year")
save_dca_df(dca_test_1,  "dca_test_1year")
save_dca_df(dca_test_2,  "dca_test_2year")
save_dca_df(dca_test_3,  "dca_test_3year")

cat("Calibration and DCA analyses finished.\n")

##############################
## 20. Nomogram construction and visualization
##############################
library(survival)
library(rms)
library(ggplot2)

dir.create(file.path(base_dir, "Nomogram"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

## -------------------------
## 20.1 Prepare training dataset
## Notes:
## - Nomogram is built exclusively on the training cohort
## - Age remains as continuous variable
## - iss_stage / treatment_type converted to categorical factors
## -------------------------
nomo_df <- risk_train[, c("OS.time", "OS", "riskscore", "Age", "iss_stage", "treatment_type"), drop = FALSE]
nomo_df <- nomo_df[complete.cases(nomo_df), , drop = FALSE]

## Variable type conversion
nomo_df$Age <- as.numeric(nomo_df$Age)
nomo_df$iss_stage <- as.factor(nomo_df$iss_stage)
nomo_df$treatment_type <- as.factor(nomo_df$treatment_type)

## Drop empty factor levels
nomo_df$iss_stage <- droplevels(nomo_df$iss_stage)
nomo_df$treatment_type <- droplevels(nomo_df$treatment_type)

## Basic sample size and event count screening
if (nrow(nomo_df) < 30) stop("Nomogram training data too small.")
if (sum(nomo_df$OS == 1, na.rm = TRUE) < 10) stop("Too few events for nomogram model.")

## -------------------------
## 20.2 Initialize datadist object (required for rms package)
## -------------------------
dd <- datadist(nomo_df)
options(datadist = "dd")

## -------------------------
## 20.3 Fit Cox proportional hazards model
## Notes:
## - x=TRUE, y=TRUE, surv=TRUE required for nomogram and calibration curves
## - time.inc set to primary time endpoint of interest (365 days for 1-year survival)
## -------------------------
nomo_fit <- cph(
   Surv(OS.time, OS) ~ riskscore + Age + iss_stage + treatment_type,
   data = nomo_df,
   x = TRUE,
   y = TRUE,
   surv = TRUE,
   time.inc = 365
)

## Print model summary statistics
print(summary(nomo_fit))

## -------------------------
## 20.4 Extract hazard ratio table for manuscript reporting
## -------------------------
coef_df <- summary(nomo_fit)

hr_table <- data.frame(
   Variable = rownames(coef_df),
   HR = exp(coef(nomo_fit)),
   stringsAsFactors = FALSE
)

ci_mat <- tryCatch(confint(nomo_fit), error = function(e) NULL)
if (!is.null(ci_mat)) {
   hr_table$Lower95CI <- exp(ci_mat[, 1])
   hr_table$Upper95CI <- exp(ci_mat[, 2])
} else {
   hr_table$Lower95CI <- NA_real_
   hr_table$Upper95CI <- NA_real_
}

write.csv(
   hr_table,
   file.path(base_dir, "Tables", "nomogram_cox_HR_table.csv"),
   row.names = FALSE
)

## -------------------------
## 20.5 Generate survival prediction function
## Used to output 1/2/3-year overall survival probability scales in nomogram
## -------------------------
surv_fun <- Survival(nomo_fit)

nom <- nomogram(
   nomo_fit,
   fun = list(
      function(x) surv_fun(365, x),
      function(x) surv_fun(730, x),
      function(x) surv_fun(1095, x)
   ),
   funlabel = c("1-year Survival Probability",
                "2-year Survival Probability",
                "3-year Survival Probability"),
   lp = FALSE,
   maxscale = 100,
   fun.at = c(0.9, 0.7, 0.5, 0.3, 0.1)
)

## -------------------------
## 20.6 Render and export nomogram figure
## Notes:
## - rms::plot.nomogram is standard for publication-quality nomograms
## - Adjust cex / lmgp / tcl / xfrac for compact, clean layout
## -------------------------
pdf(
   file.path(base_dir, "Nomogram", "nomogram_training.pdf"),
   width = 8,
   height = 6,
   family = "Helvetica"
)

par(
   mar = c(3.5, 2.5, 2.5, 1.5),
   mgp = c(2.0, 0.7, 0),
   tcl = -0.25,
   xpd = NA
)

plot(
   nom,
   xfrac = 0.42,       ## Width allocated for variable label panel on left
   cex.var = 1.05,     ## Font size for variable names
   cex.axis = 0.9,     ## Font size for axis tick labels
   lmgp = 0.25,        ## Spacing between axis labels and tick lines
   col.grid = gray(c(0.85, 0.95)),
   col.conf = "#C00000",
   conf.int = FALSE
)

title(
   main = "Nomogram for Overall Survival Prediction",
   cex.main = 1.2,
   font.main = 2
)

dev.off()

## -------------------------
## 20.7 Optional: Export linear predictor values and predicted survival probabilities
## For downstream subgroup analysis or supplementary tables
## -------------------------
nomo_df$lp <- predict(nomo_fit, type = "lp")
nomo_df$pred_1y <- surv_fun(365, nomo_df$lp)
nomo_df$pred_2y <- surv_fun(730, nomo_df$lp)
nomo_df$pred_3y <- surv_fun(1095, nomo_df$lp)

write.csv(
   nomo_df,
   file.path(base_dir, "Tables", "nomogram_training_predictions.csv"),
   row.names = FALSE
)

cat("Nomogram analysis finished.\n")




##############################
## 22. Export RSF model and full analysis workspace
##############################

dir.create(file.path(base_dir, "Model"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "Workspace"), recursive = TRUE, showWarnings = FALSE)

## 1) Save trained random survival forest model object
saveRDS(
   rsf_fit,
   file = file.path(base_dir,  "rsf_fit.rds")
)

## 2) Save training set normalization parameters for external dataset standardization
norm_info <- list(
   gene_cols = gene_cols,
   mu = mu,
   sdv = sdv,
   cutoff = cutoff
)

saveRDS(
   norm_info,
   file = file.path(base_dir,  "rsf_norm_info.rds")
)

## 3) Save complete R workspace environment
save.image(
   file = file.path(base_dir,  "RSF_full_workspace.RData")
)
load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis/RSF_full_workspace.RData")


cat("RSF model and full workspace have been saved.\n")
cat("Model file: ", file.path(base_dir, "Model", "rsf_fit.rds"), "\n")
cat("Normalization info: ", file.path(base_dir, "Model", "rsf_norm_info.rds"), "\n")
cat("Workspace file: ", file.path(base_dir, "Workspace", "RSF_full_workspace.RData"), "\n")



#Model_comparison



