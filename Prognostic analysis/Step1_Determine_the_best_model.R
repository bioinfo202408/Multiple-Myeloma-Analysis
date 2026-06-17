library(survival)
library(randomForestSRC)
library(glmnet)
library(superpc)
library(plsRcox)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(tidyr)

options(stringsAsFactors = FALSE)

## =========================
## 0) Data import & standardization
## =========================
training <- read.table("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/training_sets_symbol_clinic.csv",
                       header = TRUE, sep = ",", quote = "", check.names = FALSE)
valid <- read.table("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/valid_sets_symbol.csv",
                    header = TRUE, sep = ",", quote = "", check.names = FALSE)
test <- read.table("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/test_sets_symbol.csv",
                   header = TRUE, sep = ",", quote = "", check.names = FALSE)

# Feature columns: the first 3 columns are (Sample ID/OS.time/OS) or clinical + OS metadata, original logic excludes columns 1-6
x_cols <- colnames(training)[-c(1:6)]

# Standardize features using mean and SD calculated from training set
training_mean <- apply(training[, x_cols, drop = FALSE], 2, mean)
training_sd   <- apply(training[, x_cols, drop = FALSE], 2, sd)

training[, x_cols] <- scale(training[, x_cols, drop = FALSE])
valid[, x_cols] <- sweep(sweep(valid[, x_cols, drop = FALSE], 2, training_mean, "-"), 2, training_sd, "/")
test[, x_cols]  <- sweep(sweep(test[, x_cols, drop = FALSE], 2, training_mean, "-"), 2, training_sd, "/")

# Replace NaN values with zero
training[is.na(training)] <- 0
valid[is.na(valid)] <- 0
test[is.na(test)] <- 0

mm <- list(training = training, valid = valid, test = test)

est_data <- mm$training
val_data_list <- mm

pre_var <- colnames(est_data)[-c(1:6)]
est_dd <- est_data[, c("OS.time", "OS", pre_var)]
val_dd_list <- lapply(val_data_list, function(x) x[, c("OS.time", "OS", pre_var)])

rf_nodesize <- 5
seed <- 123

result <- data.frame()

## =========================
## 1) Helper function: Calculate C-index & append results
## =========================
get_cindex_df <- function(rs_list, model_name) {
   cc <- data.frame(
      Cindex = sapply(rs_list, function(d) {
         as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, data = d))$concordance[1])
      })
   ) %>% rownames_to_column("ID")
   cc$Model <- model_name
   cc
}

## =========================
## 2) RSF: Single training run + reuse variable importance for feature selection
## =========================
fit_rsf <- function(train_df, seed, nodesize, ntree = 1000) {
   set.seed(seed)
   rfsrc(
      Surv(OS.time, OS) ~ .,
      data = train_df,
      ntree = ntree,
      nodesize = nodesize,
      splitrule = "logrank",
      importance = TRUE,
      proximity = TRUE,
      forest = TRUE,
      seed = seed
   )
}

rs_from_rsf <- function(fit, data_list) {
   lapply(data_list, function(x) {
      cbind(x[, 1:2], RS = predict(fit, newdata = x, na.action = "na.impute")$predicted)
   })
}

select_by_rsf_vimp <- function(rsf_fit, min_imp = 0) {
   vi <- vimp(rsf_fit, importance = "permute")$importance
   rid <- names(sort(vi, decreasing = TRUE))
   rid <- rid[vi[rid] > min_imp]
   rid
}

## =========================
## 3) glmnet Cox: Elastic Net / Lasso / Ridge (standard Cox implementation)
## =========================
fit_glmnet_cox_cv <- function(train_df, feat_cols, alpha, seed, nfolds = 10) {
   x <- as.matrix(train_df[, feat_cols, drop = FALSE])
   y <- Surv(train_df$OS.time, train_df$OS)
   set.seed(seed)
   cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = nfolds)
}

rs_from_glmnet_cox <- function(cvfit, data_list) {
   lapply(data_list, function(x) {
      lp <- as.numeric(predict(cvfit,
                               type = "link",
                               newx = as.matrix(x[, -c(1,2), drop = FALSE]),
                               s = cvfit$lambda.min))
      cbind(x[, 1:2], RS = lp)
   })
}

## =========================
## 4) Single RSF model
## =========================
fit0 <- fit_rsf(est_dd, seed = seed, nodesize = rf_nodesize)
rs0  <- rs_from_rsf(fit0, val_dd_list)
result <- rbind(result, get_cindex_df(rs0, "RSF"))

# Feature screening via RSF variable importance (run once, reused for downstream combined models)
rid_rsf <- select_by_rsf_vimp(fit0, min_imp = 0)
est_dd_rsf <- est_data[, c("OS.time", "OS", rid_rsf)]
val_dd_list_rsf <- lapply(val_data_list, function(x) x[, c("OS.time", "OS", rid_rsf)])

## =========================
## 5) RSF + CoxBoost combined model
## =========================
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd_rsf$OS.time, est_dd_rsf$OS, as.matrix(est_dd_rsf[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = TRUE)
cv.res <- cv.CoxBoost(est_dd_rsf$OS.time, est_dd_rsf$OS, as.matrix(est_dd_rsf[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit_cb <- CoxBoost(est_dd_rsf$OS.time, est_dd_rsf$OS, as.matrix(est_dd_rsf[, -c(1,2)]),
                   stepno = cv.res$optimal.step, penalty = pen$penalty)

rs <- lapply(val_dd_list_rsf, function(x) {
   cbind(x[, 1:2],
         RS = as.numeric(predict(fit_cb,
                                 newdata = x[, -c(1,2), drop = FALSE],
                                 newtime = x[, 1],
                                 newstatus = x[, 2],
                                 type = "lp")))
})
result <- rbind(result, get_cindex_df(rs, "RSF + CoxBoost"))

## =========================
## 6) RSF + Elastic Net / Lasso / Ridge (standard Cox regression)
## =========================
for (alpha in seq(0.1, 0.9, 0.1)) {
   fit_en <- fit_glmnet_cox_cv(est_dd_rsf, rid_rsf, alpha = alpha, seed = seed, nfolds = 10)
   rs <- rs_from_glmnet_cox(fit_en, val_dd_list_rsf)
   result <- rbind(result, get_cindex_df(rs, paste0("RSF + Enet[α=", alpha, "]")))
}

# RSF + Lasso (alpha=1)
fit_las <- fit_glmnet_cox_cv(est_dd_rsf, rid_rsf, alpha = 1, seed = seed, nfolds = 10)
rs <- rs_from_glmnet_cox(fit_las, val_dd_list_rsf)
result <- rbind(result, get_cindex_df(rs, "RSF + Lasso"))

# RSF + Ridge (alpha=0)
fit_rid <- fit_glmnet_cox_cv(est_dd_rsf, rid_rsf, alpha = 0, seed = seed, nfolds = 10)
rs <- rs_from_glmnet_cox(fit_rid, val_dd_list_rsf)
result <- rbind(result, get_cindex_df(rs, "RSF + Ridge"))

## =========================
## 7) RSF + GBM combined model
## =========================
set.seed(seed)
fit_gbm0 <- gbm(Surv(OS.time, OS) ~ ., data = est_dd_rsf, distribution = "coxph",
                n.trees = 10000, interaction.depth = 3, n.minobsinnode = 10,
                shrinkage = 0.001, cv.folds = 10, n.cores = 6)
best <- which.min(fit_gbm0$cv.error)

set.seed(seed)
fit_gbm <- gbm(Surv(OS.time, OS) ~ ., data = est_dd_rsf, distribution = "coxph",
               n.trees = best, interaction.depth = 3, n.minobsinnode = 10,
               shrinkage = 0.001, cv.folds = 10, n.cores = 8)

rs <- lapply(val_dd_list_rsf, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_gbm, x, n.trees = best, type = "link")))
})
result <- rbind(result, get_cindex_df(rs, "RSF + GBM"))

## =========================
## 8) RSF + plsRcox combined model
## =========================
set.seed(seed)
pdf("/home/yjliu/mmProj/data_process/Human/Prognostic_models/model/cv_plsRcox_plot.pdf", width = 10, height = 8)
cv_pls <- cv.plsRcox(list(x = est_dd_rsf[, rid_rsf, drop = FALSE],
                          time = est_dd_rsf$OS.time,
                          status = est_dd_rsf$OS),
                     nt = 10, verbose = FALSE)
dev.off()
fit_pls <- plsRcox(est_dd_rsf[, rid_rsf, drop = FALSE],
                   time = est_dd_rsf$OS.time, event = est_dd_rsf$OS,
                   nt = as.numeric(cv_pls[5]))
rs <- lapply(val_dd_list_rsf, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_pls, type = "lp", newdata = x[, -c(1,2), drop = FALSE])))
})
result <- rbind(result, get_cindex_df(rs, "RSF + plsRcox"))

## =========================
## 9) RSF + SuperPC combined model
## =========================
sp_data <- list(x = t(est_dd_rsf[, -c(1,2)]),
                y = est_dd_rsf$OS.time,
                censoring.status = est_dd_rsf$OS,
                featurenames = colnames(est_dd_rsf)[-c(1,2)])
set.seed(seed)
fit_sp <- superpc.train(data = sp_data, type = "survival", s0.perc = 0.5)
cv_sp <- superpc.cv(fit_sp, sp_data, n.threshold = 20, n.fold = 10,
                    n.components = 3, min.features = 5, max.features = nrow(sp_data$x),
                    compute.fullcv = TRUE, compute.preval = TRUE)

rs <- lapply(val_dd_list_rsf, function(w) {
   sp_test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS,
                   featurenames = colnames(w)[-c(1,2)])
   th <- cv_sp$thresholds[which.max(cv_sp[["scor"]][1, ])]
   pred <- superpc.predict(fit_sp, sp_data, sp_test, threshold = th, n.components = 1)
   cbind(w[, 1:2], RS = as.numeric(pred$v.pred))
})
result <- rbind(result, get_cindex_df(rs, "RSF + SuperPC"))

## =========================
## 10) RSF + survival-SVM combined model
## =========================
fit_svm <- survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd_rsf, gamma.mu = 1)
rs <- lapply(val_dd_list_rsf, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_svm, x)$predicted))
})
result <- rbind(result, get_cindex_df(rs, "RSF + survival-SVM"))

## =========================
## 11) Single standalone models: Elastic Net / CoxBoost / plsRcox / SuperPC / GBM / survival-SVM / Lasso / Ridge
##     All Cox regression standard implementations (Lasso/Ridge use family="cox")
## =========================

## 11-1 Elastic Net (all raw features)
for (alpha in seq(0.1, 0.9, 0.1)) {
   fit_en <- fit_glmnet_cox_cv(est_dd, pre_var, alpha = alpha, seed = seed, nfolds = 10)
   rs <- rs_from_glmnet_cox(fit_en, val_dd_list)
   result <- rbind(result, get_cindex_df(rs, paste0("Enet[α=", alpha, "]")))
}

## 11-2 CoxBoost (all raw features)
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd$OS.time, est_dd$OS, as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = TRUE)
cv.res <- cv.CoxBoost(est_dd$OS.time, est_dd$OS, as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit_cb <- CoxBoost(est_dd$OS.time, est_dd$OS, as.matrix(est_dd[, -c(1,2)]),
                   stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x) {
   cbind(x[, 1:2],
         RS = as.numeric(predict(fit_cb, newdata = x[, -c(1,2), drop = FALSE],
                                 newtime = x[,1], newstatus = x[,2], type = "lp")))
})
result <- rbind(result, get_cindex_df(rs, "CoxBoost"))

## 11-3 plsRcox (all raw features)
set.seed(seed)
pdf("/home/yjliu/mmProj/data_process/Human/Prognostic_models/model/cv_plsRcox_plot2.pdf", width = 10, height = 8)
cv_pls <- cv.plsRcox(list(x = est_dd[, pre_var, drop = FALSE],
                          time = est_dd$OS.time,
                          status = est_dd$OS),
                     nt = 10, verbose = FALSE)
dev.off()
fit_pls <- plsRcox(est_dd[, pre_var, drop = FALSE], time = est_dd$OS.time, event = est_dd$OS,
                   nt = as.numeric(cv_pls[5]))
rs <- lapply(val_dd_list, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_pls, type = "lp", newdata = x[, -c(1,2), drop = FALSE])))
})
result <- rbind(result, get_cindex_df(rs, "plsRcox"))

## 11-4 SuperPC (all raw features)
sp_data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time,
                censoring.status = est_dd$OS,
                featurenames = colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit_sp <- superpc.train(data = sp_data, type = "survival", s0.perc = 0.5)
cv_sp <- superpc.cv(fit_sp, sp_data, n.threshold = 20, n.fold = 10,
                    n.components = 3, min.features = 5, max.features = nrow(sp_data$x),
                    compute.fullcv = TRUE, compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w) {
   sp_test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS,
                   featurenames = colnames(w)[-c(1,2)])
   th <- cv_sp$thresholds[which.max(cv_sp[["scor"]][1, ])]
   pred <- superpc.predict(fit_sp, sp_data, sp_test, threshold = th, n.components = 1)
   cbind(w[, 1:2], RS = as.numeric(pred$v.pred))
})
result <- rbind(result, get_cindex_df(rs, "SuperPC"))

## 11-5 GBM (all raw features)
set.seed(seed)
fit_gbm0 <- gbm(Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph",
                n.trees = 10000, interaction.depth = 3, n.minobsinnode = 10,
                shrinkage = 0.001, cv.folds = 10, n.cores = 6)
best <- which.min(fit_gbm0$cv.error)
set.seed(seed)
fit_gbm <- gbm(Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph",
               n.trees = best, interaction.depth = 3, n.minobsinnode = 10,
               shrinkage = 0.001, cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_gbm, x, n.trees = best, type = "link")))
})
result <- rbind(result, get_cindex_df(rs, "GBM"))

## 11-6 survival-SVM (all raw features)
fit_svm <- survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x) {
   cbind(x[, 1:2], RS = as.numeric(predict(fit_svm, x)$predicted))
})
result <- rbind(result, get_cindex_df(rs, "survival-SVM"))

## 11-7 Ridge regression (all raw features, standard Cox implementation)
fit_rid <- fit_glmnet_cox_cv(est_dd, pre_var, alpha = 0, seed = seed, nfolds = 10)
rs <- rs_from_glmnet_cox(fit_rid, val_dd_list)
result <- rbind(result, get_cindex_df(rs, "Ridge"))

## 11-8 Lasso regression (all raw features, standard Cox implementation)
fit_las <- fit_glmnet_cox_cv(est_dd, pre_var, alpha = 1, seed = seed, nfolds = 10)
rs <- rs_from_glmnet_cox(fit_las, val_dd_list)
result <- rbind(result, get_cindex_df(rs, "Lasso"))

## =========================
## 12) Aggregate and export output (retain original pipeline logic)
## =========================
result2 <- result %>%
   group_by(Model, ID) %>%
   summarise(Cindex = mean(Cindex), .groups = "drop") %>%
   mutate(ID = factor(ID, levels = c("training", "valid", "test"))) %>%
   arrange(Model, ID)

dd2 <- pivot_wider(result2, names_from = "ID", values_from = "Cindex") %>% as.data.frame()
dd2[, -1] <- apply(dd2[, -1, drop = FALSE], 2, as.numeric)

dd2$All <- apply(dd2[, 2:4, drop = FALSE], 1, mean)
dd2$GEO <- apply(dd2[, 3:4, drop = FALSE], 1, mean)

write.table(dd2,
            "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/output_C_index.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

save.image(file = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/model_Select.RData")
load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/model_Select.RData")
## Figure A - Heatmap of C-indices

# Sort rows by average C-index value


# # Only plot heatmap for C-indices of GEO validation cohorts
# # Prepare heatmap matrix using only valid and test dataset C-index values
# dt <- dd2[, 3:4]
# rownames(dt) <- dd2$Model
# library(circlize)
# library(ComplexHeatmap)
# 
# 
# col_ha <- HeatmapAnnotation(which = "col", Cohort = c("valid","test"),
#                             annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
#                             annotation_name_side = "left",
#                             col = list(Cohort=c("valid"="#00A087B2",
#                                                 "test"="#3C5488B2")),
#                             annotation_legend_param = list(Cohort=list(title = "Cohort",
#                                                                        title_position = "topleft",
#                                                                        title_gp = gpar(fontsize = 12, fontface = "bold"),
#                                                                        labels_rot = 0,
#                                                                        legend_height = unit(1,"cm"),
#                                                                        legend_width = unit(5,"mm"),
#                                                                        labels_gp = gpar(fontsize = 9,
#                                                                                         fontface = "bold"))
#                             )
# )
# # Row-side annotation
# row_ha <- rowAnnotation('Mean Cindex' = anno_barplot(round(rowMeans(dt), 3), bar_width = 1, add_numbers = T,
#                                                      labels = c("Mean Cindex"), height = unit(1, "mm"),
#                                                      gp = gpar(col = "white", fill = "skyblue1"), numbers_gp = gpar(fontsize = 8),
#                                                      axis_param = list(at = c(0, 0.5, 1),
#                                                                        labels = c("0", "0.5", "1")),
#                                                      width = unit(2.5, "cm")),
#                         annotation_name_side = "bottom",
#                         annotation_name_gp = gpar(fontsize = 9, fontface = "bold", angle = 90))
# 
# # Custom cell rendering function to display numeric C-index inside heatmap cells
# cell_fun <- function(j, i, x, y, width, height, fill) {
#    grid.text(
#       round(dt[i, j], 2), 
#       x, y,
#       gp = gpar(
#          fontsize = 8
#       ))
# }
# 
# # Render heatmap figure
# 
# 
# pdf("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/ComplexHeatmap.pdf", width = 10, height = 18)
# heatmap <- Heatmap(dt,name = " ",
#                    heatmap_legend_param = list(title="",title_position = "topleft", labels_rot = 0,
#                                                legend_height = unit(8,"cm"),
#                                                legend_width = unit(5,"mm"),
#                                                labels_gp = gpar(fontsize = 15, fontface = "bold")),
#                    border = TRUE,
#                    column_split = c("valid","test"),
#                    column_gap = unit(3, "mm"),
#                    show_column_names = F,
#                    show_row_names = T,
#                    col = colorRamp2(c(0.64,0.7,0.76), c("#4DBBD5B2", "white", "#E64B35B2")), 
#                    column_title ="", 
#                    column_title_side = "top",
#                    row_title_side = "left",
#                    row_title_rot = 90, 
#                    column_title_gp = gpar(fontsize = 12, fontface = "bold",col = "black"), 
#                    cluster_columns =F,
#                    cluster_rows = F,
#                    column_order = c(colnames(dt)),
#                    show_row_dend = F, 
#                    cell_fun = cell_fun,
#                    top_annotation = col_ha,
#                    right_annotation = row_ha
# )
# print(heatmap)
# dev.off()

dd2 <- dd2[order(dd2$GEO, decreasing = T),]
## ===== 1) Data preparation: extract three cohorts, ordered as training / test / valid =====
library(circlize)
library(ComplexHeatmap)
library(grid)

dt <- dd2[, 2:4]
rownames(dt) <- dd2$Model
dt <- as.matrix(dt)

# Reorder columns to specified sequence (ensure column names match)
dt <- dt[, c("training", "test", "valid"), drop = FALSE]
cohort_vec <- colnames(dt)

## ===== 2) Sort rows in descending order of average C-index =====
dt <- dt[order(rowMeans(dt), decreasing = TRUE), , drop = FALSE]

## ===== 3) Top column annotation =====
col_ha <- HeatmapAnnotation(
   which = "col",
   Cohort = cohort_vec,
   annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
   annotation_name_side = "left",
   col = list(
      Cohort = c(
         "training" = "#F4A3A3",
         "valid"    = "#00A087B2",
         "test"     = "#3C5488B2"
      )
   )
)

## ===== 4) Right row annotation: barplot for mean C-index per model =====
row_ha <- rowAnnotation(
   "Mean Cindex" = anno_barplot(
      round(rowMeans(dt), 3),
      bar_width = 1,
      add_numbers = TRUE,
      gp = gpar(col = "white", fill = "skyblue1"),
      numbers_gp = gpar(fontsize = 8),
      axis_param = list(at = c(0, 0.5, 1), labels = c("0", "0.5", "1")),
      width = unit(2.5, "cm")
   ),
   annotation_name_side = "bottom",
   annotation_name_gp = gpar(fontsize = 9, fontface = "bold", angle = 90)
)

## ===== 5) Function to print C-index values inside heatmap cells =====
cell_fun <- function(j, i, x, y, width, height, fill) {
   grid.text(sprintf("%.2f", dt[i, j]), x, y, gp = gpar(fontsize = 8))
}

## ===== 6) Export heatmap PDF =====
pdf("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/ComplexHeatmap.pdf", width = 10, height = 18)

ht <- Heatmap(
   dt,
   name = " ",
   border = TRUE,
   column_split = cohort_vec,
   column_gap = unit(3, "mm"),
   show_column_names = FALSE,
   show_row_names = TRUE,
   col = colorRamp2(c(0.65, 0.7, 0.9), c("#4DBBD5B2", "white", "#E64B35B2")),
   cluster_columns = FALSE,
   cluster_rows = FALSE,
   show_row_dend = FALSE,
   cell_fun = cell_fun,
   top_annotation = col_ha,
   right_annotation = row_ha
)

draw(ht)
dev.off()

dd2 <- dd2[order(dd2$GEO, decreasing = TRUE), ]

## ===== 1) Data preparation: extract three cohorts, ordered as training / test / valid =====
library(circlize)
library(ComplexHeatmap)
library(grid)

dt <- dd2[, 2:4]
rownames(dt) <- dd2$Model
dt <- as.matrix(dt)

# Reorder columns to target sequence
dt <- dt[, c("training", "test", "valid"), drop = FALSE]
cohort_vec <- colnames(dt)

## ===== 2) Sort rows by average C-index in descending order =====
dt <- dt[order(rowMeans(dt), decreasing = TRUE), , drop = FALSE]

## ===== 3) Column top annotation settings =====
col_ha <- HeatmapAnnotation(
   which = "col",
   Cohort = cohort_vec,
   annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
   annotation_name_side = "left",
   simple_anno_size = unit(0.7, "cm"),
   col = list(
      Cohort = c(
         "training" = "#F4A3A3",
         "valid"    = "#00A087B2",
         "test"     = "#3C5488B2"
      )
   ),
   gap = unit(2, "mm")
)

## ===== 4) Right row annotation: mean C-index bar chart =====
row_ha <- rowAnnotation(
   "Mean Cindex" = anno_barplot(
      round(rowMeans(dt), 3),
      bar_width = 0.9,
      add_numbers = TRUE,
      gp = gpar(col = NA, fill = "skyblue1"),
      border = FALSE,
      numbers_gp = gpar(fontsize = 10, fontface = "bold"),
      axis_param = list(
         at = c(0, 0.5, 1),
         labels = c("0", "0.5", "1"),
         gp = gpar(fontsize = 10, fontface = "bold")
      ),
      width = unit(3.2, "cm")
   ),
   annotation_name_side = "bottom",
   annotation_name_gp = gpar(fontsize = 12, fontface = "bold", rot = 90),
   gap = unit(2, "mm")
)

## ===== 5) Render C-index numeric labels within heatmap cells =====
cell_fun <- function(j, i, x, y, width, height, fill) {
   grid.text(
      sprintf("%.2f", dt[i, j]),
      x, y,
      gp = gpar(fontsize = 10, fontface = "bold", col = "black")
   )
}

## ===== 6) Export publication-quality heatmap PDF =====
pdf(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/ComplexHeatmap_publication.pdf",
   width = 11,
   height = 20
)

ht <- Heatmap(
   dt,
   name = "C-index",
   border = TRUE,
   rect_gp = gpar(col = "grey70", lwd = 1),
   column_split = cohort_vec,
   column_gap = unit(4, "mm"),
   row_gap = unit(1.5, "mm"),
   
   show_column_names = FALSE,
   show_row_names = TRUE,
   row_names_side = "left",
   row_names_gp = gpar(fontsize = 11, fontface = "bold"),
   
   col = colorRamp2(
      c(0.65, 0.7, 0.9),
      c("#4DBBD5B2", "white", "#E64B35B2")
   ),
   
   cluster_columns = FALSE,
   cluster_rows = FALSE,
   show_row_dend = FALSE,
   show_column_dend = FALSE,
   
   cell_fun = cell_fun,
   top_annotation = col_ha,
   right_annotation = row_ha,
   
   heatmap_legend_param = list(
      title = "C-index",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10, fontface = "bold"),
      border = "grey40",
      legend_height = unit(4.5, "cm"),
      grid_width = unit(0.6, "cm")
   )
)

draw(
   ht,
   heatmap_legend_side = "right",
   annotation_legend_side = "right",
   padding = unit(c(8, 8, 8, 8), "mm")
)

dev.off()