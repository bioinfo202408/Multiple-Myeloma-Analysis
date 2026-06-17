# ========== Load required packages ==========
library(xgboost)
library(pROC)
library(dplyr)
library(vroom)
library(ggplot2)
library(tidyr)
library(PRROC)
library(caret)

# ========== Data Loading ==========
train_data <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data_1_0.csv", 
   delim = ",", 
)
train_data <- as.data.frame(train_data)
rownames(train_data) <- make.unique(train_data[[1]])
train_data <- train_data[,-1]
# Ensure group is factor, 1=tumor, 0=health
train_data$group <- factor(train_data$group, levels = c(0, 1))

test_data <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/test_data_1_0.csv", 
   delim = ",", 
)
test_data <- as.data.frame(test_data)
rownames(test_data) <- make.unique(test_data[[1]])
test_data <- test_data[,-1]
test_data$group <- factor(test_data$group, levels = c(0, 1))

external_data1 <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/external1_1_0.csv", 
   delim = ",", 
)
external_data1 <- as.data.frame(external_data1)
rownames(external_data1) <- make.unique(external_data1[[1]])
external_data1 <- external_data1[,-1]
external_data1$group <- factor(external_data1$group, levels = c(0, 1))

external_data2 <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/external2_1_0.csv", 
   delim = ",", 
)
external_data2 <- as.data.frame(external_data2)
rownames(external_data2) <- make.unique(external_data2[[1]])
external_data2 <- external_data2[,-1]
external_data2$group <- factor(external_data2$group, levels = c(0, 1))

# ========== Prepare Data for XGBoost ==========
prepare_xgboost_data <- function(df, label_col = "group") {
   feature_cols <- setdiff(colnames(df), label_col)
   x <- as.matrix(df[, feature_cols])
   y <- as.numeric(df[[label_col]]) - 1  # Convert label to 0/1 (0=health, 1=tumor)
   return(list(x = x, y = y))
}

# Prepare training dataset
train_prepared <- prepare_xgboost_data(train_data)
dtrain <- xgb.DMatrix(data = train_prepared$x, label = train_prepared$y)

# Prepare test and external validation datasets
test_prepared <- prepare_xgboost_data(test_data)
dtest <- xgb.DMatrix(data = test_prepared$x, label = test_prepared$y)

external_prepared1 <- prepare_xgboost_data(external_data1)
dexternal1 <- xgb.DMatrix(data = external_prepared1$x, label = external_prepared1$y)

external_prepared2 <- prepare_xgboost_data(external_data2)
dexternal2 <- xgb.DMatrix(data = external_prepared2$x, label = external_prepared2$y)

# ========== Set XGBoost Hyperparameters ==========
set.seed(123)

xgb_params <- list(
   objective = "binary:logistic",
   eval_metric = "logloss",
   eta = 0.05,
   max_depth = 6,
   min_child_weight = 1,
   subsample = 0.8,
   colsample_bytree = 0.8,
   lambda = 1,
   alpha = 0,
   seed = 123
)

# ========== 1. 10-fold cross validation on training set (find optimal iteration rounds) ==========
cv_result <- xgb.cv(
   params = xgb_params,
   data = dtrain,
   nrounds = 1000,
   nfold = 10,
   early_stopping_rounds = 50,
   verbose = 100,
   showsd = TRUE,
   stratified = TRUE,
   print_every_n = 100
)

# Extract optimal iteration rounds
best_iter <- cv_result$best_iteration
cat("✅ Optimal rounds from 10-fold CV =", best_iter, "\n")

# ========== 2. Train final model with optimal iteration rounds ==========
xgb_model_full <- xgb.train(
   params = xgb_params,
   data = dtrain,
   nrounds = best_iter,
   verbose = 100
)

cat("✅ XGBoost model trained using optimal rounds from 10-fold CV\n")

# ========== Validate model on test set and calculate AUC ==========
pred_prob_test <- predict(xgb_model_full, dtest)
roc_test <- roc(test_data$group, pred_prob_test)
auc_test <- auc(roc_test)
cat("Test set AUC =", round(auc_test, 4), "\n")

# ========== Validate model on external validation set 1 and calculate AUC ==========
pred_prob_external <- predict(xgb_model_full, dexternal1)
roc_external <- roc(external_data1$group, pred_prob_external)
auc_external <- auc(roc_external)
cat("External validation set 1 AUC =", round(auc_external, 4), "\n")

# ========== Validate model on external validation set 2 and calculate AUC ==========
pred_prob_external2 <- predict(xgb_model_full, dexternal2)
roc_external2 <- roc(external_data2$group, pred_prob_external2)
auc_external2 <- auc(roc_external2)
cat("External validation set 2 AUC =", round(auc_external2, 4), "\n")

# ========== Plot and save ROC curves (3 datasets included) ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/XGBoost_ROC_all.pdf", width = 6, height = 6)
plot(roc_test, col = "darkgreen", lwd = 2, main = "XGBoost ROC: Test vs External1 vs External2")
plot(roc_external, col = "blue", lwd = 2, add = TRUE)
plot(roc_external2, col = "red", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(
          paste0("Test (AUC=", round(auc_test, 4), ")"),
          paste0("External1 (AUC=", round(auc_external, 4), ")"),
          paste0("External2 (AUC=", round(auc_external2, 4), ")")
       ), 
       col = c("darkgreen","blue","red"), lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
cat("✅ ROC curves saved\n")

# ========== Extract and save feature importance ==========
feature_importance <- xgb.importance(model = xgb_model_full)
importance_df <- data.frame(
   Feature = feature_importance$Feature,
   Importance = feature_importance$Gain
) %>% arrange(desc(Importance))
write.csv(importance_df, "/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/XGBoost_feature_importance_full.csv", row.names = FALSE)
cat("✅ Feature importance saved\n")

# ========== Plot histogram of prediction probabilities (3 datasets included) ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/Probability_distribution_all.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
hist(pred_prob_test, breaks = 20, main = "Prediction Probability Distribution (Test Set)", xlab = "Predicted Probability", col = "lightblue")
hist(pred_prob_external, breaks = 20, main = "Prediction Probability Distribution (External Set 1)", xlab = "Predicted Probability", col = "lightgreen")
hist(pred_prob_external2, breaks = 20, main = "Prediction Probability Distribution (External Set 2)", xlab = "Predicted Probability", col = "lightpink")
dev.off()

# ========== Search optimal threshold (based on test set) ==========
find_optimal_threshold <- function(probabilities, true_labels) {
   roc_obj <- roc(true_labels, probabilities)
   optimal <- coords(roc_obj, "best", best.method = "youden")
   return(optimal$threshold)
}
optimal_threshold_test <- find_optimal_threshold(pred_prob_test, test_data$group)
#optimal_threshold_test <- 0.5
cat("Optimal threshold from test set:", round(optimal_threshold_test, 4), "\n")

# ========== Calculate confusion matrices for 3 datasets ==========
# Test set
pred_class_test <- ifelse(pred_prob_test >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_test <- confusionMatrix(pred_class_test, test_data$group, positive = "1")
print("Test Set:"); print(cm_test)

# External validation set 1
pred_class_external1 <- ifelse(pred_prob_external >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external1 <- confusionMatrix(pred_class_external1, external_data1$group, positive = "1")
print("External Validation Set 1:"); print(cm_external1)

# External validation set 2
pred_class_external2 <- ifelse(pred_prob_external2 >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external2 <- confusionMatrix(pred_class_external2, external_data2$group, positive = "1")
print("External Validation Set 2:"); print(cm_external2)

# ========== Function to compute MCC ==========
calculate_mcc <- function(confusion_matrix) {
   cm <- confusion_matrix$table
   TP <- cm["1", "1"]; TN <- cm["0", "0"]; FP <- cm["1", "0"]; FN <- cm["0", "1"]
   numerator <- (TP * TN - FP * FN)
   denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
   ifelse(denominator == 0, 0, numerator / denominator)
}

# ========== Calculate MCC and extract performance metrics for 3 datasets ==========
get_metrics <- function(cm, dataset_name, mcc) {
   data.frame(
      Dataset = dataset_name,
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      F1_Score = cm$byClass["F1"],
      MCC = mcc
   )
}

metrics_list <- list(
   Test = get_metrics(cm_test, "Test", calculate_mcc(cm_test)),
   External1 = get_metrics(cm_external1, "External1", calculate_mcc(cm_external1)),
   External2 = get_metrics(cm_external2, "External2", calculate_mcc(cm_external2))
)
metrics_combined <- do.call(rbind, metrics_list)

# Save performance metrics
write.table(metrics_combined,
            file = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/Performance_all.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========== Plot bar chart of model performance metrics (3 datasets included) ==========
metric_df_long <- pivot_longer(metrics_combined, 
                               cols = c(Accuracy, Sensitivity, Specificity, F1_Score, MCC),
                               names_to = "Metric", values_to = "Value")

pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/XGBoost_Model_metrics_barplot_all.pdf", width = 10, height = 5)
ggplot(metric_df_long, aes(x = Metric, y = Value, fill = Dataset)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = sprintf("%.3f", Value)), 
             position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
   scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
   ylim(0, 1) +
   labs(y = "Score", title = "XGBoost Model Performance Metrics", fill = "Dataset") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

# ========== Plot Precision-Recall (PR) curves (3 datasets included) ==========
test_labels <- as.numeric(as.character(test_data$group))
external1_labels <- as.numeric(as.character(external_data1$group))
external2_labels <- as.numeric(as.character(external_data2$group))

# Calculate PR curves
pr_test <- pr.curve(scores.class0 = pred_prob_test, weights.class0 = test_labels, curve = TRUE)
pr_external1 <- pr.curve(scores.class0 = pred_prob_external, weights.class0 = external1_labels, curve = TRUE)
pr_external2 <- pr.curve(scores.class0 = pred_prob_external2, weights.class0 = external2_labels, curve = TRUE)

# Save PR curves
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/XGBoost_PR_curve_all.pdf", width = 6, height = 6)
plot(pr_test, color = "darkgreen", lwd = 2, main = "XGBoost PR Curve: Test vs External1 vs External2")
plot(pr_external1, color = "blue", lwd = 2, add = TRUE)
plot(pr_external2, color = "red", lwd = 2, add = TRUE)
legend("bottomleft", 
       legend = c(
          paste0("Test (AUPRC=", round(pr_test$auc.integral, 4), ")"),
          paste0("External1 (AUPRC=", round(pr_external1$auc.integral, 4), ")"),
          paste0("External2 (AUPRC=", round(pr_external2$auc.integral, 4), ")")
       ), 
       col = c("darkgreen","blue","red"), lwd = 2)
abline(h = mean(test_labels), lty = 2, col = "gray")
dev.off()

cat("✅ All analysis finished, results saved\n")

# ========== Save all model related files ==========
model_package_dir <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/XGBoost/model_package/"
dir.create(model_package_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Save XGBoost model file
xgb.save(xgb_model_full, file.path(model_package_dir, "xgboost_model.model"))

# 2. Save key preprocessing information
preprocessing_info <- list(
   feature_cols = setdiff(colnames(train_data), "group"),
   label_levels = levels(train_data$group),
   optimal_threshold = optimal_threshold_test,
   class_labels = c("Health", "Tumor")
)
saveRDS(preprocessing_info, file.path(model_package_dir, "preprocessing_info.rds"))

# 3. Save model metadata
model_metadata <- list(
   training_date = Sys.time(),
   r_version = R.version.string,
   package_versions = list(
      xgboost = packageVersion("xgboost"),
      pROC = packageVersion("pROC"),
      PRROC = packageVersion("PRROC")
   ),
   model_params = xgb_params,
   best_iteration = best_iter,
   performance = list(
      test = list(auc = auc_test, mcc = calculate_mcc(cm_test)),
      external1 = list(auc = auc_external, mcc = calculate_mcc(cm_external1)),
      external2 = list(auc = auc_external2, mcc = calculate_mcc(cm_external2))
   ),
   feature_importance = importance_df
)
saveRDS(model_metadata, file.path(model_package_dir, "model_metadata.rds"))

# 4. Save universal prediction function
predict_xgboost_model <- function(new_data, model_dir = model_package_dir, threshold = NULL) {
   # Load preprocessing dependency info
   preproc <- readRDS(file.path(model_dir, "preprocessing_info.rds"))
   model <- xgb.load(file.path(model_dir, "xgboost_model.model"))
   
   # Data preprocessing
   new_x <- new_data[, preproc$feature_cols, drop = FALSE]
   if (!is.matrix(new_x)) new_x <- as.matrix(new_x)
   
   # Create DMatrix object
   new_dmatrix <- xgb.DMatrix(data = new_x)
   
   # Predict tumor probability
   prob_tumor <- predict(model, new_dmatrix)
   
   # Determine classification threshold
   pred_threshold <- ifelse(is.null(threshold), preproc$optimal_threshold, threshold)
   
   # Predict sample class
   pred_class <- ifelse(prob_tumor >= pred_threshold, 1, 0)
   pred_class <- factor(pred_class, levels = preproc$label_levels, labels = preproc$class_labels)
   
   # Return prediction results
   return(list(
      tumor_probability = prob_tumor,
      predicted_class = pred_class,
      threshold_used = pred_threshold,
      feature_cols = preproc$feature_cols
   ))
}
saveRDS(predict_xgboost_model, file.path(model_package_dir, "predict_function.rds"))

cat("✅ XGBoost model package saved! Path: ", model_package_dir, "\n")