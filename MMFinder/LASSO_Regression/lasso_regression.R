# ========== Load required packages ==========
library(glmnet)
library(pROC)
library(dplyr)
library(vroom)
library(ggplot2)
library(tidyr)
library(PRROC)
library(caret)

# ========== Data loading ==========
train_data <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data.csv", 
   delim = ",", 
)
train_data <- as.data.frame(train_data)
rownames(train_data) <- make.unique(train_data[[1]])
train_data <- train_data[,-1]
train_data$group <- as.factor(train_data$group)

test_data <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/test_data.csv", 
   delim = ",", 
)
test_data <- as.data.frame(test_data)
rownames(test_data) <- make.unique(test_data[[1]])
test_data <- test_data[,-1]
test_data$group <- as.factor(test_data$group)

external_data1 <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/external1.csv", 
   delim = ",", 
)
external_data1 <- as.data.frame(external_data1)
rownames(external_data1) <- make.unique(external_data1[[1]])
external_data1 <- external_data1[,-1]
external_data1$group <- as.factor(external_data1$group)

external_data2 <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/external2.csv", 
   delim = ",", 
)
external_data2 <- as.data.frame(external_data2)
rownames(external_data2) <- make.unique(external_data2[[1]])
external_data2 <- external_data2[,-1]
external_data2$group <- as.factor(external_data2$group)

# ========== Ensure factor level consistency ==========
# Logistic regression requires tumor as the positive class
train_data$group <- factor(train_data$group, levels = c("health", "tumor"))
test_data$group <- factor(test_data$group, levels = c("health", "tumor"))
external_data1$group <- factor(external_data1$group, levels = c("health", "tumor"))
external_data2$group <- factor(external_data2$group, levels = c("health", "tumor"))

# ========== Set random seed ==========
set.seed(123)

# ========== 1. Basic data check ==========
cat("Basic data check...\n")
cat("Training set dimensions:", dim(train_data), "\n")
cat("Test set dimensions:", dim(test_data), "\n")
cat("External validation set 1 dimensions:", dim(external_data1), "\n")
cat("External validation set 2 dimensions:", dim(external_data2), "\n")

# Check class distribution
cat("\nClass distribution:\n")
cat("Training set - health:", sum(train_data$group == "health"), 
    ", tumor:", sum(train_data$group == "tumor"), "\n")
cat("Test set - health:", sum(test_data$group == "health"), 
    ", tumor:", sum(test_data$group == "tumor"), "\n")

# ========== 2. Use caret for 10-fold cross-validation to tune regularization parameters ==========
cat("\nStart 10-fold cross-validation tuning...\n")

# Prepare training data
X_train <- train_data[, -which(names(train_data) == "group")]
y_train <- train_data$group

# Set cross-validation control parameters
ctrl <- trainControl(
   method = "cv",
   number = 10,
   classProbs = TRUE,
   summaryFunction = twoClassSummary,
   savePredictions = "final",
   verboseIter = TRUE,
   selectionFunction = "best"
)

# Define tuning grid, regularization parameters
tune_grid <- expand.grid(
   alpha = 0,  # 0:Ridge, 0.5:Elastic Net, 1:Lasso
   lambda = c(0.001,0.01,0.1,1)
)
# tune_grid <- expand.grid(
#    alpha = c(0, 0.5, 1),  # 0:Ridge, 0.5:Elastic Net, 1:Lasso
#    lambda = 10^seq(-3, 1, length = 10)
# )
# Use glmnet for cross-validation tuning
cv_model <- train(
   x = X_train,
   y = y_train,
   method = "glmnet",
   family = "binomial",
   trControl = ctrl,
   tuneGrid = tune_grid,
   metric = "ROC",
   preProc = NULL  # No preprocessing, because the data have already been standardized
)

# Extract optimal parameters
best_params <- cv_model$bestTune
cat("✅ Optimal parameters from 10-fold CV:\n")
print(best_params)

# ========== 3. Train final model with optimal parameters ==========
cat("Start training the final logistic regression model...\n")

# Train the full model using optimal parameters
logit_model_final <- glmnet(
   x = as.matrix(X_train),
   y = y_train,
   family = "binomial",
   alpha = 1,
   lambda = 0.001,
   standardize = FALSE  # No standardization, because the data have already been standardized
)

cat("✅ Logistic regression model training completed\n")


# pred_prob_train <- predict(logit_model_final, 
#                            newx = as.matrix(X_train), 
#                            type = "response")

# ========== 4. Validate the model on the test set ==========
X_test <- test_data[, -which(names(test_data) == "group")]
y_test <- test_data$group

pred_prob_test <- predict(logit_model_final, 
                          newx = as.matrix(X_test), 
                          type = "response")
roc_test <- roc(y_test, pred_prob_test, levels = c("health", "tumor"))
auc_test <- auc(roc_test)
cat("Test set AUC =", round(auc_test, 4), "\n")

# ========== 5. Validate the model on external validation set 1 ==========
X_ext1 <- external_data1[, -which(names(external_data1) == "group")]
y_ext1 <- external_data1$group

pred_prob_external1 <- predict(logit_model_final, 
                               newx = as.matrix(X_ext1), 
                               type = "response")
roc_external1 <- roc(y_ext1, pred_prob_external1, levels = c("health", "tumor"))
auc_external1 <- auc(roc_external1)
cat("External validation set 1 AUC =", round(auc_external1, 4), "\n")

# ========== 6. Validate the model on external validation set 2 ==========
X_ext2 <- external_data2[, -which(names(external_data2) == "group")]
y_ext2 <- external_data2$group

pred_prob_external2 <- predict(logit_model_final, 
                               newx = as.matrix(X_ext2), 
                               type = "response")
roc_external2 <- roc(y_ext2, pred_prob_external2, levels = c("health", "tumor"))
auc_external2 <- auc(roc_external2)
cat("External validation set 2 AUC =", round(auc_external2, 4), "\n")

# ========== 7. Plot and save ROC curves ==========
output_dir <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/lasso_regression/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf(paste0(output_dir, "lasso_regression_ROC_all.pdf"), width = 6, height = 6)
plot(roc_test, col = "darkgreen", lwd = 2, main = "lasso_regression Regression ROC: Test vs External1 vs External2")
plot(roc_external1, col = "blue", lwd = 2, add = TRUE)
plot(roc_external2, col = "red", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(
          paste0("Test (AUC=", round(auc_test, 4), ")"),
          paste0("External1 (AUC=", round(auc_external1, 4), ")"),
          paste0("External2 (AUC=", round(auc_external2, 4), ")")
       ), 
       col = c("darkgreen","blue","red"), lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
cat("✅ ROC curve saved\n")

# ========== 8. Extract and save feature importance ==========
# Extract model coefficients
coef_matrix <- as.matrix(coef(logit_model_final))
coef_df <- data.frame(
   Feature = rownames(coef_matrix),
   Coefficient = coef_matrix[, 1],
   Abs_Coefficient = abs(coef_matrix[, 1])
) %>% arrange(desc(Abs_Coefficient))

# Remove intercept term
coef_df <- coef_df[coef_df$Feature != "(Intercept)", ]

# Calculate OR values, odds ratios
coef_df$Odds_Ratio <- exp(coef_df$Coefficient)

# Save feature importance
importance_df <- coef_df
write.csv(importance_df, paste0(output_dir, "lasso_regression_feature_importance_full.csv"), row.names = FALSE)
cat("✅ Feature importance saved\n")

# ========== 9. Plot coefficient figure ==========
pdf(paste0(output_dir, "lasso_regression_coefficient_plot.pdf"), width = 10, height = 6)
# Select the top 20 most important features
top_n <- min(20, nrow(importance_df))
top_features <- importance_df[1:top_n, ]

ggplot(top_features, aes(x = reorder(Feature, Odds_Ratio), y = Odds_Ratio)) +
   geom_bar(stat = "identity", fill = ifelse(top_features$Odds_Ratio > 1, "#d95f02", "#1b9e77")) +
   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
   coord_flip() +
   labs(x = "Feature", y = "Odds Ratio", 
        title = "Top Feature Odds Ratios (lasso_regression Regression)",
        subtitle = paste0("Odds Ratio > 1 indicates positive association with tumor")) +
   theme_bw() +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()

# ========== 10. Plot probability distribution histograms ==========
pdf(paste0(output_dir, "lasso_regression_Probability_distribution_all.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))
hist(pred_prob_test, breaks = 20, main = "Test Set Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightblue")
hist(pred_prob_external1, breaks = 20, main = "External Validation Set 1 Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightgreen")
hist(pred_prob_external2, breaks = 20, main = "External Validation Set 2 Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightpink")
dev.off()

# ========== 11. Find optimal threshold ==========
find_optimal_threshold <- function(probabilities, true_labels) {
   roc_obj <- roc(true_labels, probabilities, levels = c("health", "tumor"))
   optimal <- coords(roc_obj, "best", best.method = "youden")
   return(optimal$threshold)
}

optimal_threshold_test <- find_optimal_threshold(pred_prob_test, y_test)
cat("Optimal threshold for test set:", round(optimal_threshold_test, 4), "\n")

# ========== 12. Calculate confusion matrix ==========
# Test set
pred_class_test <- ifelse(pred_prob_test >= optimal_threshold_test, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_test <- confusionMatrix(pred_class_test, y_test, positive = "tumor")
print("Confusion matrix for test set:"); print(cm_test$table)

# External validation set 1, using the optimal threshold obtained from the test set
pred_class_external1 <- ifelse(pred_prob_external1 >= optimal_threshold_test, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_external1 <- confusionMatrix(pred_class_external1, y_ext1, positive = "tumor")
print("Confusion matrix for external validation set 1:"); print(cm_external1)

# External validation set 2, using the optimal threshold obtained from the test set
pred_class_external2 <- ifelse(pred_prob_external2 >= optimal_threshold_test, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_external2 <- confusionMatrix(pred_class_external2, y_ext2, positive = "tumor")
print("Confusion matrix for external validation set 2:"); print(cm_external2$table)

# ========== 13. Calculate MCC function ==========
calculate_mcc <- function(confusion_matrix) {
   cm <- confusion_matrix$table
   TP <- cm["tumor", "tumor"]
   TN <- cm["health", "health"]
   FP <- cm["tumor", "health"]
   FN <- cm["health", "tumor"]
   numerator <- (TP * TN - FP * FN)
   denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
   ifelse(denominator == 0, 0, numerator / denominator)
}

# ========== 14. Calculate performance metrics ==========
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
            file = paste0(output_dir, "lasso_regression_Performance_all.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========== 15. Plot performance metrics barplot ==========
metric_df_long <- pivot_longer(metrics_combined, 
                               cols = c(Accuracy, Sensitivity, Specificity, F1_Score, MCC),
                               names_to = "Metric", values_to = "Value")

pdf(paste0(output_dir, "lasso_regression_Model_metrics_barplot_all.pdf"), width = 10, height = 5)
ggplot(metric_df_long, aes(x = Metric, y = Value, fill = Dataset)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = sprintf("%.3f", Value)), 
             position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
   scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
   ylim(0, 1) +
   labs(y = "Score", title = "lasso_regression Regression Model Performance Metrics", fill = "Dataset") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

# ========== 16. Plot PR curves ==========
# Convert labels to numeric values, health=0, tumor=1
test_labels_numeric <- ifelse(y_test == "tumor", 1, 0)
external1_labels_numeric <- ifelse(y_ext1 == "tumor", 1, 0)
external2_labels_numeric <- ifelse(y_ext2 == "tumor", 1, 0)

# Calculate PR curves
pr_test <- pr.curve(scores.class0 = pred_prob_test, weights.class0 = test_labels_numeric, curve = TRUE)
pr_external1 <- pr.curve(scores.class0 = pred_prob_external1, weights.class0 = external1_labels_numeric, curve = TRUE)
pr_external2 <- pr.curve(scores.class0 = pred_prob_external2, weights.class0 = external2_labels_numeric, curve = TRUE)

# Save PR curves
pdf(paste0(output_dir, "lasso_regression_PR_curve_all.pdf"), width = 6, height = 6)
plot(pr_test, color = "darkgreen", lwd = 2, main = "lasso_regression Regression PR Curve: Test vs External1 vs External2")
plot(pr_external1, color = "blue", lwd = 2, add = TRUE)
plot(pr_external2, color = "red", lwd = 2, add = TRUE)
legend("bottomleft", 
       legend = c(
          paste0("Test (AUPRC=", round(pr_test$auc.integral, 4), ")"),
          paste0("External1 (AUPRC=", round(pr_external1$auc.integral, 4), ")"),
          paste0("External2 (AUPRC=", round(pr_external2$auc.integral, 4), ")")
       ), 
       col = c("darkgreen","blue","red"), lwd = 2)
abline(h = mean(test_labels_numeric), lty = 2, col = "gray")
dev.off()

cat("✅ All analyses completed, and results have been saved\n")

# ========== 17. Save model-related files ==========
model_package_dir <- paste0(output_dir, "model_package/")
dir.create(model_package_dir, recursive = TRUE, showWarnings = FALSE)

# Save logistic regression model
saveRDS(logit_model_final, file.path(model_package_dir, "lasso_regression_model.rds"))

# Save key preprocessing information
preprocessing_info <- list(
   feature_cols = colnames(X_train),
   label_levels = levels(y_train),
   optimal_threshold = optimal_threshold_test,
   class_labels = c("health" = "Health", "tumor" = "Tumor")
)
saveRDS(preprocessing_info, file.path(model_package_dir, "preprocessing_info.rds"))

# Save model metadata
model_metadata <- list(
   training_date = Sys.time(),
   r_version = R.version.string,
   package_versions = list(
      glmnet = packageVersion("glmnet"),
      pROC = packageVersion("pROC"),
      PRROC = packageVersion("PRROC")
   ),
   model_params = list(
      alpha = 0,
      lambda = 0.01
   ),
   performance = list(
      test = list(auc = auc_test, mcc = calculate_mcc(cm_test)),
      external1 = list(auc = auc_external1, mcc = calculate_mcc(cm_external1)),
      external2 = list(auc = auc_external2, mcc = calculate_mcc(cm_external2))
   ),
   feature_importance = importance_df
)
saveRDS(model_metadata, file.path(model_package_dir, "model_metadata.rds"))

# ========== 18. Save general prediction function ==========
predict_lasso_regression_model <- function(new_data, model_dir = model_package_dir, threshold = NULL) {
   # Load dependency information
   preproc <- readRDS(file.path(model_dir, "preprocessing_info.rds"))
   model <- readRDS(file.path(model_dir, "lasso_regression_model.rds"))
   
   # Ensure feature columns are consistent
   missing_features <- setdiff(preproc$feature_cols, colnames(new_data))
   if (length(missing_features) > 0) {
      stop(paste("Missing features:", paste(missing_features, collapse = ", ")))
   }
   
   # Select features
   new_x <- new_data[, preproc$feature_cols, drop = FALSE]
   
   # Predict probabilities
   pred <- predict(model, 
                   newx = as.matrix(new_x), 
                   type = "response")
   
   # Determine threshold
   pred_threshold <- ifelse(is.null(threshold), preproc$optimal_threshold, threshold)
   
   # Predict classes
   pred_class <- ifelse(pred >= pred_threshold, "tumor", "health")
   pred_class <- factor(pred_class, levels = preproc$label_levels)
   
   # Return results
   return(list(
      tumor_probability = as.vector(pred),
      predicted_class = pred_class,
      threshold_used = pred_threshold,
      feature_cols = preproc$feature_cols
   ))
}
saveRDS(predict_lasso_regression_model, file.path(model_package_dir, "predict_function.rds"))

cat("✅ Model package saved successfully! Path:", model_package_dir, "\n")
cat("✅ Logistic regression modeling workflow fully completed\n")