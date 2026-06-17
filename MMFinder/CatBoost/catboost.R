# ========== Load required packages ==========
library(catboost)
library(pROC)
library(dplyr)
library(vroom)
library(ggplot2)
library(tidyr)
library(PRROC)
library(caret)

# ========== Load data ==========
train_data <- vroom(
   "/home/yjliu/mmProj/data_process/Human/Machine_Learning/training_data_1_0.csv", 
   delim = ",", 
)
train_data <- as.data.frame(train_data)
rownames(train_data) <- make.unique(train_data[[1]])
train_data <- train_data[,-1]

# Ensure that group is a factor, where 1 = tumor and 0 = health
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

# ========== Prepare data for CatBoost ==========
prepare_catboost_data <- function(df, label_col = "group") {
   feature_cols <- setdiff(colnames(df), label_col)
   x <- as.matrix(df[, feature_cols])
   y <- as.numeric(df[[label_col]]) - 1  # Convert to 0/1, where 0 = health and 1 = tumor
   return(list(x = x, y = y))
}

# Prepare training data
train_prepared <- prepare_catboost_data(train_data)
train_pool <- catboost.load_pool(data = train_prepared$x, label = train_prepared$y)

# Prepare test and external validation data
test_prepared <- prepare_catboost_data(test_data)
test_pool <- catboost.load_pool(data = test_prepared$x, label = test_prepared$y)

external_prepared1 <- prepare_catboost_data(external_data1)
external_pool1 <- catboost.load_pool(data = external_prepared1$x, label = external_prepared1$y)

external_prepared2 <- prepare_catboost_data(external_data2)
external_pool2 <- catboost.load_pool(data = external_prepared2$x, label = external_prepared2$y)

# ========== Set CatBoost parameters ==========
set.seed(123)

catboost_params <- list(
   iterations = 1000,
   learning_rate = 0.05,
   depth = 6,
   l2_leaf_reg = 3,
   border_count = 32,
   loss_function = 'Logloss',
   eval_metric = 'AUC',
   random_seed = 123,
   od_type = "Iter",
   od_wait = 50,
   verbose = 100
)
 
# ========== 1. Perform 10-fold cross-validation on the training set to identify the optimal number of iterations ==========
cv_result <- catboost.cv(
   pool = train_pool,               # Use only the training set for CV
   params = catboost_params,        # Base parameters
   fold_count = 10,                 # 10-fold CV
   partition_random_seed = 123      # Ensure reproducible CV partitioning
)

# Extract the optimal number of iterations based on the minimum CV validation Logloss
best_iter <- which.min(cv_result$test.Logloss.mean)
cat("Optimal number of iterations from 10-fold CV =", best_iter, "\n")

# ========== 2. Train the final model using the optimal number of iterations ==========
# Update parameters: fix the optimal number of iterations and disable early stopping
catboost_params_opt <- catboost_params
catboost_params_opt$iterations <- best_iter  # Replace with the optimal CV iteration
catboost_params_opt$od_type <- NULL          # Remove early stopping parameter
catboost_params_opt$od_wait <- NULL

# Train the optimal model using the full training set
catboost_model_full <- catboost.train(
   learn_pool = train_pool,
   params = catboost_params_opt
)

cat("CatBoost model based on 10-fold CV has been trained successfully.\n")

# ========== 2. Validate the model on the test set and calculate AUC ==========
pred_prob_test <- catboost.predict(catboost_model_full, test_pool, prediction_type = "Probability")
roc_test <- roc(test_data$group, pred_prob_test)
auc_test <- auc(roc_test)
cat("Test set AUC =", round(auc_test, 4), "\n")

# ========== 3. Validate the model on independent test set 1 and calculate AUC ==========
pred_prob_external <- catboost.predict(catboost_model_full, external_pool1, prediction_type = "Probability")
roc_external <- roc(external_data1$group, pred_prob_external)
auc_external <- auc(roc_external)
cat("External validation set 1 AUC =", round(auc_external, 4), "\n")

# ========== 3. Validate the model on independent test set 2 and calculate AUC ==========
pred_prob_external2 <- catboost.predict(catboost_model_full, external_pool2, prediction_type = "Probability")
roc_external2 <- roc(external_data2$group, pred_prob_external2)
auc_external2 <- auc(roc_external2)
cat("External validation set 2 AUC =", round(auc_external2, 4), "\n")

# ========== 4. Plot and save ROC curves for all three datasets ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/CatBoost_ROC_all.pdf", width = 7, height = 6)
plot(roc_test, col = "darkgreen", lwd = 2, main = "CatBoost ROC: Test vs External1 vs External2")
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
cat("ROC curve has been saved.\n")

# ========== 5. Extract and save feature importance without modification ==========
feature_importance <- catboost.get_feature_importance(catboost_model_full, pool = train_pool)
importance_df <- data.frame(
   Feature = setdiff(colnames(train_data), "group"),
   Importance = as.numeric(feature_importance)
) %>% arrange(desc(Importance))
write.csv(importance_df, "//home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/CatBoost_feature_importance_full.csv", row.names = FALSE)
cat("Feature importance has been saved.\n")

# ========== Plot probability distribution histograms for all three datasets ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/Probability_distribution_all.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))  # Layout with 1 row and 3 columns
hist(pred_prob_test, breaks = 20, main = "Prediction Probability Distribution: Test Set", xlab = "Prediction Probability", col = "lightblue")
hist(pred_prob_external, breaks = 20, main = "Prediction Probability Distribution: External Set 1", xlab = "Prediction Probability", col = "lightgreen")
hist(pred_prob_external2, breaks = 20, main = "Prediction Probability Distribution: External Set 2", xlab = "Prediction Probability", col = "lightpink")
dev.off()

# ========== Find the optimal threshold based on the test set ==========
find_optimal_threshold <- function(probabilities, true_labels) {
   roc_obj <- roc(true_labels, probabilities)
   optimal <- coords(roc_obj, "best", best.method = "youden")
   return(optimal$threshold)
}

optimal_threshold_test <- find_optimal_threshold(pred_prob_test, test_data$group)
cat("Optimal threshold for the test set:", round(optimal_threshold_test, 4), "\n")

# ========== Calculate confusion matrices for all three datasets ==========
# Test set
pred_class_test <- ifelse(pred_prob_test >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_test <- confusionMatrix(pred_class_test, test_data$group, positive = "1")
print("Test set:"); print(cm_test)

# External validation set 1
pred_class_external1 <- ifelse(pred_prob_external >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external1 <- confusionMatrix(pred_class_external1, external_data1$group, positive = "1")
print("External validation set 1:"); print(cm_external1)

# External validation set 2
pred_class_external2 <- ifelse(pred_prob_external2 >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external2 <- confusionMatrix(pred_class_external2, external_data2$group, positive = "1")
print("External validation set 2:"); print(cm_external2)


# ========== Function to calculate MCC without modification ==========
calculate_mcc <- function(confusion_matrix) {
   cm <- confusion_matrix$table
   TP <- cm["1", "1"]; TN <- cm["0", "0"]; FP <- cm["1", "0"]; FN <- cm["0", "1"]
   numerator <- (TP * TN - FP * FN)
   denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
   ifelse(denominator == 0, 0, numerator / denominator)
}

# ========== Calculate MCC and extract performance metrics for all three datasets ==========
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
            file = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/Performance_all.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========== Plot performance metric barplot for all three datasets ==========
metric_df_long <- pivot_longer(metrics_combined, 
                               cols = c(Accuracy, Sensitivity, Specificity, F1_Score, MCC),
                               names_to = "Metric", values_to = "Value")

pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/CatBoost_Model_metrics_barplot_all.pdf", width = 10, height = 5)
ggplot(metric_df_long, aes(x = Metric, y = Value, fill = Dataset)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = sprintf("%.3f", Value)), 
             position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
   scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +  # Color scheme for the three datasets
   ylim(0, 1) +
   labs(y = "Score", title = "CatBoost Model Performance Metrics", fill = "Dataset") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

# ========== Plot PR curves for all three datasets ==========
# Convert labels to numeric values consistently
test_labels <- as.numeric(as.character(test_data$group))
external1_labels <- as.numeric(as.character(external_data1$group))
external2_labels <- as.numeric(as.character(external_data2$group))

# Calculate PR curves
pr_test <- pr.curve(scores.class0 = pred_prob_test, weights.class0 = test_labels, curve = TRUE)
pr_external1 <- pr.curve(scores.class0 = pred_prob_external, weights.class0 = external1_labels, curve = TRUE)
pr_external2 <- pr.curve(scores.class0 = pred_prob_external2, weights.class0 = external2_labels, curve = TRUE)

# Save PR curves
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/CatBoost_PR_curve_all.pdf", width = 7, height = 6)
plot(pr_test, color = "darkgreen", lwd = 2, main = "CatBoost PR Curve: Test vs External1 vs External2")
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

cat("All analyses have been completed and results have been saved.\n")

# ========== Save model-related files ==========
model_package_dir <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/CatBoost/model_package/"
dir.create(model_package_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Save the CatBoost model file
catboost.save_model(catboost_model_full, file.path(model_package_dir, "catboost_model.cbm"))

# 2. Save key preprocessing information, including feature columns, label levels, optimal threshold, etc.
preprocessing_info <- list(
   feature_cols = setdiff(colnames(train_data), "group"),  # Feature column names
   label_levels = levels(train_data$group),                 # Label factor levels, where 0 = health and 1 = tumor
   optimal_threshold = optimal_threshold_test,              # Optimal threshold from the test set
   class_labels = c("Health", "Tumor")                      # Class label mapping
)

saveRDS(preprocessing_info, file.path(model_package_dir, "preprocessing_info.rds"))

# 3. Save model metadata, including training information and performance metrics
model_metadata <- list(
   training_date = Sys.time(),
   r_version = R.version.string,
   package_versions = list(
      catboost = packageVersion("catboost"),
      pROC = packageVersion("pROC"),
      PRROC = packageVersion("PRROC")
   ),
   model_params = catboost_params_opt,                      # Final training parameters
   best_iteration = best_iter,                              # Optimal iteration from CV
   performance = list(
      test = list(auc = auc_test, mcc = calculate_mcc(cm_test)),
      external1 = list(auc = auc_external, mcc = calculate_mcc(cm_external1)),
      external2 = list(auc = auc_external2, mcc = calculate_mcc(cm_external2))
   ),
   feature_importance = importance_df                       # Feature importance
)

saveRDS(model_metadata, file.path(model_package_dir, "model_metadata.rds"))

# 4. Save a general prediction function for new data
predict_catboost_model <- function(new_data, model_dir = model_package_dir, threshold = NULL) {
   # Load dependency information
   preproc <- readRDS(file.path(model_dir, "preprocessing_info.rds"))
   model <- catboost.load_model(file.path(model_dir, "catboost_model.cbm"))
   
   # Data preprocessing: extract feature columns and convert to matrix
   new_x <- new_data[, preproc$feature_cols, drop = FALSE]
   if (!is.matrix(new_x)) new_x <- as.matrix(new_x)
   
   # Create CatBoost data pool
   new_pool <- catboost.load_pool(data = new_x)
   
   # Predict probabilities
   prob_tumor <- catboost.predict(model, new_pool, prediction_type = "Probability")
   
   # Determine the threshold using either the optimal threshold or a custom threshold
   pred_threshold <- ifelse(is.null(threshold), preproc$optimal_threshold, threshold)
   
   # Predict classes and map them to the original factor levels
   pred_class <- ifelse(prob_tumor >= pred_threshold, 1, 0)
   pred_class <- factor(pred_class, levels = preproc$label_levels, labels = preproc$class_labels)
   
   # Return results
   return(list(
      tumor_probability = prob_tumor,
      predicted_class = pred_class,
      threshold_used = pred_threshold,
      feature_cols = preproc$feature_cols
   ))
}

saveRDS(predict_catboost_model, file.path(model_package_dir, "predict_function.rds"))

cat("Model package has been saved. Path:", model_package_dir, "\n")