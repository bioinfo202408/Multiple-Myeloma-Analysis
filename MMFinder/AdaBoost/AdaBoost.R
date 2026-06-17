# ========== Load required packages ==========
library(adabag)
library(pROC)
library(dplyr)
library(vroom)
library(ggplot2)
library(tidyr)
library(PRROC)
library(caret)

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



# ========== Prepare data for AdaBoost ==========
# AdaBoost directly uses data frames and does not require special processing
# Ensure that group is a factor
train_data$group <- factor(train_data$group, levels = c(0, 1))
test_data$group <- factor(test_data$group, levels = c(0, 1))
external_data1$group <- factor(external_data1$group, levels = c(0, 1))
external_data2$group <- factor(external_data2$group, levels = c(0, 1))

# ========== Set AdaBoost parameters ==========
set.seed(123)

ada_params <- list(
   mfinal = 50,             # Number of iterations / boosting rounds
   maxdepth = 6,            # Maximum tree depth
   coeflearn = "Breiman"    # Weight update method
)

# ========== Train AdaBoost model ==========
ada_model_full <- boosting(
   group ~ .,               # Formula: group depends on all other variables
   data = train_data,       # Training data
   mfinal = ada_params$mfinal,
   maxdepth = ada_params$maxdepth,
   coeflearn = ada_params$coeflearn
)

cat("AdaBoost model training completed.\n")

# ========== 2. Validate the model on the test set and calculate AUC ==========
pred_prob_test <- predict.boosting(ada_model_full, newdata = test_data)$prob[,2]
roc_test <- roc(test_data$group, pred_prob_test)
auc_test <- auc(roc_test)
cat("Test set AUC =", round(auc_test, 4), "\n")

# ========== 3. Validate the model on external validation set 1 and calculate AUC ==========
pred_prob_external1 <- predict.boosting(ada_model_full, newdata = external_data1)$prob[,2]
roc_external1 <- roc(external_data1$group, pred_prob_external1)
auc_external1 <- auc(roc_external1)
cat("External validation set 1 AUC =", round(auc_external1, 4), "\n")

# ========== 3. Validate the model on external validation set 2 and calculate AUC ==========
pred_prob_external2 <- predict.boosting(ada_model_full, newdata = external_data2)$prob[,2]
roc_external2 <- roc(external_data2$group, pred_prob_external2)
auc_external2 <- auc(roc_external2)
cat("External validation set 2 AUC =", round(auc_external2, 4), "\n")

# ========== 4. Plot and save ROC curves for all three datasets ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/AdaBoost_ROC_all.pdf", width = 7, height = 6)
plot(roc_test, col = "darkgreen", lwd = 2, main = "AdaBoost ROC: Test vs External1 vs External2")
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
cat("ROC curve has been saved.\n")

# ========== 5. Extract and save feature importance ==========
feature_importance <- ada_model_full$importance
importance_df <- data.frame(
   Feature = names(feature_importance),
   Importance = as.numeric(feature_importance)
) %>% arrange(desc(Importance))

write.csv(importance_df, "/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/AdaBoost_feature_importance_full.csv", row.names = FALSE)
cat("Feature importance has been saved.\n")

# ========== Plot probability distribution histograms for all three datasets ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/Probability_distribution_all.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))  # Layout with 1 row and 3 columns
hist(pred_prob_test, breaks = 20, main = "Prediction Probability Distribution: Test Set", xlab = "Prediction Probability", col = "lightblue")
hist(pred_prob_external1, breaks = 20, main = "Prediction Probability Distribution: External Set 1", xlab = "Prediction Probability", col = "lightgreen")
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
pred_class_external1 <- ifelse(pred_prob_external1 >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external1 <- confusionMatrix(pred_class_external1, external_data1$group, positive = "1")
print("External validation set 1:"); print(cm_external1)

# External validation set 2
pred_class_external2 <- ifelse(pred_prob_external2 >= optimal_threshold_test, 1, 0) %>% factor(levels = c(0, 1))
cm_external2 <- confusionMatrix(pred_class_external2, external_data2$group, positive = "1")
print("External validation set 2:"); print(cm_external2)

# ========== Function to calculate MCC ==========
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
            file = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/Performance_all.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========== Plot performance metric barplot for all three datasets ==========
metric_df_long <- pivot_longer(metrics_combined, 
                               cols = c(Accuracy, Sensitivity, Specificity, F1_Score, MCC),
                               names_to = "Metric", values_to = "Value")

pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/AdaBoost_Model_metrics_barplot_all.pdf", width = 10, height = 5)
ggplot(metric_df_long, aes(x = Metric, y = Value, fill = Dataset)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = sprintf("%.3f", Value)), 
             position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
   scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +  # Color scheme for the three datasets
   ylim(0, 1) +
   labs(y = "Score", title = "AdaBoost Model Performance Metrics", fill = "Dataset") +
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
pr_external1 <- pr.curve(scores.class0 = pred_prob_external1, weights.class0 = external1_labels, curve = TRUE)
pr_external2 <- pr.curve(scores.class0 = pred_prob_external2, weights.class0 = external2_labels, curve = TRUE)

# Save PR curves
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/AdaBoost_PR_curve_all.pdf", width = 7, height = 6)
plot(pr_test, color = "darkgreen", lwd = 2, main = "AdaBoost PR Curve: Test vs External1 vs External2")
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
model_package_dir <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/AdaBoost/model_package/"
dir.create(model_package_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Save the AdaBoost model file
saveRDS(ada_model_full, file.path(model_package_dir, "adaboost_model.rds"))

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
      adabag = packageVersion("adabag"),
      pROC = packageVersion("pROC"),
      PRROC = packageVersion("PRROC")
   ),
   model_params = ada_params,                              # Final training parameters
   performance = list(
      test = list(auc = auc_test, mcc = calculate_mcc(cm_test)),
      external1 = list(auc = auc_external1, mcc = calculate_mcc(cm_external1)),
      external2 = list(auc = auc_external2, mcc = calculate_mcc(cm_external2))
   ),
   feature_importance = importance_df                      # Feature importance
)

saveRDS(model_metadata, file.path(model_package_dir, "model_metadata.rds"))

# 4. Save a general prediction function for new data
predict_adaboost_model <- function(new_data, model_dir = model_package_dir, threshold = NULL) {
   # Load dependency information
   preproc <- readRDS(file.path(model_dir, "preprocessing_info.rds"))
   model <- readRDS(file.path(model_dir, "adaboost_model.rds"))
   
   # Data preprocessing: ensure consistent column names
   new_data_processed <- new_data[, preproc$feature_cols, drop = FALSE]
   
   # Predict probabilities
   pred_result <- predict.boosting(model, newdata = new_data_processed)
   prob_tumor <- pred_result$prob[,2]  # The second column is the tumor probability
   
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

saveRDS(predict_adaboost_model, file.path(model_package_dir, "predict_function.rds"))

cat("AdaBoost model package has been saved. Path:", model_package_dir, "\n")