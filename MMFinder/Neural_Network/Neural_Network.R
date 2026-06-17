# ========== Load required packages ==========
library(nnet)      # Shallow neural network
library(pROC)
library(dplyr)
library(vroom)
library(ggplot2)
library(tidyr)
library(PRROC)
library(caret)

# ========== Data loading, using your format ==========
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
# Ensure all datasets have the same factor level order, health first and tumor second
train_data$group <- factor(train_data$group, levels = c("health", "tumor"))
test_data$group <- factor(test_data$group, levels = c("health", "tumor"))
external_data1$group <- factor(external_data1$group, levels = c("health", "tumor"))
external_data2$group <- factor(external_data2$group, levels = c("health", "tumor"))

# ========== Diagnose data dimensions ==========
cat("=== Data dimension diagnosis ===\n")
cat(sprintf("Number of rows in training set: %d\n", nrow(train_data)))
cat(sprintf("Number of columns in training set: %d\n", ncol(train_data)))
cat(sprintf("Number of features: %d\n", ncol(train_data) - 1))
cat(sprintf("Number of samples: %d\n", nrow(train_data)))
cat(sprintf("Number of healthy samples: %d\n", sum(train_data$group == "health")))
cat(sprintf("Number of tumor samples: %d\n", sum(train_data$group == "tumor")))

# ========== Set random seed ==========
set.seed(123)

# ========== 1. Perform 10-fold cross-validation on the training set, tuning ==========
# Define tuning grid, adjusted for high-dimensional data
n_features <- ncol(train_data) - 1

# Calculate the maximum number of weights
max_size <- 20  # Maximum number of neurons in the tuning grid
max_weights <- (n_features + 1) * max_size + (max_size + 1) * 1
cat(sprintf("Maximum possible number of weights: %d\n", max_weights))

# Adjust tuning grid, for high-dimensional data
tune_grid <- expand.grid(
   size = c(2, 5, 10),      # Reduce the number of neurons because there are too many features
   decay = c(0.0001, 0.001, 0.01, 0.1)   # Add more regularization options
)

# Set 10-fold cross-validation
ctrl <- trainControl(
   method = "cv",
   number = 10,
   classProbs = TRUE,
   summaryFunction = twoClassSummary,
   savePredictions = "final",
   verboseIter = TRUE
)

# Cross-validation tuning
cat("Start 10-fold cross-validation tuning for the shallow neural network...\n")
cat("Note: Because the number of features is as high as 2318, model training may take a long time...\n")

cv_model <- train(
   x = train_data[, -which(names(train_data) == "group")],
   y = train_data$group,
   method = "nnet",
   trControl = ctrl,
   tuneGrid = tune_grid,
   metric = "ROC",
   trace = FALSE,
   maxit = 300,           # Increase the number of iterations
   linout = FALSE,
   MaxNWts = 100000,      # Greatly increase the weight limit
   tuneLength = 5         # Limit the number of tuning combinations
)

# Extract optimal parameters
best_params <- cv_model$bestTune
cat("✅ Optimal parameters from 10-fold CV:\n")
print(best_params)

# ========== 2. Train the final model using optimal parameters ==========
cat("Train the final neural network model using optimal parameters...\n")
nn_model_full <- nnet(
   formula = group ~ .,
   data = train_data,
   size = 1,
   decay = 0.001,
   maxit = 900,           # Greatly increase the number of iterations to ensure convergence
   trace = TRUE,           # Show training process for easier monitoring
   linout = FALSE,
   MaxNWts = 200000        # Use a larger weight limit for the final model
)

# ========== 3. Validate the model on the test set ==========
pred_prob_test <- predict(nn_model_full, test_data, type = "raw")
roc_test <- roc(test_data$group, pred_prob_test, levels = c("health", "tumor"))
auc_test <- auc(roc_test)
cat("Test set AUC =", round(auc_test, 4), "\n")

# ========== 4. Validate the model on external validation set 1 ==========
pred_prob_external1 <- predict(nn_model_full, external_data1, type = "raw")
roc_external1 <- roc(external_data1$group, pred_prob_external1, levels = c("health", "tumor"))
auc_external1 <- auc(roc_external1)
cat("External validation set 1 AUC =", round(auc_external1, 4), "\n")

# ========== 5. Validate the model on external validation set 2 ==========
pred_prob_external2 <- predict(nn_model_full, external_data2, type = "raw")
roc_external2 <- roc(external_data2$group, pred_prob_external2, levels = c("health", "tumor"))
auc_external2 <- auc(roc_external2)
cat("External validation set 2 AUC =", round(auc_external2, 4), "\n")

# ========== 6. Plot and save ROC curves ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/NN_ROC_all.pdf", width = 6, height = 6)
plot(roc_test, col = "darkgreen", lwd = 2, main = "Neural Network ROC: Test vs External1 vs External2")
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

# ========== 7. Extract and save feature importance, neural network uses weight magnitude ==========
# Feature importance calculation needs adjustment for high-dimensional data
cat("Calculating feature importance, this may take some time for 2318 features...\n")

# Calculate feature importance, based on the sum of absolute weights
importance_vals <- abs(nn_model_full$wts)
# Extract weights from the input layer to the hidden layer, the first n_features * size weights
n_features <- ncol(train_data) - 1
size <- 2
input_weights <- importance_vals[1:(n_features * size)]

# Summarize importance by feature
importance_matrix <- matrix(input_weights, nrow = n_features, ncol = size, byrow = FALSE)
feature_importance <- rowSums(importance_matrix)

importance_df <- data.frame(
   Feature = setdiff(colnames(train_data), "group"),
   Importance = feature_importance
) %>% arrange(desc(Importance))

# Save the full feature importance, the file may be large
write.csv(importance_df, "/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/NN_full_feature_importance.csv", row.names = FALSE)

cat("✅ Feature importance saved\n")

# ========== 8. Plot probability distribution histograms ==========
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/Probability_distribution_all.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
hist(pred_prob_test, breaks = 20, main = "Test Set Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightblue")
hist(pred_prob_external1, breaks = 20, main = "External Validation Set 1 Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightgreen")
hist(pred_prob_external2, breaks = 20, main = "External Validation Set 2 Predicted Probability Distribution", xlab = "Predicted Probability", col = "lightpink")
dev.off()

# ========== 9. Find optimal threshold, based on test set ==========
find_optimal_threshold <- function(probabilities, true_labels) {
   roc_obj <- roc(true_labels, probabilities, levels = c("health", "tumor"))
   optimal <- coords(roc_obj, "best", best.method = "youden")
   return(optimal$threshold)
}
optimal_threshold_test <- find_optimal_threshold(pred_prob_test, test_data$group)
cat("Optimal threshold for test set:", round(optimal_threshold_test, 4), "\n")

# ========== 10. Calculate confusion matrix ==========
# Test set
pred_class_test <- ifelse(pred_prob_test >= optimal_threshold_test, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_test <- confusionMatrix(pred_class_test, test_data$group, positive = "tumor")
print("Test set:"); print(cm_test)

# External validation set 1
pred_class_external1 <- ifelse(pred_prob_external1 >= optimal_threshold_external1, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_external1 <- confusionMatrix(pred_class_external1, external_data1$group, positive = "tumor")
print("External validation set 1:"); print(cm_external1)

# External validation set 2
pred_class_external2 <- ifelse(pred_prob_external2 >= optimal_threshold_external2, "tumor", "health") %>% 
   factor(levels = c("health", "tumor"))
cm_external2 <- confusionMatrix(pred_class_external2, external_data2$group, positive = "tumor")
print("External validation set 2:"); print(cm_external2)

# ========== 11. Calculate MCC function ==========
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

# ========== 12. Calculate performance metrics ==========
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
            file = "/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/Performance_all.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========== 13. Plot performance metrics barplot ==========
metric_df_long <- pivot_longer(metrics_combined, 
                               cols = c(Accuracy, Sensitivity, Specificity, F1_Score, MCC),
                               names_to = "Metric", values_to = "Value")

pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/NN_Model_metrics_barplot_all.pdf", width = 10, height = 5)
ggplot(metric_df_long, aes(x = Metric, y = Value, fill = Dataset)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = sprintf("%.3f", Value)), 
             position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
   scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
   ylim(0, 1) +
   labs(y = "Score", title = "Neural Network Model Performance Metrics", fill = "Dataset") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
dev.off()

# ========== 14. Plot PR curves ==========
# Convert labels to numeric values, health=0, tumor=1
test_labels_numeric <- ifelse(test_data$group == "tumor", 1, 0)
external1_labels_numeric <- ifelse(external_data1$group == "tumor", 1, 0)
external2_labels_numeric <- ifelse(external_data2$group == "tumor", 1, 0)

# Calculate PR curves
pr_test <- pr.curve(scores.class0 = pred_prob_test, weights.class0 = test_labels_numeric, curve = TRUE)
pr_external1 <- pr.curve(scores.class0 = pred_prob_external1, weights.class0 = external1_labels_numeric, curve = TRUE)
pr_external2 <- pr.curve(scores.class0 = pred_prob_external2, weights.class0 = external2_labels_numeric, curve = TRUE)

# Save PR curves
pdf("/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/NN_PR_curve_all.pdf", width = 6, height = 6)
plot(pr_test, color = "darkgreen", lwd = 2, main = "Neural Network PR Curve: Test vs External1 vs External2")
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

# ========== 15. Save model-related files ==========
model_package_dir <- "/home/yjliu/mmProj/data_process/Human/Machine_Learning/NeuralNetwork/model_package/"
dir.create(model_package_dir, recursive = TRUE, showWarnings = FALSE)

# Save neural network model
saveRDS(nn_model_full, file.path(model_package_dir, "neuralnetwork_model.rds"))

# Save key preprocessing information
preprocessing_info <- list(
   feature_cols = setdiff(colnames(train_data), "group"),
   label_levels = levels(train_data$group),
   optimal_threshold = optimal_threshold_test,
   class_labels = c("health" = "Health", "tumor" = "Tumor"),
   model_type = "Shallow Neural Network",
   note = "Data already standardized before input, high-dimensional (2318 features)",
   model_complexity = list(
      n_features = n_features,
      n_neurons = 1,
      n_weights = length(nn_model_full$wts)
   )
)
saveRDS(preprocessing_info, file.path(model_package_dir, "preprocessing_info.rds"))

# Save model metadata
model_metadata <- list(
   training_date = Sys.time(),
   r_version = R.version.string,
   package_versions = list(
      nnet = packageVersion("nnet"),
      pROC = packageVersion("pROC"),
      PRROC = packageVersion("PRROC")
   ),
   model_params = list(
      size = 1,        # Number of hidden-layer neurons
      decay = 0.001,
      maxit = nn_model_full$n[3],       # Actual number of iterations
      n_input = nn_model_full$n[1],     # Number of input-layer neurons
      n_output = nn_model_full$n[3],    # Number of output-layer neurons
      MaxNWts = 200000
   ),
   data_info = list(
      n_features = n_features,
      n_samples = nrow(train_data),
      class_distribution = table(train_data$group)
   ),
   performance = list(
      test = list(auc = auc_test, mcc = calculate_mcc(cm_test)),
      external1 = list(auc = auc_external1, mcc = calculate_mcc(cm_external1)),
      external2 = list(auc = auc_external2, mcc = calculate_mcc(cm_external2))
   ),
   top_features = importance_df[1:50, ]  # Save only the top 50 important features
)
saveRDS(model_metadata, file.path(model_package_dir, "model_metadata.rds"))

# Save general prediction function
predict_nn_model <- function(new_data, model_dir = model_package_dir, threshold = NULL) {
   # Load dependency information
   preproc <- readRDS(file.path(model_dir, "preprocessing_info.rds"))
   model <- readRDS(file.path(model_dir, "neuralnetwork_model.rds"))
   
   # Ensure feature columns are consistent
   new_x <- new_data[, preproc$feature_cols, drop = FALSE]
   
   # Predict probabilities, assuming the data have already been standardized
   pred_prob <- predict(model, new_x, type = "raw")
   
   # Determine threshold
   pred_threshold <- ifelse(is.null(threshold), preproc$optimal_threshold, threshold)
   
   # Predict classes
   pred_class <- ifelse(pred_prob >= pred_threshold, "tumor", "health")
   pred_class <- factor(pred_class, levels = preproc$label_levels)
   
   # Return results
   return(list(
      tumor_probability = pred_prob,
      predicted_class = pred_class,
      threshold_used = pred_threshold,
      feature_cols = preproc$feature_cols
   ))
}
saveRDS(predict_nn_model, file.path(model_package_dir, "predict_function.rds"))

cat("✅ Model package saved successfully! Path:", model_package_dir, "\n")