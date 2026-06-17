train_data <- vroom("/home/yjliu/mmProj/homo/procdata/mRNA/Homo_mRNA_TPM_removeBatchEffect_2133training_data_newCorName_FilterJS_AvgExpres_Proportion.csv", 
                 delim = ",", col_names = TRUE)
train_data <- as.data.frame(train_data)
rownames(train_data) <- make.unique(train_data[[1]])
train_data <- train_data[, -1]
train_df_t <- as.data.frame(t(train_data))

test_data <- vroom("/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect__236test_data_FilterJS_AvgExpres_Proportion.csv", 
                 delim = ",", col_names = TRUE)
test_data <- as.data.frame(test_data)
rownames(test_data) <- make.unique(test_data[[1]])
test_data <- test_data[, -1]
test_data_t <- as.data.frame(t(test_data))

external1 <- vroom(
   "/home/yjliu/mmProj/homo/procdata/homo_mRNA_TPM_external_data1_FilterJS_AvgExpres_Proportion.csv", 
   delim = ",", 
)
external1 <- as.data.frame(external1)
rownames(external1) <- make.unique(external1[[1]])
external1 <- external1[,-1]
external1_t <- as.data.frame(t(external1))

external2 <- vroom(
   "/home/yjliu/mmProj/homo/procdata/homo_mRNA_TPM_external_data1_FilterJS_AvgExpres_Proportion.csv", 
   delim = ",", 
)
external2 <- as.data.frame(external2)
rownames(external2) <- make.unique(external2[[1]])
external2 <- external2[,-1]
external2_t <- as.data.frame(t(external2))

zscore_fit_transform <- function(train_df_t) {
   feature_cols <- colnames(train_df_t)
   train_mean <- sapply(train_df_t[, feature_cols], mean)
   train_sd   <- sapply(train_df_t[, feature_cols], sd)
   # 避免 sd = 0
   train_sd[train_sd == 0] <- 1
   
   train_scaled <- train_df_t
   train_scaled[, feature_cols] <- scale(train_df_t[, feature_cols], center = train_mean, scale = train_sd)
   
   list(train_scaled = train_scaled, mean = train_mean, sd = train_sd)
}

zscore_apply <- function(df, mean_vals, sd_vals) {
   feature_cols <- colnames(df)
   sd_vals[sd_vals == 0] <- 1
   
   df_scaled <- df
   df_scaled[, feature_cols] <- scale(df[, feature_cols], center = mean_vals, scale = sd_vals)
   return(df_scaled)
}

zfit <- zscore_fit_transform(train_df_t)
save(zfit,file = "/home/yjliu/mmProj/Machine_learning/mRNA_zfit.Rdata")
train_data_scaled <- zfit$train_scaled
test_scaled <- zscore_apply(test_data_t, zfit$mean, zfit$sd)
external1_scaled <- zscore_apply(external1_t, zfit$mean, zfit$sd)
external2_scaled <- zscore_apply(external2_t, zfit$mean, zfit$sd)

train_data_scaled_t <- as.data.frame(t(train_data_scaled))
train_data_scaled_t <- cbind(gene_name = rownames(train_data_scaled_t), train_data_scaled_t)
summary_file <- "/home/yjliu/mmProj/homo/procdata/mRNA/Homo_mRNA_TPM_removeBatchEffect_2133training_data_newCorName_FilterJS_AvgExpres_Proportion_scaled.csv"
fwrite(train_data_scaled_t, file = summary_file, row.names = FALSE, sep = ",")

test_scaled_t <- as.data.frame(t(test_scaled))
summary_file1 <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect__236test_data_FilterJS_AvgExpres_Proportion_sclaed.csv"
fwrite(test_scaled_t, file = summary_file1, row.names = TRUE, sep = ",")

external1_scaled_t <- as.data.frame(t(external1_scaled))
summary_file1 <- "/home/yjliu/mmProj/homo/procdata/homo_mRNA_TPM_external_data1_FilterJS_AvgExpres_Proportion_sclaed.csv"
fwrite(external1_scaled_t, file = summary_file1, row.names = TRUE, sep = ",")

external2_scaled_t <- as.data.frame(t(external2_scaled))
summary_file2 <- "/home/yjliu/mmProj/homo/procdata/homo_mRNA_TPM_external_data1_FilterJS_AvgExpres_Proportion_sclaed.csv"
fwrite(external2_scaled_t, file = summary_file2, row.names = TRUE, sep = ",")
