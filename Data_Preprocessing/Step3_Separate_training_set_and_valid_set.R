library(caret)
library(vroom)
sample_info <- read.csv("/home/yjliu/mmProj/homo_metadata_training_test_data.csv")
group <- sample_info$group

expression_data <- vroom("/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect.csv", 
                         delim = ",", col_names = TRUE)
expression_data <- as.data.frame(expression_data)
rownames(expression_data) <- make.unique(expression_data[[1]])
expression_data <- expression_data[, -1]  # 移除第一列
stopifnot(all(colnames(expression_data) == sample_info$Run))

expression_data_t <- as.data.frame(t(expression_data))
expression_data_t$group <- group

set.seed(123)
train_index <- createDataPartition(expression_data_t$group, p = 0.9, list = FALSE) 
train_90 <- expression_data_t[train_index, ]
test_10  <- expression_data_t[-train_index, ]

train_90_metadata <- sample_info[train_index, ]
test_10_metadata  <- sample_info[-train_index, ]
output_file <- "/home/yjliu/mmProj/homo_metadata_training_test_data_2133training_data.csv"
fwrite(train_90_metadata, file = output_file, row.names = TRUE, sep = ",")
output_file <- "/home/yjliu/mmProj/homo_metadata_training_test_data_236test_data.csv"
fwrite(test_10_metadata, file = output_file, row.names = TRUE, sep = ",")

train_90 <- train_90[, -ncol(train_90)]
train_90_t <- as.data.frame(t(train_90))

test_10 <- test_10[, -ncol(test_10)]
test_10_t <- as.data.frame(t(test_10))

output_file <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect_2133training_data.csv"
fwrite(train_90_t, file = output_file, row.names = TRUE, sep = ",")

output_file <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect_236test_data.csv"
fwrite(test_10_t, file = output_file, row.names = TRUE, sep = ",")
