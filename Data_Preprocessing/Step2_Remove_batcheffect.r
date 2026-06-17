# Step 1: Read expression data
library(vroom)
library(limma)
library(data.table)

expression_data <- vroom("/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM.csv", 
                         delim = ",", col_names = TRUE)
expression_data <- as.data.frame(expression_data)
rownames(expression_data) <- make.unique(expression_data[[1]])
expression_data <- expression_data[, -1]  # Remove the first column
expression_data <- expression_data[, -1]  # Remove GeneName column, if present
expression_data <- expression_data[rowSums(expression_data) > 0, ]

# Step 3: Convert to matrix and perform log2 transformation
expr_matrix <- as.matrix(expression_data)
expr_matrix_log <- log2(expr_matrix + 1)  # Apply log2 transformation to the data

# Step 4: Read sample information
sample_info <- read.csv("/home/yjliu/mmProj/homo_metadata_training_test_data.csv")
stopifnot(all(colnames(expr_matrix_log) == sample_info$Run))
sample_info$group <- as.factor(sample_info$group)   # Convert group information to factor
sample_info$GSE <- as.factor(sample_info$GSE)       # Convert batch information to factor

# Step 5: Construct the design matrix
design <- model.matrix(~ sample_info$group)  # Construct the design matrix containing the group variable

# Step 6: Remove batch effects using limma
expr_matrix_clean <- removeBatchEffect(expr_matrix_log, 
                                       batch = sample_info$GSE,  # Specify batch information
                                       design = design)          # Specify the design matrix to preserve group information

# Convert to data frame and save the result
limma_df <- as.data.frame(expr_matrix_clean)
limma_df <- cbind(gene_name = rownames(limma_df), limma_df)

output_file <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect.csv"
fwrite(limma_df , file = output_file, row.names = FALSE, sep = ",")