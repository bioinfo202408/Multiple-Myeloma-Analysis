# Clear the environment and initialize settings
rm(list = ls())
gc()
options(scipen = 20)

library(vroom)

expression_data <- vroom("/home/yjliu/mmProj/data_process/Human/TPM/mRNA/homo_mRNA_TPM_removeBatchEffect.csv", 
                         delim = ",", col_names = TRUE)

expression_data <- as.data.frame(expression_data)
rownames(expression_data) <- make.unique(expression_data[[1]])
expression_data <- expression_data[, -1]

metadata <- read.csv("/home/yjliu/mmProj/homo_metadata_training_test_data_2133training_data.csv", header = TRUE)
matched_indices <- match(colnames(expression_data), metadata$Run)

if (any(is.na(matched_indices))) {
   warning("Some column names were not matched. Please check whether the Run column in metadata is consistent with the column names in expression_data.")
}

new_expression_data <- expression_data
colnames(new_expression_data) <- ifelse(
   !is.na(matched_indices),
   paste0(metadata$Run[matched_indices], "_", metadata$Group[matched_indices]),
   colnames(new_expression_data)
)

write.csv(new_expression_data, file="/home/yjliu/mmProj/data_process/Human/TPM/mRNA/homo_mRNA_TPM_removeBatchEffect_newCorName.csv", row.names=TRUE)


# Load required packages
library(data.table)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# Read expression data
expression_file <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect_2133training_data_newCorName.csv"
expression_data <- fread(expression_file)
setnames(expression_data, 1, "Gene")  # The first column contains gene IDs

# Extract sample names and group information
all_samples <- colnames(expression_data)[-1]
sample_groups <- ifelse(grepl("_health·$", all_samples), "health", "tumor")
group_info <- data.frame(Sample = all_samples, Group = sample_groups, stringsAsFactors = FALSE)

# Read JSscore data
js_file <- "/home/yjliu/mmProj/JSscore/mRNA_JSscore.txt"  
js_scores_raw <- fread(js_file)
setnames(js_scores_raw, old = c("V1", "V2"), new = c("Gene", "JSScores"))
js_scores_map <- js_scores_raw[, .(Gene, JSScores)]  # Convert to mapping format

# Define gene filtering function supporting multiple parameter ranges
filter_ncRNA <- function(expr_data, avg_expression_threshold, min_sample_proportion, js_threshold) {
   numeric_matrix <- as.matrix(expr_data[, -1, with = FALSE])
   rownames(numeric_matrix) <- expr_data[[1]]
   
   avg_expression_values <- rowMeans(numeric_matrix, na.rm = TRUE)
   num_samples <- ncol(numeric_matrix)
   valid_in_samples <- rowSums(numeric_matrix > 0, na.rm = TRUE)
   
   js_scores_vec <- setNames(js_scores_map$JSScores, js_scores_map$Gene)
   js_scores_matched <- js_scores_vec[rownames(numeric_matrix)]
   
   keep <- (avg_expression_values >= avg_expression_threshold) &
      (valid_in_samples >= min_sample_proportion * num_samples) &
      (!is.na(js_scores_matched) & js_scores_matched >= js_threshold)
   
   rownames(numeric_matrix)[keep]
}

# Define within-group correlation evaluation function
evaluate_subgroup_quality <- function(filtered_expression_matrix, group_info) {
   unique_groups <- unique(group_info$Group)
   group_similarities <- list()
   for (g in unique_groups) {
      group_samples <- group_info$Sample[group_info$Group == g]
      sub_data <- filtered_expression_matrix[, group_samples, drop = FALSE]
      cormat <- cor(sub_data, method = "pearson", use = "pairwise.complete.obs")
      cluster_quality <- mean(abs(cormat[upper.tri(cormat)]), na.rm = TRUE)
      group_similarities[[g]] <- cluster_quality
   }
   mean(unlist(group_similarities), na.rm = TRUE)
}

# Parameter ranges based on the median with an appropriate range and step size
avg_expression_threshold_range <- seq(0, 0.02, by = 0.005)
min_sample_proportion_range <- seq(0, 0.1, by = 0.02)
js_threshold_range <- seq(0.5, 0.7, by = 0.02)

# Set up the parallel environment
num_cores <- 50  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, c("filter_ncRNA", "evaluate_subgroup_quality", "expression_data", "group_info", "js_scores_map"),
              envir = environment())

clusterEvalQ(cl, {
   library(data.table)
   library(dplyr)
})

# Grid search
results_list <- foreach(avg_th = avg_expression_threshold_range, .combine = rbind, 
                        .packages = c("data.table", "dplyr")) %dopar% {
                           
                           inner_res_list <- list()
                           
                           for (min_prop in min_sample_proportion_range) {
                              for (js_th in js_threshold_range) {
                                 filtered_genes <- filter_ncRNA(
                                    expr_data = expression_data,
                                    avg_expression_threshold = avg_th,
                                    min_sample_proportion = min_prop,
                                    js_threshold = js_th
                                 )
                                 
                                 if (length(filtered_genes) == 0) {
                                    inner_res_list[[length(inner_res_list) + 1]] <- data.frame(
                                       Avg_Expression_Threshold = avg_th,
                                       Min_Sample_Proportion = min_prop,
                                       JS_Threshold = js_th,
                                       Num_Genes = 0,
                                       Cluster_Quality = NA,
                                       Status = "No genes selected",
                                       stringsAsFactors = FALSE
                                    )
                                    next
                                 }
                                 
                                 filtered_expression_matrix <- expression_data[Gene %in% filtered_genes, ]
                                 mat_filtered <- as.matrix(filtered_expression_matrix[, -1, with = FALSE])
                                 rownames(mat_filtered) <- filtered_expression_matrix$Gene
                                 
                                 if (nrow(mat_filtered) < 5) {
                                    inner_res_list[[length(inner_res_list) + 1]] <- data.frame(
                                       Avg_Expression_Threshold = avg_th,
                                       Min_Sample_Proportion = min_prop,
                                       JS_Threshold = js_th,
                                       Num_Genes = nrow(mat_filtered),
                                       Cluster_Quality = NA,
                                       Status = "Too few genes (<5)",
                                       stringsAsFactors = FALSE
                                    )
                                    next
                                 }
                                 
                                 cluster_quality <- evaluate_subgroup_quality(
                                    filtered_expression_matrix = mat_filtered, 
                                    group_info = group_info
                                 )
                                 
                                 inner_res_list[[length(inner_res_list) + 1]] <- data.frame(
                                    Avg_Expression_Threshold = avg_th,
                                    Min_Sample_Proportion = min_prop,
                                    JS_Threshold = js_th,
                                    Num_Genes = nrow(mat_filtered),
                                    Cluster_Quality = cluster_quality,
                                    Status = "Done",
                                    stringsAsFactors = FALSE
                                 )
                              }
                           }
                           
                           do.call(rbind, inner_res_list)
                        }

# Save results
summary_file <- "/home/yjliu/mmProj/homo/procdata/mRNA/Filter_Avg_Propo_JS_result.csv"
fwrite(results_list, file = summary_file, row.names = FALSE, sep = ",")
cat("Filtering summary has been saved to:", summary_file, "\n")

# Stop the parallel environment
stopCluster(cl)
cat("Parallel environment has been stopped.\n")