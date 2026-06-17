## Return here after running AvgExpre_Propor_JS_threshold.R

expression_file <- "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect_2133training_data_newCorName.csv"

# Read expression data
expression_data <- as.data.frame(fread(expression_file))
setnames(expression_data, 1, "Gene")  # The first column contains gene IDs
rownames(expression_data) <- expression_data[[1]]
expression_data <- expression_data[, -1]

# Read JS score data
js_file <- "/home/yjliu/mmProj/JSscore/mRNA_JSscore.txt"  
js_scores_raw <- fread(js_file)
setnames(js_scores_raw, old = c("V1", "V2"), new = c("Gene", "JSScores"))
js_scores_map <- js_scores_raw[, .(Gene, JSScores)]  # Convert to mapping format

# Match genes with JS scores
js_scores_vec <- setNames(js_scores_map$JSScores, js_scores_map$Gene)
js_scores_matched <- js_scores_vec[rownames(expression_data)]

# Filter genes with JSScores >= XX
filtered_genes <- names(js_scores_matched[js_scores_matched >= 0.55 & !is.na(js_scores_matched)])
filtered_expression_data <- expression_data[rownames(expression_data) %in% filtered_genes, , drop = FALSE]

# Calculate the mean value for each row
row_means <- rowMeans(filtered_expression_data)

filtered_data <- filtered_expression_data[row_means >= 0, ] # Enter the selected threshold

# Count the number of values greater than 0 in each row
row_positive_count <- rowSums(filtered_data > 0)

# Calculate the proportion
row_positive_ratio <- row_positive_count / 2133 # Modify 2133 according to the actual sample size

# Retain rows with a proportion greater than XX
final_filtered_data <- filtered_data[row_positive_ratio > 0.9, ] # Enter the selected threshold

write.csv(final_filtered_data, file="/home/yjliu/mmProj/homo/procdata/mRNA/Homo_mRNA_TPM_removeBatchEffect_2133training_data_newCorName_FilterJS_AvgExpres_Proportion.csv", row.names=T)