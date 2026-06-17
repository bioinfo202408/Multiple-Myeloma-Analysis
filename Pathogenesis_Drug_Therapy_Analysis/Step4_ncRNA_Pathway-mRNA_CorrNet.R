## Data preprocessing
library(vroom)
library(org.Hs.eg.db)
library(data.table)
expData <- vroom("/home/yjliu/mmProj/data_process/Human/TPM/Standardization/mRNA/homo_mRNA_TPM_GEO_GTEx_2133training_Filter_scaled.csv", 
                 delim = ",", col_names = TRUE)
expData <- as.data.frame(expData)
rownames(expData) <- make.unique(expData[[1]])
expData <- expData[, -1]

EnsembelIDs <- rownames(expData)
EnsembelIDs <- sub("\\.[0-9]+$", "", EnsembelIDs)  # Standard method to remove transcript version suffix
rownames(expData) <- EnsembelIDs

# Map Ensembl ID to official gene symbol
gene_map <- mapIds(
   org.Hs.eg.db,
   keys = EnsembelIDs,
   column = "SYMBOL",
   keytype = "ENSEMBL",
   multiVals = "first"
)
# Optional: use gene symbol as row names
# Use make.unique() to resolve duplicated row names from gene mapping results
gene_map[is.na(gene_map)] <- "Unknown"
rownames(expData) <- make.unique(gene_map)

Pathway_mRNA <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MMrelate_Pathway_GeneSymbol.csv")
Pathway_mRNA <- Pathway_mRNA$symbol

expData_matched <- expData[rownames(expData) %in% Pathway_mRNA,]
#expData_matched <- expData_matched[match(Pathway_mRNA, rownames(expData_matched)) ,]

output_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/434PathwayRelate_mRNA_TumorHealth.csv"
fwrite(expData_matched, file = output_file, row.names = TRUE, sep = ",")

# miRNA processing
expData_miRNA <- vroom("/home/yjliu/mmProj/data_process/Human/TPM/Standardization/miRNA/homo_miRNA_TPM_GEO_GTEx_2133training_Filter_scaled.csv", 
                       delim = ",", col_names = TRUE)
expData_miRNA <- as.data.frame(expData_miRNA)
rownames(expData_miRNA) <- make.unique(expData_miRNA[[1]])
expData_miRNA <- expData_miRNA[, -1]

miRNA <- read.table("/home/yjliu/mmProj/data_process/Human/Feature_select/miRNA/tuning_groupkfold/ncRNA_all_regions_pf0.6_20251119_093042/stable_features_freq_ge9.txt")
miRNA <- miRNA$V1
expData_miRNA_matched <- expData_miRNA[rownames(expData_miRNA) %in% miRNA,]

output_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/207_Feature-miRNA_TumorHealth.csv"
fwrite(expData_miRNA_matched, file = output_file, row.names = TRUE, sep = ",")

# eRNA processing
expData_eRNA <- vroom("/home/yjliu/mmProj/data_process/Human/TPM/Standardization/eRNA/homo_eRNA_TPM_GEO_GTEx_2133training_Filter_scaled.csv", 
                      delim = ",", col_names = TRUE)
expData_eRNA <- as.data.frame(expData_eRNA)
rownames(expData_eRNA) <- make.unique(expData_eRNA[[1]])
expData_eRNA <- expData_eRNA[, -1]

eRNA <- read.table("/home/yjliu/mmProj/data_process/Human/Feature_select/eRNA/tuning_groupkfold/ncRNA_all_regions_pf0.85_20251125_085247/stable_features_freq_ge9.txt")
eRNA <- eRNA$V1
expData_eRNA_matched <- expData_eRNA[rownames(expData_eRNA) %in% eRNA,]

output_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/457_Feature-eRNA_TumorHealth.csv"
fwrite(expData_eRNA_matched, file = output_file, row.names = TRUE, sep = ",")


# lncRNA processing
expData_lncRNA <- vroom("/home/yjliu/mmProj/data_process/Human/TPM/Standardization/lncRNA/homo_lncRNA_TPM_GEO_GTEx_2133training_Filter_scaled.csv", 
                        delim = ",", col_names = TRUE)
expData_lncRNA <- as.data.frame(expData_lncRNA)
rownames(expData_lncRNA) <- make.unique(expData_lncRNA[[1]])
expData_lncRNA <- expData_lncRNA[, -1]

lncRNA <- read.table("/home/yjliu/mmProj/data_process/Human/Feature_select/lncRNA/tuning_groupkfold/ncRNA_all_regions_pf0.9_20251124_213912/stable_features_freq_ge9.txt")
lncRNA <- lncRNA$V1
expData_lncRNA_matched <- expData_lncRNA[rownames(expData_lncRNA) %in% lncRNA,]

output_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/498_Feature-lncRNA_TumorHealth.csv"
fwrite(expData_lncRNA_matched, file = output_file, row.names = TRUE, sep = ",")

# Merge all ncRNA types together
ncRNA <- rbind(expData_miRNA_matched , expData_lncRNA_matched , expData_eRNA_matched)
output_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/1162_Feature-ncRNA_TumorHealth.csv"
fwrite(ncRNA, file = output_file, row.names = TRUE, sep = ",")

## End of data preprocessing pipeline


# Input file paths
nc_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/1162_Feature-ncRNA_TumorHealth.csv"
mRNA_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/434PathwayRelate_mRNA_TumorHealth.csv"
out_file <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/ncRNA_Pathway-mRNA_CorrNet.csv"

# Parameter settings
# -----------------------------
strip_version <- FALSE   # Whether to remove Ensembl transcript version suffix
remove_self <- FALSE     # Whether to exclude self-interaction edges

# -----------------------------
# Core correlation network construction function (matrix-based efficient calculation)
# -----------------------------
build_ncRNA_targets_network <- function(ncRNA_expr, all_expr, remove_self = TRUE) {
   
   ncRNA_expr <- as.matrix(ncRNA_expr)
   all_expr   <- as.matrix(all_expr)
   
   # 1) Check consistency of sample column names
   if (!identical(colnames(ncRNA_expr), colnames(all_expr))) {
      stop("Sample columns are not identical (or not in the same order) between ncRNA_expr and all_expr.")
   }
   
   # 2) Remove zero-variance genes (cannot compute valid correlation coefficient)
   ncRNA_expr <- ncRNA_expr[apply(ncRNA_expr, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
   all_expr   <- all_expr[apply(all_expr,   1, sd, na.rm = TRUE) > 0, , drop = FALSE]
   
   n <- ncol(ncRNA_expr)
   if (n < 4) stop("Need at least 4 samples for reliable Pearson correlation significance test.")
   
   # 3) Calculate Pearson correlation matrix: rows = ncRNA, columns = mRNA targets
   r_mat <- cor(t(ncRNA_expr), t(all_expr),
                method = "pearson", use = "pairwise.complete.obs")
   
   # 4) Vectorized calculation of two-tailed P values from correlation coefficients
   df <- n - 2
   r_vec <- as.vector(r_mat)
   t_vec <- r_vec * sqrt(df / pmax(1 - r_vec^2, .Machine$double.eps))
   p_vec <- 2 * pt(-abs(t_vec), df = df)
   
   # Reshape flat P-value vector back to matrix matching correlation matrix dimension
   p_mat <- matrix(p_vec,
                   nrow = nrow(r_mat), ncol = ncol(r_mat),
                   byrow = FALSE, dimnames = dimnames(r_mat))
   
   # 5) Row-wise FDR correction (Benjamini-Hochberg) for each ncRNA's target correlations
   fdr_mat <- t(apply(p_mat, 1, p.adjust, method = "BH"))
   
   # 6) Convert full correlation matrix into long-format edge table containing all pairs
   edges <- data.frame(
      ncRNA  = rep(rownames(r_mat), times = ncol(r_mat)),
      target = rep(colnames(r_mat), each  = nrow(r_mat)),
      r      = r_vec,
      p      = as.vector(p_mat),
      FDR    = as.vector(fdr_mat),
      stringsAsFactors = FALSE
   )
   
   # Optional filter: remove self-matching edges if ncRNA identifiers overlap with mRNA names
   if (remove_self) edges <- edges[edges$ncRNA != edges$target, ]
   
   return(edges)
}


# -----------------------------
# Load input expression matrices
# -----------------------------
nc_expr <- read.table(nc_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
all_expr <- read.table(mRNA_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# all_expr <- as.data.frame(all_expr)
# rownames(all_expr) <- make.unique(all_expr[[1]])
# all_expr <- all_expr[, -1]

# Optional: strip Ensembl version suffix from gene row names
if (strip_version) {
   rownames(nc_expr)  <- gsub("\\.\\d+$", "", rownames(nc_expr))
   rownames(all_expr) <- gsub("\\.\\d+$", "", rownames(all_expr))
}

# -----------------------------
# Run correlation network construction
# -----------------------------
edges <- build_ncRNA_targets_network(nc_expr, all_expr, remove_self = remove_self)

# -----------------------------
# Export raw full correlation network table
# -----------------------------
#write.table(edges, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(edges, file = out_file, row.names = FALSE)
cat("Correlation network construction finished!\n")
cat("Total ncRNA input rows:", nrow(nc_expr), " | Valid ncRNA after zero-variance filter:", length(unique(edges$ncRNA)), "\n")
cat("Total mRNA target input rows:", nrow(all_expr), " | Valid mRNA targets after zero-variance filter:", length(unique(edges$target)), "\n")
cat("Total correlation edges in output:", nrow(edges), "\n")
cat("Output file path:", out_file, "\n")

### Filter co-expression network by ncRNA type-specific thresholds
library(dplyr)
library(stringr)

# Load raw full correlation network
edges <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/ncRNA_Pathway-mRNA_CorrNet.csv")

# Classify ncRNA types and apply correlation & significance filters separately
edges_filtered <- edges %>%
   mutate(
      type = case_when(
         str_detect(ncRNA, regex("^MI", ignore_case = TRUE)) ~ "miRNA",
         str_detect(ncRNA, regex("^hsa-mir", ignore_case = TRUE)) ~ "miRNA",
         str_detect(ncRNA, regex("^ENSG", ignore_case = TRUE)) ~ "lncRNA",
         str_detect(ncRNA, regex("^NONH", ignore_case = TRUE)) ~ "lncRNA",
         TRUE ~ "eRNA"
      )
   ) %>%
   filter(
      # miRNA: significant negative correlation only
      (type == "miRNA" & r < -0.2 & FDR < 0.05) |
         # lncRNA: significant positive/negative moderate correlation
         (type == "lncRNA" & abs(r) > 0.3 & FDR < 0.05) |
         # eRNA: significant positive correlation only
         (type == "eRNA" & r > 0.3 & FDR < 0.05)
   ) %>%
   mutate(
      FDR_safe = ifelse(FDR == 0, 1e-300, FDR),  # Avoid computational error from log10(0)
      weight = abs(r) * (-log10(FDR_safe))  # Composite edge weight combining correlation strength and significance
   ) %>%
   select(-FDR_safe, -type)  # Drop temporary auxiliary columns

# Export filtered weighted regulatory network
write.csv(edges_filtered, 
          file = "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/ncRNA_Pathway-mRNA_filtered_WeightNetwork.csv",
          row.names = FALSE)

cat("Filtering completed!\n")
cat("Retained significant regulatory edges after filtering:", nrow(edges_filtered), "\n")


## Convert miRNA MIMAT accession IDs to standard miRNA names
library(miRBaseConverter)
# Map MIMAT accession numbers to official mature miRNA names
converted <- miRNA_AccessionToName(edges_filtered$ncRNA)
edges_filtered$miRNA_name <- converted$TargetName

## Extract sub-network targeting 12 core pathogenic/drug-relevant genes
library(dplyr)
library(tidyr)
# List of 12 core multiple myeloma target genes
genes12 <- c(
   "ABCC3", "ABCC9", "AR", "CASR", "CCND1", "FGF2",
   "FGFR4", "FOXO1", "HGF", "NTF3", "TGM2", "TNFRSF11B"
)

# 1) Subset network edges where mRNA target belongs to the 12 core gene list
edges_12gene <- edges_filtered %>%
   filter(target %in% genes12)

write.csv(edges_12gene, file = "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/ncRNA_Pathway-mRNA_Corr/ncRNA_PathwayDrugTarget-12mRNA_filtered_WeightNetwork.csv", row.names = FALSE)