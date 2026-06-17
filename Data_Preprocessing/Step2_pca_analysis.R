library(ggplot2)
library(vroom)
library(FactoMineR)
library(factoextra)

# Read and preprocess expression data
ExpData <- vroom(
   "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect.csv", 
   delim = ","
)
ExpData <- as.data.frame(ExpData)
rownames(ExpData) <- make.unique(ExpData[[1]])
ExpData <- ExpData[,-1]
ExpData <- ExpData[rowSums(ExpData) > 0, ]

# Read and match metadata
metadata <- read.csv("/home/yjliu/mmProj/homo_metadata_training_test_data.csv", header = TRUE)
metadata <- metadata[match(colnames(ExpData), metadata$Run), ]

# PCA analysis
pca_res <- prcomp(t(ExpData), scale. = TRUE)
scores <- as.data.frame(pca_res$x)
scores$Group <- factor(metadata$group)

# Manually define color mapping: health = green, tumor = red
group_colors <- c("health" = "#5fbf7a", "tumor" = "#ED7675")

# Create PCA plotting function
create_pca_plot <- function(scores, pc_x, pc_y, pca_res, group_var = "Group") {
   var_x <- round(summary(pca_res)$importance[2, pc_x] * 100, 1)
   var_y <- round(summary(pca_res)$importance[2, pc_y] * 100, 1)
   
   ggplot(scores, aes(x = .data[[paste0("PC", pc_x)]], 
                      y = .data[[paste0("PC", pc_y)]], 
                      color = .data[[group_var]])) +
      geom_point(size = 3) +
      theme_bw() +
      scale_color_manual(values = group_colors) +  # Manually define colors
      xlab(paste0("PC", pc_x, " (", var_x, "%)")) +
      ylab(paste0("PC", pc_y, " (", var_y, "%)")) 
}

# Generate PCA plots
pdf("/home/yjliu/mmProj/homo/procdata/mRNA_PCA_pc1pc2.pdf", width = 6, height = 5)
print(create_pca_plot(scores, 1, 2, pca_res))
dev.off()

pdf("/home/yjliu/mmProj/homo/procdata/mRNA_PCA_pc1pc3.pdf", width = 6, height = 5)
print(create_pca_plot(scores, 1, 3, pca_res))
dev.off()

pdf("/home/yjliu/mmProj/homo/procdata/mRNA_PCA_pc2pc3.pdf", width = 6, height = 5)
print(create_pca_plot(scores, 2, 3, pca_res))
dev.off()