library(vroom)
enhExpDataMap <- vroom(
   "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_removeBatchEffect_2133training_data.csv", 
   delim = ",", 
   col_names = TRUE
)
enhExpData <- enhExpDataMap
enhExpData <- as.data.frame(enhExpData)
rownames(enhExpData) <- make.unique(enhExpData[[1]])
enhExpData <- enhExpData[,-1]

# 2. Read the metadata file
Group <- read.csv("/home/yjliu/mmProj/homo_metadata_training_test_data_2133training_data.csv", header = T)
rownames(Group) <- Group$Run
Group <- Group[colnames(enhExpData), drop = FALSE]
tumor_samples  <- rownames(Group)[ Group$group == "tumor"  ]
health_samples <- rownames(Group)[ Group$group == "health" ]

# 3. Split the expression matrix according to sample names
enhExp_tumor  <- enhExpData[, tumor_samples]
enhExp_health <- enhExpData[, health_samples]

# 4. Optional: check dimensions
dim(enhExp_tumor)   # Genes x samples in the tumor group
dim(enhExp_health)  # Genes x samples in the healthy group

# 1. Calculate the average expression of each gene in the tumor group
mean_tumor <- rowMeans(enhExp_tumor)

# 2. Calculate the average expression of each gene in the healthy group
mean_health <- rowMeans(enhExp_health)

# 3. Combine into a new matrix with columns named "tumor" and "health"
avgExp <- cbind(
   tumor  = mean_tumor,
   health = mean_health
)
avgExp[avgExp < 0] <- 0

# View the result
head(avgExp)

write.csv(avgExp,
          file = "/home/yjliu/mmProj/homo/procdata/Homo_mRNA_TPM_2group.csv",
          row.names = TRUE,
          quote = FALSE)