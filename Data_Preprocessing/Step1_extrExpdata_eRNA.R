# Three eRNA annotation sources

## Wrapper function: read StringTie gene_abund.tab files for one annotation source
readStringtieData <- function(pathPattern) {
   files <- Sys.glob(pathPattern)
   if (length(files) == 0) {
      stop(paste("No files found:", pathPattern))
   }
   
   ## Read the first file
   expData <- read.delim(file = files[1], header = TRUE, stringsAsFactors = FALSE)
   expData <- expData[, c(1, 2, 9)]  # gene_id, gene_name, TPM column
   
   ## Extract sample ID
   posIndex <- regexpr('stringtie.*gene_abund', files[1])
   sampleid <- substring(files[1], posIndex + 10, posIndex + attr(posIndex, 'match.length') - 12)
   
   ## Retain unique genes
   uniqIndex <- which(!duplicated(as.character(expData[,1])))
   uniqEnsemids <- as.character(expData[uniqIndex, 1])
   expData <- expData[uniqIndex, ]
   
   ## Initialize column names
   columnNames <- c("gene_id", "gene_name", sampleid)
   
   ## Read the remaining files
   if (length(files) > 1) {
      for (filename in files[2:length(files)]) {
         if (file.info(filename)$size > 0) {  # Check file size
            sampleMat <- read.delim(file = filename, header = TRUE, stringsAsFactors = FALSE)
            matchIndex <- match(uniqEnsemids, sampleMat[, 1])
            expData <- cbind(expData, sampleMat[matchIndex, 9])
            
            posIndex <- regexpr('stringtie.*gene_abund', filename)
            sampleid <- substring(filename, posIndex + 10, posIndex + attr(posIndex, 'match.length') - 11)
            
            cat("Reading sample:", sampleid, "\n")
            columnNames <- append(columnNames, sampleid)
         }
      }
   }
   
   colnames(expData) <- columnNames
   return(expData)
}

## ---- Read three annotation sources ---- ##
expData1 <- readStringtieData("/home/yjliu/mmProj/homo/procdata/GSE*/eRNA/EnhancerAtlas/stringtie/*/gene_abund.tab")
expData2 <- readStringtieData("/home/yjliu/mmProj/homo/procdata/GSE*/eRNA/Ensembl/stringtie/*/gene_abund.tab")
expData3 <- readStringtieData("/home/yjliu/mmProj/homo/procdata/GSE*/eRNA/FANTOM5/stringtie/*/gene_abund.tab")

## ---- Automatically align sample column order ---- ##
alignSamples <- function(df, targetNames) {
   # targetNames includes gene_id and gene_name columns
   df <- df[, c("gene_id", "gene_name", targetNames[ !(targetNames %in% c("gene_id", "gene_name")) ])]
   return(df)
}

## Get reference column names with unified sample names
refCols <- colnames(expData1)

expData2 <- alignSamples(expData2, refCols)
expData3 <- alignSamples(expData3, refCols)

## ---- Merge datasets ---- ##
expData_all <- rbind(expData1, expData2, expData3)
colnames(expData_all) <- gsub("/", "", colnames(expData_all))

## ---- Export merged expression matrix ---- ##
write.csv(expData_all, 
          file="/home/yjliu/mmProj/homo/procdata/Homo_eRNA_TPM.csv", 
          row.names=FALSE, 
          quote=FALSE)