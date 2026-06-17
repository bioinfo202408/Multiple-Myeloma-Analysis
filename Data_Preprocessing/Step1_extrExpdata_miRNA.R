# miRNA annotations from two sources

## A function to read gene_abund.tab files from a specific annotation source and return an expression matrix
readStringtieData <- function(pathPattern) {
   files <- Sys.glob(pathPattern)
   
   # Read the first file
   expData <- read.delim(file=files[1], header=TRUE, stringsAsFactors=FALSE)
   expData <- expData[, c(1,2,9)]  # gene_id, gene_name, TPM column
   
   # Extract sample ID
   posIndex <- regexpr('stringtie.*gene_abund', files[1])
   sampleid <- substring(files[1], posIndex+10, posIndex+attr(posIndex,'match.length')-12)
   
   # Retain unique genes
   uniqIndex <- which(!duplicated(as.character(expData[,1])))
   uniqEnsemids <- as.character(expData[uniqIndex, 1])
   expData <- expData[uniqIndex, ]
   
   columnNames <- c("gene_id", "gene_name", sampleid)
   
   # Read remaining files
   if (length(files) > 1) {
      for (filename in files[2:length(files)]) {
         if (file.info(filename)$size > 0) {  # Check file size
            sampleMat <- read.delim(file=filename, header=TRUE, stringsAsFactors=FALSE)
            matchIndex <- match(uniqEnsemids, sampleMat[,1])
            expData <- cbind(expData, sampleMat[matchIndex,9])
            
            posIndex <- regexpr('stringtie.*gene_abund', filename)
            sampleid <- substring(filename, posIndex+10, posIndex+attr(posIndex,'match.length')-11)
            cat("Reading sample:", sampleid, "\n")
            
            columnNames <- append(columnNames, sampleid)
         }
      }
   }
   
   colnames(expData) <- columnNames
   return(expData)
}

## Read miRBase annotation
expData1 <- readStringtieData("/home/yjliu/mmProj/homo/procdata/GSE*/miRNA/miRBase/stringtie/*/gene_abund.tab")

## Read MirGeneDB annotation
expData2 <- readStringtieData("/home/yjliu/mmProj/homo/procdata/GSE*/miRNA/MirGeneDB/stringtie/*/gene_abund.tab")

## Ensure column names are identical (excluding gene_id/gene_name)
if (!identical(colnames(expData1), colnames(expData2))) {
   stop("Sample column names differ between the two annotations. Please check before performing rbind.")
}

## Merge the two matrices
expData_all <- rbind(expData1, expData2)
colnames(expData_all) <- gsub("/", "", colnames(expData_all))

## ---- Export ---- ##
write.csv(expData_all, 
          file="/home/yjliu/mmProj/homo/procdata/Homo_miRNA_TPM.csv", 
          row.names=FALSE, 
          quote=FALSE)