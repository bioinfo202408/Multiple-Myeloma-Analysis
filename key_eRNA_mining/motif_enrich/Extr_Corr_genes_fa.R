library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)

Corr_genes <- read.table("/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_result/Ens223489_Corr_genes.txt", header=FALSE, stringsAsFactors=FALSE)
IDs <- Corr_genes$V1 

geneIDs <- transId(
   id = IDs,
   transTo = "entrez", org = "human", keepNA = FALSE
)
geneIDs <- na.omit(geneIDs)
entrezIDs <- geneIDs$entrezid

tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
valid_entrezIDs <- entrezIDs[entrezIDs %in% names(tx_by_gene)]
transcriptCoordsByGene.GRangesList <- tx_by_gene[valid_entrezIDs]

filtered_tx <- keepStandardChromosomes(transcriptCoordsByGene.GRangesList, pruning.mode="coarse")

#filtered_tx <- endoapply(filtered_tx, function(gr) gr[start(gr) > 2000])
filtered_tx <- filtered_tx[elementNROWS(filtered_tx) > 0]

promoter.seqs <- getPromoterSeq(filtered_tx, Hsapiens, upstream=5000, downstream=2000)

promoter.seqs <- unlist(promoter.seqs)
writeXStringSet(promoter.seqs, "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/Genes_promoter.fa")

