# Drug Target Enrichment Analysis

# First extract genes from pathways
GO <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy//Enrich_Analysis/MM_related_GO_filtered.csv")

KEGG <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy//Enrich_Analysis/MM_related_KEGG_filtered.csv")

HallMark <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy//Enrich_Analysis/MM_related_HallMark_filtered.csv")

# 1️⃣ Define a function to split entries in the geneID column
extract_genes <- function(gene_column) {
   genes <- unlist(strsplit(as.character(gene_column), "/"))
   return(genes)
}

# 2️⃣ Extract gene lists from three enrichment datasets
genes_GO <- extract_genes(GO$geneID)
genes_KEGG <- extract_genes(KEGG$geneID)
genes_HallMark <- extract_genes(HallMark$geneID)

# 3️⃣ Merge all genes and retain unique union set
all_genes <- unique(c(genes_GO, genes_KEGG, genes_HallMark))

suppressPackageStartupMessages({
   library(dplyr)
   library(readr)
   library(vroom)
   library(AnnotationDbi)
   library(org.Hs.eg.db)
})

dir.create("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Fisher_DrugEnrichment", 
           recursive = TRUE, showWarnings = FALSE)

# # ===================== 1. Input Files =====================
# # 1) Feature gene list
# key_gene_file <- "/home/yjliu/mmProj/data_process/Human/TPM/Standardization/mRNA/cluster1_genes..txt"

# 2) Integrated target genes of multiple myeloma therapeutic drugs
mm_target_file <- "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/MM_target_genes_integrated.tsv"

# 3) Gene expression matrix (all genes used as background population)
expr_file <- "/home/yjliu/mmProj/data_process/Human/TPM/mRNA/homo_mRNA_TPM_GEO_GTEx_removeBatchEffect_keep_negative.csv"

# Output directory for all results
outdir <- "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Fisher_DrugEnrichment/"

# ===================== 2. Load pathway-derived feature genes =====================
key_mRNA <- all_genes
key_mRNA <- as.data.frame(key_mRNA)
colnames(key_mRNA)[1] <- "symbol"
if (!"symbol" %in% colnames(key_mRNA)) {
   stop("Symbol column not found in feature gene file.")
}
write.csv(key_mRNA,"/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MMrelate_Pathway_GeneSymbol.csv",row.names = F)
key_genes <- key_mRNA |>
   transmute(symbol = toupper(trimws(as.character(symbol)))) |>
   filter(!is.na(symbol), symbol != "") |>
   distinct() |>
   pull(symbol)

# ===================== 3. Load integrated MM drug target genes =====================
mm_target_df <- read_tsv(mm_target_file, show_col_types = FALSE)

# Automatically detect target gene column: prioritize target_gene, then gene_symbol
target_col <- NULL
if ("target_gene" %in% colnames(mm_target_df)) {
   target_col <- "target_gene"
} else if ("gene_symbol" %in% colnames(mm_target_df)) {
   target_col <- "gene_symbol"
} else {
   stop("Neither target_gene nor gene_symbol column found in MM target gene file.")
}

mm_target_genes <- mm_target_df |>
   transmute(symbol = toupper(trimws(as.character(.data[[target_col]])))) |>
   filter(!is.na(symbol), symbol != "") |>
   distinct() |>
   pull(symbol)
library(vroom)
# ===================== 4. Build whole-genome background gene set =====================
expr <- vroom(expr_file, delim = ",", col_names = TRUE, progress = FALSE)
expr <- as.data.frame(expr)

# The first column is defined as Ensembl gene ID by default
rownames(expr) <- make.unique(as.character(expr[[1]]))
expr <- expr[, -1, drop = FALSE]

# Strip Ensembl version suffix
ensembl_ids <- rownames(expr)
ensembl_ids <- sub("\\.[0-9]+$", "", ensembl_ids)

# Map Ensembl IDs to official gene symbols
gene_map <- mapIds(
   org.Hs.eg.db,
   keys = ensembl_ids,
   column = "SYMBOL",
   keytype = "ENSEMBL",
   multiVals = "first"
)

bg_genes <- as.character(gene_map)
bg_genes <- toupper(trimws(bg_genes))
bg_genes <- bg_genes[!is.na(bg_genes) & bg_genes != ""]
bg_genes <- unique(bg_genes)

# ===================== 5. Filter all gene sets to retain only genes present in background =====================
key_genes_use <- intersect(key_genes, bg_genes)
mm_targets_use <- intersect(mm_target_genes, bg_genes)
overlap_genes <- intersect(key_genes_use, mm_targets_use)

# Export overlapping gene list for supplementary table
write_tsv(
   tibble(overlap_genes = sort(overlap_genes)),
   file.path(outdir, "overlap_key_mRNA_MM_targets_integrated.tsv")
)

# ===================== 6. Two-tailed Fisher's Exact Test for enrichment =====================
N <- 19421         # Total number of background human genes
K <- length(mm_targets_use)   # Count of validated MM drug target genes
n <- length(key_genes_use)    # Count of pathway feature genes
k <- length(overlap_genes)    # Overlap count between feature genes and drug targets

if (N == 0 || K == 0 || n == 0) {
   stop("Background gene set, MM target gene set or feature gene set is empty; Fisher test cannot be executed.")
}

# Construct 2x2 contingency table
fisher_mat <- matrix(
   c(
      k,          n - k,
      K - k,      N - K - n + k
   ),
   nrow = 2,
   byrow = TRUE
)

fisher_res <- fisher.test(fisher_mat, alternative = "greater")

# Compile enrichment statistics table
enrich_res <- tibble(
   N_background = N,
   K_MM_targets = K,
   n_key_genes = n,
   k_overlap = k,
   GeneRatio = paste0(k, "/", n),
   BgRatio = paste0(K, "/", N),
   OddsRatio = unname(fisher_res$estimate),
   Pvalue = fisher_res$p.value
)

write_tsv(
   enrich_res,
   file.path(outdir, "MM_target_enrichment_result_integrated.tsv")
)

# Optional: Export filtered gene sets for manual inspection
write_tsv(tibble(key_genes_use = sort(key_genes_use)),
          file.path(outdir, "key_genes_in_background.tsv"))

write_tsv(tibble(mm_targets_use = sort(mm_targets_use)),
          file.path(outdir, "MM_targets_in_background.tsv"))

# ===================== 7. Print enrichment summary to console =====================
cat("===== MM Drug Target Enrichment Summary (Integrated Database Targets) =====\n")
cat("Total background genes: ", N, "\n")
cat("Validated MM drug target genes: ", K, "\n")
cat("Pathway feature genes: ", n, "\n")
cat("Overlapping genes between feature and drug targets: ", k, "\n")
cat("Fisher's exact P value: ", fisher_res$p.value, "\n")
cat("Odds Ratio: ", unname(fisher_res$estimate), "\n")