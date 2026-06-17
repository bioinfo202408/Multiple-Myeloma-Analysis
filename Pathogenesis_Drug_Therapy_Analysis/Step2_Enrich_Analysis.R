library(ggplot2)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
library(readr)
library(stringr)

Genes <- read.table(file="/home/yjliu/mmProj/data_process/Human/Feature_select/mRNA/tuning_groupkfold/mRNA_all_regions_pf0.6_20251124_192149/stable_features_freq_ge9_genekitr.csv",header=T,sep=",",stringsAsFactors=FALSE)

geneIDs <- unique(Genes$symbol)

ego_all <- enrichGO(
   gene          = geneIDs,
   OrgDb         = org.Hs.eg.db,
   keyType       = "SYMBOL",
   ont           = "ALL",          # "BP"/"CC"/"MF"/"ALL"
   pAdjustMethod = "fdr",
   pvalueCutoff  = 0.5,
   readable      = TRUE
)
ego_all <-setReadable(ego_all,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = 'ENTREZID')
go_result = ego_all@result
go_result <- go_result[order(go_result$pvalue),]
#write.csv(go_result,'/home/yjliu/mmProj/data_process/Human/Pathogenesis/Enrich_Analysis/Feature_mRNA_1155.csv')


# KEGG enrichment analysis
ensembl_ids <- unique(Genes$input_id)
geneIDs2 <- transId(
   id = ensembl_ids,
   transTo = "entrez", org = "human", keepNA = FALSE
)
geneIDs2 <- geneIDs2[,2]
KEGG_enrich <- clusterProfiler::enrichKEGG(gene = geneIDs2,
                                           organism = 'hsa',      
                                           keyType = "kegg",      # Input gene ID type
                                           pAdjustMethod = "fdr",  # P-value correction method
                                           pvalueCutoff = 0.5   # P-value threshold
)
# Convert ENTREZ ID to Gene Symbol
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = 'org.Hs.eg.db',
                         keyType = 'ENTREZID')
KEGG_enrich_result = KEGG_enrich@result
#write.csv(KEGG_enrich_result,'/home/yjliu/mmProj/data_process/Human/Pathogenesis/Enrich_Analysis/Feature_1155mRNA_KEGG.csv')


# MSigDB Hallmark gene set enrichment
# Retrieve human Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

library(dplyr)
hallmark_df <- hallmark %>%
   select(gs_name, gene_symbol) %>%
   rename(term = gs_name, gene = gene_symbol)

hallmark_enrich <- enricher(
   gene = geneIDs,   # Differentially expressed gene symbols
   TERM2GENE = hallmark_df,
   pAdjustMethod = "fdr",
   pvalueCutoff = 0.5
)

# Inspect enrichment results
hallmark_enrich_result <- hallmark_enrich@result
#write.csv(hallmark_enrich_result,'/home/yjliu/mmProj/data_process/Human/Pathogenesis/Enrich_Analysis/Feature_1155mRNA_HallMark.csv')


### Plot enrichment figures
######
GO <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/Feature_1155mRNA_GO.csv")
# Filter 133 multiple myeloma-related pathways
GO_Filter <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_GO_filtered.csv")
GOid <- GO_Filter$ID

GO_matched <- GO[ GO$X %in% GOid ,]
# Preserve the original column order from feature_samples if needed
GO_matched <- GO_matched[ match(GOid, GO_matched$X) , ]
write.csv(GO_matched,file =  "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_GO_filtered_full_info.csv")

## Select top 10 most significant pathways for visualization:
# List of target pathway names to extract
target_paths <- c(
   "Wnt signaling pathway",
   "BMP signaling pathway",
   "cell-substrate adhesion",
   "collagen-containing extracellular matrix",
   "bone development",
   "bone mineralization",
   "vascular process in circulatory system",
   "extracellular matrix organization",
   "regulation of vascular permeability",
   "cell-cell signaling by wnt"
)

# Filter matching rows with dplyr
selected_paths <- GO_matched %>%
   filter(Description %in% target_paths)
write.csv(selected_paths,file =  "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_GO_filtered_full_info_Top10Pathway.csv")
# Load required packages
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

# ===================== 2. Data preprocessing =====================
# 2.1 Parse GeneRatio into numeric value (numerator / denominator)
selected_paths <- selected_paths %>%
   mutate(
      GeneRatio_numeric = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,1]) / 
         as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,2]),
      # Calculate -log10(pvalue) for color gradient mapping
      neg_log10_p = -log10(pvalue),
      # Order pathways by GeneRatio for vertical plot arrangement
      Description = factor(Description, levels = rev(unique(Description[order(GeneRatio_numeric)])))
   )

# ===================== 3. Composite plot construction =====================
# Color palette (light green to dark green, mapped to -log10(pvalue))
color_scale <- scale_color_gradient(low = "#e6f7ef", high = "#006633")
fill_scale <- scale_fill_gradient(low = "#e6f7ef", high = "#006633")

p <- ggplot(selected_paths, aes(y = Description)) +
   # ------------- Layer 1: Horizontal bar plot for GeneRatio visualization -------------
geom_col(
   aes(x = GeneRatio_numeric), 
   fill = "#2e8b57",  # Fixed bar green color (adjustable)
   width = 0.8,       # Bar width
   alpha = 0.8        # Bar transparency
) +
   # ------------- Layer 2: Scatter circles (size = Count, color = -log10(pvalue)) -------------
geom_point(
   aes(x = GeneRatio_numeric, size = Count, color = neg_log10_p, fill = neg_log10_p),
   shape = 21,  # Filled circle shape, supports independent fill/color mapping
   stroke = 0.5 # Circle border thickness
) +
   # ------------- Publication-standard plot aesthetic adjustments -------------
# Axis configuration
scale_x_continuous(
   expand = c(0.02, 0),  # Remove blank margins on both sides of x-axis
   labels = percent_format(accuracy = 1)  # Format GeneRatio as percentage for readability
) +
   scale_size(range = c(3, 10)) +  # Range of circle sizes, adjustable based on Count distribution
   color_scale + fill_scale +     # Circle color gradient (light green → dark green)
   # Text & theme settings
   geom_text(
      aes(x = GeneRatio_numeric, label = Description), 
      hjust = -0.05,  # Place text labels to the right of bars
      fontface = "bold",  # Bold font
      size = 3.5,         # Font size
      color = "black"     # Black text color
   ) +
   labs(
      x = "Gene Ratio", 
      y = "",  # Hide y-axis label (pathway names are displayed as text annotations)
      size = "Count", 
      color = "-log10(P-value)",
      fill = "-log10(P-value)"
   ) +
   theme_bw() +  # Clean white background theme for manuscripts
   theme(
      # Axis text formatting
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_blank(),  # Hide default y-axis labels to avoid duplication
      axis.ticks.y = element_blank(), # Remove y-axis tick marks
      # Legend & title formatting
      legend.title = element_text(size = 10, fontface = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 12, fontface = "bold"),
      # Simplified grid lines for professional appearance
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      # Plot margin reserve space on right for long pathway labels
      plot.margin = margin(10, 50, 10, 20)
   )


# ===================== 4. Export high-resolution publication-ready figure =====================
ggsave(
   "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/Top10_GO_Pathway.pdf",  # Save as vector PDF without resolution loss
   plot = p,
   width = 6,  # Plot width (inches)
   height = 4,  # Plot height (inches)
   dpi = 300,   # Render resolution
   device = "pdf"
)


### KEGG pathway visualization
KEGG <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/Feature_1155mRNA_KEGG.csv")
KEGG_Filter <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_KEGG_filtered.csv")
KEGGid <- KEGG_Filter$Unnamed..0

KEGG_matched <- KEGG[ KEGG$X %in% KEGGid ,]
# Preserve original column order if required
KEGG_matched <- KEGG_matched[ match(KEGGid, KEGG_matched$X) , ]
write.csv(KEGG_matched,file =  "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_KEGG_filtered_full_info.csv")


selected_paths <- read.csv("/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/MM_related_KEGG_filtered_full_info_Top10Pathway.csv")
# ===================== 2. Data preprocessing =====================
# 2.1 Parse GeneRatio into numeric value (numerator / denominator)
selected_paths <- selected_paths %>%
   mutate(
      GeneRatio_numeric = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,1]) / 
         as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,2]),
      # Calculate -log10(pvalue) for color mapping
      neg_log10_p = -log10(pvalue),
      # Sort pathways vertically by GeneRatio magnitude
      Description = factor(Description, levels = rev(unique(Description[order(GeneRatio_numeric)])))
   )

# ===================== 3. Composite plot construction =====================
# Blue gradient palette mapped to -log10(pvalue)
color_scale <- scale_color_gradient(low = "#BED9ED", high = "#0382E2")
fill_scale <- scale_fill_gradient(low = "#BED9ED", high = "#0382E2")

p <- ggplot(selected_paths, aes(y = Description)) +
   # ------------- Layer 1: Horizontal bar chart for GeneRatio -------------
geom_col(
   aes(x = GeneRatio_numeric), 
   fill = "#76ABE5",  # Fixed bar blue tone (adjustable)
   width = 0.8,       # Bar width
   alpha = 0.8        # Bar transparency
) +
   # ------------- Layer 2: Scatter circles (size = Count, color = -log10(pvalue)) -------------
geom_point(
   aes(x = GeneRatio_numeric, size = Count, color = neg_log10_p, fill = neg_log10_p),
   shape = 21,  # Filled circle shape with independent fill/color aesthetics
   stroke = 0.5 # Circle border width
) +
   # ------------- Publication-standard aesthetic formatting -------------
# Axis settings
scale_x_continuous(
   expand = c(0.02, 0),  # Remove empty margins on x-axis
   labels = percent_format(accuracy = 1)  # Display GeneRatio as percentage values
) +
   scale_size(range = c(3, 10)) +  # Adjustable circle size range based on Count values
   color_scale + fill_scale +     # Blue gradient for circle fill/outline
   # Text & theme formatting
   geom_text(
      aes(x = GeneRatio_numeric, label = Description), 
      hjust = -0.05,  # Place pathway labels to the right of bars
      fontface = "bold",  # Bold text
      size = 3.5,         # Font size
      color = "black"     # Black label color
   ) +
   labs(
      x = "Gene Ratio", 
      y = "",  # Hide redundant y-axis labels
      size = "Count", 
      color = "-log10(P-value)",
      fill = "-log10(P-value)"
   ) +
   theme_bw() +  # Standard white background theme for scientific manuscripts
   theme(
      # Axis text style
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_blank(),  # Remove duplicate y-axis tick labels
      axis.ticks.y = element_blank(), # Remove y-axis tick marks
      # Legend and title styling
      legend.title = element_text(size = 10, fontface = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 12, fontface = "bold"),
      # Simplified grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      # Reserve right margin space for long pathway labels
      plot.margin = margin(10, 50, 10, 20)
   )


# ===================== 4. Export high-resolution publication figure =====================
ggsave(
   "/home/yjliu/mmProj/data_process/Human/Pathogenesis_and_Drug_Therapy/Enrich_Analysis/Top10_KEGG_Pathway.pdf",  # Vector PDF output for lossless publication
   plot = p,
   width = 6,  # Plot width (inch)
   height = 4,  # Plot height (inch)
   dpi = 300,   # Rendering resolution
   device = "pdf"
)