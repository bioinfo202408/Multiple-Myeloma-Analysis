##Figure3 A
## ================================
## 0. Load packages
## ================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(concaveman)
library(grid)

## ================================
## 1. Load data and extract MM cohort
## ================================
load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF_0205/GSE232988/seu.Rdata")

seu_subset <- subset(seu, subset = cohort == "MM")
seu_subset$cohort    <- droplevels(seu_subset$cohort)
seu_subset$clonotype <- droplevels(seu_subset$clonotype)
seu_subset$batch     <- droplevels(seu_subset$batch)
seu_subset$cellType  <- droplevels(seu_subset$cellType)

## Check if reduction exists
if (!"X_umap" %in% Reductions(seu_subset)) {
   stop("Reduction 'X_umap' does not exist. Please check the dimension reduction names in the Seurat object: ", 
        paste(Reductions(seu_subset), collapse = ", "))
}

## Check if cellType metadata column exists
if (!"cellType" %in% colnames(seu_subset@meta.data)) {
   stop("Column 'cellType' does not exist in meta.data.")
}

## ================================
## 2. Map 32 cell subtypes to 10 major cell types
## ================================
major_map <- c(
   ## -------- T / NK --------
   "CD8 Cytotoxic T Cells"         = "T/NK Cells",
   "MAIT Cytotoxic T Cells"        = "T/NK Cells",
   "NKT Cells"                     = "T/NK Cells",
   "CD4 Memory T Cells"            = "T/NK Cells",
   "CD8 Memory T Cells"            = "T/NK Cells",
   "GD Cytotoxic T Cells"          = "T/NK Cells",
   "CD8 Activated T Cells"         = "T/NK Cells",
   "NK Cells"                      = "T/NK Cells",
   "CD56 Bright NK Cells"          = "T/NK Cells",
   "CD4 Naive T Cells"             = "T/NK Cells",
   "CD4 Activated Memory T Cells"  = "T/NK Cells",
   "Tregs"                         = "T/NK Cells",
   "CD4 Cytotoxic T Cells"         = "T/NK Cells",
   "CD8 Naive T Cells"             = "T/NK Cells",
   
   ## -------- Platelets --------
   "Platelets"                     = "Platelets",
   
   ## -------- Monocytic --------
   "Macrophages"                   = "Monocytic Cells",
   "CD14 Monocytes"                = "Monocytic Cells",
   "TAMs"                          = "Monocytic Cells",
   
   ## -------- Plasma --------
   "Plasma Cells"                  = "Plasma Cells",
   
   ## -------- pDC --------
   "pDCs"                          = "pDCs",
   
   ## -------- mDC --------
   "CD14 DCs"                      = "mDCs",
   "CLEC9A DCs"                    = "mDCs",
   "CD1C DCs"                      = "mDCs",
   
   ## -------- Neutrophils --------
   "Immature Neutrophils"          = "Neutrophils",
   
   ## -------- HSPCs --------
   "HSPCs"                         = "HSPCs",
   
   ## -------- B Cells --------
   "Memory B Cells"                = "B Cells",
   "Naive B Cells"                 = "B Cells",
   "PrePro B Cells"                = "B Cells",
   
   ## -------- Erythrocytes --------
   "Erythrocytes"                  = "Erythrocytes",
   "Erythroblasts"                 = "Erythrocytes",
   "Erythroid Progenitors"         = "Erythrocytes"
   
   ## "Inconsistent" entries are excluded from the 10 major groups and processed separately later
)

## Add new major cell type annotation
seu_subset$majorCellType <- unname(major_map[as.character(seu_subset$cellType)])

## Mark unmapped cell types as "Unassigned"
seu_subset$majorCellType[is.na(seu_subset$majorCellType)] <- "Unassigned"

## Print mapping results
cat("===== Cell type mapping table =====\n")
print(table(seu_subset$cellType, seu_subset$majorCellType))

cat("\n===== Cell counts per majorCellType =====\n")
print(table(seu_subset$majorCellType))

## ================================
## 3. Extract UMAP coordinates and prepare plotting data
## ================================
umap_df <- as.data.frame(Embeddings(seu_subset, reduction = "X_umap"))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

plot_df <- umap_df %>%
   tibble::rownames_to_column("cell") %>%
   mutate(
      cellType = seu_subset$cellType[match(cell, colnames(seu_subset))],
      majorCellType = seu_subset$majorCellType[match(cell, colnames(seu_subset))]
   )

## Main plot only displays the 10 major types; Unassigned cells are excluded from contour outlines by default
plot_df_main <- plot_df %>%
   filter(majorCellType != "Unassigned")

plot_df_bg <- plot_df %>%
   filter(majorCellType == "Unassigned")

## Set display order for major cell types
major_levels <- c(
   "T/NK Cells", "Platelets", "Monocytic Cells", "Plasma Cells", "pDCs",
   "mDCs", "Neutrophils", "HSPCs", "B Cells", "Erythrocytes"
)

plot_df_main$majorCellType <- factor(plot_df_main$majorCellType, levels = major_levels)

## ================================
## 4. Color palette setup (clean, publication-friendly)
## ================================
major_cols <- c(
   "T/NK Cells"       = "#D95F02",
   "Platelets"        = "#E7298A",
   "Monocytic Cells"  = "#8C510A",
   "Plasma Cells"     = "#5E72B5",
   "pDCs"             = "#D73027",
   "mDCs"             = "#7F6D1D",
   "Neutrophils"      = "#66A61E",
   "HSPCs"            = "#5CC8C2",
   "B Cells"          = "#C51B7D",
   "Erythrocytes"     = "#F1A6C8"
)

## ================================
## 5. Generate dashed concave hull contours for each major cell type
## ================================
get_hull <- function(df, group_name, concavity = 2) {
   df_sub <- df %>%
      dplyr::filter(majorCellType == group_name) %>%
      dplyr::select(UMAP1, UMAP2)
   
   ## Skip contour generation when cell count is too low
   if (nrow(df_sub) < 10) return(NULL)
   
   ## concaveman requires matrix input
   pts <- as.matrix(df_sub[, c("UMAP1", "UMAP2")])
   
   hull <- concaveman::concaveman(pts, concavity = concavity)
   
   ## Convert back to data.frame
   hull <- as.data.frame(hull)
   colnames(hull) <- c("UMAP1", "UMAP2")
   hull$majorCellType <- group_name
   
   return(hull)
}

hull_list <- lapply(major_levels, function(x) get_hull(plot_df_main, x, concavity = 2))
hull_list <- hull_list[!sapply(hull_list, is.null)]
hull_df <- do.call(rbind, hull_list)
## ================================
## 6. Calculate label positions
## Use median coordinates for robustness, combined with ggrepel to avoid label overlap
## ================================
label_df <- plot_df_main %>%
   group_by(majorCellType) %>%
   summarise(
      UMAP1 = median(UMAP1),
      UMAP2 = median(UMAP2),
      .groups = "drop"
   )

## ================================
## 7. Add UMAP axis arrows at bottom-left (high-impact journal style)
## ================================
xr <- range(plot_df$UMAP1, na.rm = TRUE)
yr <- range(plot_df$UMAP2, na.rm = TRUE)

x_span <- diff(xr)
y_span <- diff(yr)

arrow_x0 <- xr[1] + 0.04 * x_span
arrow_y0 <- yr[1] + 0.06 * y_span

arrow_x1 <- arrow_x0 + 0.12 * x_span
arrow_y1 <- arrow_y0 + 0.12 * y_span

## ================================
## 8. Plot generation
## ================================
p <- ggplot() +
   ## Background layer: unassigned cells, light gray
   geom_point(
      data = plot_df_bg,
      aes(x = UMAP1, y = UMAP2),
      color = "grey85",
      size = 0.15,
      alpha = 0.45
   ) +
   ## Main layer: 10 major cell populations
   geom_point(
      data = plot_df_main,
      aes(x = UMAP1, y = UMAP2, color = majorCellType),
      size = 0.22,
      alpha = 0.75,
      stroke = 0
   ) +
   ## Population labels
   geom_label_repel(
      data = label_df,
      aes(x = UMAP1, y = UMAP2, label = majorCellType),
      size = 4.8,
      fontface = "bold",
      family = "sans",
      color = "black",
      fill = alpha("white", 0.85),
      label.size = 0,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = NA,
      seed = 123
   ) +
   ## Custom color mapping
   scale_color_manual(values = major_cols, drop = FALSE) +
   ## Fixed coordinate ratio
   coord_fixed() +
   ## Classic white background theme
   theme_classic(base_size = 16, base_family = "sans") +
   theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin = margin(15, 20, 15, 20)
   ) +
   ## Bottom-left UMAP axis arrows
   annotate(
      "segment",
      x = arrow_x0, xend = arrow_x1,
      y = arrow_y0, yend = arrow_y0,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.8,
      color = "black"
   ) +
   annotate(
      "segment",
      x = arrow_x0, xend = arrow_x0,
      y = arrow_y0, yend = arrow_y1,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.8,
      color = "black"
   ) +
   annotate(
      "text",
      x = (arrow_x0 + arrow_x1) / 2,
      y = arrow_y0 - 0.05 * y_span,
      label = "UMAP1",
      size = 6,
      family = "sans"
   ) +
   annotate(
      "text",
      x = arrow_x0 - 0.04 * x_span,
      y = (arrow_y0 + arrow_y1) / 2,
      label = "UMAP2",
      size = 6,
      angle = 90,
      family = "sans"
   ) +
   ggtitle("Major Cell Types in MM")

print(p)

## ================================
## 9. Save figure (cairo_pdf recommended for publication)
## ================================
out_pdf <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/1A_MM_MajorCellType_X_umap_publication.pdf"

ggsave(
   filename = out_pdf,
   plot = p,
   width = 9,
   height = 8,
   units = "in",
   dpi = 300,
   bg = "white"
)

cat("\nPDF saved to path:\n", out_pdf, "\n")

###Figure3 B NMF
seu_q <- subset(seu_subset, subset = cellType == "Plasma Cells")
seu_q$cellType <- droplevels(seu_q$cellType)

# 3. Extract raw counts matrix
Plasma_matrix <- GetAssayData(seu_q, assay = "originalexp", slot = "counts")

# 4. Load feature gene list
FeatureGene <- read.table(
   "/home/yjliu/mmProj/data_process/Human/Feature_select/mRNA/tuning_groupkfold/mRNA_all_regions_pf0.6_20251124_192149/stable_features_freq_ge9_mRNA_lncRNA.txt",
   header = FALSE,
   stringsAsFactors = FALSE
)
FeatureGene <- FeatureGene$V1
FeatureGene <- sub("\\.[0-9]+$", "", FeatureGene)

# 5. Convert Ensembl IDs to Gene Symbols
gene_map <- mapIds(
   org.Hs.eg.db,
   keys = FeatureGene,
   column = "SYMBOL",
   keytype = "ENSEMBL",
   multiVals = "first"
)

# Remove NA entries and duplicate gene symbols
gene_map <- unique(na.omit(gene_map))

# 6. Subset single-cell expression matrix to selected feature genes
scRNA <- Plasma_matrix[rownames(Plasma_matrix) %in% gene_map, ]

# ----------------------------
# 7. Extract metadata
# ----------------------------
meta_q <- seu_q@meta.data

# Verify cell names as rownames of metadata
head(rownames(meta_q))
head(colnames(scRNA))

# 8. Retain cells shared by expression matrix and metadata
common_cells <- intersect(colnames(scRNA), rownames(meta_q))

# Subset and strictly align cell order
scRNA_sub <- scRNA[, common_cells]
meta_sub <- meta_q[common_cells, , drop = FALSE]

# Double-check row-column order consistency
stopifnot(identical(colnames(scRNA_sub), rownames(meta_sub)))

# ----------------------------
# 9. Construct new Seurat object
# ----------------------------
scRNA_seurat <- CreateSeuratObject(
   counts = scRNA_sub,
   meta.data = meta_sub,
   project = "MM_Plasma_FeatureGene"
)

# 10. Optional: Original dimension reductions cannot be directly inherited for downstream analysis
# Skip this block if not required for subsequent analysis

# 11. Inspect new Seurat object
scRNA_seurat
dim(scRNA_seurat)
head(scRNA_seurat@meta.data)

# ----------------------------
# 12. Save Seurat object
# ----------------------------
save(scRNA_seurat, file = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/GSE232988/seu_MM_Plasma_FeatureGene.Rdata")

load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/GSE232988/seu_MM_Plasma_FeatureGene.Rdata")

scRNA_seurat <- subset(
   scRNA_seurat,
   subset = !is.na(nCount_RNA) & nCount_RNA > 0 & !is.na(nFeature_RNA) & nFeature_RNA > 0
)

library(sciNMF)

RunNMF <- function(object, group.by, dir.output = NULL, k.range = 3:8, samples = NULL, project = "NMF",
                   normalization.method = "SCT", min.cell = 10, variable.features.n = 1000,
                   do.scale = FALSE, do.center = TRUE,
                   ncore = 1, seed = 123,
                   rm.MT = FALSE, rm.RP.S.L = FALSE, rm.HSP = FALSE,
                   loss = "mse", max.iter = 5000, method = "scd", ...) {
   
   # Check SeuratObject package version
   flag_v <- packageVersion("SeuratObject") >= "5"
   
   if (any(is.na(object@meta.data[, group.by]))) {
      warning("The ", group.by, " column contains NA values; corresponding cells will be removed!")
      idx_cell <- !is.na(object@meta.data[, group.by])
      object <- object[, idx_cell]
   }
   
   if (is.null(samples)) {
      samples <- unique(object@meta.data[, group.by])
   }
   
   if (flag_v) {
      genes <- rownames(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"))
   } else {
      genes <- rownames(Seurat::GetAssayData(object, assay = "RNA", slot = "counts"))
   }
   
   # Filter mitochondrial, ribosomal, HSP genes if enabled
   if (rm.MT) {
      genes <- grep("^MT-", genes, invert = TRUE, value = TRUE)
   }
   if (rm.RP.S.L) {
      genes <- grep("^RP[SL]", genes, invert = TRUE, value = TRUE)
   }
   if (rm.HSP) {
      genes <- grep("^HSP", genes, invert = TRUE, value = TRUE)
   }
   
   if (flag_v) {
      clean_counts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
   } else {
      clean_counts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")[genes, , drop = FALSE]
   }
   
   # Sequential execution (no parallelization)
   ls_res <- vector("list", length(samples))
   names(ls_res) <- as.character(samples)
   
   for (sam in samples) {
      message("Start sample ", sam, " -- Current time:", as.character(Sys.time()))
      idx_cell <- object@meta.data[, group.by] == sam
      
      if (sum(idx_cell) < min.cell) {
         message("Sample ", sam, " only contains ", sum(idx_cell),
                 " cells (less than ", min.cell, "), skip sample\n")
         ls_res[[as.character(sam)]] <- NULL
         next
      }
      
      idx_0_gene <- Matrix::rowSums(clean_counts[, idx_cell, drop = FALSE]) == 0
      
      srt <- Seurat::CreateSeuratObject(
         counts = clean_counts[!idx_0_gene, idx_cell, drop = FALSE],
         meta.data = object@meta.data[idx_cell, , drop = FALSE]
      )
      
      if (normalization.method == "SCT") {
         srt <- Seurat::SCTransform(
            srt,
            verbose = FALSE,
            do.scale = do.scale,
            do.center = do.center,
            variable.features.n = variable.features.n
         )
         if (flag_v) {
            data <- Seurat::GetAssayData(srt, assay = "SCT", layer = "scale.data")
         } else {
            data <- Seurat::GetAssayData(srt, assay = "SCT", slot = "scale.data")
         }
      } else if (normalization.method == "LogNormalize") {
         srt <- Seurat::NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
         srt <- Seurat::FindVariableFeatures(srt, selection.method = "vst", nfeatures = variable.features.n)
         srt <- Seurat::ScaleData(
            srt,
            features = Seurat::VariableFeatures(srt),
            do.scale = do.scale,
            do.center = do.center
         )
         if (flag_v) {
            data <- Seurat::GetAssayData(srt, assay = "RNA", layer = "scale.data")
         } else {
            data <- Seurat::GetAssayData(srt, assay = "RNA", slot = "scale.data")
         }
      } else {
         stop("Invalid normalization.method; must be either SCT or LogNormalize")
      }
      
      # Critical step: Prevent dimension reduction to vector after subsetting
      vf <- Seurat::VariableFeatures(srt)
      data <- data[vf, , drop = FALSE]
      data[data < 0] <- 0
      data <- data[apply(data, 1, var) > 0, , drop = FALSE]
      
      # Safety filter: Skip sample if usable features/cells are insufficient to avoid invalid NMF results
      if (nrow(data) < 2 || ncol(data) < 2) {
         message("Sample ", sam, " has insufficient usable features/cells post-filtering, skip sample\n")
         ls_res[[as.character(sam)]] <- NULL
         next
      }
      
      ls_WH <- lapply(k.range, function(k) {
         set.seed(seed)
         res_nmf <- NNLM::nnmf(data, k = k, loss = loss, max.iter = max.iter, method = method, ...)
         H <- res_nmf$H
         W <- res_nmf$W
         rownames(H) <- colnames(W) <- paste0(project, "_", sam, "_K", k, "_P", 1:k)
         list(H = H, W = W)
      })
      
      all_W <- lapply(ls_WH, `[[`, "W") |> do.call(what = cbind, args = _)
      all_H <- lapply(ls_WH, `[[`, "H") |> do.call(what = rbind, args = _)
      WHs <- list(W = all_W, H = all_H)
      
      if (!is.null(dir.output)) {
         if (!file.exists(dir.output)) {
            dir.create(dir.output, recursive = TRUE)
         }
         saveRDS(
            WHs,
            paste0(dir.output, "/", project, "_", sam, "_hvg", variable.features.n,
                   "_k", k.range[1], "to", tail(k.range, 1), ".rds")
         )
      }
      
      message("Sample ", sam, " completed!")
      ls_res[[as.character(sam)]] <- WHs
   }
   
   ls_res <- ls_res[!sapply(ls_res, is.null)]
   message("All samples processed!")
   return(ls_res)
}


ls_WH <- RunNMF(scRNA_seurat,
                group.by="batch",
                dir.output="/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/rds800/",
                project="MM",
                k.range=3:8,
                variable.features.n=800,
                min.cell=10,
                normalization.method="LogNormalize",
                seed=777,
                ncore=1
)


ls_ph <- list.files('/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/rds800/', full.names = TRUE)
head(ls_ph,3)
ls_WH <- lapply(ls_ph, readRDS)


pdf("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result//IQRPlot.pdf", width = 15, height = 10)
p <- IQRPlot(ls_WH[1:4], IQR.cut = 0.1, median.cut = 0.02, grid = TRUE)
print(p)
dev.off()

RobustProgram <- function(WH.list, top = 50, IQR.cut = 0.1, median.cut = 0, intra.min = 35, intra.rep = 1,
                          inter.filter = TRUE, inter.min = 10, inter.rep = 1, intra.max = 10) {
   WH.list <- WH.list[!sapply(WH.list, is.null)]
   # QC filtering based on IQR and median gene usage
   ls_pat_pg <- lapply(WH.list, function(WH) {
      # Normalize H matrix by column to calculate IQR for each factorization rank k
      Ks <- sub(".*_K([0-9]+)_P[0-9]+$", "\\1", rownames(WH$H))
      H <- split(data.frame(WH$H, check.names = FALSE), Ks) %>%
         lapply(function(sub_H) {
            H_ratio <- apply(sub_H, 2, function(me) {
               s <- sum(me, na.rm = TRUE)
               if (is.na(s) || s == 0) {
                  return(rep(NA_real_, length(me)))  # Alternatively rep(0, length(me)) per your requirement
               }
               me / s
            })
            return(H_ratio)
         }) %>%
         {
            names(.) <- NULL
            do.call(what = rbind, args = .)
         }
      
      mat_quat <- apply(H, 1, function(x) {
         if (all(is.na(x) | is.nan(x))) {
            return(c(`0%`=NA, `25%`=NA, `50%`=NA, `75%`=NA, `100%`=NA))
         }
         quantile(x, na.rm = TRUE)
      })
      idx_median <- mat_quat["50%", ] >= median.cut
      idx_IQR <- mat_quat["75%", ] - mat_quat["25%", ] >= IQR.cut
      
      W_filter <- WH$W[, idx_median & idx_IQR]
      
      pgs <- lapply(setNames(nm = colnames(W_filter)), function(pg) {
         pg <- W_filter[, pg]
         head(sort(pg, decreasing = TRUE), n = top)
      })
      # Intra-sample filtering via overlap matrix
      mat_ovlp <- OverlapMat(lapply(pgs, names))
      # mat_ovlp includes self-overlap, threshold > intra.rep
      idx_pg <- apply(mat_ovlp, 1, function(pgop) {
         sum(pgop >= intra.min)
      }) > intra.rep
      if (sum(idx_pg) == 0) {
         # Skip sample if no programs pass intra-sample overlap threshold
         return(NULL)
      }
      pgs <- pgs[idx_pg]
   })
   
   ls_pat_pg <- ls_pat_pg[!sapply(ls_pat_pg, is.null)]
   
   # Inter-sample filtering
   mat_ovlp_all <- OverlapMat(lapply(unlist(ls_pat_pg, recursive = FALSE), names))
   ls_pat_names <- lapply(ls_pat_pg, names)
   colnames(mat_ovlp_all) <- rownames(mat_ovlp_all) <- unlist(ls_pat_names, use.names = FALSE)
   
   ls_keep_pg <- lapply(1:length(ls_pat_pg), function(pat) {
      idx_pgs <- match(names(ls_pat_pg[[pat]]), colnames(mat_ovlp_all))
      sub_mat_ovlp <- mat_ovlp_all[idx_pgs, -idx_pgs, drop = FALSE]
      # sub_mat_ovlp excludes self-overlap, threshold >= inter.rep
      idx_inter <- apply(sub_mat_ovlp, 1, function(pgop) {
         sum(pgop >= inter.min) >= inter.rep
      })
      
      # Remove programs failing inter-sample overlap threshold
      pgs <- names(idx_inter)[idx_inter]
      # Skip sample if no programs pass inter-sample overlap threshold
      if (length(pgs) == 0) {
         return(NULL)
      }
      pgs <- sort(apply(sub_mat_ovlp[pgs, , drop = FALSE], 1, max), decreasing = TRUE)
      
      if (length(pgs) > 1) {
         keep_pgs <- names(pgs)[1]
         for (pg_test in names(pgs)[-1]) {
            if (max(mat_ovlp_all[keep_pgs, pg_test]) <= intra.max) {
               keep_pgs <- c(keep_pgs, pg_test)
            }
         }
      } else {
         keep_pgs <- names(pgs)
      }
      
      return(lapply(ls_pat_pg[[pat]][keep_pgs], names))
   }) %>% unlist(recursive = FALSE)
   
   return(ls_keep_pg[!sapply(ls_keep_pg, is.null)])
}

ls_RP <- RobustProgram(WH.list = ls_WH, IQR.cut = 0.1, median.cut = 0.02)
length(ls_RP)

ovlp <- OverlapMat(ls_RP)

pdf("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/ClusterMetricsPlot_gene50.pdf", width = 8, height = 6)
p1 <- ClusterMetricsPlot(mat.ovlp = ovlp, 
                         distance.clustering = 'Intersection', 
                         max.intersect = 50, 
                         method.clustering = "ward.D2")
print(p1)
dev.off()

res_cluster <- ClusterPG(mat.ovlp = ovlp, cut.num = 4)
table(res_cluster)
head(res_cluster,3)


Sample <- strsplit(names(res_cluster), '_') %>% 
   sapply(function(s){
      paste0(s[1:2],collapse = '_')
   })

df_anno <- data.frame(Sample = Sample)
rownames(df_anno) = names(res_cluster)
all_sam <- unique(df_anno$Sample)
col_sam <- setNames(as.character(paletteer::paletteer_d("ggsci::default_igv"))[1:length(all_sam)], all_sam)
ls_col_anno <- list(Sample = col_sam)

head(df_anno,3)
lapply(ls_col_anno, head, 3)

ls_res <- MetaProgram(WH.list = ls_WH, cluster.result = res_cluster,
                      color.mp = rev(as.character(paletteer::paletteer_d(`"RColorBrewer::Set1"`))),
                      show.rownames = FALSE,
                      key = "MM",
                      min.size.MP = 31, keep.rep.gene = TRUE,
                      annotation = df_anno, color.annotation = ls_col_anno) # Ignore annotation parameters if no extra metadata available

length(ls_res$MetaProgram)
options(repr.plot.width = 12, repr.plot.height = 11.5)

pdf('/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/HeatMap_gene50_4.pdf', width = 12, height = 11.5)
print(ls_res$HeatMap)
dev.off()

ls_MP <- ls_res$MetaProgram
sapply(ls_MP,head,50)


library(Seurat)
library(ggplot2)
library(patchwork)

DefaultAssay(scRNA_seurat) <- "RNA"

# 1. Extract three core meta-programs
if (all(c("MM1", "MM2", "MM3") %in% names(ls_MP))) {
   ls_MP_use <- ls_MP[c("MM1", "MM2", "MM3")]
} else {
   ls_MP_use <- ls_MP[1:3]
   names(ls_MP_use) <- c("MM1", "MM2", "MM3")
}

# 2. Retain only genes present in the Seurat object
ls_MP_use <- lapply(ls_MP_use, function(x) {
   intersect(x, rownames(scRNA_seurat))
})

print(sapply(ls_MP_use, length))

# 3. Run normalization if RNA assay lacks data layer
if (!"data" %in% Layers(scRNA_seurat[["RNA"]])) {
   scRNA_seurat <- NormalizeData(
      scRNA_seurat,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
   )
}

# 4. Calculate module scores for each meta-program using AddModuleScore
scRNA_seurat <- AddModuleScore(
   object = scRNA_seurat,
   features = ls_MP_use,
   assay = "RNA",
   name = "MM_Score",
   ctrl = 25,
   seed = 123
)

# AddModuleScore auto-generates MM_Score1, MM_Score2, MM_Score3
old_score_cols <- paste0("MM_Score", 1:3)
new_score_cols <- paste0(names(ls_MP_use), "_Score")

colnames(scRNA_seurat@meta.data)[
   match(old_score_cols, colnames(scRNA_seurat@meta.data))
] <- new_score_cols

score_cols <- new_score_cols
print(score_cols)

# 5. Assign dominant program label per cell based on maximum module score
score_mat <- scRNA_seurat@meta.data[, score_cols, drop = FALSE]

dominant_program <- apply(score_mat, 1, function(x) {
   names(x)[which.max(x)]
})

dominant_program <- sub("_Score$", "", dominant_program)

scRNA_seurat$dominant_program <- factor(
   dominant_program,
   levels = c("MM1", "MM2", "MM3")
)

print(table(scRNA_seurat$dominant_program))

# 6. Recompute PCA / UMAP / t-SNE using only meta-program marker genes
program_genes <- unique(unlist(ls_MP_use))
program_genes <- intersect(program_genes, rownames(scRNA_seurat))

scRNA_seurat <- ScaleData(
   scRNA_seurat,
   features = program_genes,
   verbose = FALSE
)

npcs_use <- min(30, length(program_genes) - 1)

scRNA_seurat <- RunPCA(
   scRNA_seurat,
   features = program_genes,
   npcs = npcs_use,
   reduction.name = "pca_program_gene",
   reduction.key = "PCPG_",
   verbose = FALSE
)

dims_umap <- 1:min(20, npcs_use)
dims_tsne <- 1:min(10, npcs_use)

scRNA_seurat <- RunUMAP(
   scRNA_seurat,
   reduction = "pca_program_gene",
   dims = dims_umap,
   reduction.name = "umap_program_gene",
   reduction.key = "UMAPPG_",
   n.neighbors = 30,
   min.dist = 0.1,
   seed.use = 123
)

scRNA_seurat <- RunTSNE(
   scRNA_seurat,
   reduction = "pca_program_gene",
   dims = dims_tsne,
   reduction.name = "tsne_program_gene",
   reduction.key = "TSNEPG_",
   seed.use = 123,
   check_duplicates = FALSE
)

# 7. Color palette for meta-programs
program_colors <- c(
   "MM1" = "#8EBD80",
   "MM2" = "#00A1D4",
   "MM3" = "#D73027"
)

# 8. UMAP FeaturePlot showing module score distributions for three programs
p_umap_score <- FeaturePlot(
   scRNA_seurat,
   features = score_cols,
   reduction = "umap_program_gene",
   ncol = 3,
   cols = c("#D7E4EF", "#7AA1BC", "#00A1D4"),
   pt.size = 0.4
) &
   theme_classic(base_size = 13) &
   theme(
      axis.line = element_line(linewidth = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
   )

ggsave(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/AddModuleScore_ProgramGene_UMAP_FeaturePlot.pdf",
   p_umap_score,
   width = 17,
   height = 5
)

# 9. t-SNE FeaturePlot showing module score distributions for three programs
p_tsne_score <- FeaturePlot(
   scRNA_seurat,
   features = score_cols,
   reduction = "tsne_program_gene",
   ncol = 3,
   cols = c("#D7E4EF", "#7AA1BC", "#00A1D4"),
   pt.size = 0.5
) &
   theme_classic(base_size = 13) &
   theme(
      axis.line = element_line(linewidth = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
   )

ggsave(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/AddModuleScore_ProgramGene_tSNE_FeaturePlot.pdf",
   p_tsne_score,
   width = 15,
   height = 5
)

# 10. UMAP DimPlot colored by dominant meta-program per cell
p_umap_dom <- DimPlot(
   scRNA_seurat,
   reduction = "umap_program_gene",
   group.by = "dominant_program",
   cols = program_colors,
   pt.size = 0.6,
   alpha = 0.8
) +
   theme_classic(base_size = 14) +
   labs(
      x = "Program-gene UMAP 1",
      y = "Program-gene UMAP 2",
      color = "Program"
   ) +
   theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.grid = element_blank()
   )

ggsave(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/AddModuleScore_DominantProgram_UMAP.pdf",
   p_umap_dom,
   width = 6,
   height = 5
)

# 11. t-SNE DimPlot colored by dominant meta-program per cell
p_tsne_dom <- DimPlot(
   scRNA_seurat,
   reduction = "tsne_program_gene",
   group.by = "dominant_program",
   cols = program_colors,
   pt.size = 0.6,
   alpha = 0.8
) +
   theme_classic(base_size = 14) +
   labs(
      x = "Program-gene t-SNE 1",
      y = "Program-gene t-SNE 2",
      color = "Program"
   ) +
   theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.grid = element_blank()
   )

ggsave(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/AddModuleScore_DominantProgram_tSNE.pdf",
   p_tsne_dom,
   width = 6,
   height = 5
)


## Visualize differences in cell states using marker gene panels
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

DefaultAssay(scRNA_seurat) <- "RNA"

# ============================================================
# 1. Define marker gene panels for distinct cell states
# ============================================================

MM3_proliferation <- c(
   "MAD2L1", "ASPM", "CDCA5", "ECT2", "RAD54L", "ORC1",
   "E2F7", "STIL", "GINS3", "CKAP2L", "FANCD2", "DTL"
)

MM2_adhesion <- c(
   "COL6A3", "COL4A5", "NCAM1", "HGF", "PCDH9",
   "RELN", "TJP1", "FMN1", "ADGRB3"
)

MM1_metabolism <- c(
   "PHGDH", "PSAT1", "SLC2A5", "BDH1", "CHCHD4"
)

MM1_identity <- c(
   "TNFRSF13B", "POU2F2", "AIM2", "DKK1"
)

marker_panels <- list(
   Proliferation = MM3_proliferation,
   Adhesion_ECM = MM2_adhesion,
   Metabolism = MM1_metabolism,
   Plasma_Identity = MM1_identity
)

# ============================================================
# 2. Filter marker lists to retain only genes present in the Seurat object
# ============================================================

marker_panels_use <- lapply(marker_panels, function(g) {
   intersect(g, rownames(scRNA_seurat))
})

marker_panels_use <- marker_panels_use[sapply(marker_panels_use, length) > 0]

print(sapply(marker_panels_use, length))

# Warning if any signature contains fewer than three genes
if (any(sapply(marker_panels_use, length) < 3)) {
   warning("Some marker panels contain fewer than 3 detected genes in the Seurat object. Please verify gene symbols.")
}

# ============================================================
# 3. Run normalization if RNA assay lacks data layer
# ============================================================

if (!"data" %in% Layers(scRNA_seurat[["RNA"]])) {
   scRNA_seurat <- NormalizeData(
      object = scRNA_seurat,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
   )
}

scRNA_seurat <- AddModuleScore(
   object = scRNA_seurat,
   features = marker_panels_use,
   assay = "RNA",
   name = "MarkerSig",
   ctrl = 25,
   nbin = 12,
   seed = 123
)

old_marker_score_cols <- paste0("MarkerSig", seq_along(marker_panels_use))
new_marker_score_cols <- paste0(names(marker_panels_use), "_AddModuleScore")

colnames(scRNA_seurat@meta.data)[
   match(old_marker_score_cols, colnames(scRNA_seurat@meta.data))
] <- new_marker_score_cols

marker_score_cols <- new_marker_score_cols

print(marker_score_cols)

# ============================================================
# 5. Set order and color for dominant_program
# ============================================================

scRNA_seurat$dominant_program <- factor(
   scRNA_seurat$dominant_program,
   levels = c("MM1", "MM2", "MM3")
)

program_colors <- c(
   "MM1" = "#8EBD80",
   "MM2" = "#00A1D4",
   "MM3" = "#D73027"
)

group_var <- "dominant_program"

# ============================================================
# 6. Boxplot: Signature scores grouped by different dominant programs
# ============================================================

df <- scRNA_seurat@meta.data

p_list <- lapply(marker_score_cols, function(feat) {
   
   ggplot(df, aes_string(x = group_var, y = feat, fill = group_var)) +
      geom_boxplot(
         outlier.size = 0.2,
         linewidth = 0.4,
         alpha = 0.85
      ) +
      scale_fill_manual(values = program_colors) +
      ggtitle(feat) +
      theme_classic(base_size = 13) +
      theme(
         plot.title = element_text(hjust = 0.5, face = "bold"),
         axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
         axis.text.y = element_text(colour = "black"),
         axis.title.x = element_blank(),
         axis.title.y = element_text(face = "bold"),
         legend.position = "none",
         axis.line = element_line(linewidth = 0.5)
      ) +
      ylab("AddModuleScore")
})

p_box <- patchwork::wrap_plots(p_list, ncol = 1)

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/MarkerSignature_AddModuleScore_Boxplot_by_state.pdf",
   plot = p_box,
   width = 12,
   height = 10
)

# ============================================================
# 7. UMAP FeaturePlot: Spatial distribution of AddModuleScore
# ============================================================

p_umap <- FeaturePlot(
   object = scRNA_seurat,
   features = marker_score_cols,
   reduction = "umap_program_gene",
   ncol = length(marker_score_cols),
   cols = c("#D7E4EF", "#7AA1BC", "#00A1D4"),
   min.cutoff = "q05",
   max.cutoff = "q95",
   pt.size = 0.4,
   order = TRUE
) &
   coord_fixed() &
   theme_classic(base_size = 13) &
   theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
   )

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/MarkerSignature_AddModuleScore_UMAP.pdf",
   plot = p_umap,
   width = 16
)

# ============================================================
# 8. t-SNE FeaturePlot: Spatial distribution of AddModuleScore
# ============================================================

p_tsne <- FeaturePlot(
   object = scRNA_seurat,
   features = marker_score_cols,
   reduction = "tsne_program_gene",
   ncol = length(marker_score_cols),
   cols = c("#F7FBFF", "#7AA1BC", "#137CB7"),
   min.cutoff = "q05",
   max.cutoff = "q95",
   pt.size = 0.4,
   order = TRUE
) &
   coord_fixed() &
   theme_classic(base_size = 13) &
   theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
   )

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/MarkerSignature_AddModuleScore_tSNE.pdf",
   plot = p_tsne,
   width = 16,
   height = 4
)

# ============================================================
# 9. Signature score DotPlot: Module activity across different dominant programs
# ============================================================

df_score <- scRNA_seurat@meta.data[, c(group_var, marker_score_cols), drop = FALSE]

df_long <- df_score %>%
   tidyr::pivot_longer(
      cols = all_of(marker_score_cols),
      names_to = "Signature",
      values_to = "Score"
   )

# Use median value of each signature as high-score threshold
score_cutoff <- df_long %>%
   group_by(Signature) %>%
   summarise(
      cutoff = median(Score, na.rm = TRUE),
      .groups = "drop"
   )

df_long <- df_long %>%
   left_join(score_cutoff, by = "Signature") %>%
   mutate(is_high = Score > cutoff)

df_dot <- df_long %>%
   group_by(.data[[group_var]], Signature) %>%
   summarise(
      mean_score = mean(Score, na.rm = TRUE),
      pct_high = mean(is_high, na.rm = TRUE) * 100,
      .groups = "drop"
   )

df_dot[[group_var]] <- factor(
   df_dot[[group_var]],
   levels = c("MM1", "MM2", "MM3")
)

df_dot$Signature <- factor(
   df_dot$Signature,
   levels = marker_score_cols
)

p_dot_score <- ggplot(df_dot, aes_string(x = group_var, y = "Signature")) +
   geom_point(aes(size = pct_high, color = mean_score)) +
   scale_size(
      name = "% high-score cells",
      range = c(1.5, 8)
   ) +
   scale_color_gradient(
      low = "#C6DBEF",
      high = "#2171B5",
      name = "Mean AddModuleScore"
   ) +
   theme_classic(base_size = 13) +
   theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
      axis.text.y = element_text(colour = "black", face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
   ) +
   ggtitle("Marker signature AddModuleScores across dominant programs") +
   xlab("Dominant program") +
   ylab("Marker signature")

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/MarkerSignature_AddModuleScore_DotPlot_by_state.pdf",
   plot = p_dot_score,
   width = 8,
   height = 5
)


DefaultAssay(scRNA_seurat) <- "RNA"
# Group variable
group_var <- "dominant_program"

# Set group display order
scRNA_seurat$dominant_program <- factor(
   scRNA_seurat$dominant_program,
   levels = c("MM1", "MM2", "MM3")
)

# Ensure normalized data exists in RNA assay
if (!"data" %in% Layers(scRNA_seurat[["RNA"]])) {
   scRNA_seurat <- NormalizeData(
      scRNA_seurat,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
   )
}

# Keep only marker genes present in Seurat object
marker_panels_use <- lapply(marker_panels, function(g) {
   intersect(g, rownames(scRNA_seurat))
})

marker_panels_use <- marker_panels_use[sapply(marker_panels_use, length) > 0]

# Print retained gene count per marker panel
print(sapply(marker_panels_use, length))

# Generate DotPlot
p_dot_gene <- DotPlot(
   object = scRNA_seurat,
   features = marker_panels_use,
   group.by = group_var,
   assay = "RNA",
   dot.scale = 6
) +
   scale_colour_gradientn(
      colours = c("#C6DBEF", "#6BAED6", "#2171B5"),
      name = "Average expression"
   ) +
   theme_classic(base_size = 13) +
   theme(
      axis.text.x = element_text(
         angle = 45,
         hjust = 1,
         vjust = 1,
         size = 10,
         face = "bold",
         colour = "black"
      ),
      axis.text.y = element_text(
         size = 10,
         colour = "black"
      ),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.border = element_rect(
         colour = "black",
         fill = NA,
         linewidth = 0.6
      )
   ) +
   ggtitle("Marker genes across dominant programs") +
   xlab("Dominant program") +
   ylab("Marker genes")

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/MarkerGene_DotPlot_by_state_blue_publication.pdf",
   plot = p_dot_gene,
   width = 12,
   height = 8
)



##TF regulation analysis

library(Seurat)
library(dplyr)
library(tidyr)
library(viper)
library(decoupleR)
library(dorothea)
library(pheatmap)
library(ggplot2)
library(patchwork)
library(Matrix)

DefaultAssay(scRNA_seurat) <- "RNA"

# ============================================================
# 1. Ensure normalized data exists in RNA assay
# ============================================================

if (!"data" %in% Layers(scRNA_seurat[["RNA"]])) {
   scRNA_seurat <- NormalizeData(
      object = scRNA_seurat,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
   )
}

# ============================================================
# 2. Calculate AddModuleScore using three program gene sets from ls_MP
# ============================================================

if (all(c("MM1", "MM2", "MM3") %in% names(ls_MP))) {
   ls_MP_use <- ls_MP[c("MM1", "MM2", "MM3")]
} else {
   ls_MP_use <- ls_MP[1:3]
   names(ls_MP_use) <- c("MM1", "MM2", "MM3")
}

ls_MP_use <- lapply(ls_MP_use, function(g) {
   intersect(g, rownames(scRNA_seurat))
})

ls_MP_use <- ls_MP_use[sapply(ls_MP_use, length) > 0]

print(sapply(ls_MP_use, length))

scRNA_seurat <- AddModuleScore(
   object = scRNA_seurat,
   features = ls_MP_use,
   assay = "RNA",
   name = "MM_AddModule",
   ctrl = 25,
   nbin = 12,
   seed = 123
)

old_score_cols <- paste0("MM_AddModule", seq_along(ls_MP_use))
new_score_cols <- paste0(names(ls_MP_use), "_Score")

colnames(scRNA_seurat@meta.data)[
   match(old_score_cols, colnames(scRNA_seurat@meta.data))
] <- new_score_cols

score_cols <- new_score_cols
print(score_cols)

# ============================================================
# 3. Redefine dominant_program based on AddModuleScore
# ============================================================

score_mat <- scRNA_seurat@meta.data[, score_cols, drop = FALSE]

dominant_program <- apply(score_mat, 1, function(x) {
   names(x)[which.max(x)]
})

dominant_program <- sub("_Score$", "", dominant_program)

scRNA_seurat$dominant_program <- factor(
   dominant_program,
   levels = c("MM1", "MM2", "MM3")
)

print(table(scRNA_seurat$dominant_program))

# ============================================================
# 4. Extract expression matrix for VIPER calculation
# ============================================================

expr_mat <- GetAssayData(
   object = scRNA_seurat,
   assay = "RNA",
   layer = "data"
)

# decoupleR/viper requires standard matrix format
expr_mat <- as.matrix(expr_mat)

cat("Expression matrix dimension:\n")
print(dim(expr_mat))

# ============================================================
# 5. Load DoRothEA regulon database
# ============================================================

data("dorothea_hs", package = "dorothea")

net <- dorothea_hs %>%
   dplyr::filter(confidence %in% c("A", "B", "C")) %>%
   dplyr::transmute(
      source = tf,
      target = target,
      mor = mor,
      confidence = confidence
   )

# Retain only target genes present in expression matrix
net_use <- net %>%
   dplyr::filter(target %in% rownames(expr_mat))

cat("Number of TFs in regulon:", length(unique(net_use$source)), "\n")
cat("Number of TF-target edges:", nrow(net_use), "\n")

# ============================================================
# 6. Calculate TF activity per cell with VIPER
# ============================================================

tf_act_long <- decoupleR::run_viper(
   mat = expr_mat,
   network = net_use,
   .source = "source",
   .target = "target",
   .mor = "mor",
   minsize = 5,
   verbose = FALSE
)

head(tf_act_long)

# ============================================================
# 7. Reshape to TF-by-cell matrix
# ============================================================

tf_act_mat <- tf_act_long %>%
   dplyr::select(source, condition, score) %>%
   tidyr::pivot_wider(
      names_from = condition,
      values_from = score
   ) %>%
   as.data.frame()

rownames(tf_act_mat) <- tf_act_mat$source
tf_act_mat$source <- NULL
tf_act_mat <- as.matrix(tf_act_mat)

# Ensure column order matches cells in Seurat object
tf_act_mat <- tf_act_mat[, colnames(scRNA_seurat), drop = FALSE]

cat("TF activity matrix dimension:\n")
print(dim(tf_act_mat))

# ============================================================
# 8. Create new TF assay in Seurat object
# ============================================================

scRNA_seurat[["TF"]] <- CreateAssayObject(data = tf_act_mat)

DefaultAssay(scRNA_seurat) <- "RNA"

# ============================================================
# 9. Differential analysis of TF activity grouped by dominant_program
# ============================================================

group_var <- "dominant_program"
df_meta <- scRNA_seurat@meta.data

stopifnot(group_var %in% colnames(df_meta))

tf_names <- rownames(tf_act_mat)

tf_kw <- lapply(tf_names, function(tf) {
   
   sub_df <- data.frame(
      tf_activity = as.numeric(tf_act_mat[tf, colnames(scRNA_seurat)]),
      group = df_meta[[group_var]]
   )
   
   sub_df <- sub_df[complete.cases(sub_df), , drop = FALSE]
   
   if (length(unique(sub_df$group)) < 2) {
      return(NULL)
   }
   
   kt <- kruskal.test(tf_activity ~ group, data = sub_df)
   
   data.frame(
      TF = tf,
      statistic = unname(kt$statistic),
      p_value = kt$p.value
   )
})

tf_kw <- do.call(rbind, tf_kw)

tf_kw$padj <- p.adjust(tf_kw$p_value, method = "BH")
tf_kw <- tf_kw[order(tf_kw$padj, -tf_kw$statistic), ]

head(tf_kw, 20)

write.csv(
   tf_kw,
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/TFActivity_Kruskal_by_AddModuleScore_dominant_program.csv",
   row.names = FALSE
)

# ============================================================
# 10. Compute average TF activity per dominant program
# ============================================================

state_levels <- c("MM1", "MM2", "MM3")

tf_mean_by_state <- sapply(state_levels, function(st) {
   cells_use <- rownames(df_meta)[df_meta[[group_var]] == st]
   rowMeans(tf_act_mat[, cells_use, drop = FALSE], na.rm = TRUE)
})

tf_mean_by_state <- as.matrix(tf_mean_by_state)
colnames(tf_mean_by_state) <- state_levels

# ============================================================
# 11. Heatmap for top differentially active TFs
# ============================================================

top_n <- 60
top_tfs <- tf_kw$TF[1:min(top_n, nrow(tf_kw))]

pdf(
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/TFActivity_Heatmap_top60_by_AddModuleScore_state1.pdf",
   width = 8,
   height = 10
)

pheatmap(
   tf_mean_by_state[top_tfs, , drop = FALSE],
   scale = "row",
   cluster_rows = TRUE,
   cluster_cols = FALSE,
   color = colorRampPalette(c("white", "#C6DBEF", "#2171B5"))(100),
   main = "Top differential TF activities by dominant program"
)

dev.off()

# ============================================================
# 12. Extract top activated TFs for each program state
# ============================================================

top_tf_each_state <- lapply(state_levels, function(st) {
   other_states <- setdiff(state_levels, st)
   
   logFC_like <- tf_mean_by_state[, st] - rowMeans(
      tf_mean_by_state[, other_states, drop = FALSE],
      na.rm = TRUE
   )
   
   data.frame(
      TF = names(logFC_like),
      state = st,
      activity_difference = logFC_like
   ) %>%
      arrange(desc(activity_difference)) %>%
      head(20)
})

top_tf_each_state <- do.call(rbind, top_tf_each_state)

write.csv(
   top_tf_each_state,
   "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/Top20_TF_each_AddModuleScore_state.csv",
   row.names = FALSE
)

print(top_tf_each_state)

# ============================================================
# 13. DotPlot for top TF activities grouped by program state
# ============================================================

top_tfs_dot <- unique(top_tf_each_state$TF)

df_tf_long <- as.data.frame(t(tf_act_mat[top_tfs_dot, , drop = FALSE]))
df_tf_long$cell <- rownames(df_tf_long)
df_tf_long[[group_var]] <- df_meta[rownames(df_tf_long), group_var]

df_tf_long <- df_tf_long %>%
   tidyr::pivot_longer(
      cols = all_of(top_tfs_dot),
      names_to = "TF",
      values_to = "Activity"
   )

df_tf_dot <- df_tf_long %>%
   group_by(.data[[group_var]], TF) %>%
   summarise(
      mean_activity = mean(Activity, na.rm = TRUE),
      pct_high = mean(Activity > median(Activity, na.rm = TRUE), na.rm = TRUE) * 100,
      .groups = "drop"
   )

df_tf_dot[[group_var]] <- factor(
   df_tf_dot[[group_var]],
   levels = state_levels
)

df_tf_dot$TF <- factor(
   df_tf_dot$TF,
   levels = rev(unique(top_tfs_dot))
)

p_tf_dot <- ggplot(df_tf_dot, aes_string(x = group_var, y = "TF")) +
   geom_point(aes(size = pct_high, color = mean_activity)) +
   scale_size(
      name = "% high activity cells",
      range = c(1.5, 7)
   ) +
   scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Mean TF activity"
   ) +
   theme_classic(base_size = 13) +
   theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
      axis.text.y = element_text(colour = "black", size = 8),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
   ) +
   xlab("Dominant program") +
   ylab("Transcription factors") +
   ggtitle("Top TF activities across dominant programs")

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/TopTF_Activity_DotPlot_by_AddModuleScore_state.pdf",
   plot = p_tf_dot,
   width = 5,
   height = 13
)

# ============================================================
# 14. FeaturePlot to visualize representative TF activity
# ============================================================

DefaultAssay(scRNA_seurat) <- "TF"

top_tfs_feature <- tf_kw$TF[1:min(9, nrow(tf_kw))]
top_tfs_feature <- "MYC"

# Select available dimension reduction
reduction_use <- if ("umap_program_gene" %in% names(scRNA_seurat@reductions)) {
   "umap_program_gene"
} else if ("umap" %in% names(scRNA_seurat@reductions)) {
   "umap"
} else if ("tsne_program_gene" %in% names(scRNA_seurat@reductions)) {
   "tsne_program_gene"
} else if ("tsne" %in% names(scRNA_seurat@reductions)) {
   "tsne"
} else {
   stop("No UMAP/tSNE reduction found. Please run dimensional reduction first.")
}

p_tf_feature <- FeaturePlot(
   object = scRNA_seurat,
   features = top_tfs_feature,
   reduction = reduction_use,
   ncol = 3,
   cols = c("#F7FBFF", "#6BAED6", "#08306B"),
   min.cutoff = "q05",
   max.cutoff = "q95",
   pt.size = 0.4,
   order = TRUE
) &
   theme_classic(base_size = 13) &
   theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
   )

ggsave(
   filename = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/NMF_Result/TopTF_Activity_FeaturePlot_MYC.pdf",
   plot = p_tf_feature,
   width = 4.8,
   height = 4
)


# Extract plasma cells from different disease stages
load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF_0205/GSE232988/seu.Rdata")
# 1. Subset SMM cohort
seu_SMM <- subset(seu, subset = cohort %in% "SMM")
seu_SMM$cohort <- droplevels(seu_SMM$cohort)
seu_SMM$clonotype <- droplevels(seu_SMM$clonotype)
seu_SMM$batch <- droplevels(seu_SMM$batch)
# Subset Plasma Cells
seu_SMM_PC <- subset(seu_SMM, subset = cellType == "Plasma Cells")
seu_SMM_PC$cellType <- droplevels(seu_SMM_PC$cellType)

# 2. Subset MGUS cohort
seu_MGUS <- subset(seu, subset = cohort %in% "MGUS")
seu_MGUS$cohort <- droplevels(seu_MGUS$cohort)
seu_MGUS$clonotype <- droplevels(seu_MGUS$clonotype)
seu_MGUS$batch <- droplevels(seu_MGUS$batch)
seu_MGUS_PC <- subset(seu_MGUS, subset = cellType == "Plasma Cells")
seu_MGUS_PC$cellType <- droplevels(seu_MGUS_PC$cellType)

# 3. Subset HV cohort
seu_HV <- subset(seu, subset = cohort %in% "HV")
seu_HV$cohort <- droplevels(seu_HV$cohort)
seu_HV$clonotype <- droplevels(seu_HV$clonotype)
seu_HV$batch <- droplevels(seu_HV$batch)
seu_HV_PC <- subset(seu_HV, subset = cellType == "Plasma Cells")
seu_HV_PC$cellType <- droplevels(seu_HV_PC$cellType)


## Calculate evolution of cell states across HV, MGUS, SMM and MM
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(scales)

## =========================
## 0. Output directory
## =========================
outdir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Evolution_of_Four_States_AddModuleScore"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


## =========================
## 1. Set reference gene set based on scRNA_seurat
## =========================
DefaultAssay(scRNA_seurat) <- "RNA"

ref_genes <- rownames(scRNA_seurat)
length(ref_genes)

## =========================
## 2. Define function to unify objects to reference gene set
##    Fill missing genes with zero expression
## =========================
pad_seurat_to_refgenes <- function(srt, ref_genes, assay = "RNA") {
   DefaultAssay(srt) <- assay
   
   # Extract raw counts
   cnt <- GetAssayData(srt, assay = assay, slot = "counts")
   
   # Genes existing in current object
   cur_genes <- rownames(cnt)
   
   # Missing genes
   miss_genes <- setdiff(ref_genes, cur_genes)
   
   # Generate zero matrix for missing genes
   if (length(miss_genes) > 0) {
      zero_mat <- Matrix(
         0,
         nrow = length(miss_genes),
         ncol = ncol(cnt),
         sparse = TRUE,
         dimnames = list(miss_genes, colnames(cnt))
      )
      cnt2 <- rbind(cnt, zero_mat)
   } else {
      cnt2 <- cnt
   }
   
   # Reorder rows to match reference gene list
   cnt2 <- cnt2[ref_genes, , drop = FALSE]
   
   # Rebuild new Seurat object
   srt_new <- CreateSeuratObject(
      counts = cnt2,
      meta.data = srt@meta.data
   )
   
   # Original reductions/ident are not preserved; recalculate downstream if needed
   return(srt_new)
}

## =========================
## 3. Unify HV / MGUS / SMM to the gene space of MM
## =========================
seu_HV_PC2   <- pad_seurat_to_refgenes(seu_HV_PC,   ref_genes = ref_genes, assay = "originalexp")
seu_MGUS_PC2 <- pad_seurat_to_refgenes(seu_MGUS_PC, ref_genes = ref_genes, assay = "originalexp")
seu_SMM_PC2  <- pad_seurat_to_refgenes(seu_SMM_PC,  ref_genes = ref_genes, assay = "originalexp")

# Reorder MM object to identical gene sequence
mm_counts <- GetAssayData(scRNA_seurat, assay = "RNA", layer = "counts")
mm_counts <- mm_counts[ref_genes, , drop = FALSE]
scRNA_seurat2 <- CreateSeuratObject(
   counts = mm_counts,
   meta.data = scRNA_seurat@meta.data
)

## =========================
## 4. Add cohort metadata label
## =========================
seu_HV_PC2$cohort   <- "HV"
seu_MGUS_PC2$cohort <- "MGUS"
seu_SMM_PC2$cohort  <- "SMM"
scRNA_seurat2$cohort <- "MM"

## =========================
## 5. Merge four groups of plasma cells
## =========================
all_pc <- merge(
   x = seu_HV_PC2,
   y = list(seu_MGUS_PC2, seu_SMM_PC2, scRNA_seurat2),
   add.cell.ids = c("HV", "MGUS", "SMM", "MM"),
   project = "MM_progression"
)

DefaultAssay(all_pc) <- "RNA"

# Global normalization to generate data layer
all_pc <- NormalizeData(
   object = all_pc,
   assay = "RNA",
   normalization.method = "LogNormalize",
   scale.factor = 10000,
   verbose = FALSE
)
all_pc[["RNA"]] <- JoinLayers(all_pc[["RNA"]])
expr_mat <- GetAssayData(all_pc, assay = "RNA", layer = "data")

expr_mat <- as.matrix(expr_mat)

## =========================
## 6. Score all cells using MM meta-program gene sets
##    ls_MP must exist with names MM1/MM2/MM3
## =========================
print(names(ls_MP))

ls_MP_use <- lapply(ls_MP, function(g) intersect(unique(g), rownames(expr_mat)))
ls_MP_use <- ls_MP_use[sapply(ls_MP_use, length) >= 5]
print(sapply(ls_MP_use, length))


state_score_mat <- sapply(names(ls_MP_use), function(nm) {
   Matrix::colMeans(expr_mat[ls_MP_use[[nm]], , drop = FALSE])
})

state_score_mat <- as.data.frame(state_score_mat)
colnames(state_score_mat) <- paste0(names(ls_MP_use), "_Score")
rownames(state_score_mat) <- colnames(all_pc)

all_pc <- AddMetaData(all_pc, metadata = state_score_mat)

score_cols <- colnames(state_score_mat)
print(score_cols)

## =========================
## 7. Assign dominant_program label per cell
## =========================
score_df <- all_pc@meta.data[, score_cols, drop = FALSE]

dominant_program <- apply(score_df, 1, function(x) {
   names(x)[which.max(x)]
})
dominant_program <- sub("_Score$", "", dominant_program)

all_pc$dominant_program <- factor(
   dominant_program,
   levels = names(ls_MP_use)
)

all_pc$cohort <- factor(all_pc$cohort, levels = c("HV", "MGUS", "SMM", "MM"))

table(all_pc$cohort, all_pc$dominant_program)

## =========================
## 8. Statistical test: Differences of state scores across cohorts
## =========================
df <- all_pc@meta.data

stat_res <- lapply(score_cols, function(feat) {
   sub_df <- df[, c("cohort", feat), drop = FALSE]
   sub_df <- sub_df[complete.cases(sub_df), , drop = FALSE]
   kruskal.test(sub_df[[feat]] ~ sub_df$cohort)
})
names(stat_res) <- score_cols
print(stat_res)

pairwise_res <- lapply(score_cols, function(feat) {
   sub_df <- df[, c("cohort", feat), drop = FALSE]
   sub_df <- sub_df[complete.cases(sub_df), , drop = FALSE]
   pairwise.wilcox.test(
      x = sub_df[[feat]],
      g = sub_df$cohort,
      p.adjust.method = "BH"
   )
})
names(pairwise_res) <- score_cols
print(pairwise_res)



library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

# Color palette setup
cohort_cols <- c(
   "HV"   = "#8EBD80",
   "MGUS" = "#00A1D4",
   "SMM"  = "#8E6BBE",
   "MM"   = "#D73027"
)

# Fix factor display order
df$cohort <- factor(df$cohort, levels = c("HV", "MGUS", "SMM", "MM"))

# Define pairwise comparison groups
comparisons_list <- list(
   c("HV", "MGUS"),
   c("MGUS", "SMM"),
   c("SMM", "MM"),
   c("HV", "MM")
)

pdf(file.path(outdir, "MM_state_score_by_cohort_violin_publication.pdf"),
    width = 4, height = 3.5 * length(score_cols))

p_vln <- lapply(score_cols, function(feat) {
   
   sub_df <- df[, c("cohort", feat), drop = FALSE]
   colnames(sub_df) <- c("cohort", "score")
   sub_df <- sub_df[complete.cases(sub_df), , drop = FALSE]
   
   # Dynamically generate y-axis positions for significance brackets to avoid overlap
   y_max <- max(sub_df$score, na.rm = TRUE)
   y_min <- min(sub_df$score, na.rm = TRUE)
   y_range <- y_max - y_min
   if (y_range == 0) y_range <- abs(y_max) * 0.2 + 1e-6
   
   stat_y <- c(
      y_max + 0.10 * y_range,
      y_max + 0.22 * y_range,
      y_max + 0.34 * y_range,
      y_max + 0.46 * y_range
   )
   
   ggplot(sub_df, aes(x = cohort, y = score, fill = cohort, color = cohort)) +
      geom_violin(
         width = 0.85,
         linewidth = 0.45,
         trim = TRUE,
         alpha = 0.9
      ) +
      geom_boxplot(
         width = 0.16,
         outlier.shape = NA,
         fill = "white",
         color = "black",
         linewidth = 0.4
      ) +
      stat_compare_means(
         comparisons = comparisons_list,
         method = "wilcox.test",
         label = "p.signif",
         hide.ns = FALSE,
         size = 4,
         bracket.size = 0.4,
         tip.length = 0.01,
         label.y = stat_y
      ) +
      scale_fill_manual(values = cohort_cols) +
      scale_color_manual(values = cohort_cols) +
      labs(
         title = gsub("_Score$", "", feat),
         x = NULL,
         y = "Program score"
      ) +
      theme_classic(base_size = 13) +
      theme(
         plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 14
         ),
         axis.title.y = element_text(
            size = 12,
            face = "bold",
            colour = "black"
         ),
         axis.text.x = element_text(
            size = 11,
            face = "bold",
            colour = "black",
            angle = 0,
            vjust = 0.8
         ),
         axis.text.y = element_text(
            size = 10,
            colour = "black"
         ),
         legend.position = "none",
         axis.line = element_line(
            linewidth = 0.5,
            colour = "black"
         ),
         axis.ticks = element_line(
            linewidth = 0.4,
            colour = "black"
         ),
         axis.ticks.length = unit(0.18, "cm"),
         plot.margin = margin(8, 15, 8, 8)
      ) +
      coord_cartesian(
         ylim = c(y_min, y_max + 0.58 * y_range),
         clip = "off"
      )
})

print(wrap_plots(p_vln, ncol = 1))
dev.off()

## =========================
## 11. Plot: Proportion of dominant_program across cohorts
## =========================
state_frac <- df %>%
   group_by(cohort, dominant_program) %>%
   summarise(n = n(), .groups = "drop") %>%
   group_by(cohort) %>%
   mutate(freq = n / sum(n))

write.csv(state_frac,
          file.path(outdir, "MM_state_fraction_by_cohort.csv"),
          row.names = FALSE)
library(scales)
pdf(file.path(outdir, "MM_state_fraction_by_cohort_stackedbar.pdf"),
    width = 7, height = 5)

p_bar <- ggplot(state_frac, aes(x = cohort, y = freq, fill = dominant_program)) +
   geom_bar(stat = "identity", position = "fill") +
   scale_y_continuous(labels = percent_format()) +
   theme_bw() +
   ggtitle("Fraction of MM-defined states across cohorts") +
   xlab("Cohort") +
   ylab("Fraction of cells")

print(p_bar)
dev.off()

## =========================
## 12. Plot: Line chart showing evolutionary trends
## =========================
pdf(file.path(outdir, "MM_state_fraction_by_cohort_line.pdf"),
    width = 7, height = 5)

p_line <- ggplot(state_frac,
                 aes(x = cohort, y = freq, color = dominant_program, group = dominant_program)) +
   geom_point(size = 3) +
   geom_line(linewidth = 1) +
   scale_y_continuous(labels = percent_format()) +
   theme_bw() +
   ggtitle("Evolution of MM-defined states from HV to MM") +
   xlab("Cohort") +
   ylab("Fraction of cells")

print(p_line)
dev.off()

## =========================
## 13. Plot: Heatmap of average program scores
## =========================
mean_score_by_cohort <- sapply(score_cols, function(feat) {
   tapply(df[[feat]], df$cohort, mean, na.rm = TRUE)
})

mean_score_by_cohort <- t(mean_score_by_cohort)
mean_score_by_cohort <- mean_score_by_cohort[, c("HV", "MGUS", "SMM", "MM"), drop = FALSE]

pdf(file.path(outdir, "MM_state_mean_score_heatmap.pdf"),
    width = 6, height = 4)

p10 <- pheatmap(
   mean_score_by_cohort,
   scale = "row",
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   display_numbers = TRUE,
   main = "Mean MM-state scores across cohorts"
)
print(p10)
dev.off()


## =========================
## 15. Overall contingency table chi-square test
## =========================
chisq_res <- chisq.test(table(df$cohort, df$dominant_program))
print(chisq_res)

## =========================
## 16. UMAP visualization
## =========================
if (!"umap" %in% names(all_pc@reductions)) {
   all_pc <- FindVariableFeatures(all_pc, verbose = FALSE)
   all_pc <- ScaleData(all_pc, verbose = FALSE)
   all_pc <- RunPCA(all_pc, verbose = FALSE)
   all_pc <- RunUMAP(all_pc, dims = 1:20, verbose = FALSE)
}

pdf(file.path(outdir, "MM_state_UMAP_by_cohort.pdf"),
    width = 12, height = 5)

p1 <- DimPlot(all_pc, group.by = "cohort", reduction = "umap") +
   ggtitle("Cohort")
p2 <- DimPlot(all_pc, group.by = "dominant_program", reduction = "umap") +
   ggtitle("MM-defined states")

print(p1 + p2)
dev.off()

## =========================
## 17. Save Seurat object and full workspace
## =========================
save.image(file = "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/all_workspace_of_NMF.RData")
load("/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/all_workspace_of_NMF.RData")
