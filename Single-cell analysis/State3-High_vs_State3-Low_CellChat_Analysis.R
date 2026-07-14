############################################################
## Step4_State3High_vs_State3Low_CellChat_Exhaustion_fixedMM3markers.R
## Purpose:
##   Directly test whether State3-High plasma cells show stronger
##   immune-exhaustion / immune-suppressive communication than
##   State3-Low plasma cells.
##
## Design:
##   1) Use all_pc$dominant_program to define State3 plasma cells.
##   2) Use the user-defined fixed MM3 marker gene module.
##   3) Use Seurat::AddModuleScore to score State3 activity.
##   4) Split State3 plasma cells into State3-High and State3-Low.
##   5) Map these labels back to the full Seurat object `seu`.
##   6) Run CellChat using:
##        Plasma_State3_High, Plasma_State3_Low, Plasma_Other,
##        T/NK subtypes, myeloid/DC/B/HSPC/other groups.
##   7) Compare outgoing communication from State3-High vs State3-Low,
##      focusing on T/NK exhaustion and immune-suppressive axes.
##
## How to run:
##   Run this script AFTER Step0/Step1/Step2 objects are available, or
##   let it load the default workspace below.
##
## Required objects in workspace:
##   - seu: full GSE232988 Seurat object with all cell types
##   - all_pc: plasma-cell Seurat object with dominant_program MM1/MM2/MM3
##
## Output:
##   /home/yjliu/mmProj/data_process/Human/SingleCell_NMF/
##      State3High_vs_State3Low_CellChat_Exhaustion
############################################################

## =========================
## 0. User settings
## =========================
workspace_file <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/all_workspace_of_NMF.RData"
base_dir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF"

outdir <- file.path(base_dir, "State3High_vs_State3Low_CellChat_Exhaustion")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Step1/Step2 paths are used for consistency and optional prioritization.
step1_outdir <- file.path(base_dir, "State3_TNK_Exhaustion")
step2_outdir_primary <- file.path(base_dir, "State3_vs_State12_CellChat_Communication", "cohort_MM")

## Primary biological focus. MM-only is the cleanest comparison.
cohort_focus_primary <- c("MM")
auto_expand_cohorts_if_low_cells <- TRUE
cohort_scope_candidates <- list(
  "MM_only" = c("MM"),
  "SMM_MM" = c("SMM", "MM"),
  "MGUS_SMM_MM" = c("MGUS", "SMM", "MM"),
  "All_cohorts" = c("HV", "MGUS", "SMM", "MM")
)

## AddModuleScore settings.
## IMPORTANT: The State3/MM3 module is fixed by the user-defined marker list below.
## No FindMarkers or data-driven marker re-selection is performed in this script.
MM3_specific <- c(
  "ASPM", "CDCA5", "CGAS", "CGREF1", "CHAF1B", "CKAP2L", "CLSPN",
  "DMRT2", "DTL", "DUSP14", "E2F7", "ECT2", "FAM111B", "FANCD2",
  "FOXRED2", "GINS3", "INTS7", "LIN9", "MAD2L1", "MSRA", "ORC1",
  "RAD54L", "SGO1", "SLC35F2", "STIL", "TEDC2", "TRIM59", "TSPAN5",
  "WDR62", "WDR76", "ZNF367"
)
state3_high_low_method <- "median"  ## "median" or "rank_half"

## CellChat settings.
min_cells_per_group <- 20
max_cells_per_group <- 2500
seed_use <- 123
population_size_weighted <- FALSE
use_tnk_subtypes <- TRUE

## Candidate immune-exhaustion / immune-suppressive axes.
candidate_ligands <- c(
  "CD274", "PDCD1LG2", "LGALS9", "PVR", "NECTIN2", "HLA-E", "HLA-F", "HLA-G",
  "MIF", "TGFB1", "TGFB2", "IL10", "VEGFA", "CD47", "FASLG", "TNFSF10",
  "TNFSF13", "TNFSF13B", "CXCL12", "IL6", "SPP1", "ICAM1", "VCAM1", "GALECTIN"
)

candidate_receptors <- c(
  "PDCD1", "HAVCR2", "TIGIT", "CD96", "CD226", "KLRC1", "LAG3", "CTLA4",
  "CD74", "CXCR4", "SIRPA", "TGFBR1", "TGFBR2", "IL10RA", "IL10RB",
  "IL6R", "IL6ST", "CD44", "ITGA4", "ITGB1", "ITGAL", "NECTIN", "HLA"
)

candidate_pathway_pattern <- "PD.L1|PD.L2|TIGIT|GALECTIN|MIF|TGF|CD47|CXCL|IL6|SPP1|ICAM|VCAM|NECTIN|MHC|HLA|FAS|TRAIL|APRIL|BAFF"

## =========================
## 1. Packages and helpers
## =========================
if (file.exists(workspace_file) && (!exists("seu") || !exists("all_pc"))) {
  load(workspace_file)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop(
    "The R package 'CellChat' is required. Install it first, e.g. devtools::install_github('jinworks/CellChat').",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(CellChat))

if (!exists("seu")) stop("Object `seu` was not found. Please load Step0 workspace first.", call. = FALSE)
if (!exists("all_pc")) stop("Object `all_pc` was not found. Please load Step0/Step1 plasma-cell object first.", call. = FALSE)
if (!"cellType" %in% colnames(seu@meta.data)) stop("`seu@meta.data` must contain `cellType`.", call. = FALSE)
if (!"cohort" %in% colnames(seu@meta.data)) stop("`seu@meta.data` must contain `cohort`.", call. = FALSE)
if (!"dominant_program" %in% colnames(all_pc@meta.data)) {
  stop("`all_pc@meta.data` must contain `dominant_program` with MM1/MM2/MM3 or State1/2/3 labels.", call. = FALSE)
}

cohort_raw_levels <- c("HV", "MGUS", "SMM", "MM")
cohort_plot_levels <- c("Health", "MGUS", "SMM", "MM")

cohort_to_plot <- function(x) {
  x <- as.character(x)
  x[x == "HV"] <- "Health"
  factor(x, levels = cohort_plot_levels)
}

sanitize_group_name <- function(x) {
  x <- as.character(x)
  x <- gsub("/", "_", x)
  x <- gsub("\\+", "pos", x)
  x <- gsub("-", "_", x)
  x <- gsub(" ", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "", x)
  x <- gsub("_+", "_", x)
  x
}

strip_cohort_prefix <- function(x) {
  gsub("^(HV|MGUS|SMM|MM)_", "", as.character(x))
}

ggsave_pdf <- function(filename, plot, width, height) {
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = filename, plot = plot, width = width, height = height,
      units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      filename = filename, plot = plot, width = width, height = height,
      units = "in", device = "pdf", bg = "white", limitsize = FALSE
    )
  }
}

save_base_pdf <- function(filename, width = 8, height = 8, expr) {
  if (!is.null(grDevices::dev.list())) {
    for (i in grDevices::dev.list()) grDevices::dev.off()
  }
  grDevices::pdf(file = filename, width = width, height = height, onefile = TRUE, useDingbats = FALSE)
  tryCatch(eval(expr), error = function(e) message("Plot warning: ", e$message))
  grDevices::dev.off()
  invisible(NULL)
}

theme_pub <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.line        = ggplot2::element_line(linewidth = 0.55, colour = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.text        = ggplot2::element_text(colour = "black", size = base_size - 2),
      axis.title       = ggplot2::element_text(colour = "black", face = "bold", size = base_size),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey30", size = base_size - 1),
      legend.title     = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.size  = grid::unit(0.42, "cm"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "black", linewidth = 0.45),
      strip.text       = ggplot2::element_text(face = "bold", colour = "black"),
      plot.margin      = ggplot2::margin(8, 10, 8, 8)
    )
}

get_assay_data_safe <- function(object, assay, layer = c("counts", "data", "scale.data")) {
  layer <- match.arg(layer)
  if (packageVersion("SeuratObject") >= "5.0.0") {
    Seurat::GetAssayData(object, assay = assay, layer = layer)
  } else {
    Seurat::GetAssayData(object, assay = assay, slot = layer)
  }
}

has_data_layer <- function(object, assay) {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    "data" %in% SeuratObject::Layers(object[[assay]])
  } else {
    nrow(Seurat::GetAssayData(object, assay = assay, slot = "data")) > 0
  }
}

ensure_normalized <- function(object, assay) {
  Seurat::DefaultAssay(object) <- assay
  if (packageVersion("SeuratObject") >= "5.0.0") {
    object[[assay]] <- tryCatch(SeuratObject::JoinLayers(object[[assay]]), error = function(e) object[[assay]])
  }
  if (!has_data_layer(object, assay)) {
    object <- Seurat::NormalizeData(
      object = object,
      assay = assay,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
  object
}

p_to_label <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
  )
}

first_non_empty <- function(x) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) NA_character_ else as.character(y[1])
}

make_regex_pattern <- function(x) {
  x <- unique(toupper(x))
  x <- gsub("-", ".", x)
  paste(x, collapse = "|")
}

## =========================
## 2. Define State3 plasma cells and derive AddModuleScore
## =========================
state_program_map <- c(
  "MM1" = "Plasma_State1",
  "MM2" = "Plasma_State2",
  "MM3" = "Plasma_State3",
  "State1" = "Plasma_State1",
  "State2" = "Plasma_State2",
  "State3" = "Plasma_State3"
)

all_pc$dominant_program <- as.character(all_pc$dominant_program)
state3_level <- if ("MM3" %in% unique(all_pc$dominant_program)) {
  "MM3"
} else if ("State3" %in% unique(all_pc$dominant_program)) {
  "State3"
} else {
  stop("Cannot find MM3 or State3 in all_pc$dominant_program.", call. = FALSE)
}

all_pc$plasma_state_step4 <- unname(state_program_map[as.character(all_pc$dominant_program)])

assay_pc <- if ("originalexp" %in% Seurat::Assays(all_pc)) {
  "originalexp"
} else if ("RNA" %in% Seurat::Assays(all_pc)) {
  "RNA"
} else {
  Seurat::DefaultAssay(all_pc)
}

all_pc <- ensure_normalized(all_pc, assay = assay_pc)
Seurat::DefaultAssay(all_pc) <- assay_pc

existing_score_col <- if (paste0(state3_level, "_Score") %in% colnames(all_pc@meta.data)) {
  paste0(state3_level, "_Score")
} else if ("MM3_Score" %in% colnames(all_pc@meta.data)) {
  "MM3_Score"
} else if ("State3_Score" %in% colnames(all_pc@meta.data)) {
  "State3_Score"
} else {
  NA_character_
}

state3_cells_pc <- rownames(all_pc@meta.data)[as.character(all_pc$dominant_program) == state3_level]
other_pc_cells <- rownames(all_pc@meta.data)[as.character(all_pc$dominant_program) != state3_level]

if (length(state3_cells_pc) < 20) {
  stop("Fewer than 20 State3 plasma cells found in all_pc. Please check dominant_program labels.", call. = FALSE)
}
if (length(other_pc_cells) < 20) {
  stop("Fewer than 20 non-State3 plasma cells found in all_pc. Please check dominant_program labels.", call. = FALSE)
}

## Use the fixed MM3 marker gene set provided by the user.
## This is deliberately NOT re-derived from FindMarkers, so that State3-High/Low
## remains anchored to the previously defined MM3 biological program.
marker_input <- unique(MM3_specific)
assay_features <- rownames(all_pc)

## Exact matching first.
marker_exact <- intersect(marker_input, assay_features)

## Case-insensitive fallback for feature naming differences.
if (length(marker_exact) < length(marker_input)) {
  feature_map <- stats::setNames(assay_features, toupper(assay_features))
  marker_upper <- toupper(marker_input)
  marker_ci <- unname(feature_map[marker_upper])
  marker_ci <- marker_ci[!is.na(marker_ci)]
  state3_signature_genes <- unique(c(marker_exact, marker_ci))
} else {
  state3_signature_genes <- marker_exact
}

state3_marker_table <- data.frame(
  gene_symbol = marker_input,
  present_exact = marker_input %in% assay_features,
  present_case_insensitive = toupper(marker_input) %in% toupper(assay_features),
  assay_feature = unname(stats::setNames(assay_features, toupper(assay_features))[toupper(marker_input)]),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    used_for_AddModuleScore = !is.na(.data$assay_feature) & .data$assay_feature %in% state3_signature_genes
  )

write.csv(
  state3_marker_table,
  file.path(outdir, "Step4_fixed_MM3_marker_presence_table.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene = state3_signature_genes, stringsAsFactors = FALSE),
  file.path(outdir, "Step4_State3_AddModule_signature_genes.csv"),
  row.names = FALSE
)

message("Fixed MM3 marker genes provided: ", length(marker_input))
message("Fixed MM3 marker genes used in AddModuleScore: ", length(state3_signature_genes))

if (length(state3_signature_genes) < 10) {
  stop(
    "Fewer than 10 fixed MM3 marker genes were found in rownames(all_pc). ",
    "Please check whether the assay uses gene symbols or Ensembl IDs. See Step4_fixed_MM3_marker_presence_table.csv.",
    call. = FALSE
  )
}

all_pc <- Seurat::AddModuleScore(
  object = all_pc,
  features = list(State3_Module = state3_signature_genes),
  assay = assay_pc,
  name = "Step4_State3AddModule",
  ctrl = 25,
  nbin = 12,
  seed = seed_use
)

all_pc$State3_AddModuleScore_step4 <- all_pc@meta.data[["Step4_State3AddModule1"]]
score_col_step4 <- "State3_AddModuleScore_step4"

## Split State3 cells into State3-High and State3-Low.
set.seed(seed_use)
all_pc$State3_HL_step4 <- NA_character_
state3_score_vec <- all_pc@meta.data[state3_cells_pc, score_col_step4, drop = TRUE]

if (state3_high_low_method == "rank_half" || length(unique(state3_score_vec[!is.na(state3_score_vec)])) < 3) {
  rank_vec <- rank(state3_score_vec, ties.method = "first", na.last = "keep")
  cutoff_rank <- stats::median(rank_vec, na.rm = TRUE)
  all_pc$State3_HL_step4[state3_cells_pc] <- ifelse(rank_vec > cutoff_rank, "State3-High", "State3-Low")
  split_cutoff <- cutoff_rank
  split_note <- "rank_half"
} else {
  score_cutoff <- stats::median(state3_score_vec, na.rm = TRUE)
  all_pc$State3_HL_step4[state3_cells_pc] <- ifelse(state3_score_vec >= score_cutoff, "State3-High", "State3-Low")
  split_cutoff <- score_cutoff
  split_note <- "median_score"
}

state3_hl_counts <- as.data.frame(table(all_pc$State3_HL_step4, useNA = "ifany"))
colnames(state3_hl_counts) <- c("State3_HL_step4", "n_cells")
write.csv(state3_hl_counts, file.path(outdir, "Step4_State3HighLow_cell_counts_in_all_pc.csv"), row.names = FALSE)

pc_state_export <- all_pc@meta.data %>%
  tibble::rownames_to_column("all_pc_cell") %>%
  dplyr::mutate(
    plasma_state_step4 = as.character(.data$plasma_state_step4),
    dominant_program = as.character(.data$dominant_program),
    State3_HL_step4 = as.character(.data$State3_HL_step4),
    State3_AddModuleScore_step4 = as.numeric(.data$State3_AddModuleScore_step4),
    split_note = split_note,
    split_cutoff = split_cutoff
  )
write.csv(pc_state_export, file.path(outdir, "Step4_all_pc_State3HighLow_metadata.csv"), row.names = FALSE)

## =========================
## 3. Map State3-High/Low labels back to full `seu`
## =========================
pc_group_df <- all_pc@meta.data %>%
  tibble::rownames_to_column("all_pc_cell") %>%
  dplyr::mutate(
    orig_cell = strip_cohort_prefix(.data$all_pc_cell),
    plasma_group_step4 = dplyr::case_when(
      .data$plasma_state_step4 == "Plasma_State3" & .data$State3_HL_step4 == "State3-High" ~ "Plasma_State3_High",
      .data$plasma_state_step4 == "Plasma_State3" & .data$State3_HL_step4 == "State3-Low"  ~ "Plasma_State3_Low",
      .data$plasma_state_step4 %in% c("Plasma_State1", "Plasma_State2") ~ "Plasma_Other",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(.data$plasma_group_step4))

match_orig <- intersect(pc_group_df$orig_cell, colnames(seu))
match_full <- intersect(pc_group_df$all_pc_cell, colnames(seu))
map_key <- if (length(match_orig) >= length(match_full)) "orig_cell" else "all_pc_cell"
matched_cells <- if (map_key == "orig_cell") match_orig else match_full

if (length(matched_cells) < 20) {
  stop("Very few all_pc cells matched full seu cell names. Please inspect cell-name prefixes.", call. = FALSE)
}

pc_group_df_use <- pc_group_df %>%
  dplyr::filter(.data[[map_key]] %in% matched_cells) %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(map_key)), .keep_all = TRUE)

plasma_group_vec <- stats::setNames(pc_group_df_use$plasma_group_step4, pc_group_df_use[[map_key]])
score_vec_map <- stats::setNames(pc_group_df_use$State3_AddModuleScore_step4, pc_group_df_use[[map_key]])
state_vec_map <- stats::setNames(pc_group_df_use$plasma_state_step4, pc_group_df_use[[map_key]])

seu$plasma_group_step4 <- NA_character_
seu$plasma_state_step4 <- NA_character_
seu$State3_AddModuleScore_step4 <- NA_real_
seu$plasma_group_step4[names(plasma_group_vec)] <- plasma_group_vec
seu$plasma_state_step4[names(state_vec_map)] <- state_vec_map
seu$State3_AddModuleScore_step4[names(score_vec_map)] <- as.numeric(score_vec_map)

mapping_diagnostic <- data.frame(
  map_key = map_key,
  n_matched_cells = length(matched_cells),
  n_State3High_mapped = sum(seu$plasma_group_step4 == "Plasma_State3_High", na.rm = TRUE),
  n_State3Low_mapped = sum(seu$plasma_group_step4 == "Plasma_State3_Low", na.rm = TRUE),
  n_PlasmaOther_mapped = sum(seu$plasma_group_step4 == "Plasma_Other", na.rm = TRUE)
)
write.csv(mapping_diagnostic, file.path(outdir, "Step4_cellname_mapping_diagnostic.csv"), row.names = FALSE)

## =========================
## 4. Build communication groups
## =========================
major_map <- c(
  ## T / NK
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
  ## Other populations
  "Platelets"                     = "Platelets",
  "Macrophages"                   = "Monocytic Cells",
  "CD14 Monocytes"                = "Monocytic Cells",
  "TAMs"                          = "Monocytic Cells",
  "Plasma Cells"                  = "Plasma Cells",
  "pDCs"                          = "pDCs",
  "CD14 DCs"                      = "mDCs",
  "CLEC9A DCs"                    = "mDCs",
  "CD1C DCs"                      = "mDCs",
  "Immature Neutrophils"          = "Neutrophils",
  "HSPCs"                         = "HSPCs",
  "Memory B Cells"                = "B Cells",
  "Naive B Cells"                 = "B Cells",
  "PrePro B Cells"                = "B Cells",
  "Erythrocytes"                  = "Erythrocytes",
  "Erythroblasts"                 = "Erythrocytes",
  "Erythroid Progenitors"         = "Erythrocytes"
)

seu$majorCellType_step4 <- unname(major_map[as.character(seu$cellType)])
seu$majorCellType_step4[is.na(seu$majorCellType_step4)] <- "Unassigned"

seu$tnk_subtype_step4 <- dplyr::case_when(
  seu$cellType %in% c("CD8 Cytotoxic T Cells", "CD8 Activated T Cells", "CD8 Memory T Cells", "CD8 Naive T Cells", "MAIT Cytotoxic T Cells", "GD Cytotoxic T Cells") ~ "CD8_Cytotoxic_T",
  seu$cellType %in% c("NK Cells", "CD56 Bright NK Cells", "NKT Cells") ~ "NK_NKT",
  seu$cellType %in% c("Tregs") ~ "Treg",
  seu$cellType %in% c("CD4 Memory T Cells", "CD4 Naive T Cells", "CD4 Activated Memory T Cells", "CD4 Cytotoxic T Cells") ~ "CD4_T",
  seu$majorCellType_step4 == "T/NK Cells" ~ "Other_TNK",
  TRUE ~ NA_character_
)

meta_tmp <- seu@meta.data
comm_group <- rep(NA_character_, nrow(meta_tmp))
names(comm_group) <- rownames(meta_tmp)

is_plasma <- as.character(meta_tmp$cellType) == "Plasma Cells"
comm_group[is_plasma] <- meta_tmp$plasma_group_step4[is_plasma]

if (use_tnk_subtypes) {
  is_tnk <- meta_tmp$majorCellType_step4 == "T/NK Cells"
  comm_group[is_tnk] <- meta_tmp$tnk_subtype_step4[is_tnk]
} else {
  is_tnk <- meta_tmp$majorCellType_step4 == "T/NK Cells"
  comm_group[is_tnk] <- "T_NK_Cells"
}

is_nonplasma_nontnk <- !is_plasma & meta_tmp$majorCellType_step4 != "T/NK Cells"
comm_group[is_nonplasma_nontnk] <- sanitize_group_name(meta_tmp$majorCellType_step4[is_nonplasma_nontnk])
comm_group[is.na(comm_group) | comm_group == ""] <- "Unassigned"

seu$comm_group_step4 <- factor(comm_group)
seu$cohort_plot_step4 <- cohort_to_plot(seu$cohort)

cell_counts_all <- as.data.frame(table(seu$cohort, seu$comm_group_step4)) %>%
  dplyr::rename(cohort = Var1, comm_group_step4 = Var2, n_cells = Freq)
write.csv(cell_counts_all, file.path(outdir, "Step4_cell_counts_by_cohort_and_group_all_cells.csv"), row.names = FALSE)

## =========================
## 5. Select cohort scope for CellChat
## =========================
required_source_groups <- c("Plasma_State3_High", "Plasma_State3_Low")
immune_target_candidates <- c(
  "CD8_Cytotoxic_T", "NK_NKT", "Treg", "CD4_T", "Other_TNK", "T_NK_Cells",
  "Monocytic_Cells", "pDCs", "mDCs", "B_Cells", "Neutrophils"
)

scope_diagnostic <- dplyr::bind_rows(lapply(names(cohort_scope_candidates), function(scope_name) {
  cohorts_use <- cohort_scope_candidates[[scope_name]]
  cells_scope <- colnames(seu)[
    as.character(seu$cohort) %in% cohorts_use &
      as.character(seu$comm_group_step4) != "Unassigned" &
      !is.na(seu$comm_group_step4)
  ]
  counts <- table(as.character(seu$comm_group_step4[cells_scope]))
  n_high <- ifelse("Plasma_State3_High" %in% names(counts), counts[["Plasma_State3_High"]], 0)
  n_low <- ifelse("Plasma_State3_Low" %in% names(counts), counts[["Plasma_State3_Low"]], 0)
  n_immune_groups <- sum(names(counts) %in% immune_target_candidates & as.numeric(counts) >= min_cells_per_group)
  data.frame(
    scope_name = scope_name,
    cohorts = paste(cohorts_use, collapse = ";"),
    n_cells_total = length(cells_scope),
    n_State3High = as.integer(n_high),
    n_State3Low = as.integer(n_low),
    n_immune_groups_passing = as.integer(n_immune_groups),
    pass = n_high >= min_cells_per_group & n_low >= min_cells_per_group & n_immune_groups >= 1,
    stringsAsFactors = FALSE
  )
}))
write.csv(scope_diagnostic, file.path(outdir, "Step4_CellChat_scope_selection_diagnostic.csv"), row.names = FALSE)

if (auto_expand_cohorts_if_low_cells) {
  pass_scope <- scope_diagnostic %>%
    dplyr::filter(.data$pass) %>%
    dplyr::slice_head(n = 1)
  if (nrow(pass_scope) == 0) {
    stop("No cohort scope has enough State3-High, State3-Low, and immune target cells for CellChat.", call. = FALSE)
  }
  selected_scope_name <- pass_scope$scope_name[1]
  cohort_focus <- cohort_scope_candidates[[selected_scope_name]]
} else {
  selected_scope_name <- "manual_primary"
  cohort_focus <- cohort_focus_primary
}

message("Selected Step4 CellChat scope: ", selected_scope_name, " / cohorts: ", paste(cohort_focus, collapse = ", "))

cells_use <- colnames(seu)[
  as.character(seu$cohort) %in% cohort_focus &
    as.character(seu$comm_group_step4) != "Unassigned" &
    !is.na(seu$comm_group_step4)
]

seu_comm <- subset(seu, cells = cells_use)

assay_use <- if ("originalexp" %in% Seurat::Assays(seu_comm)) {
  "originalexp"
} else if ("RNA" %in% Seurat::Assays(seu_comm)) {
  "RNA"
} else {
  Seurat::DefaultAssay(seu_comm)
}

seu_comm <- ensure_normalized(seu_comm, assay = assay_use)
Seurat::DefaultAssay(seu_comm) <- assay_use

## Keep groups with enough cells.
group_counts <- table(seu_comm$comm_group_step4)
groups_keep <- names(group_counts)[group_counts >= min_cells_per_group]
seu_comm <- subset(seu_comm, subset = comm_group_step4 %in% groups_keep)
seu_comm$comm_group_step4 <- droplevels(seu_comm$comm_group_step4)

if (!all(required_source_groups %in% levels(seu_comm$comm_group_step4))) {
  stop("Plasma_State3_High and/or Plasma_State3_Low was removed by min_cells_per_group filtering.", call. = FALSE)
}

## Downsample to keep CellChat tractable.
set.seed(seed_use)
cells_down <- unlist(lapply(split(colnames(seu_comm), seu_comm$comm_group_step4), function(x) {
  if (length(x) > max_cells_per_group) sample(x, max_cells_per_group) else x
}), use.names = FALSE)
seu_comm <- subset(seu_comm, cells = cells_down)
seu_comm$comm_group_step4 <- droplevels(seu_comm$comm_group_step4)

cell_counts_used <- as.data.frame(table(seu_comm$comm_group_step4)) %>%
  dplyr::rename(comm_group_step4 = Var1, n_cells = Freq)
write.csv(cell_counts_used, file.path(outdir, "Step4_CellChat_cell_counts_used.csv"), row.names = FALSE)

immune_targets_use <- intersect(immune_target_candidates, levels(seu_comm$comm_group_step4))
write.csv(
  data.frame(immune_target = immune_targets_use, stringsAsFactors = FALSE),
  file.path(outdir, "Step4_immune_targets_used.csv"),
  row.names = FALSE
)

## =========================
## 6. Overview plots before CellChat
## =========================
pc_score_plot_df <- all_pc@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::filter(.data$plasma_state_step4 == "Plasma_State3") %>%
  dplyr::mutate(
    State3_HL_step4 = factor(.data$State3_HL_step4, levels = c("State3-Low", "State3-High"))
  )

p_score_violin <- ggplot2::ggplot(
  pc_score_plot_df,
  ggplot2::aes(x = .data$State3_HL_step4, y = .data$State3_AddModuleScore_step4, fill = .data$State3_HL_step4)
) +
  ggplot2::geom_violin(width = 0.85, trim = FALSE, alpha = 0.72, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.45) +
  ggplot2::geom_jitter(width = 0.12, size = 0.28, alpha = 0.28) +
  ggplot2::scale_fill_manual(values = c("State3-Low" = "#4575B4", "State3-High" = "#D73027"), drop = FALSE) +
  ggplot2::labs(
    title = "State3 AddModuleScore stratification",
    subtitle = paste0("Split: ", split_note, "; State3 signature genes: ", ifelse(is.null(state3_signature_genes), "existing score", length(state3_signature_genes))),
    x = NULL,
    y = "State3 AddModuleScore"
  ) +
  theme_pub(13) +
  ggplot2::theme(legend.position = "none")

ggsave_pdf(file.path(outdir, "A_State3_AddModuleScore_HighLow_definition.pdf"), p_score_violin, width = 5.0, height = 4.6)

reduction_use <- if ("X_umap" %in% Seurat::Reductions(seu_comm)) {
  "X_umap"
} else if ("umap" %in% Seurat::Reductions(seu_comm)) {
  "umap"
} else {
  NA_character_
}

if (!is.na(reduction_use)) {
  p_umap <- Seurat::DimPlot(
    seu_comm,
    reduction = reduction_use,
    group.by = "comm_group_step4",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.22
  ) +
    theme_pub(12) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::ggtitle("Step4 CellChat communication groups")

  ggsave_pdf(file.path(outdir, "B_CellChat_groups_UMAP_State3HighLow.pdf"), p_umap, width = 9.5, height = 7.5)
}

## =========================
## 7. Run CellChat
## =========================
options(future.globals.maxSize = 20 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

expr_mat <- get_assay_data_safe(seu_comm, assay = assay_use, layer = "data")
meta_input <- seu_comm@meta.data %>%
  dplyr::select(dplyr::all_of(c("comm_group_step4", "cohort", "cohort_plot_step4", "cellType", "majorCellType_step4", "plasma_group_step4"))) %>%
  as.data.frame()
meta_input$comm_group_step4 <- droplevels(factor(meta_input$comm_group_step4))
stopifnot(identical(colnames(expr_mat), rownames(meta_input)))

cellchat <- CellChat::createCellChat(
  object = expr_mat,
  meta = meta_input,
  group.by = "comm_group_step4"
)

data("CellChatDB.human", package = "CellChat")
cellchat@DB <- CellChatDB.human

cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  raw.use = TRUE,
  population.size = population_size_weighted
)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_per_group)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)
cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file.path(outdir, "Step4_CellChat_object_State3High_vs_State3Low.rds"))

## Overall CellChat network plots.
group_size <- as.numeric(table(cellchat@idents))
names(group_size) <- names(table(cellchat@idents))

save_base_pdf(
  file.path(outdir, "C_CellChat_overall_network_count_circle.pdf"),
  width = 8.5,
  height = 8.5,
  expr = {
    CellChat::netVisual_circle(
      cellchat@net$count,
      vertex.weight = group_size,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = "Number of interactions"
    )
  }
)

save_base_pdf(
  file.path(outdir, "D_CellChat_overall_network_weight_circle.pdf"),
  width = 8.5,
  height = 8.5,
  expr = {
    CellChat::netVisual_circle(
      cellchat@net$weight,
      vertex.weight = group_size,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = "Interaction strength"
    )
  }
)

## =========================
## 8. Extract and compare communication: State3-High vs State3-Low
## =========================
comm_df <- CellChat::subsetCommunication(cellchat)
comm_df <- as.data.frame(comm_df)

if (!"prob" %in% colnames(comm_df)) {
  stop("subsetCommunication(cellchat) did not return a `prob` column. Please inspect your CellChat version.", call. = FALSE)
}

if (!"interaction_name_2" %in% colnames(comm_df)) {
  if (all(c("ligand", "receptor") %in% colnames(comm_df))) {
    comm_df$interaction_name_2 <- paste(comm_df$ligand, comm_df$receptor, sep = " - ")
  } else if ("interaction_name" %in% colnames(comm_df)) {
    comm_df$interaction_name_2 <- comm_df$interaction_name
  } else {
    comm_df$interaction_name_2 <- paste0("LR_", seq_len(nrow(comm_df)))
  }
}

if (!"pathway_name" %in% colnames(comm_df)) comm_df$pathway_name <- "Unknown_pathway"
if (!"ligand" %in% colnames(comm_df)) comm_df$ligand <- NA_character_
if (!"receptor" %in% colnames(comm_df)) comm_df$receptor <- NA_character_
if (!"pval" %in% colnames(comm_df)) comm_df$pval <- NA_real_

comm_df <- comm_df %>%
  dplyr::mutate(
    source = as.character(.data$source),
    target = as.character(.data$target),
    pathway_name = as.character(.data$pathway_name),
    interaction_name_2 = as.character(.data$interaction_name_2),
    ligand = as.character(.data$ligand),
    receptor = as.character(.data$receptor)
  )

write.csv(comm_df, file.path(outdir, "Step4_CellChat_all_significant_LR_interactions.csv"), row.names = FALSE)

source_groups <- c("Plasma_State3_High", "Plasma_State3_Low")
immune_targets_use <- intersect(immune_target_candidates, unique(comm_df$target))
if (length(immune_targets_use) == 0) {
  stop("No immune target group was found in CellChat communication table.", call. = FALSE)
}

## 8.1 Total outgoing communication by immune target.
out_target <- comm_df %>%
  dplyr::filter(.data$source %in% source_groups, .data$target %in% immune_targets_use) %>%
  dplyr::group_by(.data$source, .data$target) %>%
  dplyr::summarise(
    prob_sum = sum(.data$prob, na.rm = TRUE),
    n_LR = dplyr::n_distinct(.data$interaction_name_2),
    .groups = "drop"
  )

out_target_full <- tidyr::expand_grid(
  source = source_groups,
  target = immune_targets_use
) %>%
  dplyr::left_join(out_target, by = c("source", "target")) %>%
  dplyr::mutate(
    prob_sum = tidyr::replace_na(.data$prob_sum, 0),
    n_LR = tidyr::replace_na(.data$n_LR, 0L)
  )

out_target_diff <- out_target_full %>%
  tidyr::pivot_wider(names_from = source, values_from = c(prob_sum, n_LR), values_fill = 0) %>%
  dplyr::mutate(
    High_prob = .data[["prob_sum_Plasma_State3_High"]],
    Low_prob = .data[["prob_sum_Plasma_State3_Low"]],
    High_minus_Low = .data$High_prob - .data$Low_prob,
    High_ratio_vs_Low = (.data$High_prob + 1e-8) / (.data$Low_prob + 1e-8),
    High_nLR = .data[["n_LR_Plasma_State3_High"]],
    Low_nLR = .data[["n_LR_Plasma_State3_Low"]]
  ) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low))

write.csv(out_target_diff, file.path(outdir, "Step4_State3High_vs_Low_outgoing_total_by_immune_target.csv"), row.names = FALSE)

## 8.2 Pathway-level comparison.
out_pathway <- comm_df %>%
  dplyr::filter(.data$source %in% source_groups, .data$target %in% immune_targets_use) %>%
  dplyr::group_by(.data$source, .data$target, .data$pathway_name) %>%
  dplyr::summarise(
    prob_sum = sum(.data$prob, na.rm = TRUE),
    n_LR = dplyr::n_distinct(.data$interaction_name_2),
    .groups = "drop"
  )

pathway_key <- out_pathway %>%
  dplyr::distinct(.data$target, .data$pathway_name) %>%
  dplyr::mutate(pathway_id = dplyr::row_number())

out_pathway_full <- tidyr::expand_grid(
  source = source_groups,
  pathway_id = pathway_key$pathway_id
) %>%
  dplyr::left_join(pathway_key, by = "pathway_id") %>%
  dplyr::select(-dplyr::all_of("pathway_id")) %>%
  dplyr::left_join(out_pathway, by = c("source", "target", "pathway_name")) %>%
  dplyr::mutate(
    prob_sum = tidyr::replace_na(.data$prob_sum, 0),
    n_LR = tidyr::replace_na(.data$n_LR, 0L)
  )

out_pathway_diff <- out_pathway_full %>%
  tidyr::pivot_wider(names_from = source, values_from = c(prob_sum, n_LR), values_fill = 0) %>%
  dplyr::mutate(
    High_prob = .data[["prob_sum_Plasma_State3_High"]],
    Low_prob = .data[["prob_sum_Plasma_State3_Low"]],
    High_minus_Low = .data$High_prob - .data$Low_prob,
    High_ratio_vs_Low = (.data$High_prob + 1e-8) / (.data$Low_prob + 1e-8),
    High_nLR = .data[["n_LR_Plasma_State3_High"]],
    Low_nLR = .data[["n_LR_Plasma_State3_Low"]]
  ) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low), dplyr::desc(.data$High_prob))

write.csv(out_pathway_diff, file.path(outdir, "Step4_State3High_vs_Low_outgoing_pathway_by_immune_target.csv"), row.names = FALSE)

## 8.3 Ligand-receptor-level comparison.
lr_cols <- c("source", "target", "pathway_name", "interaction_name_2", "ligand", "receptor")
out_lr <- comm_df %>%
  dplyr::filter(.data$source %in% source_groups, .data$target %in% immune_targets_use) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(lr_cols))) %>%
  dplyr::summarise(
    prob_sum = sum(.data$prob, na.rm = TRUE),
    min_pval = suppressWarnings(min(.data$pval, na.rm = TRUE)),
    .groups = "drop"
  )

lr_key <- out_lr %>%
  dplyr::group_by(.data$target, .data$pathway_name, .data$interaction_name_2) %>%
  dplyr::summarise(
    ligand = first_non_empty(.data$ligand),
    receptor = first_non_empty(.data$receptor),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lr_id = dplyr::row_number())

out_lr_full <- tidyr::expand_grid(
  source = source_groups,
  lr_id = lr_key$lr_id
) %>%
  dplyr::left_join(lr_key, by = "lr_id") %>%
  dplyr::select(-dplyr::all_of("lr_id")) %>%
  dplyr::left_join(
    out_lr %>%
      dplyr::select(dplyr::all_of(c("source", "target", "pathway_name", "interaction_name_2", "prob_sum", "min_pval"))),
    by = c("source", "target", "pathway_name", "interaction_name_2")
  ) %>%
  dplyr::mutate(prob_sum = tidyr::replace_na(.data$prob_sum, 0))

ligand_pattern <- make_regex_pattern(candidate_ligands)
receptor_pattern <- make_regex_pattern(candidate_receptors)

out_lr_diff <- out_lr_full %>%
  tidyr::pivot_wider(names_from = source, values_from = prob_sum, values_fill = 0) %>%
  dplyr::mutate(
    High_prob = .data[["Plasma_State3_High"]],
    Low_prob = .data[["Plasma_State3_Low"]],
    High_minus_Low = .data$High_prob - .data$Low_prob,
    High_ratio_vs_Low = (.data$High_prob + 1e-8) / (.data$Low_prob + 1e-8),
    ligand_upper = toupper(.data$ligand),
    receptor_upper = toupper(.data$receptor),
    pathway_upper = toupper(.data$pathway_name),
    is_exhaustion_axis =
      grepl(ligand_pattern, .data$ligand_upper) |
      grepl(receptor_pattern, .data$receptor_upper) |
      grepl(candidate_pathway_pattern, .data$pathway_upper),
    High_biased = .data$High_minus_Low > 0,
    LR_label = paste0(.data$ligand, " → ", .data$receptor)
  ) %>%
  dplyr::arrange(dplyr::desc(.data$is_exhaustion_axis), dplyr::desc(.data$High_minus_Low), dplyr::desc(.data$High_prob))

write.csv(out_lr_diff, file.path(outdir, "Step4_State3High_vs_Low_outgoing_LR_pairs_to_immune_targets.csv"), row.names = FALSE)

high_biased_exhaustion_lr <- out_lr_diff %>%
  dplyr::filter(.data$is_exhaustion_axis, .data$High_prob > 0, .data$High_minus_Low > 0) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low), dplyr::desc(.data$High_prob))
write.csv(high_biased_exhaustion_lr, file.path(outdir, "Step4_HighBiased_exhaustion_LR_axes.csv"), row.names = FALSE)

exhaustion_summary_by_target <- out_lr_diff %>%
  dplyr::group_by(.data$target) %>%
  dplyr::summarise(
    High_total_immune_strength = sum(.data$High_prob, na.rm = TRUE),
    Low_total_immune_strength = sum(.data$Low_prob, na.rm = TRUE),
    High_exhaustion_strength = sum(.data$High_prob[.data$is_exhaustion_axis], na.rm = TRUE),
    Low_exhaustion_strength = sum(.data$Low_prob[.data$is_exhaustion_axis], na.rm = TRUE),
    n_high_biased_exhaustion_LR = sum(.data$is_exhaustion_axis & .data$High_minus_Low > 0 & .data$High_prob > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    High_minus_Low_exhaustion = .data$High_exhaustion_strength - .data$Low_exhaustion_strength,
    High_ratio_vs_Low_exhaustion = (.data$High_exhaustion_strength + 1e-8) / (.data$Low_exhaustion_strength + 1e-8)
  ) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low_exhaustion))
write.csv(exhaustion_summary_by_target, file.path(outdir, "Step4_exhaustion_communication_summary_by_target.csv"), row.names = FALSE)

overall_exhaustion_summary <- out_lr_diff %>%
  dplyr::summarise(
    High_total_immune_strength = sum(.data$High_prob, na.rm = TRUE),
    Low_total_immune_strength = sum(.data$Low_prob, na.rm = TRUE),
    High_exhaustion_strength = sum(.data$High_prob[.data$is_exhaustion_axis], na.rm = TRUE),
    Low_exhaustion_strength = sum(.data$Low_prob[.data$is_exhaustion_axis], na.rm = TRUE),
    n_high_biased_exhaustion_LR = sum(.data$is_exhaustion_axis & .data$High_minus_Low > 0 & .data$High_prob > 0, na.rm = TRUE),
    n_exhaustion_LR_total = sum(.data$is_exhaustion_axis, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    High_minus_Low_exhaustion = .data$High_exhaustion_strength - .data$Low_exhaustion_strength,
    High_ratio_vs_Low_exhaustion = (.data$High_exhaustion_strength + 1e-8) / (.data$Low_exhaustion_strength + 1e-8)
  )
write.csv(overall_exhaustion_summary, file.path(outdir, "Step4_overall_exhaustion_communication_summary.csv"), row.names = FALSE)

## Optional: bring Step2 LR-pair results forward, if present, to identify overlap with prior State3-vs-State1/2 axes.
step2_lr_file <- file.path(step2_outdir_primary, "State3_vs_State12_outgoing_LR_pairs.csv")
if (file.exists(step2_lr_file)) {
  step2_lr <- read.csv(step2_lr_file, stringsAsFactors = FALSE)
  if ("interaction_name_2" %in% colnames(step2_lr)) {
    step2_overlap <- out_lr_diff %>%
      dplyr::mutate(Step2_State3_vs_State12_candidate = .data$interaction_name_2 %in% step2_lr$interaction_name_2) %>%
      dplyr::arrange(dplyr::desc(.data$Step2_State3_vs_State12_candidate), dplyr::desc(.data$is_exhaustion_axis), dplyr::desc(.data$High_minus_Low))
    write.csv(step2_overlap, file.path(outdir, "Step4_LR_pairs_with_Step2_overlap_annotation.csv"), row.names = FALSE)
  }
}

## =========================
## 9. Publication-style figures
## =========================
source_cols <- c("Plasma_State3_Low" = "#4575B4", "Plasma_State3_High" = "#D73027")
hilo_cols <- c("State3-Low" = "#4575B4", "State3-High" = "#D73027")

## 9.1 Total outgoing communication to immune targets.
out_target_long <- out_target_full %>%
  dplyr::mutate(
    source_label = dplyr::recode(.data$source, "Plasma_State3_High" = "State3-High", "Plasma_State3_Low" = "State3-Low"),
    source_label = factor(.data$source_label, levels = c("State3-Low", "State3-High")),
    target = factor(.data$target, levels = out_target_diff$target[order(out_target_diff$High_minus_Low)])
  )

p_target_bar <- ggplot2::ggplot(
  out_target_long,
  ggplot2::aes(x = .data$target, y = .data$prob_sum, fill = .data$source_label)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.66, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = hilo_cols, name = NULL) +
  ggplot2::labs(
    title = "State3-High plasma cells show stronger outgoing communication",
    subtitle = paste0("Cohorts: ", paste(cohort_focus, collapse = ", "), "; CellChat source comparison"),
    x = "Immune / microenvironment target",
    y = "Outgoing communication strength"
  ) +
  theme_pub(13) +
  ggplot2::theme(legend.position = "top")

ggsave_pdf(file.path(outdir, "E_State3High_vs_Low_total_outgoing_to_immune_targets.pdf"), p_target_bar, width = 7.5, height = 5.8)

## 9.2 Exhaustion communication summary by target.
exhaustion_long <- exhaustion_summary_by_target %>%
  dplyr::select(dplyr::all_of(c("target", "High_exhaustion_strength", "Low_exhaustion_strength", "High_minus_Low_exhaustion"))) %>%
  tidyr::pivot_longer(
    cols = c(High_exhaustion_strength, Low_exhaustion_strength),
    names_to = "source_label",
    values_to = "exhaustion_strength"
  ) %>%
  dplyr::mutate(
    source_label = dplyr::recode(
      .data$source_label,
      "High_exhaustion_strength" = "State3-High",
      "Low_exhaustion_strength" = "State3-Low"
    ),
    source_label = factor(.data$source_label, levels = c("State3-Low", "State3-High")),
    target = factor(.data$target, levels = exhaustion_summary_by_target$target[order(exhaustion_summary_by_target$High_minus_Low_exhaustion)])
  )

p_exh_bar <- ggplot2::ggplot(
  exhaustion_long,
  ggplot2::aes(x = .data$target, y = .data$exhaustion_strength, fill = .data$source_label)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.66, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = hilo_cols, name = NULL) +
  ggplot2::labs(
    title = "Immune-exhaustion communication is enriched in State3-High",
    subtitle = "Candidate checkpoint, TIGIT, galectin, MIF, TGF, CD47, CXCL12 and adhesion axes",
    x = "Receiver cell group",
    y = "Exhaustion-associated communication strength"
  ) +
  theme_pub(13) +
  ggplot2::theme(legend.position = "top")

ggsave_pdf(file.path(outdir, "F_State3High_vs_Low_exhaustion_communication_by_target.pdf"), p_exh_bar, width = 7.6, height = 5.8)

## 9.3 High-biased exhaustion LR dotplot.
top_exh_lr_plot <- high_biased_exhaustion_lr %>%
  dplyr::filter(.data$High_prob > 0) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low), dplyr::desc(.data$High_prob)) %>%
  dplyr::slice_head(n = 35)

if (nrow(top_exh_lr_plot) > 0) {
  top_exh_lr_plot <- top_exh_lr_plot %>%
    dplyr::mutate(
      LR_label_unique = paste0(.data$LR_label, " | ", .data$pathway_name),
      LR_label_unique = factor(.data$LR_label_unique, levels = rev(unique(.data$LR_label_unique))),
      target = factor(.data$target, levels = immune_targets_use)
    )

  p_lr_dot <- ggplot2::ggplot(top_exh_lr_plot, ggplot2::aes(x = .data$target, y = .data$LR_label_unique)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$High_prob, color = .data$High_minus_Low), alpha = 0.92) +
    ggplot2::scale_size(range = c(1.4, 7.8), name = "State3-High strength") +
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "High - Low"
    ) +
    ggplot2::labs(
      title = "State3-High-biased exhaustion ligand-receptor axes",
      subtitle = "Top immune-suppressive LR pairs with higher communication in State3-High plasma cells",
      x = "Receiver cell group",
      y = "Ligand → receptor | pathway"
    ) +
    theme_pub(11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
      legend.position = "right"
    )

  ggsave_pdf(file.path(outdir, "G_Top_State3High_biased_exhaustion_LR_pairs.pdf"), p_lr_dot, width = 10.2, height = 10.5)
}

## 9.4 Pathway dotplot.
top_pathways <- out_pathway_diff %>%
  dplyr::filter(.data$High_prob > 0) %>%
  dplyr::arrange(dplyr::desc(.data$High_minus_Low), dplyr::desc(.data$High_prob)) %>%
  dplyr::slice_head(n = 30)

if (nrow(top_pathways) > 0) {
  top_pathways <- top_pathways %>%
    dplyr::mutate(
      pathway_target = paste0(.data$pathway_name, " → ", .data$target),
      pathway_target = factor(.data$pathway_target, levels = rev(unique(.data$pathway_target)))
    )

  p_pathway_dot <- ggplot2::ggplot(top_pathways, ggplot2::aes(x = .data$target, y = .data$pathway_target)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$High_prob, color = .data$High_minus_Low), alpha = 0.92) +
    ggplot2::scale_size(range = c(1.2, 7.5), name = "State3-High strength") +
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "High - Low"
    ) +
    ggplot2::labs(
      title = "State3-High-biased communication pathways",
      subtitle = "Outgoing pathways from State3-High plasma cells to immune targets",
      x = "Receiver cell group",
      y = "Pathway → target"
    ) +
    theme_pub(11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
      legend.position = "right"
    )

  ggsave_pdf(file.path(outdir, "H_Top_State3High_biased_pathways_to_immune_targets.pdf"), p_pathway_dot, width = 9.6, height = 9.2)
}

## 9.5 Heatmap-style tile of High - Low exhaustion communication by pathway and target.
exh_pathway_heat <- out_pathway_diff %>%
  dplyr::mutate(
    pathway_upper = toupper(.data$pathway_name),
    is_candidate_pathway = grepl(candidate_pathway_pattern, .data$pathway_upper)
  ) %>%
  dplyr::filter(.data$is_candidate_pathway | .data$High_minus_Low > 0) %>%
  dplyr::group_by(.data$pathway_name) %>%
  dplyr::mutate(max_abs_diff = max(abs(.data$High_minus_Low), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(.data$max_abs_diff)) %>%
  dplyr::filter(.data$pathway_name %in% unique(.data$pathway_name)[seq_len(min(20, length(unique(.data$pathway_name))))]) %>%
  dplyr::mutate(
    pathway_name = factor(.data$pathway_name, levels = rev(unique(.data$pathway_name))),
    target = factor(.data$target, levels = immune_targets_use)
  )

if (nrow(exh_pathway_heat) > 0) {
  p_heat <- ggplot2::ggplot(exh_pathway_heat, ggplot2::aes(x = .data$target, y = .data$pathway_name, fill = .data$High_minus_Low)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$High_prob > 0, sprintf("%.2g", .data$High_prob), "")), size = 2.8) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "High - Low"
    ) +
    ggplot2::labs(
      title = "State3-High communication bias across immune pathways",
      subtitle = "Tile color: High - Low; label: State3-High strength",
      x = "Receiver cell group",
      y = "CellChat pathway"
    ) +
    theme_pub(11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5)
    )

  ggsave_pdf(file.path(outdir, "I_State3High_minus_Low_pathway_heatmap.pdf"), p_heat, width = 8.8, height = 7.6)
}

## 9.6 CellChat bubble plot for direct inspection.
save_base_pdf(
  file.path(outdir, "J_CellChat_bubble_State3HighLow_to_immune_targets.pdf"),
  width = 12,
  height = 8,
  expr = {
    CellChat::netVisual_bubble(
      cellchat,
      sources.use = source_groups,
      targets.use = immune_targets_use,
      remove.isolate = FALSE,
      angle.x = 45
    )
  }
)

## 9.7 Combined summary figure.
combined_panels <- p_score_violin + p_target_bar + p_exh_bar + patchwork::plot_layout(ncol = 3, widths = c(0.85, 1.25, 1.25))
ggsave_pdf(file.path(outdir, "Figure_Step4_State3High_ExhaustionCommunication_summary.pdf"), combined_panels, width = 17, height = 5.5)

## =========================
## 10. Concise result text table for manuscript writing
## =========================
result_summary_text <- data.frame(
  item = c(
    "selected_scope",
    "cohorts_used",
    "State3_high_cells_used",
    "State3_low_cells_used",
    "immune_targets_used",
    "overall_High_exhaustion_strength",
    "overall_Low_exhaustion_strength",
    "overall_High_minus_Low_exhaustion",
    "overall_High_ratio_vs_Low_exhaustion",
    "n_high_biased_exhaustion_LR"
  ),
  value = c(
    selected_scope_name,
    paste(cohort_focus, collapse = ";"),
    as.character(cell_counts_used$n_cells[cell_counts_used$comm_group_step4 == "Plasma_State3_High"]),
    as.character(cell_counts_used$n_cells[cell_counts_used$comm_group_step4 == "Plasma_State3_Low"]),
    paste(immune_targets_use, collapse = ";"),
    as.character(overall_exhaustion_summary$High_exhaustion_strength[1]),
    as.character(overall_exhaustion_summary$Low_exhaustion_strength[1]),
    as.character(overall_exhaustion_summary$High_minus_Low_exhaustion[1]),
    as.character(overall_exhaustion_summary$High_ratio_vs_Low_exhaustion[1]),
    as.character(overall_exhaustion_summary$n_high_biased_exhaustion_LR[1])
  ),
  stringsAsFactors = FALSE
)
write.csv(result_summary_text, file.path(outdir, "Step4_result_summary_for_manuscript.csv"), row.names = FALSE)

message("Step4 completed successfully.")
message("Output directory: ", outdir)
message("Key files:")
message("  - Step4_HighBiased_exhaustion_LR_axes.csv")
message("  - Step4_exhaustion_communication_summary_by_target.csv")
message("  - Figure_Step4_State3High_ExhaustionCommunication_summary.pdf")
message("  - G_Top_State3High_biased_exhaustion_LR_pairs.pdf")
