#!/usr/bin/env Rscript
# ============================================================================
# State 2 / State 3 independent validation in GSE122231 / SWOG S0120
#
# Logic inherited from GSE117847:
#   1. Fixed State2 / State3 gene programs.
#   2. Collapse duplicate gene symbols.
#   3. Rank-based within-sample program score.
#   4. Export gene coverage, sample-level scores, statistics, and figures.
#
# Key difference from GSE117847:
#   GSE122231 is a baseline AMG / MGUS / AMM cohort, not paired P-SMM -> MM.
#   Therefore:
#      - no paired Wilcoxon test
#      - no P-SMM -> MM delta table
#      - main validation is progression-risk / time-to-progression analysis
# ============================================================================

suppressPackageStartupMessages({
   library(GEOquery)
   library(Biobase)
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(pROC)
   library(pheatmap)
   library(survival)
   library(scales)
   library(grid)
})

# -----------------------------------------------------------------------------
# 0. User settings
# -----------------------------------------------------------------------------

family_soft_path <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE122231/GSE122231_family.soft"

outdir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/State2State3_Analysis/GSE122231_State2State3_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260703)

# Optional:
# If you have a separate clinical metadata table containing progression status
# and time-to-progression, put the path here.
# It should contain gsm_id or patient_id, plus event/time columns.
clinical_csv_path <- NULL
# clinical_csv_path <- "/path/to/GSE122231_clinical_metadata.csv"

# Optional manual overrides.
# Leave as NULL first. If automatic detection fails, inspect:
#   Table_GSE122231_metadata_dictionary.csv
# then fill these.
manual_patient_col   <- NULL
manual_diagnosis_col <- NULL
manual_event_col     <- NULL
manual_time_col      <- NULL

# -----------------------------------------------------------------------------
# 1. Fixed State2 / State3 programs
# -----------------------------------------------------------------------------

State2_genes <- c(
   "MCC", "CADPS2", "SFMBT2", "COBLL1", "SAMD12", "KIAA1217",
   "COL6A3", "ADGRB3", "COL4A5", "EPHA6", "NEU3", "NBEA",
   "RASSF6", "STARD9", "CNTN5", "NCAM1", "SYT1", "FMN1",
   "ATP10B", "PCDH9", "RELN", "MEF2C", "KHDRBS2", "DDX31",
   "RAPGEF5", "TMTC2", "HOMER1", "NDNF", "AJAP1", "NEB",
   "ESRRG", "MAP2", "TTC28", "PRKG1", "DCC", "SP4"
)

State3_genes <- c(
   "MAD2L1", "CLSPN", "ASPM", "SGO1", "CDCA5", "WDR76", "TEDC2", "ZNF367",
   "DTL", "FAM111B", "NUP155", "FANCD2", "CKAP2L", "ECT2", "WDR62", "RAD54L",
   "CGAS", "LIN9", "INTS7", "BDH1", "STIL", "E2F7", "FAM133A", "CHAF1B",
   "TSPAN5", "SLC35F2", "TRIM59", "GINS3", "ORC1", "CGREF1", "DMRT2",
   "FOXRED2", "DUSP14", "MSRA"
)

State2_genes <- unique(toupper(State2_genes))
State3_genes <- unique(toupper(State3_genes))

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

collapse_duplicate_genes <- function(mat) {
   mat <- as.matrix(mat)
   storage.mode(mat) <- "numeric"
   
   rn <- toupper(trimws(rownames(mat)))
   keep <- !is.na(rn) & rn != "" & rn != "---"
   
   mat <- mat[keep, , drop = FALSE]
   rn <- rn[keep]
   
   collapsed <- rowsum(mat, group = rn, reorder = FALSE, na.rm = TRUE)
   n_per_gene <- as.numeric(table(rn)[rownames(collapsed)])
   collapsed <- sweep(collapsed, 1, n_per_gene, "/")
   
   collapsed
}

rank_program_score <- function(expr, signature) {
   signature <- intersect(toupper(signature), rownames(expr))
   
   if (length(signature) < 5) {
      stop(
         "Fewer than 5 signature genes are detected. ",
         "Please check probe-to-gene-symbol mapping."
      )
   }
   
   rank_mat <- apply(expr, 2, function(z) {
      rank(z, ties.method = "average", na.last = "keep") / (sum(!is.na(z)) + 1)
   })
   
   if (is.null(dim(rank_mat))) {
      rank_mat <- matrix(rank_mat, ncol = 1)
   }
   
   rownames(rank_mat) <- rownames(expr)
   colnames(rank_mat) <- colnames(expr)
   
   score <- colMeans(rank_mat[signature, , drop = FALSE], na.rm = TRUE)
   attr(score, "detected_genes") <- signature
   
   score
}

cliffs_delta <- function(x, y) {
   x <- x[is.finite(x)]
   y <- y[is.finite(y)]
   mean(sign(outer(x, y, "-")))
}

bootstrap_cliffs_delta <- function(x, y, B = 3000, seed = 20260703) {
   set.seed(seed)
   boot <- replicate(B, {
      cliffs_delta(
         sample(x, length(x), replace = TRUE),
         sample(y, length(y), replace = TRUE)
      )
   })
   unname(stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE))
}

fmt_p <- function(p) {
   ifelse(
      is.na(p),
      "NA",
      ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = 3)))
   )
}

pub_theme <- function(base_size = 13) {
   theme_classic(base_size = base_size, base_family = "sans") +
      theme(
         axis.title = element_text(face = "bold", colour = "black"),
         axis.text = element_text(colour = "black"),
         axis.line = element_line(linewidth = 0.55, colour = "black"),
         plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
         plot.subtitle = element_text(colour = "grey30", size = base_size - 1),
         plot.margin = margin(8, 10, 8, 8),
         legend.title = element_text(face = "bold"),
         legend.background = element_blank()
      )
}

save_plot_dual <- function(plot, filename, width, height) {
   ggsave(
      file.path(outdir, paste0(filename, ".pdf")),
      plot,
      width = width,
      height = height,
      units = "in",
      device = "pdf",
      bg = "white"
   )
   
   ggsave(
      file.path(outdir, paste0(filename, ".png")),
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
   )
}

clean_colname <- function(x) {
   x <- tolower(x)
   x <- gsub("[^a-z0-9]+", "_", x)
   x <- gsub("^_|_$", "", x)
   make.unique(x)
}

extract_first_symbol <- function(x) {
   x <- as.character(x)
   x[is.na(x)] <- ""
   
   # GPL570 often uses "GENE1 /// GENE2".
   x <- sapply(strsplit(x, "///|//|;|,"), function(z) {
      z <- trimws(z)
      z <- z[!is.na(z) & z != "" & z != "---"]
      if (length(z) == 0) NA_character_ else z[1]
   })
   
   toupper(trimws(x))
}

numeric_from_text <- function(x) {
   x <- as.character(x)
   x <- gsub(",", "", x)
   suppressWarnings(as.numeric(gsub(".*?(-?[0-9]+\\.?[0-9]*).*", "\\1", x)))
}

parse_binary_event <- function(x) {
   sx <- tolower(trimws(as.character(x)))
   
   event <- rep(NA_integer_, length(sx))
   
   neg <- grepl(
      "no progression|not progress|non.?progress|without progression|no cmm|no clinical myeloma|censored|censor|stable|^0$|^no$|^false$",
      sx,
      ignore.case = TRUE
   )
   
   pos <- grepl(
      "progressed|progression|progress|clinical myeloma|cmm|requiring therapy|required therapy|symptomatic|^1$|^yes$|^true$",
      sx,
      ignore.case = TRUE
   )
   
   event[neg] <- 0L
   event[pos & !neg] <- 1L
   
   if (all(is.na(event))) {
      suppressWarnings(num <- as.numeric(sx))
      if (any(num %in% c(0, 1), na.rm = TRUE)) {
         event <- ifelse(is.na(num), NA_integer_, ifelse(num > 0, 1L, 0L))
      }
   }
   
   event
}

find_event_col <- function(pd) {
   name_hits <- grep(
      "progress|cmm|clinical_myeloma|therapy|event|status|relapse|outcome",
      colnames(pd),
      ignore.case = TRUE,
      value = TRUE
   )
   
   value_hits <- names(pd)[vapply(pd, function(x) {
      any(grepl("progress|cmm|clinical myeloma|requiring therapy|censored", as.character(x), ignore.case = TRUE), na.rm = TRUE)
   }, logical(1))]
   
   candidates <- unique(c(name_hits, value_hits))
   
   for (cc in candidates) {
      ev <- parse_binary_event(pd[[cc]])
      if (length(unique(na.omit(ev))) == 2) {
         return(cc)
      }
   }
   
   NULL
}

find_time_col <- function(pd) {
   candidates <- grep(
      "time|ttp|progression_time|time_to|follow|followup|follow_up|months|days|censor",
      colnames(pd),
      ignore.case = TRUE,
      value = TRUE
   )
   
   for (cc in candidates) {
      vv <- numeric_from_text(pd[[cc]])
      if (sum(is.finite(vv)) >= 10 && length(unique(vv[is.finite(vv)])) > 5) {
         return(cc)
      }
   }
   
   NULL
}

find_diagnosis_col <- function(pd) {
   name_hits <- grep(
      "diagnosis|disease|condition|group|classification|amm|mgus|smm",
      colnames(pd),
      ignore.case = TRUE,
      value = TRUE
   )
   
   value_hits <- names(pd)[vapply(pd, function(x) {
      any(grepl("\\bMGUS\\b|\\bAMM\\b|\\bSMM\\b|smolder", as.character(x), ignore.case = TRUE), na.rm = TRUE)
   }, logical(1))]
   
   candidates <- unique(c(name_hits, value_hits))
   
   for (cc in candidates) {
      vv <- as.character(pd[[cc]])
      if (any(grepl("\\bMGUS\\b|\\bAMM\\b|\\bSMM\\b|smolder", vv, ignore.case = TRUE), na.rm = TRUE)) {
         return(cc)
      }
   }
   
   NULL
}

standardize_diagnosis <- function(x) {
   sx <- toupper(trimws(as.character(x)))
   
   out <- rep(NA_character_, length(sx))
   out[grepl("MGUS", sx)] <- "MGUS"
   out[grepl("AMM|ASYMPTOMATIC", sx)] <- "AMM"
   out[grepl("SMM|SMOLDER", sx)] <- "SMM"
   
   out[is.na(out) & sx != "" & sx != "NA"] <- sx[is.na(out) & sx != "" & sx != "NA"]
   out
}

# -----------------------------------------------------------------------------
# 3. Load GSE122231 from SOFT or existing object
# -----------------------------------------------------------------------------

if (!exists("gse", inherits = TRUE)) {
   if (!file.exists(family_soft_path)) {
      stop("Cannot find GSE122231 family.soft file: ", family_soft_path)
   }
   
   message("Loading GSE122231 from: ", family_soft_path)
   gse <- GEOquery::getGEO(filename = family_soft_path)
}

# This function handles both:
#   1. ExpressionSet / ExpressionSet list from series matrix
#   2. GSE object from family.soft
extract_geo_expression <- function(gse_obj) {
   
   if (inherits(gse_obj, "ExpressionSet")) {
      expr <- Biobase::exprs(gse_obj)
      pd <- Biobase::pData(gse_obj)
      fd <- Biobase::fData(gse_obj)
      
      if (!"gsm_id" %in% colnames(pd)) {
         if ("geo_accession" %in% colnames(pd)) {
            pd$gsm_id <- pd$geo_accession
         } else {
            pd$gsm_id <- colnames(expr)
         }
      }
      
      return(list(expr = expr, pd = pd, fd = fd))
   }
   
   if (is.list(gse_obj) && length(gse_obj) > 0 && inherits(gse_obj[[1]], "ExpressionSet")) {
      eset <- gse_obj[[1]]
      
      expr <- Biobase::exprs(eset)
      pd <- Biobase::pData(eset)
      fd <- Biobase::fData(eset)
      
      if (!"gsm_id" %in% colnames(pd)) {
         if ("geo_accession" %in% colnames(pd)) {
            pd$gsm_id <- pd$geo_accession
         } else {
            pd$gsm_id <- colnames(expr)
         }
      }
      
      return(list(expr = expr, pd = pd, fd = fd))
   }
   
   if (inherits(gse_obj, "GSE")) {
      message("Detected GEOquery GSE object. Building expression matrix from GSM sample tables.")
      
      gsm_list <- GEOquery::GSMList(gse_obj)
      gpl_list <- GEOquery::GPLList(gse_obj)
      
      if (length(gsm_list) == 0) {
         stop("No GSM entries found in the GSE object.")
      }
      
      sample_ids <- names(gsm_list)
      
      sample_tables <- lapply(gsm_list, function(gsm) {
         tab <- GEOquery::Table(gsm)
         
         if (!all(c("ID_REF", "VALUE") %in% colnames(tab))) {
            stop("GSM table does not contain ID_REF and VALUE columns.")
         }
         
         tab <- tab[, c("ID_REF", "VALUE")]
         tab$VALUE <- suppressWarnings(as.numeric(tab$VALUE))
         tab
      })
      
      common_ids <- Reduce(intersect, lapply(sample_tables, function(z) as.character(z$ID_REF)))
      
      if (length(common_ids) < 1000) {
         stop("Too few common probe IDs across GSM sample tables.")
      }
      
      expr <- sapply(sample_tables, function(tab) {
         tab$VALUE[match(common_ids, tab$ID_REF)]
      })
      
      expr <- as.matrix(expr)
      rownames(expr) <- common_ids
      colnames(expr) <- sample_ids
      storage.mode(expr) <- "numeric"
      
      # GSM metadata
      parse_one_gsm <- function(gsm, gsm_id) {
         mm <- GEOquery::Meta(gsm)
         
         base <- list(
            gsm_id = gsm_id,
            title = paste(mm$title %||% NA_character_, collapse = "; "),
            geo_accession = paste(mm$geo_accession %||% gsm_id, collapse = "; "),
            source_name_ch1 = paste(mm$source_name_ch1 %||% NA_character_, collapse = "; ")
         )
         
         chars <- mm$characteristics_ch1
         char_df <- list()
         
         if (!is.null(chars) && length(chars) > 0) {
            for (ii in seq_along(chars)) {
               z <- chars[[ii]]
               
               if (grepl(":", z)) {
                  key <- sub(":.*$", "", z)
                  val <- sub("^[^:]+:\\s*", "", z)
               } else {
                  key <- paste0("characteristics_ch1_", ii)
                  val <- z
               }
               
               key <- clean_colname(key)
               char_df[[key]] <- val
            }
         }
         
         as.data.frame(c(base, char_df), stringsAsFactors = FALSE)
      }
      
      `%||%` <- function(a, b) {
         if (is.null(a) || length(a) == 0) b else a
      }
      
      pd <- dplyr::bind_rows(
         Map(parse_one_gsm, gsm_list, sample_ids)
      )
      
      # GPL annotation
      if (length(gpl_list) == 0) {
         fd <- data.frame(ID = rownames(expr), stringsAsFactors = FALSE)
      } else {
         fd <- GEOquery::Table(gpl_list[[1]])
      }
      
      return(list(expr = expr, pd = pd, fd = fd))
   }
   
   stop("Unsupported object returned by getGEO(). Please use family.soft or series matrix.")
}

geo_obj <- extract_geo_expression(gse)
expr_probe <- geo_obj$expr
pd <- geo_obj$pd
fd <- geo_obj$fd

pd <- as.data.frame(pd, stringsAsFactors = FALSE, check.names = FALSE)
colnames(pd) <- clean_colname(colnames(pd))

if (!"gsm_id" %in% colnames(pd)) {
   if ("geo_accession" %in% colnames(pd)) {
      pd$gsm_id <- pd$geo_accession
   } else {
      pd$gsm_id <- colnames(expr_probe)
   }
}

# Align sample order
pd$gsm_id <- as.character(pd$gsm_id)
expr_probe <- expr_probe[, pd$gsm_id, drop = FALSE]

# Save metadata dictionary for inspection
meta_dictionary <- tibble(
   column = colnames(pd),
   example_1 = vapply(pd, function(x) as.character(x[which(!is.na(x))[1]]), character(1)),
   non_missing = vapply(pd, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1)),
   n_unique = vapply(pd, function(x) length(unique(na.omit(as.character(x)))), integer(1))
)

write.csv(
   meta_dictionary,
   file.path(outdir, "Table_GSE122231_metadata_dictionary.csv"),
   row.names = FALSE
)

message("Expression matrix from GEO: ", nrow(expr_probe), " probes x ", ncol(expr_probe), " samples")

# -----------------------------------------------------------------------------
# 4. Probe-to-gene-symbol mapping
# -----------------------------------------------------------------------------

fd <- as.data.frame(fd, stringsAsFactors = FALSE, check.names = FALSE)
colnames(fd) <- clean_colname(colnames(fd))

id_col <- intersect(c("id", "id_ref", "probe_id", "probeset_id"), colnames(fd))[1]

if (is.na(id_col)) {
   if (!is.null(rownames(fd)) && all(rownames(expr_probe) %in% rownames(fd))) {
      fd$id <- rownames(fd)
      id_col <- "id"
   } else {
      stop("Cannot identify probe ID column in feature/platform annotation.")
   }
}

symbol_col <- intersect(
   c(
      "gene_symbol",
      "genesymbol",
      "symbol",
      "gene_assignment",
      "gene_title",
      "gene_name"
   ),
   colnames(fd)
)[1]

if (is.na(symbol_col)) {
   stop(
      "Cannot find gene-symbol column in platform annotation. ",
      "Please inspect fData/Table(GPL570)."
   )
}

probe_to_symbol <- fd[, c(id_col, symbol_col)]
colnames(probe_to_symbol) <- c("probe_id", "gene_symbol")
probe_to_symbol$probe_id <- as.character(probe_to_symbol$probe_id)
probe_to_symbol$gene_symbol <- extract_first_symbol(probe_to_symbol$gene_symbol)

gene_symbols <- probe_to_symbol$gene_symbol[match(rownames(expr_probe), probe_to_symbol$probe_id)]

expr_gene <- expr_probe
rownames(expr_gene) <- gene_symbols

# Auto log2 transform if expression looks unlogged
expr_range <- range(expr_gene, na.rm = TRUE)

if (is.finite(expr_range[2]) && expr_range[2] > 100) {
   message("Expression values appear unlogged. Applying log2(x + 1).")
   expr_gene <- log2(expr_gene + 1)
}

expr122231 <- collapse_duplicate_genes(expr_gene)

saveRDS(
   expr122231,
   file.path(outdir, "GSE122231_gene_symbol_expression_collapsed.rds")
)

message("Collapsed expression matrix: ", nrow(expr122231), " genes x ", ncol(expr122231), " samples")

# -----------------------------------------------------------------------------
# 5. Build sample_info for independent validation
# -----------------------------------------------------------------------------

if (!is.null(clinical_csv_path) && file.exists(clinical_csv_path)) {
   clin <- read.csv(clinical_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
   colnames(clin) <- clean_colname(colnames(clin))
   
   if ("gsm_id" %in% colnames(clin)) {
      pd <- pd %>% left_join(clin, by = "gsm_id")
   } else if ("patient_id" %in% colnames(clin)) {
      # Patient IDs are assigned below if needed; for now keep clinical columns.
      pd <- pd %>% left_join(clin, by = "patient_id")
   } else {
      warning("clinical_csv_path was provided but contains neither gsm_id nor patient_id. Skipping merge.")
   }
}

patient_col <- manual_patient_col

if (is.null(patient_col)) {
   patient_col <- intersect(c("patient_id", "subject_id", "title", "description"), colnames(pd))[1]
}

if (is.na(patient_col) || is.null(patient_col)) {
   patient_id <- pd$gsm_id
} else {
   patient_id <- as.character(pd[[patient_col]])
   patient_id[is.na(patient_id) | patient_id == ""] <- pd$gsm_id[is.na(patient_id) | patient_id == ""]
}

diagnosis_col <- manual_diagnosis_col
if (is.null(diagnosis_col)) {
   diagnosis_col <- find_diagnosis_col(pd)
}

event_col <- manual_event_col
if (is.null(event_col)) {
   event_col <- find_event_col(pd)
}

time_col <- manual_time_col
if (is.null(time_col)) {
   time_col <- find_time_col(pd)
}

diagnosis <- if (!is.null(diagnosis_col)) standardize_diagnosis(pd[[diagnosis_col]]) else NA_character_
event <- if (!is.null(event_col)) parse_binary_event(pd[[event_col]]) else rep(NA_integer_, nrow(pd))
time_to_progression <- if (!is.null(time_col)) numeric_from_text(pd[[time_col]]) else rep(NA_real_, nrow(pd))

sample_info <- tibble(
   gsm_id = pd$gsm_id,
   patient_id = patient_id,
   diagnosis = diagnosis,
   event = event,
   time_to_progression = time_to_progression
)

sample_info <- sample_info %>%
   mutate(
      progression_group = case_when(
         event == 1 ~ "Progressed",
         event == 0 ~ "Non-progressed",
         TRUE ~ NA_character_
      )
   )

write.csv(
   sample_info,
   file.path(outdir, "GSE122231_metadata_used.csv"),
   row.names = FALSE
)

message("Auto-detected columns:")
message("  patient_col   = ", ifelse(is.null(patient_col), "NULL", patient_col))
message("  diagnosis_col = ", ifelse(is.null(diagnosis_col), "NULL", diagnosis_col))
message("  event_col     = ", ifelse(is.null(event_col), "NULL", event_col))
message("  time_col      = ", ifelse(is.null(time_col), "NULL", time_col))

print(table(sample_info$diagnosis, useNA = "ifany"))
print(table(sample_info$progression_group, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 6. Score State2 and State3
# -----------------------------------------------------------------------------

expr122231 <- expr122231[, sample_info$gsm_id, drop = FALSE]

state2_score <- rank_program_score(expr122231, State2_genes)
state3_score <- rank_program_score(expr122231, State3_genes)

coverage <- bind_rows(
   tibble(
      program = "State 2",
      genes_requested = length(State2_genes),
      genes_detected = length(attr(state2_score, "detected_genes")),
      genes_missing = paste(setdiff(State2_genes, attr(state2_score, "detected_genes")), collapse = "; ")
   ),
   tibble(
      program = "State 3",
      genes_requested = length(State3_genes),
      genes_detected = length(attr(state3_score, "detected_genes")),
      genes_missing = paste(setdiff(State3_genes, attr(state3_score, "detected_genes")), collapse = "; ")
   )
)

write.csv(
   coverage,
   file.path(outdir, "Table_GSE122231_gene_coverage.csv"),
   row.names = FALSE
)

score_df <- sample_info %>%
   mutate(
      State2_rank = unname(state2_score[gsm_id]),
      State3_rank = unname(state3_score[gsm_id]),
      State2_z = as.numeric(scale(State2_rank)),
      State3_z = as.numeric(scale(State3_rank)),
      State3_minus_State2 = State3_rank - State2_rank,
      State2_plus_State3 = State2_rank + State3_rank
   )

write.csv(
   score_df,
   file.path(outdir, "Table_GSE122231_patient_level_program_scores.csv"),
   row.names = FALSE
)

print(coverage)

##
# -----------------------------------------------------------------------------
# 7. GSE122231 validation available from SOFT metadata:
#    MGUS versus AMM cross-sectional disease-stage validation
# -----------------------------------------------------------------------------

# 当前 SOFT 文件没有 progression event 和 time-to-progression。
# 所以主分析先做 MGUS vs AMM。
# MM 和 WM 样本数太少，排除在主统计外。

stage_df <- score_df %>%
   filter(diagnosis %in% c("MGUS", "AMM")) %>%
   mutate(
      diagnosis = factor(diagnosis, levels = c("MGUS", "AMM"))
   )

message("Samples used for MGUS vs AMM validation:")
print(table(stage_df$diagnosis, useNA = "ifany"))

if (nrow(stage_df) < 20 || length(unique(stage_df$diagnosis)) < 2) {
   stop("Too few MGUS/AMM samples for cross-sectional validation.")
}

run_stage_wilcox <- function(score_name) {
   x_mgus <- stage_df %>%
      filter(diagnosis == "MGUS") %>%
      pull(.data[[score_name]])
   
   x_amm <- stage_df %>%
      filter(diagnosis == "AMM") %>%
      pull(.data[[score_name]])
   
   wt <- wilcox.test(
      x_amm,
      x_mgus,
      alternative = "two.sided",
      exact = FALSE
   )
   
   cd <- cliffs_delta(x_amm, x_mgus)
   cd_ci <- bootstrap_cliffs_delta(x_amm, x_mgus)
   
   tibble(
      analysis = "GSE122231 MGUS versus AMM",
      program = gsub("_rank", "", score_name),
      n_MGUS = length(x_mgus),
      n_AMM = length(x_amm),
      median_MGUS = median(x_mgus, na.rm = TRUE),
      median_AMM = median(x_amm, na.rm = TRUE),
      delta_median_AMM_minus_MGUS = median(x_amm, na.rm = TRUE) - median(x_mgus, na.rm = TRUE),
      statistic_W = unname(wt$statistic),
      p_value = wt$p.value,
      effect = "Cliff's delta, AMM minus MGUS",
      effect_estimate = cd,
      ci_low = cd_ci[1],
      ci_high = cd_ci[2]
   )
}

stage_stats <- bind_rows(
   run_stage_wilcox("State2_rank"),
   run_stage_wilcox("State3_rank"),
   run_stage_wilcox("State3_minus_State2"),
   run_stage_wilcox("State2_plus_State3")
) %>%
   mutate(p_BH = p.adjust(p_value, method = "BH"))

write.csv(
   stage_stats,
   file.path(outdir, "Table_GSE122231_MGUS_vs_AMM_statistics.csv"),
   row.names = FALSE
)

print(stage_stats)

# -----------------------------------------------------------------------------
# 8. Exploratory ROC: AMM versus MGUS
# -----------------------------------------------------------------------------

make_stage_roc <- function(score_name, label) {
   rr <- pROC::roc(
      response = factor(stage_df$diagnosis, levels = c("MGUS", "AMM")),
      predictor = stage_df[[score_name]],
      levels = c("MGUS", "AMM"),
      direction = "<",
      ci = TRUE,
      quiet = TRUE
   )
   
   tibble(
      analysis = "AMM versus MGUS",
      program = label,
      AUC = as.numeric(pROC::auc(rr)),
      CI_low = as.numeric(pROC::ci.auc(rr)[1]),
      CI_high = as.numeric(pROC::ci.auc(rr)[3])
   )
}

stage_roc_tbl <- bind_rows(
   make_stage_roc("State2_rank", "State 2"),
   make_stage_roc("State3_rank", "State 3"),
   make_stage_roc("State3_minus_State2", "State 3 minus State 2"),
   make_stage_roc("State2_plus_State3", "State 2 plus State 3")
)

write.csv(
   stage_roc_tbl,
   file.path(outdir, "Table_GSE122231_MGUS_vs_AMM_ROC.csv"),
   row.names = FALSE
)

print(stage_roc_tbl)

# -----------------------------------------------------------------------------
# 9.1 Boxplot / violin plot: MGUS versus AMM
# -----------------------------------------------------------------------------

stage_cols <- c("MGUS" = "#4C78A8", "AMM" = "#E45756")

plot_stage_box <- function(score_var, title) {
   plot_dat <- stage_df %>%
      filter(is.finite(.data[[score_var]]))
   
   y_max <- max(plot_dat[[score_var]], na.rm = TRUE)
   y_min <- min(plot_dat[[score_var]], na.rm = TRUE)
   y_pad <- max((y_max - y_min) * 0.18, 0.025)
   
   wt <- wilcox.test(
      plot_dat[[score_var]][plot_dat$diagnosis == "AMM"],
      plot_dat[[score_var]][plot_dat$diagnosis == "MGUS"],
      exact = FALSE
   )
   
   ggplot(
      plot_dat,
      aes(
         x = diagnosis,
         y = .data[[score_var]],
         fill = diagnosis,
         colour = diagnosis
      )
   ) +
      geom_violin(trim = FALSE, alpha = 0.22, linewidth = 0.75, width = 0.88) +
      geom_boxplot(width = 0.27, alpha = 0.42, outlier.shape = NA, linewidth = 0.65) +
      geom_jitter(width = 0.075, size = 2.2, alpha = 0.9, show.legend = FALSE) +
      annotate(
         "text",
         x = 1.5,
         y = y_max + y_pad,
         label = fmt_p(wt$p.value),
         fontface = "bold",
         size = 4
      ) +
      scale_fill_manual(values = stage_cols) +
      scale_colour_manual(values = stage_cols) +
      coord_cartesian(
         ylim = c(y_min - y_pad * 0.1, y_max + y_pad * 1.4),
         clip = "off"
      ) +
      labs(
         title = title,
         x = NULL,
         y = "Rank-based program score"
      ) +
      pub_theme() +
      theme(legend.position = "none")
}

p_stage2 <- plot_stage_box(
   "State2_rank",
   "GSE122231: State 2 score in MGUS versus AMM"
)

p_stage3 <- plot_stage_box(
   "State3_rank",
   "GSE122231: State 3 score in MGUS versus AMM"
)

p_stage_box <- p_stage2 + p_stage3

save_plot_dual(
   p_stage_box,
   "Figure_GSE122231_State2_State3_MGUS_vs_AMM_boxplots",
   width = 11,
   height = 5
)



##ROC曲线
# -----------------------------------------------------------------------------
# 9.2 ROC plots: AMM versus MGUS
# -----------------------------------------------------------------------------

roc_stage2 <- pROC::roc(
   response = factor(stage_df$diagnosis, levels = c("MGUS", "AMM")),
   predictor = stage_df$State2_rank,
   levels = c("MGUS", "AMM"),
   direction = "<",
   ci = TRUE,
   quiet = TRUE
)

roc_stage3 <- pROC::roc(
   response = factor(stage_df$diagnosis, levels = c("MGUS", "AMM")),
   predictor = stage_df$State3_rank,
   levels = c("MGUS", "AMM"),
   direction = "<",
   ci = TRUE,
   quiet = TRUE
)

plot_roc <- function(roc_obj, title) {
   auc_val <- as.numeric(pROC::auc(roc_obj))
   ci_val <- as.numeric(pROC::ci.auc(roc_obj))
   
   roc_dat <- tibble(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
   )
   
   ggplot(roc_dat, aes(x = fpr, y = tpr)) +
      geom_abline(linetype = "dashed", colour = "grey70", linewidth = 0.7) +
      geom_path(linewidth = 1.15) +
      annotate(
         "text",
         x = 0.55,
         y = 0.15,
         hjust = 0,
         fontface = "bold",
         size = 4,
         label = paste0(
            "AUC = ", sprintf("%.3f", auc_val),
            "\n95% CI ",
            sprintf("%.3f", ci_val[1]),
            "–",
            sprintf("%.3f", ci_val[3])
         )
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
         title = title,
         x = "False-positive rate",
         y = "True-positive rate"
      ) +
      pub_theme()
}

p_stage_roc <- plot_roc(roc_stage2, "State 2: AMM versus MGUS") +
   plot_roc(roc_stage3, "State 3: AMM versus MGUS")

save_plot_dual(
   p_stage_roc,
   "Figure_GSE122231_State2_State3_MGUS_vs_AMM_ROC",
   width = 11,
   height = 5
)




##
# -----------------------------------------------------------------------------
# 9.3 State-space plot
# -----------------------------------------------------------------------------

p_stage_space <- score_df %>%
   mutate(
      diagnosis = factor(diagnosis, levels = c("MGUS", "AMM", "MM", "WM"))
   ) %>%
   ggplot(aes(x = State2_z, y = State3_z)) +
   geom_point(
      aes(colour = diagnosis),
      size = 2.7,
      alpha = 0.9
   ) +
   labs(
      title = "GSE122231 State-space validation",
      subtitle = "Baseline samples positioned by State 2 and State 3 program activity",
      x = "State 2 score, z-scaled",
      y = "State 3 score, z-scaled",
      colour = "Diagnosis"
   ) +
   pub_theme()

save_plot_dual(
   p_stage_space,
   "Figure_GSE122231_State2_State3_state_space_by_diagnosis",
   width = 7,
   height = 6
)

# -----------------------------------------------------------------------------
# 9.4 Program score heatmap
# -----------------------------------------------------------------------------

score_order <- score_df %>%
   mutate(
      diagnosis = factor(diagnosis, levels = c("MGUS", "AMM", "MM", "WM"))
   ) %>%
   arrange(diagnosis, State2_rank, State3_rank) %>%
   pull(gsm_id)

program_mat <- rbind(
   `State 2` = score_df$State2_rank[match(score_order, score_df$gsm_id)],
   `State 3` = score_df$State3_rank[match(score_order, score_df$gsm_id)]
)

program_z <- t(scale(t(program_mat)))
colnames(program_z) <- score_order

ann_col <- score_df %>%
   select(gsm_id, diagnosis) %>%
   arrange(match(gsm_id, score_order)) %>%
   as.data.frame()

rownames(ann_col) <- ann_col$gsm_id
ann_col$gsm_id <- NULL

pdf(
   file.path(outdir, "Figure_GSE122231_State2_State3_program_score_heatmap.pdf"),
   width = 12,
   height = 3.2,
   useDingbats = FALSE
)

pheatmap::pheatmap(
   program_z,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   annotation_col = ann_col,
   color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
   breaks = seq(-2.5, 2.5, length.out = 102),
   border_color = NA,
   fontsize_row = 12,
   fontsize_col = 6,
   main = "State 2 and State 3 activity across GSE122231 baseline samples"
)

dev.off()


# # -----------------------------------------------------------------------------
# # 7. Cross-sectional progression comparison
# # -----------------------------------------------------------------------------
# 
# analysis_df <- score_df %>%
#    filter(!is.na(event)) %>%
#    mutate(
#       progression_group = factor(
#          progression_group,
#          levels = c("Non-progressed", "Progressed")
#       )
#    )
# 
# if (nrow(analysis_df) >= 10 && length(unique(analysis_df$event)) == 2) {
#    
#    run_progression_wilcox <- function(score_name) {
#       x0 <- analysis_df %>%
#          filter(event == 0) %>%
#          pull(.data[[score_name]])
#       
#       x1 <- analysis_df %>%
#          filter(event == 1) %>%
#          pull(.data[[score_name]])
#       
#       wt <- wilcox.test(x1, x0, alternative = "two.sided", exact = FALSE)
#       cd <- cliffs_delta(x1, x0)
#       cd_ci <- bootstrap_cliffs_delta(x1, x0)
#       
#       tibble(
#          analysis = "Baseline progressed versus non-progressed",
#          program = gsub("_rank", "", score_name),
#          n_non_progressed = length(x0),
#          n_progressed = length(x1),
#          statistic_W = unname(wt$statistic),
#          p_value = wt$p.value,
#          effect = "Cliff's delta, progressed minus non-progressed",
#          effect_estimate = cd,
#          ci_low = cd_ci[1],
#          ci_high = cd_ci[2]
#       )
#    }
#    
#    wilcox_tbl <- bind_rows(
#       run_progression_wilcox("State2_rank"),
#       run_progression_wilcox("State3_rank"),
#       run_progression_wilcox("State3_minus_State2"),
#       run_progression_wilcox("State2_plus_State3")
#    ) %>%
#       mutate(p_BH = p.adjust(p_value, method = "BH"))
#    
#    write.csv(
#       wilcox_tbl,
#       file.path(outdir, "Table_GSE122231_progression_group_statistics.csv"),
#       row.names = FALSE
#    )
#    
#    print(wilcox_tbl)
#    
#    # ROC for ever-progressed status
#    make_roc <- function(score_name, label) {
#       rr <- pROC::roc(
#          response = factor(analysis_df$event, levels = c(0, 1)),
#          predictor = analysis_df[[score_name]],
#          levels = c("0", "1"),
#          direction = "<",
#          ci = TRUE,
#          quiet = TRUE
#       )
#       
#       tibble(
#          program = label,
#          AUC = as.numeric(pROC::auc(rr)),
#          CI_low = as.numeric(pROC::ci.auc(rr)[1]),
#          CI_high = as.numeric(pROC::ci.auc(rr)[3])
#       )
#    }
#    
#    roc_tbl <- bind_rows(
#       make_roc("State2_rank", "State 2"),
#       make_roc("State3_rank", "State 3"),
#       make_roc("State3_minus_State2", "State 3 minus State 2"),
#       make_roc("State2_plus_State3", "State 2 plus State 3")
#    )
#    
#    write.csv(
#       roc_tbl,
#       file.path(outdir, "Table_GSE122231_progression_ROC.csv"),
#       row.names = FALSE
#    )
#    
#    print(roc_tbl)
#    
# } else {
#    warning(
#       "No usable binary progression event detected. ",
#       "Scores were computed, but progression statistics were skipped. ",
#       "Inspect Table_GSE122231_metadata_dictionary.csv and set manual_event_col."
#    )
# }
# 
# # -----------------------------------------------------------------------------
# # 8. Cox regression: time-to-progression validation
# # -----------------------------------------------------------------------------
# 
# surv_df <- score_df %>%
#    filter(
#       !is.na(event),
#       is.finite(time_to_progression),
#       time_to_progression > 0
#    )
# 
# if (nrow(surv_df) >= 10 && length(unique(surv_df$event)) == 2) {
#    
#    run_cox <- function(formula_text, model_name) {
#       fit <- survival::coxph(as.formula(formula_text), data = surv_df)
#       ss <- summary(fit)
#       
#       out <- as.data.frame(ss$coefficients)
#       out$term <- rownames(out)
#       
#       ci <- as.data.frame(ss$conf.int)
#       ci$term <- rownames(ci)
#       
#       out <- out %>%
#          left_join(ci, by = "term") %>%
#          transmute(
#             model = model_name,
#             term = term,
#             beta = coef,
#             HR = `exp(coef)`,
#             HR_low = `lower .95`,
#             HR_high = `upper .95`,
#             p_value = `Pr(>|z|)`,
#             concordance = ss$concordance[1]
#          )
#       
#       out
#    }
#    
#    cox_tbl <- bind_rows(
#       run_cox("Surv(time_to_progression, event) ~ State2_z", "State2_z only"),
#       run_cox("Surv(time_to_progression, event) ~ State3_z", "State3_z only"),
#       run_cox("Surv(time_to_progression, event) ~ State2_z + State3_z", "State2_z + State3_z"),
#       run_cox("Surv(time_to_progression, event) ~ State3_minus_State2", "State3_minus_State2"),
#       run_cox("Surv(time_to_progression, event) ~ State2_plus_State3", "State2_plus_State3")
#    ) %>%
#       mutate(p_BH = p.adjust(p_value, method = "BH"))
#    
#    # Optional diagnosis-adjusted model if diagnosis is available
#    if (sum(!is.na(surv_df$diagnosis)) >= 10 && length(unique(na.omit(surv_df$diagnosis))) > 1) {
#       surv_df$diagnosis <- factor(surv_df$diagnosis)
#       
#       cox_adj_tbl <- bind_rows(
#          run_cox("Surv(time_to_progression, event) ~ State2_z + diagnosis", "State2_z + diagnosis"),
#          run_cox("Surv(time_to_progression, event) ~ State3_z + diagnosis", "State3_z + diagnosis"),
#          run_cox("Surv(time_to_progression, event) ~ State2_z + State3_z + diagnosis", "State2_z + State3_z + diagnosis")
#       )
#       
#       cox_tbl <- bind_rows(cox_tbl, cox_adj_tbl)
#    }
#    
#    write.csv(
#       cox_tbl,
#       file.path(outdir, "Table_GSE122231_Cox_time_to_progression.csv"),
#       row.names = FALSE
#    )
#    
#    print(cox_tbl)
#    
# } else {
#    warning(
#       "No usable time-to-progression column detected. ",
#       "Cox models were skipped. ",
#       "Inspect Table_GSE122231_metadata_dictionary.csv and set manual_time_col."
#    )
# }
# 
# # -----------------------------------------------------------------------------
# # 9. Figures
# # -----------------------------------------------------------------------------
# 
# # 9.1 Progressed vs non-progressed box/violin plots
# if (exists("analysis_df") && nrow(analysis_df) >= 10 && length(unique(analysis_df$event)) == 2) {
#    
#    group_cols <- c("Non-progressed" = "#4C78A8", "Progressed" = "#E45756")
#    
#    plot_progression_box <- function(score_var, title) {
#       plot_dat <- analysis_df %>%
#          filter(is.finite(.data[[score_var]]))
#       
#       y_max <- max(plot_dat[[score_var]], na.rm = TRUE)
#       y_min <- min(plot_dat[[score_var]], na.rm = TRUE)
#       y_pad <- max((y_max - y_min) * 0.18, 0.025)
#       
#       wt <- wilcox.test(
#          plot_dat[[score_var]][plot_dat$event == 1],
#          plot_dat[[score_var]][plot_dat$event == 0],
#          exact = FALSE
#       )
#       
#       ggplot(
#          plot_dat,
#          aes(x = progression_group, y = .data[[score_var]],
#              fill = progression_group, colour = progression_group)
#       ) +
#          geom_violin(trim = FALSE, alpha = 0.22, linewidth = 0.75, width = 0.88) +
#          geom_boxplot(width = 0.27, alpha = 0.42, outlier.shape = NA, linewidth = 0.65) +
#          geom_jitter(width = 0.075, size = 2.2, alpha = 0.9, show.legend = FALSE) +
#          annotate("text", x = 1.5, y = y_max + y_pad, label = fmt_p(wt$p.value), fontface = "bold", size = 4) +
#          scale_fill_manual(values = group_cols) +
#          scale_colour_manual(values = group_cols) +
#          coord_cartesian(ylim = c(y_min - y_pad * 0.1, y_max + y_pad * 1.4), clip = "off") +
#          labs(title = title, x = NULL, y = "Rank-based program score") +
#          pub_theme() +
#          theme(legend.position = "none")
#    }
#    
#    p_box2 <- plot_progression_box("State2_rank", "State 2 baseline score by later progression")
#    p_box3 <- plot_progression_box("State3_rank", "State 3 baseline score by later progression")
#    
#    p_box <- p_box2 + p_box3
#    
#    save_plot_dual(
#       p_box,
#       "Figure_GSE122231_State2_State3_progression_boxplots",
#       width = 11,
#       height = 5
#    )
#    
#    # ROC figure
#    roc_state2 <- pROC::roc(
#       response = factor(analysis_df$event, levels = c(0, 1)),
#       predictor = analysis_df$State2_rank,
#       levels = c("0", "1"),
#       direction = "<",
#       ci = TRUE,
#       quiet = TRUE
#    )
#    
#    roc_state3 <- pROC::roc(
#       response = factor(analysis_df$event, levels = c(0, 1)),
#       predictor = analysis_df$State3_rank,
#       levels = c("0", "1"),
#       direction = "<",
#       ci = TRUE,
#       quiet = TRUE
#    )
#    
#    plot_roc <- function(roc_obj, title) {
#       auc_val <- as.numeric(pROC::auc(roc_obj))
#       ci_val <- as.numeric(pROC::ci.auc(roc_obj))
#       
#       roc_dat <- tibble(
#          fpr = 1 - roc_obj$specificities,
#          tpr = roc_obj$sensitivities
#       )
#       
#       ggplot(roc_dat, aes(x = fpr, y = tpr)) +
#          geom_abline(linetype = "dashed", colour = "grey70", linewidth = 0.7) +
#          geom_path(linewidth = 1.15) +
#          annotate(
#             "text",
#             x = 0.55,
#             y = 0.15,
#             hjust = 0,
#             fontface = "bold",
#             size = 4,
#             label = paste0(
#                "AUC = ", sprintf("%.3f", auc_val),
#                "\n95% CI ", sprintf("%.3f", ci_val[1]), "–", sprintf("%.3f", ci_val[3])
#             )
#          ) +
#          coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
#          scale_x_continuous(labels = percent_format(accuracy = 1)) +
#          scale_y_continuous(labels = percent_format(accuracy = 1)) +
#          labs(title = title, x = "False-positive rate", y = "True-positive rate") +
#          pub_theme()
#    }
#    
#    p_roc <- plot_roc(roc_state2, "State 2: progression ROC") +
#       plot_roc(roc_state3, "State 3: progression ROC")
#    
#    save_plot_dual(
#       p_roc,
#       "Figure_GSE122231_State2_State3_progression_ROC",
#       width = 11,
#       height = 5
#    )
# }
# 
# # 9.2 KM curves using median split
# survfit_to_df <- function(fit) {
#    ss <- summary(fit)
#    
#    tibble(
#       time = ss$time,
#       surv = ss$surv,
#       lower = ss$lower,
#       upper = ss$upper,
#       strata = if (!is.null(ss$strata)) names(ss$strata) else "All"
#    )
# }
# 
# plot_km_median <- function(score_var, title) {
#    dat <- surv_df %>%
#       filter(is.finite(.data[[score_var]])) %>%
#       mutate(
#          score_group = ifelse(
#             .data[[score_var]] >= median(.data[[score_var]], na.rm = TRUE),
#             "High",
#             "Low"
#          ),
#          score_group = factor(score_group, levels = c("Low", "High"))
#       )
#    
#    fit <- survival::survfit(Surv(time_to_progression, event) ~ score_group, data = dat)
#    lr <- survival::survdiff(Surv(time_to_progression, event) ~ score_group, data = dat)
#    p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
#    
#    km_df <- survfit_to_df(fit) %>%
#       mutate(
#          score_group = sub("^score_group=", "", strata)
#       )
#    
#    ggplot(km_df, aes(x = time, y = surv, colour = score_group)) +
#       geom_step(linewidth = 1.05) +
#       coord_cartesian(ylim = c(0, 1), expand = FALSE) +
#       scale_y_continuous(labels = percent_format(accuracy = 1)) +
#       labs(
#          title = title,
#          subtitle = paste0("Median split; log-rank ", fmt_p(p)),
#          x = "Time to progression",
#          y = "Progression-free probability",
#          colour = NULL
#       ) +
#       pub_theme()
# }
# 
# if (exists("surv_df") && nrow(surv_df) >= 10 && length(unique(surv_df$event)) == 2) {
#    p_km2 <- plot_km_median("State2_rank", "State 2 and time to progression")
#    p_km3 <- plot_km_median("State3_rank", "State 3 and time to progression")
#    
#    p_km <- p_km2 + p_km3
#    
#    save_plot_dual(
#       p_km,
#       "Figure_GSE122231_State2_State3_KM_median_split",
#       width = 11,
#       height = 5
#    )
# }
# 
# # 9.3 State-space plot
# p_space <- score_df %>%
#    ggplot(aes(x = State2_z, y = State3_z)) +
#    geom_point(
#       aes(colour = progression_group, shape = diagnosis),
#       size = 2.7,
#       alpha = 0.9
#    ) +
#    labs(
#       title = "GSE122231 State-space validation",
#       subtitle = "Baseline samples positioned by State 2 and State 3 program activity",
#       x = "State 2 score, z-scaled",
#       y = "State 3 score, z-scaled",
#       colour = "Progression",
#       shape = "Diagnosis"
#    ) +
#    pub_theme()
# 
# save_plot_dual(
#    p_space,
#    "Figure_GSE122231_State2_State3_state_space",
#    width = 7,
#    height = 6
# )
# 
# # 9.4 Program score heatmap
# score_order <- score_df %>%
#    arrange(event, diagnosis, time_to_progression) %>%
#    pull(gsm_id)
# 
# program_mat <- rbind(
#    `State 2` = score_df$State2_rank[match(score_order, score_df$gsm_id)],
#    `State 3` = score_df$State3_rank[match(score_order, score_df$gsm_id)]
# )
# 
# program_z <- t(scale(t(program_mat)))
# colnames(program_z) <- score_order
# 
# ann_col <- score_df %>%
#    select(gsm_id, diagnosis, progression_group) %>%
#    arrange(match(gsm_id, score_order)) %>%
#    as.data.frame()
# 
# rownames(ann_col) <- ann_col$gsm_id
# ann_col$gsm_id <- NULL
# 
# pdf(
#    file.path(outdir, "Figure_GSE122231_State2_State3_program_score_heatmap.pdf"),
#    width = 12,
#    height = 3.2,
#    useDingbats = FALSE
# )
# 
# pheatmap::pheatmap(
#    program_z,
#    cluster_rows = FALSE,
#    cluster_cols = FALSE,
#    annotation_col = ann_col,
#    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
#    breaks = seq(-2.5, 2.5, length.out = 102),
#    border_color = NA,
#    fontsize_row = 12,
#    fontsize_col = 6,
#    main = "State 2 and State 3 activity across GSE122231 baseline samples"
# )
# 
# dev.off()
# 
# # -----------------------------------------------------------------------------
# # 10. Optional marker heatmaps
# # -----------------------------------------------------------------------------
# 
# make_marker_heatmap <- function(signature, filename, title, width = 10, height = 8) {
#    genes <- intersect(toupper(signature), rownames(expr122231))
#    samples <- score_order
#    
#    if (length(genes) < 2 || length(samples) < 2) {
#       warning("Too few genes or samples for heatmap: ", filename)
#       return(invisible(NULL))
#    }
#    
#    mat <- expr122231[genes, samples, drop = FALSE]
#    
#    keep_gene <- apply(mat, 1, function(x) {
#       all(is.finite(x)) && sd(x, na.rm = TRUE) > 0
#    })
#    
#    mat <- mat[keep_gene, , drop = FALSE]
#    
#    if (nrow(mat) < 2) {
#       warning("Too few variable genes after filtering: ", filename)
#       return(invisible(NULL))
#    }
#    
#    mat_z <- t(scale(t(mat)))
#    mat_z[mat_z > 2] <- 2
#    mat_z[mat_z < -2] <- -2
#    
#    ann <- score_df %>%
#       filter(gsm_id %in% samples) %>%
#       arrange(match(gsm_id, samples)) %>%
#       select(gsm_id, diagnosis, progression_group) %>%
#       as.data.frame()
#    
#    rownames(ann) <- ann$gsm_id
#    ann$gsm_id <- NULL
#    ann <- ann[samples, , drop = FALSE]
#    
#    out_file <- file.path(outdir, filename)
#    
#    p <- pheatmap::pheatmap(
#       mat_z,
#       cluster_rows = TRUE,
#       cluster_cols = FALSE,
#       annotation_col = ann,
#       color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
#       border_color = NA,
#       fontsize_row = 6.5,
#       fontsize_col = 5.5,
#       main = title,
#       silent = TRUE
#    )
#    
#    pdf(out_file, width = width, height = height, useDingbats = FALSE)
#    grid::grid.newpage()
#    grid::grid.draw(p$gtable)
#    dev.off()
#    
#    message("Saved: ", out_file)
# }
# 
# make_marker_heatmap(
#    State2_genes,
#    "FigureS_GSE122231_State2_marker_heatmap.pdf",
#    "State 2 marker expression in GSE122231",
#    width = 12,
#    height = 8
# )
# 
# make_marker_heatmap(
#    State3_genes,
#    "FigureS_GSE122231_State3_marker_heatmap.pdf",
#    "State 3 marker expression in GSE122231",
#    width = 12,
#    height = 8
# )
# 
# # -----------------------------------------------------------------------------
# # 11. Session info
# # -----------------------------------------------------------------------------
# 
# writeLines(
#    capture.output(sessionInfo()),
#    file.path(outdir, "sessionInfo_GSE122231_State2State3.txt")
# )
# 
# message("Done. Results saved to: ", outdir)
# 

