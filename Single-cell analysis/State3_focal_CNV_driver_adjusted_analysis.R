#!/usr/bin/env Rscript

## =============================================================================
## Stand-alone State3 focal CNV landscape and CNV-adjusted signature analysis
##
## Analyses
##   1. Reconstruct cnv_delta and gene_info automatically when this script is
##      started in a fresh R session.
##   2. Compare State3 versus State2 for:
##        - gain-like 1q21 signal (CKS1B/MCL1 region);
##        - gain-like 8q24 signal (MYC region);
##        - loss-like 17p13 signal (TP53 region).
##   3. Test sample-paired focal-event differences.
##   4. Model State3 signature and State3 identity after adjustment for global
##      inferred CNV burden, cell cycle, library size, and sample identity.
##   5. Generate publication-style PDFs using grDevices::pdf(), not cairo_pdf.
##
## Interpretation
##   inferCNV produces RNA-derived relative copy-number signals. Terms such as
##   "gain-like" and "loss-like" are intentionally used here. DNA/FISH/WES/WGS
##   validation is required before calling a genomic amplification or deletion.
## =============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260713)

## =============================================================================
## 0. Packages
## =============================================================================

required_pkgs <- c(
  "infercnv",
  "matrixStats",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall CRAN packages with install.packages() and Bioconductor packages ",
    "with BiocManager::install()."
  )
}

suppressPackageStartupMessages({
  library(infercnv)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

## =============================================================================
## 1. Paths and analysis parameters
## =============================================================================

## Keep this path consistent with State3_inferCNV_malignancy_analysis.R.
out_dir <- paste0(
  "/home/yjliu/mmProj/data_process/Human/",
  "SingleCell_NMF/NMF_Result/State3_inferCNV"
)

infercnv_rds_file <- file.path(
  out_dir,
  "State3_inferCNV_final_object.rds"
)

metric_csv_file <- file.path(
  out_dir,
  "Table_cell_level_inferCNV_metrics.csv"
)

analysis_out_dir <- file.path(
  out_dir,
  "State3_focal_CNV_analysis"
)

dir.create(
  analysis_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

## Minimum numbers used for stable summaries.
minimum_genes_per_region <- 3L
minimum_cells_per_sample_state <- 5L
minimum_paired_samples <- 3L
bootstrap_iterations <- 2000L

state_levels <- c("State1", "State2", "State3")

state_colors <- c(
  "State1" = "#79A96B",
  "State2" = "#1D9BC2",
  "State3" = "#C83E3A"
)

## Broad windows are used because inferCNV is smoothed across neighboring genes.
## Coordinates are GRCh38.
driver_regions <- tibble::tribble(
  ~event_id,           ~event_label,          ~chr,    ~start,    ~end,      ~direction, ~direction_label, ~anchor_genes,
  "gain_1q21",         "1q21 gain-like",       "chr1",  145e6,     160e6,      1,          "Gain-like",       "CKS1B;MCL1",
  "gain_MYC_8q24",     "8q24 MYC gain-like",   "chr8",  125e6,     132e6,      1,          "Gain-like",       "MYC",
  "loss_TP53_17p13",   "17p13 TP53 loss-like", "chr17", 5e6,       10e6,      -1,          "Loss-like",       "TP53"
)

driver_regions$event_label <- factor(
  driver_regions$event_label,
  levels = driver_regions$event_label
)

## =============================================================================
## 2. Helper functions
## =============================================================================

message2 <- function(...) {
  message(paste0(...))
}

save_pdf <- function(plot_object, filename, width, height) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    onefile = TRUE,
    useDingbats = FALSE,
    paper = "special"
  )
  print(plot_object)
  grDevices::dev.off()
  invisible(filename)
}

theme_publication <- function(base_size = 13) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        colour = "#20262D"
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        colour = "#20262D"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 2
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        colour = "#58616A"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = "#20262D"
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F3F5F6",
        colour = NA
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.key = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

format_p <- function(p_value) {
  if (!is.finite(p_value)) {
    return("P = NA")
  }
  if (p_value < 1e-4) {
    return("P < 1×10\u207b\u2074")
  }
  paste0(
    "P = ",
    formatC(
      p_value,
      format = "g",
      digits = 3
    )
  )
}

safe_zscore <- function(x) {
  x <- as.numeric(x)
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(x_sd) || x_sd == 0) {
    return(rep(0, length(x)))
  }

  (x - x_mean) / x_sd
}

safe_paired_wilcox <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2L) {
    return(
      list(
        p_value = NA_real_,
        estimate = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  fit <- tryCatch(
    suppressWarnings(
      stats::wilcox.test(
        x = x,
        y = y,
        paired = TRUE,
        exact = FALSE,
        conf.int = TRUE,
        conf.level = 0.95
      )
    ),
    error = function(e) NULL
  )

  estimate <- if (!is.null(fit) && length(fit$estimate) > 0L) {
    unname(fit$estimate)
  } else {
    stats::median(x - y, na.rm = TRUE)
  }

  conf_int <- if (!is.null(fit) && length(fit$conf.int) == 2L) {
    as.numeric(fit$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }

  list(
    p_value = if (!is.null(fit)) fit$p.value else NA_real_,
    estimate = estimate,
    conf_low = conf_int[1],
    conf_high = conf_int[2]
  )
}

bootstrap_median_difference_ci <- function(
    difference,
    iterations = 2000L,
    seed = 20260713
) {
  difference <- difference[is.finite(difference)]

  if (length(difference) < 2L) {
    return(c(NA_real_, NA_real_))
  }

  set.seed(seed)

  bootstrap_values <- replicate(
    iterations,
    stats::median(
      sample(
        difference,
        size = length(difference),
        replace = TRUE
      ),
      na.rm = TRUE
    )
  )

  stats::quantile(
    bootstrap_values,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
}

extract_model_coefficients <- function(model_object, model_type) {
  if (is.null(model_object)) {
    return(tibble::tibble())
  }

  coefficient_matrix <- summary(model_object)$coefficients

  coefficient_table <- as.data.frame(
    coefficient_matrix,
    stringsAsFactors = FALSE
  ) %>%
    tibble::rownames_to_column("term")

  if (inherits(model_object, "glm")) {
    coefficient_table <- coefficient_table %>%
      dplyr::transmute(
        term = term,
        estimate = exp(Estimate),
        conf_low = exp(
          Estimate - 1.96 * `Std. Error`
        ),
        conf_high = exp(
          Estimate + 1.96 * `Std. Error`
        ),
        p_value = `Pr(>|z|)`,
        model_type = model_type
      )
  } else {
    coefficient_table <- coefficient_table %>%
      dplyr::transmute(
        term = term,
        estimate = Estimate,
        conf_low = Estimate - 1.96 * `Std. Error`,
        conf_high = Estimate + 1.96 * `Std. Error`,
        p_value = `Pr(>|t|)`,
        model_type = model_type
      )
  }

  coefficient_table
}

make_placeholder_plot <- function(title_text, body_text) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.58,
      label = title_text,
      fontface = "bold",
      size = 5
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.42,
      label = body_text,
      size = 4,
      colour = "#58616A"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

## =============================================================================
## 3. Load or reconstruct metric_data, cnv_delta, and gene_info
## =============================================================================

## metric_data ---------------------------------------------------------------
if (!exists("metric_data", inherits = FALSE)) {
  if (!file.exists(metric_csv_file)) {
    stop(
      "metric_data is absent from the current R session and the saved table ",
      "cannot be found:\n  ",
      metric_csv_file,
      "\nRun State3_inferCNV_malignancy_analysis.R first or update out_dir."
    )
  }

  metric_data <- utils::read.csv(
    metric_csv_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

required_metric_columns <- c(
  "cell",
  "sample_id",
  "reference_status",
  "state_label",
  "cnv_burden",
  "State3_signature_score",
  "S.Score",
  "G2M.Score",
  "nCount_RNA"
)

missing_metric_columns <- setdiff(
  required_metric_columns,
  colnames(metric_data)
)

if (length(missing_metric_columns) > 0L) {
  stop(
    "metric_data lacks required column(s): ",
    paste(missing_metric_columns, collapse = ", ")
  )
}

metric_data <- metric_data %>%
  dplyr::mutate(
    cell = as.character(cell),
    sample_id = as.character(sample_id),
    reference_status = as.character(reference_status),
    state_label = factor(
      as.character(state_label),
      levels = state_levels
    )
  )

## inferCNV object ------------------------------------------------------------
if (
  !exists("cnv_delta", inherits = FALSE) ||
    !exists("gene_info", inherits = FALSE)
) {
  if (!file.exists(infercnv_rds_file)) {
    stop(
      "cnv_delta/gene_info are absent from the current R session and the ",
      "saved inferCNV object cannot be found:\n  ",
      infercnv_rds_file,
      "\nRun State3_inferCNV_malignancy_analysis.R first or update out_dir."
    )
  }

  message2(
    "Loading saved inferCNV object and reconstructing cnv_delta/gene_info."
  )

  infercnv_final <- readRDS(
    infercnv_rds_file
  )

  cnv_expression <- infercnv_final@expr.data

  if (
    is.null(cnv_expression) ||
      nrow(cnv_expression) == 0L ||
      ncol(cnv_expression) == 0L
  ) {
    stop("The saved inferCNV object contains an empty expr.data matrix.")
  }

  reference_cells <- metric_data %>%
    dplyr::filter(
      reference_status == "HV reference"
    ) %>%
    dplyr::pull(cell) %>%
    intersect(colnames(cnv_expression))

  if (length(reference_cells) < 3L) {
    stop(
      "Fewer than three healthy-reference cells overlap the inferCNV matrix. ",
      "Check cell identifiers in metric_data and infercnv_final@expr.data."
    )
  }

  reference_center <- matrixStats::rowMedians(
    as.matrix(
      cnv_expression[
        ,
        reference_cells,
        drop = FALSE
      ]
    ),
    na.rm = TRUE
  )

  cnv_delta <- sweep(
    cnv_expression,
    MARGIN = 1,
    STATS = reference_center,
    FUN = "-"
  )

  final_gene_order <- tryCatch(
    as.data.frame(
      infercnv_final@gene_order,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (
    is.null(final_gene_order) ||
      nrow(final_gene_order) == 0L
  ) {
    stop(
      "The saved inferCNV object does not contain usable gene-order ",
      "information."
    )
  }

  if (
    is.null(rownames(final_gene_order)) ||
      any(!nzchar(rownames(final_gene_order))) ||
      all(grepl("^[0-9]+$", rownames(final_gene_order)))
  ) {
    if (nrow(final_gene_order) != nrow(cnv_expression)) {
      stop(
        "Gene-order row names are unavailable and the number of rows does not ",
        "match infercnv_final@expr.data."
      )
    }
    rownames(final_gene_order) <- rownames(cnv_expression)
  }

  if (!all(rownames(cnv_expression) %in% rownames(final_gene_order))) {
    if (nrow(final_gene_order) == nrow(cnv_expression)) {
      rownames(final_gene_order) <- rownames(cnv_expression)
    } else {
      stop(
        "The genes in infercnv_final@expr.data cannot be matched to the ",
        "saved gene-order table."
      )
    }
  }

  gene_info <- final_gene_order[
    rownames(cnv_expression),
    ,
    drop = FALSE
  ]

  chromosome_candidates <- c(
    "chr",
    "chromosome",
    "chrom",
    "seqnames"
  )
  start_candidates <- c(
    "start",
    "gene_start",
    "tx_start"
  )
  stop_candidates <- c(
    "stop",
    "end",
    "gene_end",
    "tx_end"
  )

  chromosome_column <- chromosome_candidates[
    chromosome_candidates %in% colnames(gene_info)
  ][1]

  start_column <- start_candidates[
    start_candidates %in% colnames(gene_info)
  ][1]

  stop_column <- stop_candidates[
    stop_candidates %in% colnames(gene_info)
  ][1]

  if (
    is.na(chromosome_column) ||
      is.na(start_column)
  ) {
    stop(
      "Cannot identify chromosome/start columns in inferCNV gene_order. ",
      "Available columns: ",
      paste(colnames(gene_info), collapse = ", ")
    )
  }

  gene_info$gene <- rownames(gene_info)
  gene_info$chr <- as.character(
    gene_info[[chromosome_column]]
  )
  gene_info$chr <- paste0(
    "chr",
    sub("^chr", "", gene_info$chr)
  )
  gene_info$start <- as.numeric(
    gene_info[[start_column]]
  )

  if (is.na(stop_column)) {
    gene_info$stop <- gene_info$start
  } else {
    gene_info$stop <- as.numeric(
      gene_info[[stop_column]]
    )
  }
} else {
  message2(
    "Using cnv_delta and gene_info already present in the R session."
  )

  if (!"gene" %in% colnames(gene_info)) {
    gene_info$gene <- rownames(gene_info)
  }

  gene_info$chr <- paste0(
    "chr",
    sub("^chr", "", as.character(gene_info$chr))
  )
  gene_info$start <- as.numeric(gene_info$start)

  if (!"stop" %in% colnames(gene_info)) {
    gene_info$stop <- gene_info$start
  }
  gene_info$stop <- as.numeric(gene_info$stop)
}

common_cells <- intersect(
  colnames(cnv_delta),
  metric_data$cell
)

if (length(common_cells) == 0L) {
  stop("No cells overlap between cnv_delta and metric_data.")
}

cnv_delta <- cnv_delta[
  ,
  common_cells,
  drop = FALSE
]

metric_data <- metric_data %>%
  dplyr::filter(
    cell %in% common_cells
  )

## =============================================================================
## 4. Calculate signed, absolute, and event-direction focal CNV scores
## =============================================================================

region_gene_manifest <- list()
focal_score_list <- list()

for (i in seq_len(nrow(driver_regions))) {
  region_row <- driver_regions[i, , drop = FALSE]

  genes_in_region <- gene_info %>%
    dplyr::filter(
      chr == as.character(region_row$chr),
      is.finite(start),
      is.finite(stop),
      start <= as.numeric(region_row$end),
      stop >= as.numeric(region_row$start)
    ) %>%
    dplyr::pull(gene) %>%
    unique() %>%
    intersect(rownames(cnv_delta))

  region_gene_manifest[[i]] <- tibble::tibble(
    event_id = as.character(region_row$event_id),
    event_label = as.character(region_row$event_label),
    gene = genes_in_region
  )

  if (length(genes_in_region) < minimum_genes_per_region) {
    warning(
      "Region ", as.character(region_row$event_label),
      " contains only ", length(genes_in_region),
      " retained inferCNV genes and will be excluded."
    )
    next
  }

  region_matrix <- cnv_delta[
    genes_in_region,
    common_cells,
    drop = FALSE
  ]

  signed_score <- colMeans(
    region_matrix,
    na.rm = TRUE
  )

  absolute_score <- colMeans(
    abs(region_matrix),
    na.rm = TRUE
  )

  event_direction_score <- signed_score *
    as.numeric(region_row$direction)

  focal_score_list[[i]] <- tibble::tibble(
    cell = names(signed_score),
    event_id = as.character(region_row$event_id),
    event_label = as.character(region_row$event_label),
    direction_label = as.character(region_row$direction_label),
    n_genes = length(genes_in_region),
    signed_focal_cnv = as.numeric(signed_score),
    absolute_focal_deviation = as.numeric(absolute_score),
    event_direction_score = as.numeric(event_direction_score)
  )
}

region_gene_manifest <- dplyr::bind_rows(
  region_gene_manifest
)

focal_cnv_scores <- dplyr::bind_rows(
  focal_score_list
)

if (nrow(focal_cnv_scores) == 0L) {
  stop(
    "None of the focal regions contained enough genes. ",
    "Inspect Table_focal_region_gene_manifest.csv and coordinate annotation."
  )
}

focal_cnv_scores <- focal_cnv_scores %>%
  dplyr::left_join(
    metric_data %>%
      dplyr::select(
        cell,
        sample_id,
        reference_status,
        state_label,
        cnv_burden,
        State3_signature_score,
        S.Score,
        G2M.Score,
        nCount_RNA
      ),
    by = "cell"
  ) %>%
  dplyr::mutate(
    event_label = factor(
      event_label,
      levels = as.character(driver_regions$event_label)
    ),
    state_label = factor(
      state_label,
      levels = state_levels
    )
  )

utils::write.csv(
  region_gene_manifest,
  file = file.path(
    analysis_out_dir,
    "Table_focal_region_gene_manifest.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  focal_cnv_scores,
  file = file.path(
    analysis_out_dir,
    "Table_cell_level_focal_CNV_scores.csv"
  ),
  row.names = FALSE
)

## =============================================================================
## 5. Sample-state summaries and paired State3 versus State2 tests
## =============================================================================

focal_sample_summary <- focal_cnv_scores %>%
  dplyr::filter(
    reference_status == "MM query",
    state_label %in% c("State2", "State3")
  ) %>%
  dplyr::group_by(
    sample_id,
    state_label,
    event_id,
    event_label
  ) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_signed_focal_cnv = stats::median(
      signed_focal_cnv,
      na.rm = TRUE
    ),
    median_absolute_focal_deviation = stats::median(
      absolute_focal_deviation,
      na.rm = TRUE
    ),
    median_event_direction_score = stats::median(
      event_direction_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_cells >= minimum_cells_per_sample_state
  )

focal_sample_wide <- focal_sample_summary %>%
  dplyr::select(
    sample_id,
    event_id,
    event_label,
    state_label,
    median_event_direction_score
  ) %>%
  tidyr::pivot_wider(
    names_from = state_label,
    values_from = median_event_direction_score
  )

if (!"State2" %in% colnames(focal_sample_wide)) {
  focal_sample_wide$State2 <- NA_real_
}
if (!"State3" %in% colnames(focal_sample_wide)) {
  focal_sample_wide$State3 <- NA_real_
}

paired_test_list <- lapply(
  split(
    focal_sample_wide,
    focal_sample_wide$event_id
  ),
  function(one_event) {
    keep <- is.finite(one_event$State2) &
      is.finite(one_event$State3)

    one_event <- one_event[
      keep,
      ,
      drop = FALSE
    ]

    difference <- one_event$State3 -
      one_event$State2

    wilcox_result <- safe_paired_wilcox(
      one_event$State3,
      one_event$State2
    )

    bootstrap_ci <- bootstrap_median_difference_ci(
      difference,
      iterations = bootstrap_iterations
    )

    tibble::tibble(
      event_id = one_event$event_id[1],
      event_label = as.character(
        one_event$event_label[1]
      ),
      n_paired_samples = nrow(one_event),
      median_State2 = if (nrow(one_event) > 0L) {
        stats::median(
          one_event$State2,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      median_State3 = if (nrow(one_event) > 0L) {
        stats::median(
          one_event$State3,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      median_paired_difference = if (nrow(one_event) > 0L) {
        stats::median(
          difference,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      hodges_lehmann_shift = wilcox_result$estimate,
      wilcox_conf_low = wilcox_result$conf_low,
      wilcox_conf_high = wilcox_result$conf_high,
      bootstrap_conf_low = bootstrap_ci[1],
      bootstrap_conf_high = bootstrap_ci[2],
      p_value = wilcox_result$p_value
    )
  }
)

focal_paired_tests <- dplyr::bind_rows(
  paired_test_list
) %>%
  dplyr::mutate(
    padj = stats::p.adjust(
      p_value,
      method = "BH"
    ),
    event_label = factor(
      event_label,
      levels = as.character(driver_regions$event_label)
    )
  )

utils::write.csv(
  focal_sample_summary,
  file = file.path(
    analysis_out_dir,
    "Table_sample_state_focal_CNV_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  focal_paired_tests,
  file = file.path(
    analysis_out_dir,
    "Table_State3_vs_State2_paired_focal_CNV_tests.csv"
  ),
  row.names = FALSE
)

## =============================================================================
## 6. Regional signed-CNV profiles using sample-state aggregation
## =============================================================================

query_metadata <- metric_data %>%
  dplyr::filter(
    reference_status == "MM query",
    state_label %in% c("State2", "State3")
  ) %>%
  dplyr::select(
    cell,
    sample_id,
    state_label
  )

profile_list <- list()
anchor_list <- list()

for (i in seq_len(nrow(driver_regions))) {
  region_row <- driver_regions[i, , drop = FALSE]

  region_gene_info <- gene_info %>%
    dplyr::filter(
      chr == as.character(region_row$chr),
      is.finite(start),
      is.finite(stop),
      start <= as.numeric(region_row$end),
      stop >= as.numeric(region_row$start),
      gene %in% rownames(cnv_delta)
    ) %>%
    dplyr::arrange(start) %>%
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )

  if (nrow(region_gene_info) < minimum_genes_per_region) {
    next
  }

  region_cells <- intersect(
    query_metadata$cell,
    colnames(cnv_delta)
  )

  region_cell_meta <- query_metadata[
    match(region_cells, query_metadata$cell),
    ,
    drop = FALSE
  ]

  region_matrix <- as.matrix(
    cnv_delta[
      region_gene_info$gene,
      region_cells,
      drop = FALSE
    ]
  )

  group_id <- paste(
    region_cell_meta$sample_id,
    region_cell_meta$state_label,
    sep = "|||"
  )

  split_indices <- split(
    seq_along(region_cells),
    group_id
  )

  sample_state_matrix <- vapply(
    split_indices,
    function(index_vector) {
      rowMeans(
        region_matrix[
          ,
          index_vector,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    },
    numeric(nrow(region_matrix))
  )

  if (is.null(dim(sample_state_matrix))) {
    sample_state_matrix <- matrix(
      sample_state_matrix,
      ncol = 1L
    )
    colnames(sample_state_matrix) <- names(split_indices)
  }

  sample_state_long <- as.data.frame(
    sample_state_matrix,
    stringsAsFactors = FALSE
  ) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "sample_state",
      values_to = "signed_relative_cnv"
    ) %>%
    tidyr::separate(
      sample_state,
      into = c("sample_id", "state_label"),
      sep = "\\|\\|\\|",
      remove = TRUE,
      extra = "merge"
    ) %>%
    dplyr::left_join(
      region_gene_info %>%
        dplyr::select(
          gene,
          start,
          stop
        ),
      by = "gene"
    ) %>%
    dplyr::mutate(
      event_id = as.character(region_row$event_id),
      event_label = as.character(region_row$event_label),
      position_mb = (
        as.numeric(start) +
          as.numeric(stop)
      ) / 2e6
    )

  profile_summary <- sample_state_long %>%
    dplyr::group_by(
      event_id,
      event_label,
      state_label,
      gene,
      position_mb
    ) %>%
    dplyr::summarise(
      median_signed_cnv = stats::median(
        signed_relative_cnv,
        na.rm = TRUE
      ),
      q25_signed_cnv = stats::quantile(
        signed_relative_cnv,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      ),
      q75_signed_cnv = stats::quantile(
        signed_relative_cnv,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      ),
      n_sample_states = sum(
        is.finite(signed_relative_cnv)
      ),
      .groups = "drop"
    )

  profile_list[[i]] <- profile_summary

  anchor_gene_vector <- strsplit(
    as.character(region_row$anchor_genes),
    split = ";",
    fixed = TRUE
  )[[1]]

  anchor_data <- region_gene_info %>%
    dplyr::filter(
      gene %in% anchor_gene_vector
    ) %>%
    dplyr::transmute(
      event_id = as.character(region_row$event_id),
      event_label = as.character(region_row$event_label),
      gene = gene,
      position_mb = (
        as.numeric(start) +
          as.numeric(stop)
      ) / 2e6
    )

  anchor_list[[i]] <- anchor_data
}

regional_profile_summary <- dplyr::bind_rows(
  profile_list
) %>%
  dplyr::mutate(
    event_label = factor(
      event_label,
      levels = as.character(driver_regions$event_label)
    ),
    state_label = factor(
      state_label,
      levels = c("State2", "State3")
    )
  )

anchor_data <- dplyr::bind_rows(
  anchor_list
) %>%
  dplyr::mutate(
    event_label = factor(
      event_label,
      levels = as.character(driver_regions$event_label)
    )
  )

utils::write.csv(
  regional_profile_summary,
  file = file.path(
    analysis_out_dir,
    "Table_regional_signed_CNV_profiles.csv"
  ),
  row.names = FALSE
)

## =============================================================================
## 7. State3 signature and State3 identity adjusted for global CNV burden
## =============================================================================

focal_wide <- focal_cnv_scores %>%
  dplyr::filter(
    reference_status == "MM query"
  ) %>%
  dplyr::select(
    cell,
    event_id,
    event_direction_score
  ) %>%
  tidyr::pivot_wider(
    names_from = event_id,
    values_from = event_direction_score
  )

model_data <- metric_data %>%
  dplyr::filter(
    reference_status == "MM query",
    state_label %in% c("State2", "State3")
  ) %>%
  dplyr::left_join(
    focal_wide,
    by = "cell"
  ) %>%
  dplyr::mutate(
    is_State3 = as.integer(
      state_label == "State3"
    ),
    log_library_size = log10(
      as.numeric(nCount_RNA) + 1
    )
  )

required_model_variables <- c(
  "State3_signature_score",
  "cnv_burden",
  "gain_1q21",
  "gain_MYC_8q24",
  "loss_TP53_17p13",
  "S.Score",
  "G2M.Score",
  "log_library_size"
)

model_complete <- stats::complete.cases(
  model_data[
    ,
    required_model_variables,
    drop = FALSE
  ]
)

model_data <- model_data[
  model_complete,
  ,
  drop = FALSE
]

## Retain samples containing both State2 and State3 to make the sample fixed
## effect represent a within-sample comparison.
eligible_samples <- model_data %>%
  dplyr::group_by(
    sample_id
  ) %>%
  dplyr::summarise(
    n_State2 = sum(
      state_label == "State2"
    ),
    n_State3 = sum(
      state_label == "State3"
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_State2 >= minimum_cells_per_sample_state,
    n_State3 >= minimum_cells_per_sample_state
  ) %>%
  dplyr::pull(sample_id)

model_data <- model_data %>%
  dplyr::filter(
    sample_id %in% eligible_samples
  ) %>%
  dplyr::mutate(
    sample_id = factor(sample_id),
    global_cnv_z = safe_zscore(cnv_burden),
    focal_1q21_z = safe_zscore(gain_1q21),
    focal_MYC_z = safe_zscore(gain_MYC_8q24),
    focal_TP53loss_z = safe_zscore(loss_TP53_17p13),
    S_score_z = safe_zscore(S.Score),
    G2M_score_z = safe_zscore(G2M.Score),
    library_size_z = safe_zscore(log_library_size)
  )

signature_model <- NULL
identity_model <- NULL

if (
  nrow(model_data) >= 50L &&
    length(unique(model_data$sample_id)) >= 2L
) {
  signature_model <- tryCatch(
    stats::lm(
      State3_signature_score ~
        global_cnv_z +
        focal_1q21_z +
        focal_MYC_z +
        focal_TP53loss_z +
        S_score_z +
        G2M_score_z +
        library_size_z +
        sample_id,
      data = model_data
    ),
    error = function(e) {
      warning(
        "The CNV-adjusted State3 signature model failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  identity_model <- tryCatch(
    suppressWarnings(
      stats::glm(
        is_State3 ~
          global_cnv_z +
          focal_1q21_z +
          focal_MYC_z +
          focal_TP53loss_z +
          S_score_z +
          G2M_score_z +
          library_size_z +
          sample_id,
        data = model_data,
        family = stats::binomial()
      )
    ),
    error = function(e) {
      warning(
        "The CNV-adjusted State3 identity model failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
} else {
  warning(
    "Insufficient complete within-sample State2/State3 data for adjusted ",
    "models. Required: at least 50 cells from at least two eligible samples."
  )
}

signature_coefficients <- extract_model_coefficients(
  signature_model,
  model_type = "State3 signature (standardized beta)"
)

identity_coefficients <- extract_model_coefficients(
  identity_model,
  model_type = "State3 identity (odds ratio)"
)

core_terms <- c(
  "global_cnv_z",
  "focal_1q21_z",
  "focal_MYC_z",
  "focal_TP53loss_z"
)

term_labels <- c(
  "global_cnv_z" = "Global CNV burden",
  "focal_1q21_z" = "1q21 gain-like signal",
  "focal_MYC_z" = "8q24 MYC gain-like signal",
  "focal_TP53loss_z" = "17p13 TP53 loss-like signal"
)

signature_core_coefficients <- signature_coefficients %>%
  dplyr::filter(
    term %in% core_terms
  ) %>%
  dplyr::mutate(
    label = unname(term_labels[term]),
    label = factor(
      label,
      levels = rev(unname(term_labels[core_terms]))
    )
  )

identity_core_coefficients <- identity_coefficients %>%
  dplyr::filter(
    term %in% core_terms
  ) %>%
  dplyr::mutate(
    label = unname(term_labels[term]),
    label = factor(
      label,
      levels = rev(unname(term_labels[core_terms]))
    )
  )

utils::write.csv(
  model_data,
  file = file.path(
    analysis_out_dir,
    "Table_adjusted_model_input_cells.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  signature_coefficients,
  file = file.path(
    analysis_out_dir,
    "Table_State3_signature_adjusted_model_coefficients.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  identity_coefficients,
  file = file.path(
    analysis_out_dir,
    "Table_State3_identity_adjusted_model_coefficients.csv"
  ),
  row.names = FALSE
)

signature_model_text <- c(
  "Model:",
  paste(
    "State3_signature_score ~ global_cnv_z + focal_1q21_z +",
    "focal_MYC_z + focal_TP53loss_z + S_score_z + G2M_score_z +",
    "library_size_z + sample_id"
  ),
  "",
  paste0(
    "Eligible samples: ",
    length(unique(model_data$sample_id))
  ),
  paste0(
    "Analyzed cells: ",
    nrow(model_data)
  ),
  "",
  if (is.null(signature_model)) {
    "Model was not fitted."
  } else {
    capture.output(summary(signature_model))
  }
)

writeLines(
  signature_model_text,
  con = file.path(
    analysis_out_dir,
    "Model_State3_signature_adjusted_global_and_focal_CNV.txt"
  )
)

identity_model_text <- c(
  "Model:",
  paste(
    "is_State3 ~ global_cnv_z + focal_1q21_z + focal_MYC_z +",
    "focal_TP53loss_z + S_score_z + G2M_score_z +",
    "library_size_z + sample_id"
  ),
  "",
  paste0(
    "Eligible samples: ",
    length(unique(model_data$sample_id))
  ),
  paste0(
    "Analyzed cells: ",
    nrow(model_data)
  ),
  "",
  if (is.null(identity_model)) {
    "Model was not fitted."
  } else {
    capture.output(summary(identity_model))
  }
)

writeLines(
  identity_model_text,
  con = file.path(
    analysis_out_dir,
    "Model_State3_identity_adjusted_global_and_focal_CNV.txt"
  )
)

## =============================================================================
## 8. Publication-style figures
## =============================================================================

## 8.1 Regional signed-CNV profiles -------------------------------------------
if (nrow(regional_profile_summary) > 0L) {
  anchor_y <- regional_profile_summary %>%
    dplyr::group_by(
      event_label
    ) %>%
    dplyr::summarise(
      y_position = max(
        q75_signed_cnv,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  anchor_plot_data <- anchor_data %>%
    dplyr::left_join(
      anchor_y,
      by = "event_label"
    )

  p_regional_profile <- ggplot2::ggplot(
    regional_profile_summary,
    ggplot2::aes(
      x = position_mb,
      y = median_signed_cnv,
      colour = state_label,
      fill = state_label
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "#8B9298"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = q25_signed_cnv,
        ymax = q75_signed_cnv
      ),
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_line(
      linewidth = 1.0,
      alpha = 0.95
    ) +
    ggplot2::geom_vline(
      data = anchor_plot_data,
      ggplot2::aes(
        xintercept = position_mb
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.5,
      colour = "#303840"
    ) +
    ggrepel::geom_text_repel(
      data = anchor_plot_data,
      ggplot2::aes(
        x = position_mb,
        y = y_position,
        label = gene
      ),
      inherit.aes = FALSE,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      seed = 20260713,
      size = 3.3,
      fontface = "bold",
      colour = "#303840"
    ) +
    ggplot2::facet_wrap(
      ~event_label,
      scales = "free_x",
      nrow = 1
    ) +
    ggplot2::scale_colour_manual(
      values = state_colors[c("State2", "State3")],
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = state_colors[c("State2", "State3")],
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Regional inferred CNV architecture at candidate myeloma drivers",
      subtitle = paste(
        "Lines show medians across sample-state summaries;",
        "ribbons show interquartile ranges"
      ),
      x = "Genomic position (Mb, GRCh38)",
      y = "Signed relative CNV",
      colour = NULL,
      fill = NULL
    ) +
    theme_publication() +
    ggplot2::theme(
      legend.position = "top",
      panel.spacing.x = grid::unit(1.2, "lines")
    )
} else {
  p_regional_profile <- make_placeholder_plot(
    "Regional CNV profiles",
    "No focal region passed the minimum gene threshold."
  )
}

## 8.2 Paired sample-state focal-event plot -----------------------------------
paired_plot_data <- focal_sample_summary %>%
  dplyr::filter(
    state_label %in% c("State2", "State3")
  ) %>%
  dplyr::mutate(
    state_label = factor(
      state_label,
      levels = c("State2", "State3")
    ),
    event_label = factor(
      event_label,
      levels = as.character(driver_regions$event_label)
    )
  )

paired_annotations <- focal_paired_tests %>%
  dplyr::left_join(
    paired_plot_data %>%
      dplyr::group_by(
        event_label
      ) %>%
      dplyr::summarise(
        annotation_y = max(
          median_event_direction_score,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "event_label"
  ) %>%
  dplyr::mutate(
    annotation_text = paste0(
      "paired n = ",
      n_paired_samples,
      "\n",
      format_p(p_value)
    )
  )

p_paired_focal <- ggplot2::ggplot(
  paired_plot_data,
  ggplot2::aes(
    x = state_label,
    y = median_event_direction_score,
    group = sample_id
  )
) +
  ggplot2::geom_line(
    colour = "#B8BEC3",
    linewidth = 0.5,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = state_label
    ),
    shape = 21,
    size = 2.6,
    stroke = 0.35,
    colour = "#303840"
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(
      group = 1
    ),
    fun = stats::median,
    geom = "line",
    linewidth = 1.2,
    colour = "#111820"
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(
      group = 1
    ),
    fun = stats::median,
    geom = "point",
    shape = 23,
    size = 3.8,
    fill = "white",
    colour = "#111820"
  ) +
  ggplot2::geom_text(
    data = paired_annotations,
    ggplot2::aes(
      x = 1.5,
      y = annotation_y,
      label = annotation_text
    ),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 3.2,
    colour = "#303840"
  ) +
  ggplot2::facet_wrap(
    ~event_label,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::scale_fill_manual(
    values = state_colors[c("State2", "State3")],
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Within-sample focal driver signals",
    subtitle = paste(
      "Positive values indicate gain-like 1q21/MYC or loss-like TP53",
      "signals in the expected event direction"
    ),
    x = NULL,
    y = "Median event-direction CNV score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(
      face = "bold"
    )
  )

## 8.3 Forest plot of paired State3 minus State2 shifts ------------------------
forest_data <- focal_paired_tests %>%
  dplyr::mutate(
    plot_estimate = median_paired_difference,
    plot_conf_low = bootstrap_conf_low,
    plot_conf_high = bootstrap_conf_high,
    p_label = paste0(
      format_p(p_value),
      "; FDR = ",
      ifelse(
        is.finite(padj),
        formatC(
          padj,
          format = "g",
          digits = 3
        ),
        "NA"
      )
    )
  )

p_paired_forest <- ggplot2::ggplot(
  forest_data,
  ggplot2::aes(
    x = plot_estimate,
    y = event_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "#8B9298"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = plot_conf_low,
      xmax = plot_conf_high
    ),
    height = 0.14,
    linewidth = 0.75,
    colour = "#303840"
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 3.5,
    stroke = 0.45,
    fill = state_colors[["State3"]],
    colour = "#303840"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = p_label
    ),
    hjust = -0.05,
    nudge_x = 0.01 * max(
      abs(c(
        forest_data$plot_conf_low,
        forest_data$plot_conf_high
      )),
      na.rm = TRUE
    ),
    size = 3.1,
    colour = "#303840"
  ) +
  ggplot2::labs(
    title = "Paired focal-event effect sizes",
    subtitle = paste(
      "Median State3−State2 shift with bootstrap 95% confidence interval;",
      "positive values favor State3"
    ),
    x = "State3 − State2 event-direction score",
    y = NULL
  ) +
  theme_publication() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "bold"
    )
  )

## 8.4 Adjusted State3 signature coefficients ---------------------------------
if (nrow(signature_core_coefficients) > 0L) {
  p_signature_coefficients <- ggplot2::ggplot(
    signature_core_coefficients,
    ggplot2::aes(
      x = estimate,
      y = label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.55,
      colour = "#8B9298"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = conf_low,
        xmax = conf_high
      ),
      height = 0.15,
      linewidth = 0.75,
      colour = "#303840"
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 3.5,
      stroke = 0.45,
      fill = "#D59B45",
      colour = "#303840"
    ) +
    ggplot2::labs(
      title = "Independent associations with State3 signature",
      subtitle = paste(
        "Linear model adjusted for cell cycle, library size,",
        "and sample fixed effects"
      ),
      x = "Adjusted coefficient per 1-SD increase",
      y = NULL
    ) +
    theme_publication() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        face = "bold"
      )
    )
} else {
  p_signature_coefficients <- make_placeholder_plot(
    "State3 signature model",
    "The adjusted linear model could not be fitted."
  )
}

## 8.5 Adjusted State3 identity odds ratios -----------------------------------
if (nrow(identity_core_coefficients) > 0L) {
  p_identity_coefficients <- ggplot2::ggplot(
    identity_core_coefficients,
    ggplot2::aes(
      x = estimate,
      y = label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.55,
      colour = "#8B9298"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = conf_low,
        xmax = conf_high
      ),
      height = 0.15,
      linewidth = 0.75,
      colour = "#303840"
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 3.5,
      stroke = 0.45,
      fill = state_colors[["State3"]],
      colour = "#303840"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "State3 identity after adjustment for global CNV",
      subtitle = paste(
        "Logistic model comparing State3 with State2;",
        "sample fixed effects included"
      ),
      x = "Odds ratio per 1-SD increase (log scale)",
      y = NULL
    ) +
    theme_publication() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        face = "bold"
      )
    )
} else {
  p_identity_coefficients <- make_placeholder_plot(
    "State3 identity model",
    "The adjusted logistic model could not be fitted."
  )
}

## =============================================================================
## 9. Save figures
## =============================================================================

save_pdf(
  plot_object = p_regional_profile,
  filename = file.path(
    analysis_out_dir,
    "FigA_regional_signed_focal_CNV_profiles.pdf"
  ),
  width = 14,
  height = 5.5
)

save_pdf(
  plot_object = p_paired_focal,
  filename = file.path(
    analysis_out_dir,
    "FigB_sample_paired_focal_CNV_signals.pdf"
  ),
  width = 12,
  height = 5.5
)

save_pdf(
  plot_object = p_paired_forest,
  filename = file.path(
    analysis_out_dir,
    "FigC_paired_focal_CNV_effect_sizes.pdf"
  ),
  width = 9.5,
  height = 5.5
)

save_pdf(
  plot_object = p_signature_coefficients,
  filename = file.path(
    analysis_out_dir,
    "FigD_State3_signature_adjusted_coefficients.pdf"
  ),
  width = 8.5,
  height = 5.5
)

save_pdf(
  plot_object = p_identity_coefficients,
  filename = file.path(
    analysis_out_dir,
    "FigE_State3_identity_adjusted_odds_ratios.pdf"
  ),
  width = 8.5,
  height = 5.5
)

integrated_figure <- (
  p_regional_profile
) / (
  p_paired_focal | p_paired_forest
) / (
  p_signature_coefficients | p_identity_coefficients
) +
  patchwork::plot_layout(
    heights = c(1.05, 1, 1)
  ) +
  patchwork::plot_annotation(
    title = paste(
      "State3 focal driver architecture after accounting for",
      "global inferred CNV burden"
    ),
    subtitle = paste(
      "Testing whether selective gain-like/loss-like events, rather than",
      "maximal genome-wide CNV burden, characterize State3"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 18,
        colour = "#20262D"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 11,
        colour = "#58616A"
      )
    )
  )

save_pdf(
  plot_object = integrated_figure,
  filename = file.path(
    analysis_out_dir,
    "Fig_State3_focal_CNV_and_adjusted_signature_integrated.pdf"
  ),
  width = 17,
  height = 16
)

## =============================================================================
## 10. Save reproducibility objects
## =============================================================================

analysis_results <- list(
  driver_regions = driver_regions,
  region_gene_manifest = region_gene_manifest,
  focal_cnv_scores = focal_cnv_scores,
  focal_sample_summary = focal_sample_summary,
  focal_paired_tests = focal_paired_tests,
  regional_profile_summary = regional_profile_summary,
  signature_model = signature_model,
  identity_model = identity_model,
  signature_coefficients = signature_coefficients,
  identity_coefficients = identity_coefficients
)

saveRDS(
  analysis_results,
  file = file.path(
    analysis_out_dir,
    "State3_focal_CNV_analysis_results.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    analysis_out_dir,
    "sessionInfo_State3_focal_CNV_analysis.txt"
  )
)

message2(
  "\nCompleted State3 focal CNV analysis.\n",
  "Output directory:\n  ",
  analysis_out_dir,
  "\n\nPrimary outputs:\n",
  "  Fig_State3_focal_CNV_and_adjusted_signature_integrated.pdf\n",
  "  Table_State3_vs_State2_paired_focal_CNV_tests.csv\n",
  "  Model_State3_signature_adjusted_global_and_focal_CNV.txt\n",
  "  Model_State3_identity_adjusted_global_and_focal_CNV.txt\n"
)
