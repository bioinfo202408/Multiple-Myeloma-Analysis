############################################################
## Step5_State3Risk_SomaticMutation_MAF_Oncoplot.R
## Purpose:
##   1) Convert UCSC Xena MMRF-CoMMpass WXS somatic mutation TSV to MAF-like format
##   2) Match mutation samples to State3/RSF training-set risk groups
##   3) Compare somatic mutation landscapes between High-risk and Low-risk groups
##   4) Generate publication-quality mutation waterfall plots and summary figures
##
## Input files:
##   - MMRF-COMMPASS.somaticmutation_wxs.tsv
##   - risk_train.csv from RSF subgroup analysis
##
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

if (!requireNamespace("maftools", quietly = TRUE)) {
  stop("The R package 'maftools' is required. Install it first: BiocManager::install('maftools')", call. = FALSE)
}
suppressPackageStartupMessages(library(maftools))

## =========================
## 0. User settings
## =========================

mut_file <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/somatic_mutation_analys/MMRF-COMMPASS.somaticmutation_wxs.tsv"
risk_file <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/RSF_Subgroup_Analysis/Tables/risk_train.csv"

outdir <- "/home/yjliu/mmProj/data_process/Human/SingleCell_NMF/Clinic/somatic_mutation_analys/State3_RSF_Risk_Mutation_Analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Matching strategy:
##   "auto"     : exact sample-level matching first; fallback to patient-level if needed
##   "sample"   : force mutation$sample <-> risk$Sample
##   "patient"  : force patient_id matching, e.g. MMRF_2401
match_strategy <- "auto"

## Main oncoplot setting
n_top_genes <- 30
min_mut_for_maf_compare <- 5
exome_size_mb <- 38   # approximate coding exome size for manual nonsilent mutation burden

## MM-relevant genes to keep if present; these are appended to top mutated genes
mm_driver_genes <- c(
  "TP53", "NRAS", "KRAS", "BRAF", "DIS3", "FAM46C", "TRAF3", "CYLD",
  "IRF4", "PRDM1", "CCND1", "CCND2", "CCND3", "MYC", "FGFR3", "NSD2",
  "ATM", "ATR", "RB1", "MAX", "HIST1H1E", "SP140", "EGR1", "NFKBIA"
)

## If you know exact mappings for your clinical variables, edit these labels.
gender_map <- c("1" = "Gender_1", "2" = "Gender_2")
age_map <- c("1" = "Age_1", "2" = "Age_2", "3" = "Age_3")
treatment_map <- c("0" = "Treatment_0", "1" = "Treatment_1")

## =========================
## 1. Helper functions
## =========================

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names)
}

safe_write_tsv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file = file, sep = "\t", quote = FALSE, na = "")
}

safe_saveRDS <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = file)
}

open_pdf <- function(file, width = 8, height = 6) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file = file, width = width, height = height, useDingbats = FALSE, onefile = TRUE)
}

ggsave_pdf <- function(filename, plot, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                    units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  } else {
    ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                    units = "in", device = "pdf", bg = "white", limitsize = FALSE)
  }
}

theme_pub <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(colour = "black", face = "bold", size = base_size),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey30", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.key.size = grid::unit(0.38, "cm"),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

p_to_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

extract_patient_id <- function(x) {
  stringr::str_extract(as.character(x), "MMRF_[0-9]+")
}

clean_risk_group <- function(x) {
  x <- stringr::str_squish(as.character(x))
  dplyr::case_when(
    grepl("^high", x, ignore.case = TRUE) ~ "High risk",
    grepl("^low", x, ignore.case = TRUE) ~ "Low risk",
    TRUE ~ x
  )
}

variant_type_from_alleles <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  dplyr::case_when(
    ref == "-" & alt != "-" ~ "INS",
    ref != "-" & alt == "-" ~ "DEL",
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNP",
    nchar(ref) > nchar(alt) ~ "DEL",
    nchar(ref) < nchar(alt) ~ "INS",
    TRUE ~ "ONP"
  )
}

effect_to_maf <- function(effect, ref, alt) {
  effect <- as.character(effect)
  ref <- as.character(ref)
  alt <- as.character(alt)
  dplyr::case_when(
    grepl("missense_variant", effect) ~ "Missense_Mutation",
    grepl("synonymous_variant", effect) ~ "Silent",
    grepl("stop_gained", effect) ~ "Nonsense_Mutation",
    grepl("stop_lost", effect) ~ "Nonstop_Mutation",
    grepl("start_lost|initiator_codon_variant", effect) ~ "Translation_Start_Site",
    grepl("splice_acceptor_variant|splice_donor_variant|splice_region_variant", effect) ~ "Splice_Site",
    grepl("frameshift_variant", effect) & (ref == "-" | nchar(ref) <= nchar(alt)) ~ "Frame_Shift_Ins",
    grepl("frameshift_variant", effect) & (alt == "-" | nchar(ref) > nchar(alt)) ~ "Frame_Shift_Del",
    grepl("inframe_deletion|disruptive_inframe_deletion|conservative_inframe_deletion", effect) ~ "In_Frame_Del",
    grepl("inframe_insertion|disruptive_inframe_insertion|conservative_inframe_insertion", effect) ~ "In_Frame_Ins",
    grepl("protein_altering_variant", effect) ~ "Missense_Mutation",
    grepl("3_prime_UTR_variant", effect) ~ "3'UTR",
    grepl("5_prime_UTR_variant", effect) ~ "5'UTR",
    grepl("intron_variant", effect) ~ "Intron",
    grepl("upstream_gene_variant|downstream_gene_variant|intergenic_variant", effect) ~ "IGR",
    grepl("non_coding_transcript|nc_transcript", effect) ~ "RNA",
    TRUE ~ "Unknown"
  )
}

non_silent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site",
  "Nonstop_Mutation"
)

## =========================
## 2. Read input data
## =========================

if (!file.exists(mut_file)) stop("Mutation file does not exist: ", mut_file, call. = FALSE)
if (!file.exists(risk_file)) stop("Risk file does not exist: ", risk_file, call. = FALSE)

mut_raw <- data.table::fread(mut_file, data.table = FALSE, check.names = FALSE)
risk_raw <- data.table::fread(risk_file, data.table = FALSE, check.names = FALSE)

required_mut_cols <- c("sample", "gene", "chrom", "start", "end", "ref", "alt", "effect")
missing_mut_cols <- setdiff(required_mut_cols, colnames(mut_raw))
if (length(missing_mut_cols) > 0) {
  stop("Mutation file is missing required columns: ", paste(missing_mut_cols, collapse = ", "), call. = FALSE)
}

required_risk_cols <- c("Sample", "risk_group")
missing_risk_cols <- setdiff(required_risk_cols, colnames(risk_raw))
if (length(missing_risk_cols) > 0) {
  stop("Risk file is missing required columns: ", paste(missing_risk_cols, collapse = ", "), call. = FALSE)
}

risk_df <- risk_raw %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    risk_group = clean_risk_group(risk_group),
    patient_id = extract_patient_id(Sample),
    Gender_label = dplyr::recode(as.character(Gender), !!!gender_map, .default = paste0("Gender_", as.character(Gender))),
    Age_label = dplyr::recode(as.character(Age), !!!age_map, .default = paste0("Age_", as.character(Age))),
    ISS_label = paste0("ISS_", as.character(iss_stage)),
    Treatment_label = dplyr::recode(as.character(treatment_type), !!!treatment_map, .default = paste0("Treatment_", as.character(treatment_type)))
  ) %>%
  dplyr::filter(!is.na(Sample), Sample != "") %>%
  dplyr::filter(risk_group %in% c("High risk", "Low risk")) %>%
  dplyr::arrange(Sample, dplyr::desc(riskscore)) %>%
  dplyr::distinct(Sample, .keep_all = TRUE)

safe_write_csv(risk_df, file.path(outdir, "Step5_cleaned_risk_train_table.csv"))

## =========================
## 3. Convert Xena TSV to MAF-like table
## =========================

maf_all <- mut_raw %>%
  dplyr::mutate(
    sample = as.character(sample),
    patient_id = extract_patient_id(sample),
    Hugo_Symbol = as.character(gene),
    Chromosome = gsub("^chr", "", as.character(chrom)),
    Start_Position = as.integer(start),
    End_Position = as.integer(end),
    Reference_Allele = as.character(ref),
    Tumor_Seq_Allele1 = as.character(ref),
    Tumor_Seq_Allele2 = as.character(alt),
    Variant_Classification = effect_to_maf(effect, ref, alt),
    Variant_Type = variant_type_from_alleles(ref, alt),
    Tumor_Sample_Barcode = sample,
    Protein_Change = if ("Amino_Acid_Change" %in% colnames(mut_raw)) as.character(Amino_Acid_Change) else "",
    dna_vaf = if ("dna_vaf" %in% colnames(mut_raw)) as.numeric(dna_vaf) else NA_real_,
    callers = if ("callers" %in% colnames(mut_raw)) as.character(callers) else NA_character_
  ) %>%
  dplyr::filter(!is.na(Hugo_Symbol), Hugo_Symbol != "") %>%
  dplyr::filter(!is.na(Tumor_Sample_Barcode), Tumor_Sample_Barcode != "") %>%
  dplyr::filter(!is.na(Start_Position), !is.na(End_Position)) %>%
  dplyr::select(
    Hugo_Symbol, Chromosome, Start_Position, End_Position,
    Variant_Classification, Variant_Type,
    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    Tumor_Sample_Barcode, Protein_Change,
    sample, patient_id, dna_vaf, effect, callers
  )

maf_nonsilent <- maf_all %>%
  dplyr::filter(Variant_Classification %in% non_silent_classes)

safe_write_tsv(maf_all, file.path(outdir, "Step5_MMRF_COMMPASS_wxs_converted_all_variants.maf"))
safe_write_tsv(maf_nonsilent, file.path(outdir, "Step5_MMRF_COMMPASS_wxs_converted_nonsilent.maf"))

safe_write_csv(
  maf_all %>% dplyr::count(Variant_Classification, sort = TRUE, name = "n_variants"),
  file.path(outdir, "Step5_variant_classification_counts_all_variants.csv")
)
safe_write_csv(
  maf_nonsilent %>% dplyr::count(Variant_Classification, sort = TRUE, name = "n_variants"),
  file.path(outdir, "Step5_variant_classification_counts_nonsilent.csv")
)

## =========================
## 4. Match mutation data with risk group
## =========================

sample_exact_overlap <- intersect(unique(maf_nonsilent$Tumor_Sample_Barcode), unique(risk_df$Sample))
patient_overlap <- intersect(unique(maf_nonsilent$patient_id), unique(risk_df$patient_id))

match_diagnostic <- data.frame(
  metric = c(
    "n_mut_samples_nonsilent",
    "n_risk_samples",
    "n_exact_sample_overlap",
    "n_mut_patients_nonsilent",
    "n_risk_patients",
    "n_patient_overlap"
  ),
  value = c(
    length(unique(maf_nonsilent$Tumor_Sample_Barcode)),
    length(unique(risk_df$Sample)),
    length(sample_exact_overlap),
    length(unique(maf_nonsilent$patient_id)),
    length(unique(risk_df$patient_id)),
    length(patient_overlap)
  )
)
safe_write_csv(match_diagnostic, file.path(outdir, "Step5_mutation_risk_matching_diagnostic.csv"))

if (match_strategy == "sample") {
  selected_match <- "sample"
} else if (match_strategy == "patient") {
  selected_match <- "patient"
} else {
  ## Prefer exact sample matching if it captures enough risk samples.
  selected_match <- ifelse(length(sample_exact_overlap) >= max(10, 0.5 * length(patient_overlap)), "sample", "patient")
}

message("Selected matching strategy: ", selected_match)

if (selected_match == "sample") {
  clinical_for_maf <- risk_df %>%
    dplyr::filter(Sample %in% unique(maf_nonsilent$Tumor_Sample_Barcode)) %>%
    dplyr::transmute(
      Tumor_Sample_Barcode = Sample,
      Patient_ID = patient_id,
      risk_group = risk_group,
      riskscore = riskscore,
      OS.time = OS.time,
      OS = OS,
      set = set,
      Gender = Gender_label,
      Age = Age_label,
      ISS = ISS_label,
      Treatment = Treatment_label
    ) %>%
    dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

  maf_matched <- maf_nonsilent %>%
    dplyr::filter(Tumor_Sample_Barcode %in% clinical_for_maf$Tumor_Sample_Barcode)

} else {
  ## Patient-level fallback: annotate every mutation sample whose patient is in risk_df.
  ## If one patient has multiple risk records, use the record with the highest absolute information retained by first order.
  risk_patient <- risk_df %>%
    dplyr::arrange(patient_id, dplyr::desc(riskscore)) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)

  maf_matched <- maf_nonsilent %>%
    dplyr::filter(patient_id %in% risk_patient$patient_id) %>%
    dplyr::left_join(
      risk_patient %>%
        dplyr::select(patient_id, risk_group, riskscore, OS.time, OS, set, Gender_label, Age_label, ISS_label, Treatment_label),
      by = "patient_id"
    )

  sample_clin <- maf_matched %>%
    dplyr::distinct(Tumor_Sample_Barcode, patient_id, risk_group, riskscore, OS.time, OS, set, Gender_label, Age_label, ISS_label, Treatment_label)

  clinical_for_maf <- sample_clin %>%
    dplyr::transmute(
      Tumor_Sample_Barcode = Tumor_Sample_Barcode,
      Patient_ID = patient_id,
      risk_group = risk_group,
      riskscore = riskscore,
      OS.time = OS.time,
      OS = OS,
      set = set,
      Gender = Gender_label,
      Age = Age_label,
      ISS = ISS_label,
      Treatment = Treatment_label
    ) %>%
    dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

  maf_matched <- maf_matched %>%
    dplyr::select(-dplyr::any_of(c("risk_group", "riskscore", "OS.time", "OS", "set", "Gender_label", "Age_label", "ISS_label", "Treatment_label")))
}

if (nrow(clinical_for_maf) < 4) {
  stop("Fewer than 4 matched mutation samples. Please inspect Step5_mutation_risk_matching_diagnostic.csv.", call. = FALSE)
}

risk_group_counts <- clinical_for_maf %>%
  dplyr::count(risk_group, name = "n_samples")
safe_write_csv(risk_group_counts, file.path(outdir, "Step5_matched_sample_counts_by_risk_group.csv"))
safe_write_csv(clinical_for_maf, file.path(outdir, "Step5_clinical_annotation_for_maftools.csv"))
safe_write_tsv(maf_matched, file.path(outdir, "Step5_MMRF_COMMPASS_wxs_nonsilent_matched_to_risk.maf"))

if (!all(c("High risk", "Low risk") %in% clinical_for_maf$risk_group)) {
  stop("Matched samples do not contain both High risk and Low risk groups.", call. = FALSE)
}

## =========================
## 5. Build maftools MAF objects
## =========================

maf_file_matched <- file.path(outdir, "Step5_MMRF_COMMPASS_wxs_nonsilent_matched_to_risk.maf")
clinical_file <- file.path(outdir, "Step5_clinical_annotation_for_maftools.csv")

maf_obj <- maftools::read.maf(
  maf = maf_file_matched,
  clinicalData = clinical_file,
  isTCGA = FALSE,
  verbose = FALSE
)

high_samples <- clinical_for_maf %>%
  dplyr::filter(risk_group == "High risk") %>%
  dplyr::pull(Tumor_Sample_Barcode)
low_samples <- clinical_for_maf %>%
  dplyr::filter(risk_group == "Low risk") %>%
  dplyr::pull(Tumor_Sample_Barcode)

maf_high <- maftools::subsetMaf(maf = maf_obj, tsb = high_samples, mafObj = TRUE, includeSyn = FALSE)
maf_low <- maftools::subsetMaf(maf = maf_obj, tsb = low_samples, mafObj = TRUE, includeSyn = FALSE)

safe_saveRDS(maf_obj, file.path(outdir, "Step5_maftools_matched_all_risk_groups.rds"))
safe_saveRDS(maf_high, file.path(outdir, "Step5_maftools_high_risk.rds"))
safe_saveRDS(maf_low, file.path(outdir, "Step5_maftools_low_risk.rds"))

## =========================
## 6. Select common genes for comparable oncoplots
## =========================

gene_freq_all <- maf_matched %>%
  dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  dplyr::count(Hugo_Symbol, name = "n_mut_samples") %>%
  dplyr::mutate(freq = n_mut_samples / length(unique(clinical_for_maf$Tumor_Sample_Barcode))) %>%
  dplyr::arrange(dplyr::desc(n_mut_samples), Hugo_Symbol)

genes_show <- unique(c(
  gene_freq_all %>% dplyr::slice_head(n = n_top_genes) %>% dplyr::pull(Hugo_Symbol),
  intersect(mm_driver_genes, gene_freq_all$Hugo_Symbol)
))
genes_show <- genes_show[genes_show %in% gene_freq_all$Hugo_Symbol]

safe_write_csv(gene_freq_all, file.path(outdir, "Step5_all_matched_gene_mutation_frequency.csv"))
safe_write_csv(data.frame(gene = genes_show), file.path(outdir, "Step5_genes_used_for_comparable_oncoplots.csv"))

## Color settings for maftools plots
vc_cols <- c(
  Missense_Mutation = "#4DBBD5FF",
  Nonsense_Mutation = "#E64B35FF",
  Frame_Shift_Del = "#3C5488FF",
  Frame_Shift_Ins = "#00A087FF",
  In_Frame_Del = "#F39B7FFF",
  In_Frame_Ins = "#8491B4FF",
  Splice_Site = "#7E6148FF",
  Translation_Start_Site = "#B09C85FF",
  Nonstop_Mutation = "#91D1C2FF"
)

anno_cols <- list(
  risk_group = c("High risk" = "#D73027", "Low risk" = "#4575B4"),
  Gender = c("Gender_1" = "#4DBBD5FF", "Gender_2" = "#E64B35FF"),
  Age = c("Age_1" = "#00A087FF", "Age_2" = "#F39B7FFF", "Age_3" = "#7E6148FF"),
  ISS = c("ISS_1" = "#91D1C2FF", "ISS_2" = "#8491B4FF", "ISS_3" = "#E64B35FF"),
  Treatment = c("Treatment_0" = "#B09C85FF", "Treatment_1" = "#3C5488FF")
)

clinical_features <- c("Gender", "Age", "ISS", "Treatment")

## =========================
## 7. Waterfall / oncoplot figures
## =========================

open_pdf(file.path(outdir, "A_Oncoplot_HighRisk_common_genes.pdf"), width = 12.5, height = 8.5)
maftools::oncoplot(
  maf = maf_high,
  genes = genes_show,
  colors = vc_cols,
  clinicalFeatures = clinical_features,
  annotationColor = anno_cols,
  sortByAnnotation = TRUE,
  draw_titv = FALSE,
  showTumorSampleBarcodes = FALSE,
  titleText = paste0("High-risk group (n = ", length(high_samples), ")")
)
grDevices::dev.off()

open_pdf(file.path(outdir, "B_Oncoplot_LowRisk_common_genes.pdf"), width = 12.5, height = 8.5)
maftools::oncoplot(
  maf = maf_low,
  genes = genes_show,
  colors = vc_cols,
  clinicalFeatures = clinical_features,
  annotationColor = anno_cols,
  sortByAnnotation = TRUE,
  draw_titv = FALSE,
  showTumorSampleBarcodes = FALSE,
  titleText = paste0("Low-risk group (n = ", length(low_samples), ")")
)
grDevices::dev.off()

open_pdf(file.path(outdir, "C_Oncoplot_AllMatchedRiskGroups_common_genes.pdf"), width = 14, height = 9)
maftools::oncoplot(
  maf = maf_obj,
  genes = genes_show,
  colors = vc_cols,
  clinicalFeatures = c("risk_group", clinical_features),
  annotationColor = anno_cols,
  sortByAnnotation = TRUE,
  draw_titv = FALSE,
  showTumorSampleBarcodes = FALSE,
  titleText = "Matched CoMMpass WXS cohort stratified by State3-RSF risk group"
)
grDevices::dev.off()

open_pdf(file.path(outdir, "D_MAF_summary_matched_samples.pdf"), width = 9, height = 7)
maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
grDevices::dev.off()

## =========================
## 8. Mutation burden comparison
## =========================

mutation_burden <- maf_matched %>%
  dplyr::distinct(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Hugo_Symbol, Reference_Allele, Tumor_Seq_Allele2) %>%
  dplyr::count(Tumor_Sample_Barcode, name = "n_nonsilent_mutations") %>%
  dplyr::right_join(clinical_for_maf, by = "Tumor_Sample_Barcode") %>%
  dplyr::mutate(
    n_nonsilent_mutations = tidyr::replace_na(n_nonsilent_mutations, 0L),
    nonsilent_mutations_per_mb = n_nonsilent_mutations / exome_size_mb,
    risk_group = factor(risk_group, levels = c("Low risk", "High risk"))
  )

safe_write_csv(mutation_burden, file.path(outdir, "Step5_sample_level_nonsilent_mutation_burden.csv"))

burden_test <- wilcox.test(n_nonsilent_mutations ~ risk_group, data = mutation_burden, exact = FALSE)
burden_stats <- data.frame(
  comparison = "High risk vs Low risk",
  metric = "n_nonsilent_mutations",
  p_value = burden_test$p.value,
  signif = p_to_label(burden_test$p.value),
  n_low = sum(mutation_burden$risk_group == "Low risk"),
  n_high = sum(mutation_burden$risk_group == "High risk"),
  median_low = median(mutation_burden$n_nonsilent_mutations[mutation_burden$risk_group == "Low risk"], na.rm = TRUE),
  median_high = median(mutation_burden$n_nonsilent_mutations[mutation_burden$risk_group == "High risk"], na.rm = TRUE)
)
safe_write_csv(burden_stats, file.path(outdir, "Step5_mutation_burden_wilcox_test.csv"))

p_burden <- ggplot2::ggplot(mutation_burden, ggplot2::aes(x = risk_group, y = n_nonsilent_mutations, fill = risk_group, color = risk_group)) +
  ggplot2::geom_violin(width = 0.75, trim = FALSE, alpha = 0.55, linewidth = 0.25, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.45) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 123), size = 2.1, alpha = 0.85) +
  ggplot2::scale_fill_manual(values = c("Low risk" = "#4575B4", "High risk" = "#D73027"), drop = FALSE) +
  ggplot2::scale_color_manual(values = c("Low risk" = "#4575B4", "High risk" = "#D73027"), drop = FALSE) +
  ggplot2::labs(
    title = "Nonsilent mutation burden by State3-RSF risk group",
    subtitle = paste0("Wilcoxon p = ", signif(burden_test$p.value, 3)),
    x = NULL,
    y = "Nonsilent mutations per sample"
  ) +
  theme_pub(base_size = 13) +
  ggplot2::theme(legend.position = "none")

ggsave_pdf(file.path(outdir, "E_NonsilentMutationBurden_High_vs_Low.pdf"), p_burden, width = 4.8, height = 4.5)

p_tmb <- ggplot2::ggplot(mutation_burden, ggplot2::aes(x = risk_group, y = nonsilent_mutations_per_mb, fill = risk_group, color = risk_group)) +
  ggplot2::geom_violin(width = 0.75, trim = FALSE, alpha = 0.55, linewidth = 0.25, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.45) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 123), size = 2.1, alpha = 0.85) +
  ggplot2::scale_fill_manual(values = c("Low risk" = "#4575B4", "High risk" = "#D73027"), drop = FALSE) +
  ggplot2::scale_color_manual(values = c("Low risk" = "#4575B4", "High risk" = "#D73027"), drop = FALSE) +
  ggplot2::labs(
    title = "Estimated nonsilent mutation burden per Mb",
    subtitle = paste0("Assumed exome size = ", exome_size_mb, " Mb"),
    x = NULL,
    y = "Nonsilent mutations / Mb"
  ) +
  theme_pub(base_size = 13) +
  ggplot2::theme(legend.position = "none")

ggsave_pdf(file.path(outdir, "F_EstimatedNonsilentMutationBurdenPerMb_High_vs_Low.pdf"), p_tmb, width = 4.8, height = 4.5)

## =========================
## 9. Differentially mutated genes: Fisher exact tests
## =========================

sample_info <- clinical_for_maf %>%
  dplyr::select(Tumor_Sample_Barcode, risk_group, riskscore) %>%
  dplyr::mutate(risk_group = factor(risk_group, levels = c("Low risk", "High risk")))

mut_binary <- maf_matched %>%
  dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  dplyr::mutate(mutated = 1L)

all_samples <- unique(sample_info$Tumor_Sample_Barcode)
all_genes <- sort(unique(mut_binary$Hugo_Symbol))

fisher_res <- lapply(all_genes, function(g) {
  sub <- sample_info %>%
    dplyr::left_join(
      mut_binary %>% dplyr::filter(Hugo_Symbol == g) %>% dplyr::select(Tumor_Sample_Barcode, mutated),
      by = "Tumor_Sample_Barcode"
    ) %>%
    dplyr::mutate(mutated = tidyr::replace_na(mutated, 0L))

  high_mut <- sum(sub$risk_group == "High risk" & sub$mutated == 1L, na.rm = TRUE)
  high_wt  <- sum(sub$risk_group == "High risk" & sub$mutated == 0L, na.rm = TRUE)
  low_mut  <- sum(sub$risk_group == "Low risk" & sub$mutated == 1L, na.rm = TRUE)
  low_wt   <- sum(sub$risk_group == "Low risk" & sub$mutated == 0L, na.rm = TRUE)

  mat <- matrix(c(high_mut, high_wt, low_mut, low_wt), nrow = 2, byrow = TRUE)
  ft <- suppressWarnings(fisher.test(mat))

  data.frame(
    Hugo_Symbol = g,
    high_mut = high_mut,
    high_total = high_mut + high_wt,
    high_freq = high_mut / (high_mut + high_wt),
    low_mut = low_mut,
    low_total = low_mut + low_wt,
    low_freq = low_mut / (low_mut + low_wt),
    freq_diff_high_minus_low = high_mut / (high_mut + high_wt) - low_mut / (low_mut + low_wt),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    fdr = p.adjust(p_value, method = "BH"),
    direction = dplyr::case_when(
      freq_diff_high_minus_low > 0 ~ "Higher in High risk",
      freq_diff_high_minus_low < 0 ~ "Higher in Low risk",
      TRUE ~ "No difference"
    ),
    total_mut = high_mut + low_mut,
    signif = p_to_label(fdr)
  ) %>%
  dplyr::arrange(p_value, dplyr::desc(abs(freq_diff_high_minus_low)))

safe_write_csv(fisher_res, file.path(outdir, "Step5_differentially_mutated_genes_Fisher_exact_all_genes.csv"))

## Top genes for frequency barplot: include frequent and differential genes
top_freq_genes <- fisher_res %>%
  dplyr::filter(total_mut >= 2) %>%
  dplyr::arrange(dplyr::desc(total_mut), p_value) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(Hugo_Symbol)

top_diff_genes <- fisher_res %>%
  dplyr::filter(total_mut >= 2) %>%
  dplyr::arrange(p_value, dplyr::desc(abs(freq_diff_high_minus_low))) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(Hugo_Symbol)

plot_genes_freq <- unique(c(top_diff_genes, top_freq_genes, intersect(mm_driver_genes, fisher_res$Hugo_Symbol)))
plot_genes_freq <- plot_genes_freq[plot_genes_freq %in% fisher_res$Hugo_Symbol]
plot_genes_freq <- head(plot_genes_freq, 35)

freq_long <- fisher_res %>%
  dplyr::filter(Hugo_Symbol %in% plot_genes_freq) %>%
  dplyr::select(Hugo_Symbol, high_freq, low_freq, high_mut, low_mut, p_value, fdr) %>%
  tidyr::pivot_longer(cols = c(high_freq, low_freq), names_to = "group", values_to = "frequency") %>%
  dplyr::mutate(
    risk_group = dplyr::recode(group, high_freq = "High risk", low_freq = "Low risk"),
    risk_group = factor(risk_group, levels = c("Low risk", "High risk")),
    Hugo_Symbol = factor(Hugo_Symbol, levels = rev(plot_genes_freq))
  )

p_freq <- ggplot2::ggplot(freq_long, ggplot2::aes(x = Hugo_Symbol, y = frequency, fill = risk_group)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.7, alpha = 0.95) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c("Low risk" = "#4575B4", "High risk" = "#D73027"), drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::labs(
    title = "Mutation frequency of recurrent genes by risk group",
    subtitle = "Non-silent coding mutations; genes selected by recurrence and differential frequency",
    x = NULL,
    y = "Mutation frequency",
    fill = "Risk group"
  ) +
  theme_pub(base_size = 12) +
  ggplot2::theme(legend.position = "top")

ggsave_pdf(file.path(outdir, "G_RecurrentGeneMutationFrequency_High_vs_Low.pdf"), p_freq, width = 7.2, height = 8.5)

## Volcano-like differential plot
volcano_df <- fisher_res %>%
  dplyr::filter(total_mut >= 2) %>%
  dplyr::mutate(
    neg_log10_p = -log10(p_value + 1e-300),
    label = ifelse(Hugo_Symbol %in% head(top_diff_genes, 12) | Hugo_Symbol %in% intersect(mm_driver_genes, Hugo_Symbol), Hugo_Symbol, NA_character_),
    direction = factor(direction, levels = c("Higher in Low risk", "No difference", "Higher in High risk"))
  )

p_volcano <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = freq_diff_high_minus_low, y = neg_log10_p, color = direction)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linewidth = 0.35, linetype = "dashed", color = "grey50") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(size = 2.0, alpha = 0.85) +
  ggplot2::scale_color_manual(values = c("Higher in Low risk" = "#4575B4", "No difference" = "grey70", "Higher in High risk" = "#D73027"), drop = FALSE) +
  ggplot2::labs(
    title = "Differentially mutated genes by risk group",
    subtitle = "Fisher's exact test; x-axis shows High-risk minus Low-risk mutation frequency",
    x = "Mutation frequency difference: High risk - Low risk",
    y = "-log10(P value)",
    color = NULL
  ) +
  theme_pub(base_size = 12) +
  ggplot2::theme(legend.position = "top")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = volcano_df %>% dplyr::filter(!is.na(label)),
      ggplot2::aes(label = label),
      size = 3,
      min.segment.length = 0,
      box.padding = 0.35,
      max.overlaps = 60,
      show.legend = FALSE
    )
}

ggsave_pdf(file.path(outdir, "H_DifferentiallyMutatedGenes_Volcano.pdf"), p_volcano, width = 6.4, height = 5.2)

## maftools mafCompare + forestPlot
maf_compare_res <- tryCatch({
  maftools::mafCompare(
    m1 = maf_high,
    m2 = maf_low,
    m1Name = "High risk",
    m2Name = "Low risk",
    minMut = min_mut_for_maf_compare
  )
}, error = function(e) {
  message("mafCompare failed: ", e$message)
  NULL
})

if (!is.null(maf_compare_res)) {
  safe_saveRDS(maf_compare_res, file.path(outdir, "Step5_maftools_mafCompare_high_vs_low.rds"))
  if (!is.null(maf_compare_res$results)) {
    safe_write_csv(maf_compare_res$results, file.path(outdir, "Step5_maftools_mafCompare_high_vs_low_results.csv"))
  }

  open_pdf(file.path(outdir, "I_maftools_ForestPlot_DifferentiallyMutatedGenes.pdf"), width = 7, height = 6)
  tryCatch({
    maftools::forestPlot(
      mafCompareRes = maf_compare_res,
      pVal = 0.2,
      color = c("#D73027", "#4575B4")
    )
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("forestPlot failed:", e$message), cex = 0.8)
  })
  grDevices::dev.off()
}

## =========================
## 10. Risk score by key mutation status
## =========================

key_genes <- intersect(c("TP53", "NRAS", "KRAS", "BRAF", "DIS3", "FAM46C", "TRAF3", "CYLD"), fisher_res$Hugo_Symbol)

if (length(key_genes) > 0) {
  key_status <- tidyr::expand_grid(
    Tumor_Sample_Barcode = sample_info$Tumor_Sample_Barcode,
    Hugo_Symbol = key_genes
  ) %>%
    dplyr::left_join(mut_binary %>% dplyr::filter(Hugo_Symbol %in% key_genes), by = c("Tumor_Sample_Barcode", "Hugo_Symbol")) %>%
    dplyr::mutate(mutated = tidyr::replace_na(mutated, 0L)) %>%
    dplyr::left_join(sample_info, by = "Tumor_Sample_Barcode") %>%
    dplyr::mutate(
      mutation_status = ifelse(mutated == 1L, "Mutant", "Wildtype"),
      mutation_status = factor(mutation_status, levels = c("Wildtype", "Mutant")),
      Hugo_Symbol = factor(Hugo_Symbol, levels = key_genes)
    )

  safe_write_csv(key_status, file.path(outdir, "Step5_key_driver_mutation_status_by_sample.csv"))

  key_tests <- key_status %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::group_modify(~ {
      if (length(unique(.x$mutation_status)) < 2 || sum(.x$mutation_status == "Mutant") < 2) {
        return(data.frame(p_value = NA_real_, n_mutant = sum(.x$mutation_status == "Mutant"), n_wildtype = sum(.x$mutation_status == "Wildtype")))
      }
      wt <- wilcox.test(riskscore ~ mutation_status, data = .x, exact = FALSE)
      data.frame(p_value = wt$p.value, n_mutant = sum(.x$mutation_status == "Mutant"), n_wildtype = sum(.x$mutation_status == "Wildtype"))
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH"), signif = p_to_label(fdr))

  safe_write_csv(key_tests, file.path(outdir, "Step5_riskScore_by_key_driver_mutation_status_wilcox.csv"))

  p_key <- ggplot2::ggplot(key_status, ggplot2::aes(x = mutation_status, y = riskscore, fill = mutation_status, color = mutation_status)) +
    ggplot2::geom_violin(width = 0.75, trim = FALSE, alpha = 0.50, linewidth = 0.25, color = NA) +
    ggplot2::geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 123), size = 1.4, alpha = 0.65) +
    ggplot2::facet_wrap(~ Hugo_Symbol, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = c("Wildtype" = "grey75", "Mutant" = "#D73027"), drop = FALSE) +
    ggplot2::scale_color_manual(values = c("Wildtype" = "grey50", "Mutant" = "#D73027"), drop = FALSE) +
    ggplot2::labs(
      title = "State3-RSF risk score by recurrent driver mutation status",
      x = NULL,
      y = "State3-RSF risk score"
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_rect(fill = "grey95", color = "black", linewidth = 0.35), strip.text = ggplot2::element_text(face = "bold"))

  ggsave_pdf(file.path(outdir, "J_RiskScore_by_KeyDriverMutationStatus.pdf"), p_key, width = 8.5, height = 6.4)
}

## =========================
## 11. Final integrated summary table
## =========================

summary_table <- data.frame(
  item = c(
    "Mutation input file",
    "Risk input file",
    "Selected matching strategy",
    "Matched total samples",
    "Matched High-risk samples",
    "Matched Low-risk samples",
    "Matched nonsilent variants",
    "Matched nonsilent mutated genes",
    "Oncoplot genes displayed",
    "Mutation burden Wilcoxon p"
  ),
  value = c(
    mut_file,
    risk_file,
    selected_match,
    as.character(nrow(clinical_for_maf)),
    as.character(length(high_samples)),
    as.character(length(low_samples)),
    as.character(nrow(maf_matched)),
    as.character(length(unique(maf_matched$Hugo_Symbol))),
    as.character(length(genes_show)),
    signif(burden_test$p.value, 4)
  )
)
safe_write_csv(summary_table, file.path(outdir, "Step5_final_analysis_summary.csv"))

message("\nStep5 somatic mutation analysis completed.")
message("Output directory: ", outdir)
message("Key files:")
message("  A_Oncoplot_HighRisk_common_genes.pdf")
message("  B_Oncoplot_LowRisk_common_genes.pdf")
message("  C_Oncoplot_AllMatchedRiskGroups_common_genes.pdf")
message("  E_NonsilentMutationBurden_High_vs_Low.pdf")
message("  G_RecurrentGeneMutationFrequency_High_vs_Low.pdf")
message("  H_DifferentiallyMutatedGenes_Volcano.pdf")
message("  I_maftools_ForestPlot_DifferentiallyMutatedGenes.pdf")
message("  J_RiskScore_by_KeyDriverMutationStatus.pdf")

