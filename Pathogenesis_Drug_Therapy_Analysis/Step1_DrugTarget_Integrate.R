rm(list = ls())
gc()

suppressPackageStartupMessages({
   library(httr2)
   library(jsonlite)
   library(dplyr)
   library(purrr)
   library(stringr)
   library(tidyr)
   library(tibble)
   library(readr)
})

dir.create("./Outdata/18.Drug_target_integrated", recursive = TRUE, showWarnings = FALSE)

# ===================== 1. List of multiple myeloma (MM) therapeutic drugs =====================
mm_drugs <- c(
   "bortezomib",
   "carfilzomib",
   "ixazomib",
   "lenalidomide",
   "pomalidomide",
   "thalidomide",
   "daratumumab",
   "isatuximab",
   "elotuzumab",
   "selinexor",
   "teclistamab",
   "elranatamab",
   "talquetamab",
   "belantamab mafodotin",
   "idecabtagene vicleucel",
   "ciltacabtagene autoleucel",
   "melphalan",
   "cyclophosphamide",
   "carmustine",
   "doxorubicin",
   "pamidronate",
   "zoledronic acid",
   "plerixafor"
)

# ===================== 2. Basic API & file configuration =====================
CHEMBL_BASE <- "https://www.ebi.ac.uk/chembl/api/data"
OPENTARGETS_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

# DGIdb official TSV file is recommended for supplementary data retrieval.
# Replace this path with your local file or a downloadable public TSV URL.
# The official DGIdb download page provides the full drug-gene interaction claims TSV file.
DGIDB_TSV <- "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/lw_Method/interactions.tsv"

# ===================== 3. General helper functions =====================
`%||%` <- function(a, b) {
   if (is.null(a) || length(a) == 0) b else a
}

safe_req_json <- function(req, max_tries = 3, sleep_sec = 1) {
   last_err <- NULL
   for (i in seq_len(max_tries)) {
      out <- tryCatch(
         req |> req_perform() |> resp_body_json(simplifyVector = FALSE),
         error = function(e) {
            last_err <<- e
            NULL
         }
      )
      if (!is.null(out)) return(out)
      Sys.sleep(sleep_sec * i)
   }
   warning("Request failed after retries: ", conditionMessage(last_err))
   NULL
}

clean_text <- function(x) {
   x |>
      as.character() |>
      str_squish()
}

clean_drug_name <- function(x) {
   clean_text(x) |>
      str_to_lower()
}

is_human <- function(x) {
   is.na(x) | str_to_lower(x) %in% c("homo sapiens", "human", "")
}

coalesce_chr <- function(...) {
   x <- list(...)
   for (i in seq_along(x)) {
      xi <- x[[i]]
      if (!is.null(xi) && length(xi) > 0 && !all(is.na(xi))) return(xi)
   }
   NA_character_
}

# Extract target column from dataframe/list by matching candidate column names
grab_col <- function(df, candidates) {
   if (is.null(df)) return(NULL)
   nms <- names(df)
   hit <- nms[tolower(nms) %in% tolower(candidates)]
   if (length(hit) == 0) return(NULL)
   df[[hit[1]]]
}

# ===================== 4. ChEMBL API query pipeline =====================

chembl_search_drug <- function(drug_name) {
   url <- paste0(
      CHEMBL_BASE, "/molecule/search?q=",
      URLencode(drug_name, reserved = TRUE),
      "&format=json"
   )
   
   x <- safe_req_json(
      request(url) |> req_headers(Accept = "application/json")
   )
   
   if (is.null(x) || is.null(x$molecules) || nrow(as.data.frame(x$molecules)) == 0) {
      return(tibble(
         input_drug = drug_name,
         chembl_pref_name = NA_character_,
         molecule_chembl_id = NA_character_
      ))
   }
   
   mol <- as_tibble(x$molecules)
   
   pref <- if ("pref_name" %in% names(mol)) mol$pref_name else rep(NA_character_, nrow(mol))
   idx <- which(str_to_upper(pref %||% "") == str_to_upper(drug_name))
   if (length(idx) >= 1) {
      mol <- mol[idx[1], , drop = FALSE]
   } else {
      mol <- mol[1, , drop = FALSE]
   }
   
   tibble(
      input_drug = drug_name,
      chembl_pref_name = mol$pref_name %||% NA_character_,
      molecule_chembl_id = mol$molecule_chembl_id %||% NA_character_
   )
}

chembl_get_mechanisms <- function(molecule_chembl_id) {
   if (is.na(molecule_chembl_id) || molecule_chembl_id == "") {
      return(tibble(
         action_type = NA_character_,
         mechanism_of_action = NA_character_,
         target_chembl_id = NA_character_,
         target_pref_name = NA_character_
      ))
   }
   
   url <- paste0(
      CHEMBL_BASE, "/mechanism.json?molecule_chembl_id=",
      molecule_chembl_id
   )
   
   x <- safe_req_json(
      request(url) |> req_headers(Accept = "application/json")
   )
   
   if (is.null(x) || is.null(x$mechanisms) || nrow(as.data.frame(x$mechanisms)) == 0) {
      return(tibble(
         action_type = NA_character_,
         mechanism_of_action = NA_character_,
         target_chembl_id = NA_character_,
         target_pref_name = NA_character_
      ))
   }
   
   mech <- as_tibble(x$mechanisms)
   
   tibble(
      action_type = mech$action_type %||% NA_character_,
      mechanism_of_action = mech$mechanism_of_action %||% NA_character_,
      target_chembl_id = mech$target_chembl_id %||% NA_character_,
      target_pref_name = mech$target_pref_name %||% NA_character_
   )
}

chembl_get_target_genes <- function(target_chembl_id) {
   if (is.na(target_chembl_id) || target_chembl_id == "") {
      return(tibble(
         target_gene = NA_character_,
         target_component_description = NA_character_,
         target_organism = NA_character_
      ))
   }
   
   url <- paste0(CHEMBL_BASE, "/target/", target_chembl_id, ".json")
   
   x <- safe_req_json(
      request(url) |> req_headers(Accept = "application/json")
   )
   
   if (is.null(x)) {
      return(tibble(
         target_gene = NA_character_,
         target_component_description = NA_character_,
         target_organism = NA_character_
      ))
   }
   
   comps <- x$target_components
   
   if (is.null(comps) || nrow(as.data.frame(comps)) == 0) {
      return(tibble(
         target_gene = NA_character_,
         target_component_description = x$pref_name %||% NA_character_,
         target_organism = x$organism %||% NA_character_
      ))
   }
   
   bind_rows(lapply(seq_len(nrow(as.data.frame(comps))), function(i) {
      comp <- comps[i, , drop = FALSE]
      
      syns <- comp$target_component_synonyms[[1]]
      gene_symbol <- NA_character_
      
      if (!is.null(syns) && nrow(as.data.frame(syns)) > 0) {
         syns_df <- as_tibble(syns)
         idx <- which(syns_df$syn_type == "GENE_SYMBOL")
         if (length(idx) >= 1) {
            gene_symbol <- syns_df$component_synonym[idx[1]]
         }
      }
      
      tibble(
         target_gene = gene_symbol,
         target_component_description = comp$component_description %||% NA_character_,
         target_organism = x$organism %||% NA_character_
      )
   }))
}

get_targets_from_chembl <- function(drug_name) {
   drug_tbl <- chembl_search_drug(drug_name)
   
   mech_tbl <- drug_tbl |>
      mutate(mech = map(molecule_chembl_id, chembl_get_mechanisms)) |>
      unnest(mech)
   
   target_tbl <- mech_tbl |>
      mutate(target_info = map(target_chembl_id, chembl_get_target_genes)) |>
      unnest(target_info)
   
   out <- target_tbl |>
      filter(is_human(target_organism)) |>
      transmute(
         input_drug = drug_name,
         matched_drug = chembl_pref_name,
         source = "ChEMBL",
         target_gene = na_if(clean_text(target_gene), ""),
         target_name = coalesce_chr(target_pref_name, target_component_description),
         evidence = coalesce_chr(action_type, mechanism_of_action),
         source_id = target_chembl_id,
         organism = target_organism
      ) |>
      distinct()
   
   if (nrow(out) == 0) {
      out <- tibble(
         input_drug = drug_name,
         matched_drug = NA_character_,
         source = "ChEMBL",
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = NA_character_,
         source_id = NA_character_,
         organism = NA_character_
      )
   }
   
   out
}

# ===================== 5. DGIdb (Local TSV file parsing recommended) =====================

# Notes:
# The most stable DGIdb retrieval method is reading the official downloadable TSV file.
# Column names vary across different DGIdb versions; automatic column matching is implemented below.
read_dgidb_interactions <- function(dgidb_tsv) {
   if (!file.exists(dgidb_tsv)) {
      warning("DGIdb TSV not found: ", dgidb_tsv)
      return(NULL)
   }
   
   x <- read_tsv(dgidb_tsv, show_col_types = FALSE, progress = FALSE)
   nms <- names(x)
   
   find_col <- function(cands) {
      hit <- nms[tolower(nms) %in% tolower(cands)]
      if (length(hit) == 0) return(NA_character_)
      hit[1]
   }
   
   drug_col <- find_col(c("drug_name", "drug_claim_name", "drug", "drug_concept_name"))
   gene_col <- find_col(c("gene_name", "gene_claim_name", "gene", "entrez_gene_symbol"))
   source_col <- find_col(c("interaction_claim_source", "source_db_name", "source", "sources"))
   type_col <- find_col(c("interaction_types", "interaction_type", "interaction_claim_type"))
   score_col <- find_col(c("interaction_score", "score"))
   pmid_col <- find_col(c("pmids", "pubmed_ids", "pmid"))
   approval_col <- find_col(c("drug_approved", "approved"))
   
   if (is.na(drug_col) || is.na(gene_col)) {
      warning("DGIdb TSV columns not recognized.")
      return(NULL)
   }
   
   x |>
      transmute(
         drug_name_raw = .data[[drug_col]],
         target_gene = .data[[gene_col]],
         dgidb_source = if (!is.na(source_col)) as.character(.data[[source_col]]) else NA_character_,
         dgidb_interaction_type = if (!is.na(type_col)) as.character(.data[[type_col]]) else NA_character_,
         dgidb_score = if (!is.na(score_col)) as.character(.data[[score_col]]) else NA_character_,
         dgidb_pmids = if (!is.na(pmid_col)) as.character(.data[[pmid_col]]) else NA_character_,
         dgidb_approved = if (!is.na(approval_col)) as.character(.data[[approval_col]]) else NA_character_
      ) |>
      mutate(
         drug_name_std = clean_drug_name(drug_name_raw),
         target_gene = na_if(clean_text(target_gene), "")
      ) |>
      filter(!is.na(drug_name_std), drug_name_std != "", !is.na(target_gene), target_gene != "")
}

get_targets_from_dgidb <- function(drug_name, dgidb_tbl) {
   if (is.null(dgidb_tbl)) {
      return(tibble(
         input_drug = drug_name,
         matched_drug = NA_character_,
         source = "DGIdb",
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = NA_character_,
         source_id = NA_character_,
         organism = NA_character_
      ))
   }
   
   key <- clean_drug_name(drug_name)
   
   hit_exact <- dgidb_tbl |>
      filter(drug_name_std == key)
   
   hit_fuzzy <- if (nrow(hit_exact) == 0) {
      dgidb_tbl |>
         filter(str_detect(drug_name_std, fixed(key)) | str_detect(key, fixed(drug_name_std)))
   } else {
      tibble()
   }
   
   hit <- bind_rows(hit_exact, hit_fuzzy) |>
      distinct()
   
   if (nrow(hit) == 0) {
      return(tibble(
         input_drug = drug_name,
         matched_drug = NA_character_,
         source = "DGIdb",
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = NA_character_,
         source_id = NA_character_,
         organism = NA_character_
      ))
   }
   
   hit |>
      transmute(
         input_drug = drug_name,
         matched_drug = drug_name_raw,
         source = "DGIdb",
         target_gene = target_gene,
         target_name = target_gene,
         evidence = paste(
            na.omit(c(dgidb_interaction_type, dgidb_source, dgidb_score, dgidb_pmids)),
            collapse = " | "
         ),
         source_id = NA_character_,
         organism = "Homo sapiens"
      ) |>
      distinct()
}

# ===================== 6. Open Targets GraphQL API query =====================

ot_post <- function(query, variables = list()) {
   req <- request(OPENTARGETS_URL) |>
      req_method("POST") |>
      req_headers(
         `Content-Type` = "application/json",
         Accept = "application/json"
      ) |>
      req_body_json(list(
         query = query,
         variables = variables
      ), auto_unbox = TRUE)
   
   safe_req_json(req)
}

# First search drug name to retrieve OpenTargets drug ID / ChEMBL ID
ot_search_drug <- function(drug_name) {
   query <- '
  query SearchEntity($q: String!) {
    search(queryString: $q, entityNames: ["drug"]) {
      hits {
        id
        entity
        object {
          ... on Drug {
            id
            name
          }
        }
      }
    }
  }'
   
   x <- ot_post(query, list(q = drug_name))
   hits <- x$data$search$hits
   
   if (is.null(hits) || nrow(as.data.frame(hits)) == 0) {
      return(tibble(
         input_drug = drug_name,
         drug_id = NA_character_,
         matched_drug = NA_character_
      ))
   }
   
   # Prioritize exact name match; take first hit if no exact match
   hits_tbl <- tibble(
      drug_id = map_chr(hits, ~ .x$object$id %||% .x$id %||% NA_character_),
      matched_drug = map_chr(hits, ~ .x$object$name %||% NA_character_)
   )
   
   idx <- which(str_to_upper(hits_tbl$matched_drug) == str_to_upper(drug_name))
   if (length(idx) == 0) idx <- 1
   
   tibble(
      input_drug = drug_name,
      drug_id = hits_tbl$drug_id[idx[1]],
      matched_drug = hits_tbl$matched_drug[idx[1]]
   )
}

# Retrieve mechanism of action and target genes via drug ID
ot_get_drug_targets <- function(drug_id) {
   if (is.na(drug_id) || drug_id == "") {
      return(tibble(
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = NA_character_,
         source_id = NA_character_,
         organism = NA_character_
      ))
   }
   
   query <- '
  query DrugMoA($drugId: String!) {
    drug(drugId: $drugId) {
      id
      name
      mechanismsOfAction {
        rows {
          mechanismOfAction
          actionType
          targets {
            id
            approvedSymbol
            approvedName
          }
        }
      }
    }
  }'
   
   x <- ot_post(query, list(drugId = drug_id))
   
   rows <- x$data$drug$mechanismsOfAction$rows
   if (is.null(rows) || length(rows) == 0) {
      return(tibble(
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = paste(na.omit(c(r$actionType, r$mechanismOfAction)), collapse = " | "),
         source_id = NA_character_,
         organism = "Homo sapiens"
      ))
   }
   
   out <- bind_rows(lapply(rows, function(r) {
      targets <- r$targets
      if (is.null(targets) || length(targets) == 0) {
         return(tibble(
            target_gene = NA_character_,
            target_name = NA_character_,
            evidence = paste(na.omit(c(r$actionType, r$mechanismOfAction)), collapse = " | "),
            source_id = NA_character_,
            organism = "Homo sapiens"
         ))
      }
      
      bind_rows(lapply(targets, function(tg) {
         tibble(
            target_gene = tg$approvedSymbol %||% NA_character_,
            target_name = tg$approvedName %||% NA_character_,
            evidence = paste(na.omit(c(r$actionType, r$mechanismOfAction)), collapse = " | "),
            source_id = tg$id %||% NA_character_,
            organism = "Homo sapiens"
         )
      }))
   })) |>
      distinct()
   
   out
}

get_targets_from_opentargets <- function(drug_name) {
   hit <- ot_search_drug(drug_name)
   tar <- ot_get_drug_targets(hit$drug_id[1])
   
   out <- tar |>
      transmute(
         input_drug = drug_name,
         matched_drug = hit$matched_drug[1],
         source = "OpenTargets",
         target_gene = na_if(clean_text(target_gene), ""),
         target_name = na_if(clean_text(target_name), ""),
         evidence = na_if(clean_text(evidence), ""),
         source_id = na_if(clean_text(source_id), ""),
         organism = organism
      ) |>
      distinct()
   
   if (nrow(out) == 0) {
      out <- tibble(
         input_drug = drug_name,
         matched_drug = hit$matched_drug[1] %||% NA_character_,
         source = "OpenTargets",
         target_gene = NA_character_,
         target_name = NA_character_,
         evidence = NA_character_,
         source_id = NA_character_,
         organism = NA_character_
      )
   }
   
   out
}

# ===================== 7. Run integration pipeline across three databases =====================

message("Reading DGIdb interaction table...")
dgidb_tbl <- read_dgidb_interactions(DGIDB_TSV)

message("Querying target information from ChEMBL...")
chembl_res <- map_dfr(mm_drugs, get_targets_from_chembl)

message("Querying target information from DGIdb...")
dgidb_res <- map_dfr(mm_drugs, ~ get_targets_from_dgidb(.x, dgidb_tbl))

message("Querying target information from Open Targets...")
opentargets_res <- map_dfr(mm_drugs, get_targets_from_opentargets)

all_res_raw <- bind_rows(
   chembl_res,
   dgidb_res,
   opentargets_res
)

# ===================== 8. Data cleaning & multi-source integration =====================
all_res <- all_res_raw |>
   mutate(
      input_drug = clean_text(input_drug),
      matched_drug = na_if(clean_text(matched_drug), ""),
      target_gene = na_if(clean_text(target_gene), ""),
      target_name = na_if(clean_text(target_name), ""),
      evidence = na_if(clean_text(evidence), ""),
      source_id = na_if(clean_text(source_id), ""),
      organism = na_if(clean_text(organism), "")
   ) |>
   filter(is.na(organism) | is_human(organism)) |>
   distinct()

# Retain only records with valid target gene/name annotation
all_hits <- all_res |>
   filter(!is.na(target_gene) | !is.na(target_name))

# Aggregate records from multiple databases into unified integrated table
integrated_drug_target <- all_hits |>
   mutate(target_key = dplyr::coalesce(target_gene, target_name)) |>
   group_by(input_drug, target_key) |>
   summarise(
      target_gene = dplyr::first(na.omit(target_gene), default = NA_character_),
      target_name = dplyr::first(na.omit(target_name), default = NA_character_),
      matched_drug_names = paste(sort(unique(na.omit(matched_drug))), collapse = " ; "),
      sources = paste(sort(unique(source)), collapse = " ; "),
      evidence = paste(sort(unique(na.omit(evidence))), collapse = " || "),
      source_ids = paste(sort(unique(na.omit(source_id))), collapse = " ; "),
      organism = dplyr::first(na.omit(organism), default = "Homo sapiens"),
      n_sources = n_distinct(source),
      .groups = "drop"
   ) |>
   arrange(input_drug, desc(n_sources), target_gene, target_name)

# Unique target gene list extracted from integrated table
integrated_target_genes <- integrated_drug_target |>
   filter(!is.na(target_gene), target_gene != "") |>
   distinct(target_gene) |>
   arrange(target_gene)

# Source-wise summary statistics
source_summary <- all_hits |>
   group_by(source, input_drug) |>
   summarise(
      n_targets = n_distinct(coalesce(target_gene, target_name)),
      .groups = "drop"
   ) |>
   arrange(input_drug, source)

# ===================== 9. Export output tables =====================
write_tsv(all_res, "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/MM_drug_targets_all_sources_raw.tsv")
write_tsv(integrated_drug_target, "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/MM_drug_targets_integrated.tsv")
write_tsv(integrated_target_genes, "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/MM_target_genes_integrated.tsv")
write_tsv(source_summary, "/home/yjliu/mmProj/data_process/Human/Drug_Therapy/MM_drug_targets_source_summary.tsv")

cat("===== Multi-myeloma drug-target integration pipeline completed =====\n")
cat("Total input drugs: ", length(mm_drugs), "\n")
cat("Integrated drug-target interaction pairs: ", nrow(integrated_drug_target), "\n")
cat("Unique target genes after integration: ", nrow(integrated_target_genes), "\n")
cat("Raw unfiltered database records: ", nrow(all_res), "\n")