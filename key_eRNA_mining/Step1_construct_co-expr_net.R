nc_file <- "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNAFeature_expression_tumor.csv"
all_file <- "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/TF_allFeature_expression_tumor.csv"
out_file <- "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_target_TumorGroup.csv"

-------------------
strip_version <- TRUE 
remove_self <- TRUE 

# -----------------------------
# Core function (matrix version)
# -----------------------------
build_ncRNA_targets_network <- function(ncRNA_expr, all_expr, remove_self = TRUE) {
   
   ncRNA_expr <- as.matrix(ncRNA_expr)
   all_expr   <- as.matrix(all_expr)
   
   # 1) Check sample columns
   if (!identical(colnames(ncRNA_expr), colnames(all_expr))) {
      stop("Sample columns are not identical (or not in the same order) between ncRNA_expr and all_expr.")
   }
   
   # 2) Remove zero-variance rows (cannot compute correlation)
   ncRNA_expr <- ncRNA_expr[apply(ncRNA_expr, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
   all_expr   <- all_expr[apply(all_expr,   1, sd, na.rm = TRUE) > 0, , drop = FALSE]
   
   n <- ncol(ncRNA_expr)
   if (n < 4) stop("Need at least 4 samples for Pearson correlation test.")
   
   # 3) Pearson r matrix: ncRNA x all
   r_mat <- cor(t(ncRNA_expr), t(all_expr),
                method = "pearson", use = "pairwise.complete.obs")
   
   # 4) Vectorized p-values from r
   df <- n - 2
   r_vec <- as.vector(r_mat)
   t_vec <- r_vec * sqrt(df / pmax(1 - r_vec^2, .Machine$double.eps))
   p_vec <- 2 * pt(-abs(t_vec), df = df)
   
   # reshape p-values into matrix
   p_mat <- matrix(p_vec,
                   nrow = nrow(r_mat), ncol = ncol(r_mat),
                   byrow = FALSE, dimnames = dimnames(r_mat))
   
   # 5) FDR per-ncRNA (row-wise BH)
   fdr_mat <- t(apply(p_mat, 1, p.adjust, method = "BH"))
   
   # 6) Edge list (ALL pairs)
   edges <- data.frame(
      ncRNA  = rep(rownames(r_mat), times = ncol(r_mat)),
      target = rep(colnames(r_mat), each  = nrow(r_mat)),
      r      = r_vec,
      p      = as.vector(p_mat),
      FDR    = as.vector(fdr_mat),
      stringsAsFactors = FALSE
   )
   
   # Optional: remove self edges (if ncRNA ids appear in all_expr)
   if (remove_self) edges <- edges[edges$ncRNA != edges$target, ]
   
   return(edges)
}


# -----------------------------
# Read input
# -----------------------------
nc_expr <- read.table(nc_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
all_expr <- read.table(all_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# all_expr <- as.data.frame(all_expr)
# rownames(all_expr) <- make.unique(all_expr[[1]])
# all_expr <- all_expr[, -1]

# Optional: strip Ensembl version suffix (.1, .2 ...)
if (strip_version) {
   rownames(nc_expr)  <- gsub("\\.\\d+$", "", rownames(nc_expr))
   rownames(all_expr) <- gsub("\\.\\d+$", "", rownames(all_expr))
}

# -----------------------------
# Build network
# -----------------------------
edges <- build_ncRNA_targets_network(nc_expr, all_expr, remove_self = remove_self)

# -----------------------------
# Output
# -----------------------------
write.table(edges, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(edges, file = out_file, row.names = FALSE)
cat("Done!\n")
cat("ncRNA:", nrow(nc_expr), " rows (after sd>0 filter:", length(unique(edges$ncRNA)), ")\n")
cat("targets:", nrow(all_expr), " rows (after sd>0 filter:", length(unique(edges$target)), ")\n")
cat("edges:", nrow(edges), "\n")
cat("Output:", out_file, "\n")
