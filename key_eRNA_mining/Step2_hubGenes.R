library(igraph)
eRNA_target <- read.table("/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_target_TumorGroup.txt",header = T)
eRNA_target_filt <- subset(
   eRNA_target,
   FDR < 0.05 & abs(r) > 0.3
)
output_file <- "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/eRNA_target_TumorGroup_Significant.csv"
fwrite(eRNA_target_filt, file = output_file, row.names = FALSE, sep = ",")
eps <- .Machine$double.xmin
eRNA_target_filt$weight <- abs(eRNA_target_filt$r) * 
   (-log10(pmax(eRNA_target_filt$FDR, eps)))

summary(eRNA_target_filt$weight)
head(
   eRNA_target_filt[order(eRNA_target_filt$weight, decreasing = TRUE), ],
   10
)

eRNA_summary <- aggregate(
   cbind(degree = weight, strength = weight) ~ ncRNA,
   data = eRNA_target_filt,
   FUN = function(x) c(length(x), sum(x))
)

eRNA_summary <- data.frame(
   ncRNA = eRNA_summary$ncRNA,
   degree = eRNA_summary$degree[, 1],
   strength = eRNA_summary$strength[, 2]
)

eRNA_summary <- eRNA_summary[order(eRNA_summary$strength, decreasing = TRUE), ]
head(eRNA_summary, 10)

write.csv(eRNA_summary, "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/eRNA_degree_strength_TumorGroup.csv", row.names = FALSE)


stopifnot(all(c("ncRNA","target","r","FDR","weight") %in% colnames(eRNA_target_filt)))

eRNA_target_filt <- eRNA_target_filt[is.finite(eRNA_target_filt$weight) & !is.na(eRNA_target_filt$weight), ]
g <- graph_from_data_frame(
   d = eRNA_target_filt[, c("ncRNA", "target", "weight")],
   directed = FALSE
)

E(g)$weight <- eRNA_target_filt$weight

v_strength <- strength(g, vids = V(g), weights = E(g)$weight)

E(g)$dist <- 1 / pmax(E(g)$weight, .Machine$double.xmin)

v_bet <- betweenness(g, v = V(g), directed = FALSE, weights = E(g)$dist, normalized = TRUE)

v_pr <- page_rank(g, directed = FALSE, weights = E(g)$weight)$vector

eRNA_nodes <- unique(eRNA_target_filt$ncRNA)

metrics_all <- data.frame(
   node = names(v_strength),
   strength = as.numeric(v_strength),
   betweenness = as.numeric(v_bet),
   pagerank = as.numeric(v_pr),
   stringsAsFactors = FALSE
)

metrics_eRNA <- metrics_all[metrics_all$node %in% eRNA_nodes, ]

top_pct <- 0.05

cut_strength <- quantile(metrics_eRNA$strength, 1 - top_pct, na.rm = TRUE)
cut_bet      <- quantile(metrics_eRNA$betweenness, 1 - top_pct, na.rm = TRUE)
cut_pr       <- quantile(metrics_eRNA$pagerank, 1 - top_pct, na.rm = TRUE)

hub_strict <- metrics_eRNA[
   metrics_eRNA$strength   >= cut_strength &
      metrics_eRNA$betweenness>= cut_bet &
      metrics_eRNA$pagerank   >= cut_pr, 
]

hub_strict <- hub_strict[order(hub_strict$strength, decreasing = TRUE), ]
head(hub_strict, 20)

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

metrics_eRNA$z_strength <- z(metrics_eRNA$strength)
metrics_eRNA$z_bet      <- z(metrics_eRNA$betweenness)
metrics_eRNA$z_pr       <- z(metrics_eRNA$pagerank)

metrics_eRNA$hub_score <- metrics_eRNA$z_strength + metrics_eRNA$z_bet + metrics_eRNA$z_pr

hub_soft <- metrics_eRNA[order(metrics_eRNA$hub_score, decreasing = TRUE), ]
head(hub_soft, 20)

cut_score <- quantile(hub_soft$hub_score, 1 - top_pct, na.rm = TRUE)
hub_soft_top <- hub_soft[hub_soft$hub_score >= cut_score, ]

write.csv(metrics_eRNA, "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/eRNA_strength_betweenness_pagerank.csv", row.names = FALSE)
write.csv(hub_strict,   "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/hub_eRNA_strict_intersection.csv", row.names = FALSE)
write.csv(hub_soft_top, "/home/yjliu/mmProj/data_process/Human/Key_ncRNA/eRNA_resultNew/hub_eRNA_soft_hubscore_top.csv", row.names = FALSE)

