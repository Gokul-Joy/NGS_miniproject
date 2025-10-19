

# ========== 1) Libraries ==========
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(igraph)
  library(dorothea)
  library(viper)
  library(decoupleR)   # used as fallback
})

# ========== 2) Load secretome and expression matrix ==========
secretome <- read_csv("secretome_genes.csv", show_col_types = FALSE)
expr <- read_csv("GSE109233_expression_prepared.csv", show_col_types = FALSE) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# ensure numeric
mode(expr) <- "numeric"

message("Genes in expression matrix: ", nrow(expr), "  Samples: ", ncol(expr))

# Subset (if you need later)
secretome_genes <- intersect(secretome$hgnc_symbol, rownames(expr))
message("Secretome genes present in expression: ", length(secretome_genes))

# ========== 3) Load DoRothEA regulons and filter by confidence ==========
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A","B","C"))
message("Rows in dorothea regulons (A/B/C): ", nrow(regulons))

# quick inspect (optional)
# print(names(regulons)); print(head(regulons))

# ========== 4) Convert to VIPER-style regulon list (df2regulon) ==========
# df2regulon should exist in this dorothea version
regulon_obj <- df2regulon(regulons)

# Checks
message("Regulon object class: ", class(regulon_obj))
message("Number of TFs in regulon_obj: ", length(regulon_obj))
if (length(regulon_obj) == 0) stop("regulon_obj is empty — check df2regulon/regulons input.")

# ========== 5) Check overlap between expression and regulon targets ==========
targets_all <- unique(unlist(lapply(regulon_obj, function(x) names(x$tfmode))))
message("Targets in regulon: ", length(targets_all))
message("Genes in expression: ", nrow(expr))
ov <- sum(rownames(expr) %in% targets_all)
message("Overlap (matching gene symbols): ", ov)

# If overlap is unexpectedly small, try simple heuristics to harmonize
min_overlap_frac <- 0.2
if (ov < (length(targets_all) * min_overlap_frac)) {
  message("Low overlap detected -> attempting simple harmonization (toupper).")
  # uppercase rownames and regulon targets
  rownames(expr) <- toupper(rownames(expr))
  regulon_obj <- lapply(regulon_obj, function(x) {
    names(x$tfmode)     <- toupper(names(x$tfmode))
    names(x$likelihood) <- toupper(names(x$likelihood))
    x
  })
  targets_all2 <- unique(unlist(lapply(regulon_obj, function(x) names(x$tfmode))))
  ov2 <- sum(rownames(expr) %in% targets_all2)
  message("Overlap after toupper: ", ov2)
  if (ov2 <= ov) {
    message("toupper did not help much. If overlap still low, consider mapping aliases -> HGNC.")
  } else {
    targets_all <- targets_all2
    ov <- ov2
  }
}

# ========== 6) Run VIPER (viper::viper) on regulon list ==========
viper_activity <- NULL
try({
  message("Running viper::viper(...) with regulon list (fast path).")
  viper_activity <- viper::viper(expr, regulon_obj, verbose = FALSE)
  message("viper() succeeded. Result dim (TFs x samples): ", paste(dim(viper_activity), collapse = " x "))
}, silent = TRUE)

# ========== 7) Fallback: use decoupleR::run_viper if viper() failed ==========
if (is.null(viper_activity)) {
  message("viper() failed or returned NULL. Falling back to decoupleR::run_viper().")
  # convert regulon_obj -> long network df
  network_df <- do.call(rbind, lapply(names(regulon_obj), function(tf) {
    tfm <- regulon_obj[[tf]]$tfmode
    lik <- regulon_obj[[tf]]$likelihood
    data.frame(
      source = tf,
      target = names(tfm),
      .mor = as.numeric(tfm),
      likelihood = as.numeric(lik),
      stringsAsFactors = FALSE
    )
  }))
  # sanitize columns
  network_df <- network_df %>% distinct(source, target, .mor, likelihood)
  message("Network edges (rows): ", nrow(network_df), "  Unique TFs: ", length(unique(network_df$source)))
  # ensure gene names match expr rows (network_df$target)
  ov_net <- sum(rownames(expr) %in% unique(network_df$target))
  message("Overlap with expression (network_df): ", ov_net)
  # run decoupleR (use minsize to avoid tiny regulons)
  viper_res <- tryCatch({
    decoupleR::run_viper(
      mat = expr,
      network = network_df,
      .source = "source",
      .target = "target",
      .mor = ".mor",
      .likelihood = "likelihood",
      minsize = 5,
      verbose = TRUE
    )
  }, error = function(e) {
    message("decoupleR::run_viper error: ", conditionMessage(e))
    NULL
  })
  if (!is.null(viper_res)) {
    # decoupleR returns a data.frame: rows = regulators, columns = samples (may be tibble)
    viper_activity <- as.matrix(viper_res)
    message("decoupleR::run_viper succeeded. Result dim: ", paste(dim(viper_activity), collapse = " x "))
  } else {
    stop("Both viper() and decoupleR::run_viper() failed. Inspect the earlier messages to debug.")
  }
}

# ========== 8) Inspect & save results ==========
# viper_activity should be TFs x samples numeric matrix
stopifnot(!is.null(viper_activity))
message("Final VIPER activity matrix dims: ", paste(dim(viper_activity), collapse = " x "))

# Show a small view
print(dim(viper_activity))
print(head(viper_activity[, 1:min(6, ncol(viper_activity))]))

# Save result to CSV (optional)
write.csv(as.data.frame(viper_activity), file = "viper_activity_matrix.csv", quote = FALSE)

message("Done. 'viper_activity_matrix.csv' written to working directory.")



# ==== 3) Compute TF activity via VIPER ====
tf_activity <- viper_activity

tf_activity_mean <- rowMeans(tf_activity, na.rm=TRUE)

# ==== 4) Test enrichment: TF targets vs secretome genes ====
tf_targets <- regulons %>% split(.$tf)

results <- tibble(TF=character(), TF_activity=numeric(),
                  n_targets=integer(), n_targets_in_secretome=integer(),
                  fisher_p=numeric())

all_genes_universe <- rownames(expr)

for(tf in names(tf_targets)){
  targets <- intersect(tf_targets[[tf]]$target, all_genes_universe)
  n_targets <- length(targets)
  n_in_secret <- length(intersect(targets, secretome_genes))
  # Fisher contingency table
  a <- n_in_secret
  b <- length(secretome_genes) - a
  c <- n_targets - a
  d <- length(all_genes_universe) - (a+b+c)
  mat <- matrix(c(a,b,c,d), nrow=2)
  ft <- fisher.test(mat, alternative = "greater")
  activity <- ifelse(tf %in% names(tf_activity_mean), tf_activity_mean[tf], NA)
  results <- results %>% add_row(TF=tf, TF_activity=activity,
                                 n_targets=n_targets,
                                 n_targets_in_secretome=n_in_secret,
                                 fisher_p=ft$p.value)
}

# Adjust p-values
results <- results %>% arrange(fisher_p) %>% mutate(p_adj = p.adjust(fisher_p, method="BH"))
message("Top TF candidates:")
print(head(results,10))

# ==== 5) Build Cytoscape edge & node tables ====
# Select significant TFs or top activity
sig_tfs <- results %>% filter(p_adj<0.05 | abs(TF_activity) >= quantile(abs(TF_activity),0.9, na.rm=TRUE)) %>% pull(TF)

edges <- tibble(TF=character(), Target=character())
for(tf in sig_tfs){
  targets <- intersect(tf_targets[[tf]]$target, secretome_genes)
  if(length(targets)>0){
    edges <- bind_rows(edges, tibble(TF=tf, Target=targets))
  }
}

# Attach ProteinIntensity if available
if("ProteinIntensity" %in% colnames(secretome)){
  edges <- edges %>% left_join(secretome %>% select(hgnc_symbol, ProteinIntensity),
                               by=c("Target"="hgnc_symbol"))
}

write_csv(edges, "TF_secretome_edges.csv")
message("Edges saved to TF_secretome_edges.csv (", nrow(edges), " edges)")

# Node table
tf_nodes <- tibble(name=unique(edges$TF), type="TF",
                   TF_activity=tf_activity_mean[unique(edges$TF)])
gene_nodes <- tibble(
  name = unique(edges$Target),
  type = "Gene",
  ProteinIntensity = NA,  # fill with NA
  Expression = rowMeans(expr[unique(edges$Target), , drop = FALSE])
)


colnames(secretome)






nodes <- bind_rows(tf_nodes, gene_nodes)
write_csv(nodes, "TF_secretome_nodes.csv")
message("Nodes saved to TF_secretome_nodes.csv (", nrow(nodes), " nodes)")

# ==== 6) Quick barplot of top TF activity ====
top_tf <- results %>% arrange(p_adj) %>% slice(1:10)
library(ggplot2)
ggplot(top_tf, aes(x=reorder(TF,-TF_activity), y=TF_activity)) +
  geom_col(fill="steelblue") + coord_flip() +
  labs(title="Top TF activity (mean across samples)", y="TF activity", x="TF")

