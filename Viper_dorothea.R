BiocManager::install("viper")
BiocManager::install("dorothea")


library(tidyverse)
library(dorothea)
library(viper)
library(data.table)
library(igraph)

# ==== 1) Load secretome genes and expression ====
secretome <- read_csv("secretome_genes.csv", show_col_types = FALSE)
expr <- read_csv("GSE109233_expression_prepared.csv", show_col_types = FALSE) %>%
  column_to_rownames("Gene") %>% as.matrix()

# Subset to secretome genes present in expression
secretome_genes <- intersect(secretome$hgnc_symbol, rownames(expr))
message("Secretome genes used in network: ", length(secretome_genes))

# ==== 2) Load DoRothEA regulons (human) ====
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A","B","C"))

# Prepare regulon object for VIPER
regulon_vip <- dorothea::dorothea2viper(regulons)

# ==== 3) Compute TF activity via VIPER ====
tf_activity <- viper(expr, regulon_vip, verbose=FALSE)
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
gene_nodes <- tibble(name=unique(edges$Target), type="Gene",
                     ProteinIntensity=secretome %>% filter(hgnc_symbol %in% unique(edges$Target)) %>% pull(ProteinIntensity),
                     Expression=rowMeans(expr[unique(edges$Target), , drop=FALSE]))
nodes <- bind_rows(tf_nodes, gene_nodes)
write_csv(nodes, "TF_secretome_nodes.csv")
message("Nodes saved to TF_secretome_nodes.csv (", nrow(nodes), " nodes)")

# ==== 6) Quick barplot of top TF activity ====
top_tf <- results %>% arrange(p_adj) %>% slice(1:10)
library(ggplot2)
ggplot(top_tf, aes(x=reorder(TF,-TF_activity), y=TF_activity)) +
  geom_col(fill="steelblue") + coord_flip() +
  labs(title="Top TF activity (mean across samples)", y="TF activity", x="TF")
