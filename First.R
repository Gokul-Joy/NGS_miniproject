# ==== Setup & install packages (run once) ====
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("biomaRt","tidyverse","GEOquery","limma","data.table")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) BiocManager::install(to_install, ask=FALSE)

library(biomaRt)
library(tidyverse)
library(GEOquery)
library(limma)
library(data.table)

# ==== 1) Your UniProt list ====
# Read the file (adjust the path to where your file actually is)
uniprot_list <- readLines("D:/MSC/NGS_mini/Galaxy231-[overlapping_uniprot_ids.txt].txt")

# Check what got loaded
head(uniprot_list)


# ==== 2) Map UniProt → HGNC gene symbols via biomaRt ====
# ==== 2) Map UniProt → HGNC gene symbols via biomaRt ====
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query mapping using correct UniProt attribute
mapping <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = uniprot_list,
  mart = ensembl
)

# View results
head(mapping)

# Check mapping
message("UniProt → Gene mapping:")
print(mapping)

# Save mapped secretome file
write_csv(mapping, "D:/MSC/NGS_mini/secretome_genes.csv")
message("Saved mapped secretome to secretome_genes.csv")

# ==== 3) Download GSE109233 expression data from GEO ====
message("Downloading GSE109233 from GEO...")
gset <- getGEO("GSE109233", GSEMatrix=TRUE)
if(length(gset) > 1) gset <- gset[[1]]
expr <- exprs(gset)
pheno <- pData(gset)
platform <- annotation(gset)
message("Downloaded GSE: ", "GSE109233", " platform: ", platform, " dims: ", paste(dim(expr), collapse=" x "))

# ==== 4) Normalize / log2 transform if needed ====
q <- quantile(expr, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1))
if(q[5] > 50) {
  message("Log2-transforming expression matrix...")
  expr <- log2(expr + 1)
}
expr <- normalizeBetweenArrays(expr, method = "quantile")

# ==== 5) Map Illumina probes → gene symbols if needed ====
# Check if rownames are ILMN IDs
if(any(grepl("^ILMN", rownames(expr)))) {
  message("Mapping Illumina probes to gene symbols via platform table...")
  platform_table <- Table(getGEO(platform))
  possible_cols <- c("Symbol","GENE_SYMBOL","GeneSymbol","gene_assignment")
  col_found <- intersect(possible_cols, colnames(platform_table))
  if(length(col_found) == 0) stop("Could not find gene symbol column in platform table.")
  col <- col_found[1]
  probe_map <- platform_table %>% select(ID = ID, symbol = !!sym(col)) %>% filter(symbol != "")
  
  # subset expression matrix
  common_probes <- intersect(rownames(expr), probe_map$ID)
  expr_sub <- expr[common_probes, , drop=FALSE]
  probe_map_sub <- probe_map %>% filter(ID %in% common_probes)
  
  # collapse multiple probes per gene (take probe with highest mean expression)
  probe_map_sub$mean_expr <- rowMeans(expr_sub[probe_map_sub$ID, , drop=FALSE])
  probe_map_sub <- probe_map_sub %>% group_by(symbol) %>% slice_max(order_by=mean_expr, n=1)
  expr_gene <- expr_sub[probe_map_sub$ID, , drop=FALSE]
  rownames(expr_gene) <- probe_map_sub$symbol
  expr <- expr_gene
  message("Collapsed to ", nrow(expr), " unique gene symbols.")
} else {
  message("Expression matrix already uses gene symbols as rownames.")
}

# ==== 6) Save prepared expression matrix ====
write_csv(as.data.frame(expr) %>% rownames_to_column(var="Gene"), "GSE109233_expression_prepared.csv")
message("Prepared expression matrix saved as GSE109233_expression_prepared.csv")

# ==== 7) Check overlap with secretome genes ====
secretome_genes <- mapping$hgnc_symbol
present_genes <- intersect(secretome_genes, rownames(expr))
missing_genes <- setdiff(secretome_genes, rownames(expr))
message("Secretome genes present in expression data: ", length(present_genes))
message("Secretome genes missing in expression data: ", length(missing_genes))
if(length(missing_genes) > 0) message("Missing genes: ", paste(missing_genes, collapse=", "))
