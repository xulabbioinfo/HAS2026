#!/usr/bin/env Rscript
## custom GSEA
##   Build TERM2GENE: enhancer-linked genes x liver/colon specific genes (+ random control)
##   Run clusterProfiler GSEA on a ranked gene list, plot with GseaVis
##   
## Usage:
##   
##   enhancer.xlsx: enhancer table with gene_anchor and anno columns
##   ranked_genes.tsv: tab-separated <gene> <score> with header (e.g. log2FC), score column first

suppressPackageStartupMessages({
  library(openxlsx)
  library(clusterProfiler)
  library(GseaVis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript enhancer_gsea_pipeline.R <enhancer.xlsx> <colon_specific.csv> <liver_specific.csv> <ranked_genes.tsv> [outdir]")
}
enhancer_file <- args[1]
colon_file <- args[2]
liver_file <- args[3]
rank_file <- args[4]
outdir <- if (length(args) >= 5) args[5] else "."

set.seed(123456)
nperm <- 10000
min_gssize <- 1
max_gssize <- 1e8
pvalue_cutoff <- 1
p_adjust <- "fdr"
random_size <- 1000

dir.create(outdir, showWarnings = FALSE)

target1 <- read.xlsx(enhancer_file)
colon_spe_list <- unique(read.table(colon_file, sep = ",", header = TRUE)$Description)
liver_spe_list <- na.omit(unique(read.table(liver_file, sep = ",", header = TRUE)$Description))

rank_df <- read.table(rank_file, header = TRUE, sep = "\t", row.names = 1)
id <- rank_df[[1]]
names(id) <- rownames(rank_df)
id <- sort(id, decreasing = TRUE)

one_set <- function(term, genes) {
  if (length(genes) == 0) return(data.frame(term = character(), gene = character()))
  data.frame(term = term, gene = genes)
}

build_termm2gene <- function(enhancer, suffix) {
  term2gene <- rbind(
    one_set("enhancer_liver", intersect(enhancer, liver_spe_list)),
    one_set("enhancer_colon", intersect(enhancer, colon_spe_list)),
    one_set("random1000", sample(names(id), random_size, replace = FALSE))
  )
  write.table(term2gene,
              file.path(outdir, paste0("enahncer_overlap_colon_liver", suffix, ".txt")),
              quote = FALSE, row.names = FALSE, sep = "\t")
  term2gene
}

run_gsea <- function(term2gene, suffix) {
  set.seed(123456)
  egmt <- GSEA(id, TERM2GENE = term2gene, verbose = TRUE,
               nPerm = nperm, minGSSize = min_gssize, maxGSSize = max_gssize,
               pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust)
  pdf(file.path(outdir, paste0("gsea_enhancer", suffix, ".pdf")), width = 10, height = 8)
  GseaVis::gseaNb(object = egmt, geneSetID = c("enhancer_liver", "enhancer_colon"), addPval = TRUE)
  dev.off()
}

cat("Part 1: all enhancers\n")
t1 <- build_termm2gene(unique(na.omit(target1$gene_anchor)), "")
run_gsea(t1, "")

if ("anno" %in% colnames(target1)) {
  cat("Part 2: intergenic enhancers only\n")
  t2 <- build_termm2gene(unique(na.omit(target1$gene_anchor[target1$anno == "intergenic"])), "_nointron")
  run_gsea(t2, "_nointron")
}

cat("Done\n")
