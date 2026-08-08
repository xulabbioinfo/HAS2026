#!/usr/bin/env Rscript
## Positive/negative expression analysis for markers of interest
##   Per marker: violin of marker set in positive cells + pairwise wilcox + spearman corr heatmap
##   Across markers: violin of each marker's expression among all positive groups
## 

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(psych)
  library(pheatmap)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript positive_analysis_pipeline.R <input_seurat.Rds> [outdir]")
obj_path <- args[1]
outdir <- if (length(args) >= 2) args[2] else "."

markers <- c("NR1H4", "HHEX", "PROX1", "AFP", "CDX2")
group_cols <- paste0(markers, "_exp")
positive_groups <- paste0(markers, "_exp_positive")
corr_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

dir.create(outdir, showWarnings = FALSE)

violin_theme <- function() {
  cowplot::theme_cowplot(font_family = "Arial") +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x = element_text(color = "black", size = 11),
      axis.title.x = element_text(color = "black", size = 15),
      axis.ticks.x = element_line(color = "black", lineend = "round")
    )
}

make_violin <- function(df, x_col, facet_col, xlab_text) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[["expression"]], fill = .data[[x_col]])) +
    geom_violin(scale = "width", adjust = 1) +
    facet_grid(cols = vars(.data[[facet_col]]), scales = "free") +
    violin_theme() +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category20c_d3")) +
    xlab(xlab_text) +
    ylab("")
}

analyze_gene <- function(obj, gene, outdir) {
  ident_col <- paste0(gene, "_exp")

  d <- FetchData(obj, vars = c(markers, ident_col))
  d <- d[d[[ident_col]] == paste0(gene, "_exp_positive"), ]
  d[[ident_col]] <- NULL

  dm <- melt(d, variable.name = "marker", value.name = "expression")

  p <- make_violin(dm, "marker", "marker", paste0(gene, "_exp_positive Expression Level"))
  ggsave(file.path(outdir, paste0(gene, "_exp_positive_violin.svg")), p, height = 4, width = 12)

  stats <- compare_means(expression ~ marker, data = dm, method = "wilcox.test")
  write.xlsx(stats, file.path(outdir, paste0(gene, "_exp_positive_wilcox.xlsx")))

  res <- corr.test(d[, markers, drop = FALSE], method = "spearman")
  png(file.path(outdir, paste0(gene, "_exp_positive_corr.png")))
  pheatmap(res$r, scale = "none", color = corr_colors, border = FALSE,
           cluster_row = TRUE, cluster_col = TRUE,
           display_numbers = res$stars, fontsize_number = 10)
  dev.off()
}

plot_across_groups <- function(obj, target, outdir) {
  expr <- as.data.frame(t(as.data.frame(GetAssayData(obj, assay = "RNA", slot = "data")[markers, ])))
  meta <- obj@meta.data[, group_cols, drop = FALSE]
  data4group <- cbind(expr, meta)

  parts <- lapply(group_cols, function(g) {
    d <- data4group[, c(target, g), drop = FALSE]
    colnames(d) <- c("expression", "group")
    d
  })
  plot_df <- do.call(rbind, parts)
  plot_df <- plot_df[plot_df$group %in% positive_groups, ]

  p <- make_violin(plot_df, "group", "group", paste0(target, "_exp_positive Expression Level"))
  ggsave(file.path(outdir, paste0(target, "_across_groups_violin.svg")), p, height = 4, width = 12)
}

cat("Loading:", obj_path, "\n")
obj <- readRDS(obj_path)

cat("Adding positive/negative metadata\n")
for (m in markers) {
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")[m, ]
  label <- ifelse(counts > 0, paste0(m, "_exp_positive"), paste0(m, "_exp_negative"))
  obj <- AddMetaData(obj, metadata = label, col.name = paste0(m, "_exp"))
}

for (m in markers) analyze_gene(obj, m, outdir)
for (m in markers) plot_across_groups(obj, m, outdir)

cat("Done\n")
