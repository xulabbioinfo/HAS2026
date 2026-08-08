
# ==============================================================================
# 0. Initialization & Universal Themes
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(fastCNV)
library(monocle)
library(tidyverse)
library(clusterProfiler)
library(GseaVis)
library(tools)
library(openxlsx)
library(ggsci)

future::plan("multisession", workers = 28) 
options(future.globals.maxSize = 10 * 1024^3)

base_dir <- "/home/data2/penglab"
setwd(base_dir)

# --- Define Universal Plotting Themes for Consistency ---
custom_umap_theme <- theme_minimal(base_family = "Arial") + 
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 16, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.text = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

custom_traj_theme <- theme_classic(base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position = "right"
  )

# ==============================================================================
# 1. Data Reading, Merging & Preprocessing
# ==============================================================================
samples <- c("SC-1", "SC-2", "SC-3", "SC-4", "SC-5")
objs <- list()

for (s in samples) {
  path <- file.path(base_dir, s)
  mat <- Read10X(data.dir = path)
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3, min.features = 200)
  obj$orig.ident <- s
  objs[[s]] <- obj
}

seurat_obj <- merge(x = objs[[1]], y = objs[2:5], add.cell.ids = samples, project = "Merged_Data")
rm(objs, mat, obj); gc()

# Quality Control (QC)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
HB.genes <- CaseMatch(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBQ2"), rownames(seurat_obj))
seurat_obj[["percent.HB"]] <- PercentageFeatureSet(seurat_obj, features = HB.genes)

seurat_obj <- subset(seurat_obj, subset = nCount_RNA < 15000 & nFeature_RNA > 200 & 
                       nFeature_RNA < 7500 & percent.mt < 20 & percent.HB < 5)

# Normalization & Dimensionality Reduction
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.2, 0.3), verbose = FALSE)

# ==============================================================================
# 2. Global Cell Annotation Visualization (Fig 2G & Fig S7a)
# ==============================================================================
# Visualize key lineage markers on UMAP to guide cell type annotation.
# --- Figure 2G: UMAP of all Cell Types ---
fig_2g_colors <- c("CD8+ T Cells" = "#1f78b4", "Neutrophil Cells" = "#a6cee3", 
                   "CD4+ T Cells" = "#33a02c", "Epithelial Cells" = "lightseagreen",
                   "Tumor" = "#e31a1c", "Macrophages" = "#fb9a99", 
                   "CAFs" = "#ff7f00", "Endothelial Cells" = "#b2df8a", 
                   "Plasma Cells" = "#6a3d9a", "B Cells" = "#cab2d6", "Mast Cells" = "#7E4909")

fig_2g <- DimPlot(seurat_obj, group.by = "Tumor_Group", alpha = 0.5) +
  scale_color_manual(values = fig_2g_colors) + 
  labs(title = "Global Cell Types", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme

print(fig_2g)

# --- Figure S7a: Marker Gene DotPlot ---
gene_to_check <- c("MS4A1","CD79A","CD79B","CD22","JCHAIN","MZB1","IGHG1","IGHG3",
                   "CD4","CD28","CTLA4","IL2RA","CD8A","CD8B","GZMB","NKG7",
                   "CD68","CD14","MRC1","AIF1","S100A8","S100A9","CSF3R","CXCR2",
                   "TPSAB1","TPSB2","MS4A2","CMA1","COL1A1","COL1A2","COL3A1","COL5A1",
                   "PECAM1","VWF","ENG","CDH5","EPCAM","KRT18","KRT8","MUC5AC")

fig_s7a <- DotPlot(seurat_obj, features = gene_to_check, group.by = "Tumor_Group") +
  scale_color_gradientn(colors = rev(colorspace::sequential_hcl(25, palette = "Reds"))) + 
  RotatedAxis() + 
  theme(
    text = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# ==============================================================================
# 3. Malignant Cell Identification via fastCNV (Fig S7b, S7c)
# ==============================================================================

seurat_obj$isTumor <- ifelse(seurat_obj$cnv_fraction > 0.2, "Tumor", "Epithelial Cells")

# --- Figure S7b: Aneuploid vs Diploid ---
fig_s7b <- DimPlot(seurat_obj, group.by = "CNV", alpha = 1) +
  scale_color_manual(values = c("Aneuploid" = "#e31a1c", "Diploid" = "grey")) +
  labs(title = "CNV Status", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme

# --- Figure S7c: Sample Origin ---
fig_s7c <- DimPlot(seurat_obj, group.by = "orig.ident", alpha = 1) +
  scale_color_manual(values = c("SC-1"="#e31a1c", "SC-2"="lightseagreen", 
                                "SC-3"="#6a3d9a", "SC-4"="#ff7f00", "SC-5"="#33a02c")) +
  labs(title = "Sample Identity", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme

# ==============================================================================
# 4. Epithelial Subsetting & Cell Cycle (Fig S7d series)
# ==============================================================================
tumor_obj <- subset(seurat_obj, subset = SCT_snn_res.0.2 == "3")
epi_obj <- subset(seurat_obj, subset = SCT_snn_res.0.2 %in% c("3", "4", "6"))

epi_obj <- CellCycleScoring(epi_obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
epi_obj$Pseudotime <- cds_epi$Pseudotime
# --- Figure S7d:  CNV fraction
fig_s7d_cnv <- FeaturePlot_scCustom(epi_obj, features = "cnv_fraction", reduction = "umap") + 
  labs(title = "CNV Fraction", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme +
  guides(
    color = guide_colorbar(
      title = "CNV_fraction", 
      barwidth = 1.5, 
      barheight = 10, 
      ticks = TRUE
    )
  )

# --- Figure S7d: Cell Cycle Phase ---
fig_s7d_phase <- DimPlot(epi_obj, group.by = "Phase", alpha = 1) +
  scale_color_manual(values = c("G1" = "grey", "G2M" = "#e31a1c", "S" = "#7E4909")) +
  labs(title = "Cell Cycle Phase", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme

# --- Figure S7d: Tumor vs Normal Epithelial ---
fig_s7d_tumor <- DimPlot(epi_obj, group.by = "isTumor", alpha = 0.5) +
  scale_color_manual(values = c("Epithelial Cells" = "lightseagreen", "Tumor" = "#e31a1c")) +
  labs(title = "Tumor vs Normal Epithelial", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme

# --- Figure S7d:  Pseudotime

fig_s7d_pseudotime <- FeaturePlot_scCustom(epi_obj, features = "Pseudotime", reduction = "umap") + 
  labs(title = "Pseudotime Projection", x = "UMAP 1", y = "UMAP 2") +
  custom_umap_theme +
  guides(
    color = guide_colorbar(
      title = "Pseudotime", 
      barwidth = 1.5, 
      barheight = 10, 
      ticks = TRUE
    )
  )


# ==============================================================================
# 5. Non-Zero Expression Violin Plot (Fig 2H)
# ==============================================================================
genes_2h <- c("AFP", "CDX2","NR1H4","PROX1","HHEX")
expr <- GetAssayData(tumor_obj, slot = "data")[genes_2h, ]
is_pos <- expr > 0

violin_df <- map_dfr(genes_2h, function(g) {
  pos_cells <- colnames(expr)[is_pos[g, ]]
  tibble(cell = pos_cells, gene = paste0(g, "+"), expression = expr[g, pos_cells])
})
violin_df$gene <- factor(violin_df$gene, levels = paste0(genes_2h, "+"))

fig_2h <- ggplot(violin_df, aes(x = gene, y = expression, fill = gene)) +
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("AFP+"   = "#2C7FB8",
                               "CDX2+"  = "#F03B20",
                               "HHEX+"  = "#31A354",
                               "NR1H4+" = "#756BB1",
                               "PROX1+" = "#636363")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) + ylab("Expression level (Non-zero)")

print(fig_2h)

# ==============================================================================
# 6. Trajectory Inference with Monocle 2
# ==============================================================================

data_epi <- epi_obj@assays$RNA$counts
pd_epi <- new('AnnotatedDataFrame', data = epi_obj@meta.data)
fd_epi <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(data_epi), row.names = rownames(data_epi)))

cds_epi <- newCellDataSet(
  data_epi, 
  phenoData = pd_epi, 
  featureData = fd_epi, 
  lowerDetectionLimit = 0.5, 
  expressionFamily = negbinomial.size()
)

cds_epi <- estimateSizeFactors(cds_epi)
cds_epi <- estimateDispersions(cds_epi)
cds_epi <- detectGenes(cds_epi, min_expr = 0.1)

expressed_genes_epi <- row.names(subset(fData(cds_epi), num_cells_expressed >= 10))

diff_test_res_epi <- differentialGeneTest(
  cds_epi[expressed_genes_epi, ], 
  fullModelFormulaStr = "~SCT_snn_res.0.2"
)

ordering_genes_epi <- row.names(subset(diff_test_res_epi, qval < 0.01))
cds_epi <- setOrderingFilter(cds_epi, ordering_genes_epi)
cds_epi <- reduceDimension(cds_epi, max_components = 2, method = 'DDRTree')
cds_epi <- orderCells(cds_epi)

epi_obj$Pseudotime <- pData(cds_epi)$Pseudotime

data_tumor <- tumor_obj@assays$RNA$counts
pd_tumor <- new('AnnotatedDataFrame', data = tumor_obj@meta.data)
fd_tumor <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(data_tumor), row.names = rownames(data_tumor)))

cds_tumor <- newCellDataSet(data_tumor, phenoData = pd_tumor, featureData = fd_tumor, 
                            lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

cds_tumor <- estimateSizeFactors(cds_tumor)
cds_tumor <- estimateDispersions(cds_tumor)
cds_tumor <- detectGenes(cds_tumor, min_expr = 0.1)

expressed_genes_tumor <- row.names(subset(fData(cds_tumor), num_cells_expressed >= 10))
diff_test_res_tumor <- differentialGeneTest(cds_tumor[expressed_genes_tumor, ], fullModelFormulaStr = "~Media")

ordering_genes_tumor <- row.names(subset(diff_test_res_tumor, qval < 0.01))
cds_tumor <- setOrderingFilter(cds_tumor, ordering_genes_tumor)
cds_tumor <- reduceDimension(cds_tumor, max_components = 2, method = 'DDRTree')
cds_tumor <- orderCells(cds_tumor)

tumor_obj$Pseudotime <- pData(cds_tumor)$Pseudotime

# ==============================================================================
# 7. Trajectory Visualization (Fig 2K, 2I & Module Scores)
# ==============================================================================
# --- Figure 2K: Pseudotime & State Projection ---
fig_2k_pseudotime <- plot_cell_trajectory(cds_tumor, color_by = "Pseudotime", show_branch_points = FALSE) +
  scale_color_viridis_c(option = "C", na.value = "grey80", name = "Pseudotime") +
  custom_traj_theme

fig_2k_state <- plot_cell_trajectory(cds_tumor, color_by = "State", show_branch_points = FALSE) +
  scale_color_npg(name = "State") +
  custom_traj_theme

# --- Figure 2I: Specific Marker / Module Mapping ---
# 1. Map single gene 
expr_afp <- GetAssayData(tumor_obj, slot = "data")["AFP", ]
pData(cds_tumor)$AFP_expr <- as.numeric(expr_afp)

fig_2i_marker <- plot_cell_trajectory(cds_tumor, color_by = "AFP_expr", show_branch_points = FALSE) +
  scale_color_viridis_c(option = "C", na.value = "grey80", name = "AFP Expr") +
  custom_traj_theme

# 2. Map Module Score (Gene Set)
target_genes <- c("genes") 
tumor_obj <- AddModuleScore(tumor_obj, features = list(target_genes), name = "Target_Module")
pData(cds_tumor)$Module_Score <- tumor_obj$Target_Module1

fig_2i_module <- plot_cell_trajectory(cds_tumor, color_by = "Module_Score", show_branch_points = FALSE) +
  scale_color_viridis_c(option = "C", na.value = "grey80", name = "Module Score") +
  custom_traj_theme

# ==============================================================================
# 8. GSEA Analysis & Plotting (Fig 2J)
# ==============================================================================
dir_p <- "/media/bio/新加卷/peng_gsea/" 
setwd(dir_p)
deg_file <- "deg_tumor_vs_epi_nofilter.csv"

data <- read.csv(file.path(dir_p, deg_file)) %>%
  rename(Gene_name = X, log2FoldChange = avg_log2FC, pvalue = p_val) %>%
  filter(!is.na(Gene_name) & !is.na(log2FoldChange))

data_s <- data[!duplicated(data$Gene_name), ]
alldiff <- data_s[order(data_s$log2FoldChange, decreasing = TRUE), ]
id <- setNames(alldiff$log2FoldChange, alldiff$Gene_name)

term2gene <- read.gmt('/home/data/database/MsigDB/human/symbol/h.all.v2023.1.Hs.symbols.gmt') 

set.seed(777) 
gsea_res <- GSEA(geneList = id, TERM2GENE = term2gene, minGSSize = 1, 
                 maxGSSize = 100000, pvalueCutoff = 0.05, seed = 777)

# --- Save Complete Results ---
g_all_df <- as.data.frame(gsea_res)
g_all_df$leadingEdge_count <- sapply(strsplit(as.character(g_all_df$core_enrichment), "/"), length)
g_all_df <- rename(g_all_df, leadingEdge_genes = core_enrichment)

write.xlsx(g_all_df, file.path(dir_p, "gsea_h_complete_results.xlsx"))

# --- Figure 2J: GseaVis Pathway Plots ---
target_pathways <- c("HALLMARK_COAGULATION", "HALLMARK_XENOBIOTIC_METABOLISM", 
                     "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM", 
                     "HALLMARK_ADIPOGENESIS", "HALLMARK_GLYCOLYSIS")

valid_pathways <- target_pathways[target_pathways %in% g_all_df$ID]

gsea_res_plot <- gsea_res
gsea_res_plot@result$NES <- round(gsea_res_plot@result$NES, 2)
gsea_res_plot@result$pvalue <- signif(gsea_res_plot@result$pvalue, 3)

for(p_name in valid_pathways) {
  clean_title <- toTitleCase(tolower(gsub("_", " ", gsub("HALLMARK_", "", p_name))))
  
  fig_2j <- gseaNb(object = gsea_res_plot, geneSetID = p_name, addPval = TRUE) + 
    ggtitle(clean_title) 
  
  # ggsave(fig_2j, file = paste0(p_name, "_GseaVis.pdf"), width = 6, height = 6.7)
}