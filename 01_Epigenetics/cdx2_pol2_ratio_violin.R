#!/usr/bin/env Rscript
## Integrated pipeline: last-CDX2-site Pol II post/pre ratio analysis + paired violin plot
##   Analysis: per gene, select the last CDX2 peak along transcription direction,
##     define pre/post CDX2 gene-body regions, quantify spike-in-normalized Pol II
##     bigWig mean0 signal per sample, compute log2(post/pre) per sample and
##     delta log2(post/pre) = treatment - control per gene.
##   Plot: paired control vs treatment violin + boxplot (all genes; lines for a subset)
## Usage:
##   Rscript cdx2_pol2_ratio_violin.R \
##     (--genes genes.bed | --gene_list genes.txt --gtf annot.gtf) \
##     --cdx2 peaks.bed --samples sample_info.tsv --out prefix \
##     [--gene_key auto] [--control control] [--treatment treatment] \
##     [--tss_exclude 2000] [--cdx2_buffer 0] [--min_pre 1000] [--min_post 1000] [--pseudo 0.01] \
##     [--outdir plots_dir] [--dpi 300] [--summary_stat both] \
##     [--control_color #4DAF4A] [--treatment_color #984EA3] \
##     [--decrease_color #2C7BB6] [--increase_color #D7191C] \
##     [--max_pair_lines 120] [--pair_line_alpha 0.12] [--pair_line_width 0.28] \
##     [--pair_point_alpha 0.28] [--pair_point_size 0.7] [--pair_line_seed 1] \
##     [--left_group_expand 0.55] [--right_group_expand 0.3] \
##     [--pair_ymin auto] [--pair_ymax auto]
##   sample_info.tsv columns: sample_id, condition, bw_file

suppressPackageStartupMessages({
  need_pkgs <- c("data.table", "GenomicRanges", "IRanges", "rtracklayer", "GenomeInfoDb", "ggplot2")
  missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R/Bioconductor packages: ", paste(missing_pkgs, collapse = ", "),
      "\nBiocManager::install(c('GenomicRanges','IRanges','rtracklayer','GenomeInfoDb'))\n",
      "install.packages(c('data.table','ggplot2'))\n",
      call. = FALSE
    )
  }
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(ggplot2)
})

usage <- function() {
  cat("\nUsage, Mode A: target gene BED6 already prepared\n")
  cat("  Rscript cdx2_pol2_ratio_violin.R --genes target_genes.bed --cdx2 CDX2_peaks.bed --samples sample_info.tsv --out output_prefix [options]\n\n")
  cat("Usage, Mode B: generate target gene BED6 from gene list + GTF\n")
  cat("  Rscript cdx2_pol2_ratio_violin.R --gene_list target_genes.txt --gtf annotation.gtf --cdx2 CDX2_peaks.bed --samples sample_info.tsv --out output_prefix [options]\n\n")
  cat("sample_info.tsv columns: sample_id, condition, bw_file\n")
  cat("Gene BED6 format: chr start end gene_id score strand\n\n")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(save = "no", status = 0)
  }
  if (length(args) %% 2 != 0) {
    usage()
    stop("Arguments must be supplied as --key value pairs.", call. = FALSE)
  }
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- args[[i + 1]]
    if (!grepl("^--", key)) stop("Invalid argument name: ", key, call. = FALSE)
    out[[sub("^--", "", key)]] <- val
    i <- i + 2
  }
  out
}

get_arg <- function(args, key, default = NULL, required = FALSE, type = "character") {
  val <- args[[key]]
  if (is.null(val)) {
    if (required) stop("Missing required argument --", key, call. = FALSE)
    val <- default
  }
  if (type == "numeric") val <- as.numeric(val)
  if (type == "integer") val <- as.integer(val)
  val
}

clean_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

read_bed_any <- function(file, type = c("gene", "peak")) {
  type <- match.arg(type)
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  dt <- fread(file, header = FALSE, sep = "\t", fill = TRUE, data.table = TRUE)
  if (ncol(dt) < 3) stop("BED file must have at least 3 columns: ", file, call. = FALSE)
  setnames(dt, paste0("V", seq_len(ncol(dt))))
  out <- data.table(
    chr = as.character(dt$V1),
    start0 = as.integer(dt$V2),
    end0 = as.integer(dt$V3)
  )
  out[, id := if (ncol(dt) >= 4) as.character(dt$V4) else paste0(type, "_", .I)]
  out[, score := if (ncol(dt) >= 5) as.character(dt$V5) else "0"]
  out[, strand := if (ncol(dt) >= 6) as.character(dt$V6) else "."]
  out <- out[!is.na(chr) & !is.na(start0) & !is.na(end0)]
  out <- out[end0 > start0]
  if (nrow(out) == 0) stop("No valid intervals found in: ", file, call. = FALSE)
  if (type == "gene") {
    bad <- !out$strand %in% c("+", "-")
    if (any(bad)) {
      stop(
        "Gene BED must be BED6 and strand must be + or -. Bad rows: ",
        paste(which(bad)[seq_len(min(5, sum(bad)))], collapse = ", "),
        call. = FALSE
      )
    }
  }
  out
}

read_gene_list <- function(file) {
  if (!file.exists(file)) stop("Gene list file not found: ", file, call. = FALSE)
  x <- fread(file, header = FALSE, sep = "\t", fill = TRUE, data.table = TRUE)
  if (ncol(x) < 1) stop("Empty gene list file: ", file, call. = FALSE)
  genes <- trimws(as.character(x[[1]]))
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- genes[!grepl("^#", genes)]
  genes <- genes[!tolower(genes) %in% c("gene", "genes", "gene_id", "gene_name", "symbol", "target_gene")]
  unique(genes)
}

make_genes_from_gtf <- function(gtf_file, gene_list_file, gene_key = "auto") {
  if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file, call. = FALSE)
  gene_key <- match.arg(gene_key, c("auto", "gene_name", "gene_id"))
  query <- read_gene_list(gene_list_file)
  if (length(query) == 0) stop("No valid genes found in gene list: ", gene_list_file, call. = FALSE)

  message("  - Importing GTF: ", gtf_file)
  gtf <- rtracklayer::import(gtf_file)
  meta <- as.data.table(as.data.frame(mcols(gtf)))
  if (!"gene_id" %in% colnames(meta)) {
    stop("GTF must contain gene_id attribute.", call. = FALSE)
  }
  if (!"gene_name" %in% colnames(meta)) {
    meta[, gene_name := NA_character_]
  }
  if (!"type" %in% colnames(meta)) {
    meta[, type := NA_character_]
  }

  dt <- data.table(
    chr = as.character(seqnames(gtf)),
    start0 = start(gtf) - 1L,
    end0 = end(gtf),
    strand = as.character(strand(gtf)),
    type = as.character(meta$type),
    gene_id = as.character(meta$gene_id),
    gene_name = as.character(meta$gene_name)
  )
  dt <- dt[!is.na(gene_id) & gene_id != "" & strand %in% c("+", "-")]
  dt[, gene_id_clean := clean_version(gene_id)]
  dt[is.na(gene_name) | gene_name == "", gene_name := gene_id_clean]

  if (any(dt$type == "gene", na.rm = TRUE)) {
    feat <- dt[type == "gene"]
  } else {
    message("  - No GTF feature type == 'gene' found; deriving gene bodies from transcript/exon-like features")
    keep_types <- c("transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon")
    feat <- dt[type %in% keep_types | is.na(type)]
  }
  if (nrow(feat) == 0) stop("No usable gene/transcript/exon features found in GTF.", call. = FALSE)

  query_clean <- clean_version(query)
  if (gene_key == "gene_name") {
    feat <- feat[gene_name %in% query]
  } else if (gene_key == "gene_id") {
    feat <- feat[gene_id %in% query | gene_id_clean %in% query_clean]
  } else {
    feat <- feat[gene_name %in% query | gene_id %in% query | gene_id_clean %in% query_clean]
  }

  if (nrow(feat) == 0) {
    stop("No genes from gene_list matched GTF. Check gene symbols/Ensembl IDs and GTF version.", call. = FALSE)
  }

  gene_dt <- feat[, .(
    start0 = min(start0, na.rm = TRUE),
    end0 = max(end0, na.rm = TRUE),
    n_gtf_features = .N
  ), by = .(gene_id_clean, gene_id, gene_name, chr, strand)]

  gene_dt[, width := end0 - start0]
  setorder(gene_dt, gene_id_clean, -width)
  gene_dt <- gene_dt[, .SD[1], by = gene_id_clean]

  gene_dt[, display_id := ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id_clean)]
  dup_display <- gene_dt$display_id[duplicated(gene_dt$display_id) | duplicated(gene_dt$display_id, fromLast = TRUE)]
  gene_dt[display_id %in% dup_display, display_id := paste0(display_id, "|", gene_id_clean)]

  out <- gene_dt[, .(
    chr,
    start0 = as.integer(start0),
    end0 = as.integer(end0),
    id = display_id,
    score = "0",
    strand,
    gene_id_clean,
    gene_id,
    gene_name,
    n_gtf_features
  )]
  out <- out[end0 > start0]

  matched_symbol <- unique(gene_dt$gene_name)
  matched_id <- unique(c(gene_dt$gene_id, gene_dt$gene_id_clean))
  unmatched <- query[!(query %in% matched_symbol | query %in% matched_id | clean_version(query) %in% matched_id)]

  list(genes = out, input_genes = query, unmatched = unmatched)
}

make_gr <- function(dt, strand_col = NULL) {
  st <- if (is.null(strand_col)) "*" else dt[[strand_col]]
  GRanges(
    seqnames = dt$chr,
    ranges = IRanges(start = dt$start0 + 1L, end = dt$end0),
    strand = st
  )
}

resolve_bw_paths <- function(sample_dt, sample_file) {
  sample_dir <- dirname(normalizePath(sample_file, mustWork = TRUE))
  sample_dt[, bw_file_original := bw_file]
  sample_dt[, bw_file := as.character(bw_file)]
  for (i in seq_len(nrow(sample_dt))) {
    f <- sample_dt$bw_file[i]
    if (!file.exists(f)) {
      f2 <- file.path(sample_dir, f)
      if (file.exists(f2)) sample_dt$bw_file[i] <- normalizePath(f2, mustWork = TRUE)
    } else {
      sample_dt$bw_file[i] <- normalizePath(f, mustWork = TRUE)
    }
  }
  missing <- !file.exists(sample_dt$bw_file)
  if (any(missing)) {
    stop(
      "These bigWig files do not exist:\n",
      paste(sample_dt$bw_file_original[missing], collapse = "\n"),
      call. = FALSE
    )
  }
  sample_dt
}

harmonize_chr_style <- function(dt, bw_file) {
  bw <- BigWigFile(bw_file)
  bw_chrs <- tryCatch(seqlevels(seqinfo(bw)), error = function(e) character())
  if (length(bw_chrs) == 0) return(dt)
  bed_has_chr <- mean(grepl("^chr", dt$chr)) > 0.5
  bw_has_chr <- mean(grepl("^chr", bw_chrs)) > 0.5
  if (bed_has_chr && !bw_has_chr) {
    dt[, chr := sub("^chr", "", chr)]
  }
  if (!bed_has_chr && bw_has_chr) {
    dt[, chr := paste0("chr", chr)]
  }
  dt
}

select_last_cdx2 <- function(genes, cdx2, tss_exclude = 2000) {
  genes <- copy(genes)
  cdx2 <- copy(cdx2)
  genes[, gene_i := .I]
  cdx2[, peak_i := .I]
  genes[, `:=`(
    gene_id = id,
    tss0 = ifelse(strand == "+", start0, end0),
    tes0 = ifelse(strand == "+", end0, start0)
  )]
  cdx2[, peak_center0 := floor((start0 + end0) / 2)]

  gene_gr <- make_gr(genes, strand_col = "strand")
  cdx2_gr <- make_gr(cdx2, strand_col = NULL)
  hits <- findOverlaps(gene_gr, cdx2_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(list(last = data.table(), skipped = genes[, .(gene_id, reason = "no_cdx2_peak_in_gene_body")]))
  }

  ov <- data.table(gene_i = queryHits(hits), peak_i = subjectHits(hits))
  ov <- merge(
    ov,
    genes[, .(gene_i, gene_id, chr, gene_start0 = start0, gene_end0 = end0, strand, tss0, tes0)],
    by = "gene_i",
    all.x = TRUE
  )
  ov <- merge(
    ov,
    cdx2[, .(
      peak_i,
      cdx2_id = id,
      cdx2_chr = chr,
      cdx2_start0 = start0,
      cdx2_end0 = end0,
      cdx2_center0 = peak_center0,
      cdx2_score = score
    )],
    by = "peak_i",
    all.x = TRUE
  )

  # Exclude CDX2 candidates within the promoter-proximal part of gene body.
  ov <- ov[
    (strand == "+" & cdx2_center0 >= gene_start0 + tss_exclude) |
      (strand == "-" & cdx2_center0 <= gene_end0 - tss_exclude)
  ]

  if (nrow(ov) == 0) {
    return(list(last = data.table(), skipped = genes[, .(gene_id, reason = "no_cdx2_peak_after_tss_exclusion")]))
  }

  last_plus <- ov[strand == "+", .SD[which.max(cdx2_center0)], by = gene_i]
  last_minus <- ov[strand == "-", .SD[which.min(cdx2_center0)], by = gene_i]
  last <- rbindlist(list(last_plus, last_minus), use.names = TRUE, fill = TRUE)

  genes_with_last <- unique(last$gene_id)
  skipped <- genes[!id %in% genes_with_last, .(gene_id = id, reason = "no_valid_cdx2_peak")]
  list(last = last, skipped = skipped)
}

make_pre_post_regions <- function(last_dt, tss_exclude = 2000, cdx2_buffer = 0, min_pre = 1000, min_post = 1000) {
  if (nrow(last_dt) == 0) return(list(regions = data.table(), sites = data.table(), skipped = data.table()))
  x <- copy(last_dt)

  x[, `:=`(
    pre_start0 = NA_integer_,
    pre_end0 = NA_integer_,
    post_start0 = NA_integer_,
    post_end0 = NA_integer_
  )]

  # + strand: TSS -> CDX2 -> TES; pre is upstream of CDX2, post downstream.
  x[strand == "+", `:=`(
    pre_start0 = gene_start0 + tss_exclude,
    pre_end0 = cdx2_start0 - cdx2_buffer,
    post_start0 = cdx2_end0 + cdx2_buffer,
    post_end0 = gene_end0
  )]

  # - strand: transcription runs from high coordinate (TSS) to low; pre is the
  # CDX2-to-TSS interval and post is the TES-to-CDX2 interval.
  x[strand == "-", `:=`(
    pre_start0 = cdx2_end0 + cdx2_buffer,
    pre_end0 = gene_end0 - tss_exclude,
    post_start0 = gene_start0,
    post_end0 = cdx2_start0 - cdx2_buffer
  )]

  x[, `:=`(
    pre_len = pre_end0 - pre_start0,
    post_len = post_end0 - post_start0
  )]

  x[, valid_region := TRUE]
  x[pre_len < min_pre, valid_region := FALSE]
  x[post_len < min_post, valid_region := FALSE]
  x[pre_start0 < 0 | post_start0 < 0, valid_region := FALSE]
  x[pre_end0 <= pre_start0 | post_end0 <= post_start0, valid_region := FALSE]

  skipped <- x[valid_region == FALSE, .(
    gene_id,
    reason = paste0("invalid_or_short_pre_post_region;pre_len=", pre_len, ";post_len=", post_len)
  )]

  keep <- x[valid_region == TRUE]
  if (nrow(keep) == 0) return(list(regions = data.table(), sites = x, skipped = skipped))

  pre <- keep[, .(
    chr,
    start0 = pre_start0,
    end0 = pre_end0,
    region_id = paste0(gene_id, "|pre"),
    score = "0",
    strand,
    gene_id,
    region = "pre",
    cdx2_id,
    cdx2_start0,
    cdx2_end0,
    cdx2_center0,
    gene_start0,
    gene_end0,
    tss0,
    tes0
  )]

  post <- keep[, .(
    chr,
    start0 = post_start0,
    end0 = post_end0,
    region_id = paste0(gene_id, "|post"),
    score = "0",
    strand,
    gene_id,
    region = "post",
    cdx2_id,
    cdx2_start0,
    cdx2_end0,
    cdx2_center0,
    gene_start0,
    gene_end0,
    tss0,
    tes0
  )]

  regions <- rbindlist(list(pre, post), use.names = TRUE)
  setorder(regions, chr, start0, end0)
  list(regions = regions, sites = keep, skipped = skipped)
}

bw_mean0_regions <- function(bw_file, regions_dt) {
  reg_gr <- GRanges(
    seqnames = regions_dt$chr,
    ranges = IRanges(start = regions_dt$start0 + 1L, end = regions_dt$end0),
    strand = "*"
  )
  bw_gr <- tryCatch(
    import(BigWigFile(bw_file), which = reg_gr, format = "BigWig"),
    error = function(e) {
      stop("Failed to import bigWig: ", bw_file, "\n", conditionMessage(e), call. = FALSE)
    }
  )

  result <- numeric(length(reg_gr))
  if (length(bw_gr) == 0) return(result)
  if (!"score" %in% colnames(mcols(bw_gr))) {
    stop("Imported bigWig has no score column: ", bw_file, call. = FALSE)
  }

  hits <- findOverlaps(reg_gr, bw_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(result)

  q <- queryHits(hits)
  s <- subjectHits(hits)
  inter <- pintersect(reg_gr[q], bw_gr[s], ignore.strand = TRUE)
  contrib <- width(inter) * as.numeric(mcols(bw_gr)$score[s])
  sums <- rowsum(contrib, group = q, reorder = FALSE)
  idx <- as.integer(rownames(sums))
  result[idx] <- as.numeric(sums[, 1]) / width(reg_gr)[idx]
  result
}

write_bed6 <- function(dt, file, name_col = "region_id") {
  bed <- dt[, .(chr, start0, end0, name = get(name_col), score, strand)]
  fwrite(bed, file, sep = "\t", col.names = FALSE)
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 1e-4")
  paste0("p = ", signif(p, 3))
}

safe_wilcox_paired <- function(treat, ctrl, alternative = "two.sided") {
  ok <- is.finite(treat) & is.finite(ctrl)
  treat <- treat[ok]
  ctrl <- ctrl[ok]
  if (length(treat) < 2) return(list(p.value = NA_real_, statistic = NA_real_, n = length(treat)))
  res <- tryCatch(
    wilcox.test(treat, ctrl, paired = TRUE, exact = FALSE, alternative = alternative),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(p.value = NA_real_, statistic = NA_real_, n = length(treat)))
  list(p.value = res$p.value, statistic = unname(res$statistic), n = length(treat))
}

format_value <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) >= 100) return(signif(x, 3))
  if (abs(x) >= 10) return(round(x, 2))
  round(x, 3)
}

summary_label <- function(x, summary_stat = "both") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("NA")
  labels <- character(0)
  if (summary_stat %in% c("median", "both")) {
    labels <- c(labels, paste0("Median = ", format_value(median(x, na.rm = TRUE))))
  }
  if (summary_stat %in% c("mean", "both")) {
    labels <- c(labels, paste0("Mean = ", format_value(mean(x, na.rm = TRUE))))
  }
  paste(labels, collapse = "\n")
}

make_direction <- function(delta) {
  out <- ifelse(
    delta < 0,
    "Decreased (delta < 0)",
    ifelse(delta > 0, "Increased (delta > 0)", NA_character_)
  )
  factor(out, levels = c(
    "Decreased (delta < 0)",
    "Increased (delta > 0)"
  ))
}

add_direction_caption <- function() {
  "Direction: delta = log2(post/pre)treatment - log2(post/pre)control. delta < 0 means treatment reduces downstream Pol II relative to upstream after the last intronic CDX2 site."
}

save_plot <- function(p, filename_base, width, height, dpi) {
  ggplot2::ggsave(paste0(filename_base, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  ggplot2::ggsave(paste0(filename_base, ".png"), p, width = width, height = height, dpi = dpi)
}

plot_violin <- function(ratio_delta, opt) {
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  control_col <- paste0("mean_log2_post_pre_ratio_", opt$control)
  treatment_col <- paste0("mean_log2_post_pre_ratio_", opt$treatment)
  if (!all(c(control_col, treatment_col) %in% colnames(ratio_delta))) {
    stop("Cannot find control/treatment ratio columns. Available columns: ", paste(colnames(ratio_delta), collapse = ", "), call. = FALSE)
  }

  plot_main <- data.frame(
    gene_id = as.character(ratio_delta$gene_id),
    control = as.numeric(ratio_delta[[control_col]]),
    treatment = as.numeric(ratio_delta[[treatment_col]]),
    stringsAsFactors = FALSE
  )
  plot_main$delta <- plot_main$treatment - plot_main$control
  plot_main <- plot_main[is.finite(plot_main$control) & is.finite(plot_main$treatment) & is.finite(plot_main$delta), , drop = FALSE]
  if (nrow(plot_main) == 0) stop("No finite paired control/treatment/delta values to plot.", call. = FALSE)
  plot_main$direction <- make_direction(plot_main$delta)

  condition_cols <- c(opt$control_color, opt$treatment_color)
  names(condition_cols) <- c(opt$control, opt$treatment)
  direction_cols <- c(opt$decrease_color, opt$increase_color)
  names(direction_cols) <- levels(plot_main$direction)

  paired_two <- safe_wilcox_paired(plot_main$treatment, plot_main$control, alternative = "two.sided")

  long_pair <- rbind(
    data.frame(gene_id = plot_main$gene_id, condition = opt$control, ratio = plot_main$control, delta = plot_main$delta, direction = plot_main$direction, stringsAsFactors = FALSE),
    data.frame(gene_id = plot_main$gene_id, condition = opt$treatment, ratio = plot_main$treatment, delta = plot_main$delta, direction = plot_main$direction, stringsAsFactors = FALSE)
  )
  long_pair$condition <- factor(long_pair$condition, levels = c(opt$control, opt$treatment))

  all_pair_genes <- unique(as.character(long_pair$gene_id))
  if (length(all_pair_genes) > opt$max_pair_lines && opt$max_pair_lines > 0) {
    set.seed(opt$pair_line_seed)
    show_pair_genes <- sample(all_pair_genes, opt$max_pair_lines)
  } else if (opt$max_pair_lines > 0) {
    show_pair_genes <- all_pair_genes
  } else {
    show_pair_genes <- character(0)
  }
  pair_line_df <- long_pair[long_pair$gene_id %in% show_pair_genes, , drop = FALSE]

  control_values <- plot_main$control[is.finite(plot_main$control)]
  treatment_values <- plot_main$treatment[is.finite(plot_main$treatment)]

  control_summary_text <- summary_label(control_values, opt$summary_stat)
  treatment_summary_text <- summary_label(treatment_values, opt$summary_stat)
  control_summary_inline <- gsub("\\n", "; ", control_summary_text)
  treatment_summary_inline <- gsub("\\n", "; ", treatment_summary_text)
  paired_stat_block <- paste0(
    opt$control, ": ", control_summary_inline,
    "\n", opt$treatment, ": ", treatment_summary_inline,
    "\nPaired Wilcoxon signed-rank test: ", format_p(paired_two$p.value),
    "; n = ", paired_two$n
  )

  pair_range <- range(long_pair$ratio[is.finite(long_pair$ratio)], na.rm = TRUE)
  pair_span <- diff(pair_range)
  if (!is.finite(pair_span) || pair_span == 0) pair_span <- max(abs(pair_range), 1)
  pair_ylim_auto <- c(pair_range[1] - 0.05 * pair_span, pair_range[2] + 0.05 * pair_span)
  pair_ylim <- if (is.finite(opt$pair_ymin) && is.finite(opt$pair_ymax)) {
    c(opt$pair_ymin, opt$pair_ymax)
  } else {
    pair_ylim_auto
  }
  message("Paired-plot y-axis range: ", paste(signif(pair_ylim, 5), collapse = " to "),
          if (is.finite(opt$pair_ymin)) " [manual]" else " [automatic]")

  p3 <- ggplot(long_pair, aes(x = condition, y = ratio)) +
    geom_violin(
      aes(fill = condition),
      width = 0.62,
      alpha = 0.48,
      color = "black",
      linewidth = 0.55,
      trim = TRUE
    ) +
    geom_line(
      data = pair_line_df,
      aes(group = gene_id, color = direction),
      alpha = opt$pair_line_alpha,
      linewidth = opt$pair_line_width
    ) +
    geom_point(
      data = pair_line_df,
      aes(group = gene_id, color = direction),
      alpha = opt$pair_point_alpha,
      size = opt$pair_point_size
    ) +
    geom_boxplot(
      width = 0.12,
      fill = "white",
      alpha = 0.92,
      outlier.shape = NA,
      linewidth = 0.55
    ) +
    scale_fill_manual(values = condition_cols, guide = "none") +
    scale_color_manual(values = direction_cols, drop = FALSE, na.value = "transparent", name = "Gene-wise direction") +
    scale_x_discrete(
      labels = c(opt$control, opt$treatment),
      expand = expansion(add = c(opt$left_group_expand, opt$right_group_expand))
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 0.8, linewidth = 0.8))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 9.5, face = "bold", margin = margin(t = 7)),
      axis.title.x = element_text(size = 8.5, lineheight = 1.10, margin = margin(t = 10)),
      axis.line = element_line(linewidth = 0.8),
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      axis.line.x = element_blank(),
      plot.margin = margin(8, 12, 10, 12)
    ) +
    labs(
      x = paired_stat_block,
      y = expression(log[2]~"Pol II post/pre ratio"),
      title = "Paired control vs treatment Pol II post/pre ratio",
      subtitle = paste0(
        "Trajectories shown for ", length(show_pair_genes), "/", length(all_pair_genes),
        " genes; distributions and statistics use all genes"
      ),
      caption = add_direction_caption()
    ) +
    coord_cartesian(ylim = pair_ylim, clip = "on")
  save_plot(
    p3,
    file.path(opt$outdir, paste0(opt$prefix, ".colored_control_vs_treatment_paired_boxplot_direction")),
    width = 3.5,
    height = 6.6,
    dpi = opt$dpi
  )

  p3_clean <- ggplot(long_pair, aes(x = condition, y = ratio)) +
    geom_violin(
      aes(fill = condition),
      width = 0.62,
      alpha = 0.48,
      color = "black",
      linewidth = 0.55,
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.12,
      fill = "white",
      alpha = 0.92,
      outlier.shape = NA,
      linewidth = 0.55
    ) +
    scale_fill_manual(values = condition_cols, guide = "none") +
    scale_x_discrete(
      labels = c(opt$control, opt$treatment),
      expand = expansion(add = c(opt$left_group_expand, opt$right_group_expand))
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 9.5, face = "bold", margin = margin(t = 7)),
      axis.title.x = element_text(size = 8.5, lineheight = 1.10, margin = margin(t = 10)),
      axis.line = element_line(linewidth = 0.8),
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      axis.line.x = element_blank(),
      plot.margin = margin(8, 10, 10, 10)
    ) +
    labs(
      x = paired_stat_block,
      y = expression(log[2]~"Pol II post/pre ratio"),
      title = "Control vs treatment Pol II post/pre ratio",
      subtitle = paste0("Distributions and paired statistics use all ", length(all_pair_genes), " genes"),
      caption = add_direction_caption()
    ) +
    coord_cartesian(ylim = pair_ylim, clip = "on")
  save_plot(
    p3_clean,
    file.path(opt$outdir, paste0(opt$prefix, ".colored_control_vs_treatment_clean_violin_boxplot")),
    width = 3.5,
    height = 6.6,
    dpi = opt$dpi
  )

  message("Violin plots written to: ", opt$outdir)
}

main <- function() {
  args <- parse_args()

  genes_file <- get_arg(args, "genes", default = NULL)
  gene_list_file <- get_arg(args, "gene_list", default = NULL)
  gtf_file <- get_arg(args, "gtf", default = NULL)
  gene_key <- get_arg(args, "gene_key", default = "auto")

  cdx2_file <- get_arg(args, "cdx2", required = TRUE)
  samples_file <- get_arg(args, "samples", required = TRUE)
  out_prefix <- get_arg(args, "out", required = TRUE)
  control_name <- get_arg(args, "control", default = "control")
  treatment_name <- get_arg(args, "treatment", default = "treatment")
  tss_exclude <- get_arg(args, "tss_exclude", default = 2000, type = "integer")
  cdx2_buffer <- get_arg(args, "cdx2_buffer", default = 0, type = "integer")
  min_pre <- get_arg(args, "min_pre", default = 1000, type = "integer")
  min_post <- get_arg(args, "min_post", default = 1000, type = "integer")
  pseudo <- get_arg(args, "pseudo", default = 0.01, type = "numeric")

  plot_outdir <- get_arg(args, "outdir", default = paste0(out_prefix, ".violin_plots"))
  plot_dpi <- get_arg(args, "dpi", default = 300, type = "numeric")
  summary_stat <- tolower(get_arg(args, "summary_stat", default = "both"))
  if (!summary_stat %in% c("median", "mean", "both")) stop("--summary_stat must be one of: median, mean, both", call. = FALSE)
  decrease_color <- get_arg(args, "decrease_color", default = "#2C7BB6")
  increase_color <- get_arg(args, "increase_color", default = "#D7191C")
  control_color <- get_arg(args, "control_color", default = "#4DAF4A")
  treatment_color <- get_arg(args, "treatment_color", default = "#984EA3")
  pair_line_alpha <- get_arg(args, "pair_line_alpha", default = 0.12, type = "numeric")
  pair_line_width <- get_arg(args, "pair_line_width", default = 0.28, type = "numeric")
  pair_point_alpha <- get_arg(args, "pair_point_alpha", default = 0.28, type = "numeric")
  pair_point_size <- get_arg(args, "pair_point_size", default = 0.7, type = "numeric")
  max_pair_lines <- get_arg(args, "max_pair_lines", default = 120, type = "numeric")
  pair_line_seed <- get_arg(args, "pair_line_seed", default = 1, type = "numeric")
  left_group_expand <- get_arg(args, "left_group_expand", default = 0.55, type = "numeric")
  right_group_expand <- get_arg(args, "right_group_expand", default = 0.3, type = "numeric")
  pair_ymin <- get_arg(args, "pair_ymin", default = NA_real_, type = "numeric")
  pair_ymax <- get_arg(args, "pair_ymax", default = NA_real_, type = "numeric")

  if (!is.finite(max_pair_lines) || max_pair_lines < 0) stop("--max_pair_lines must be >= 0", call. = FALSE)
  if ((is.finite(pair_ymin) || is.finite(pair_ymax)) && !(is.finite(pair_ymin) && is.finite(pair_ymax))) {
    stop("Provide both --pair_ymin and --pair_ymax, or neither for automatic limits.", call. = FALSE)
  }
  if (is.finite(pair_ymin) && is.finite(pair_ymax) && pair_ymin >= pair_ymax) {
    stop("--pair_ymin must be smaller than --pair_ymax.", call. = FALSE)
  }

  opt <- list(
    prefix = out_prefix,
    outdir = plot_outdir,
    control = control_name,
    treatment = treatment_name,
    dpi = plot_dpi,
    summary_stat = summary_stat,
    decrease_color = decrease_color,
    increase_color = increase_color,
    control_color = control_color,
    treatment_color = treatment_color,
    pair_line_alpha = pair_line_alpha,
    pair_line_width = pair_line_width,
    pair_point_alpha = pair_point_alpha,
    pair_point_size = pair_point_size,
    max_pair_lines = max_pair_lines,
    pair_line_seed = pair_line_seed,
    left_group_expand = left_group_expand,
    right_group_expand = right_group_expand,
    pair_ymin = pair_ymin,
    pair_ymax = pair_ymax
  )

  message("[1/9] Reading sample information and CDX2 BED")
  cdx2 <- read_bed_any(cdx2_file, type = "peak")

  samples <- fread(samples_file, header = TRUE, sep = "\t", data.table = TRUE)
  required_cols <- c("sample_id", "condition", "bw_file")
  if (!all(required_cols %in% colnames(samples))) {
    stop("sample_info.tsv must contain columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  samples[, sample_id := as.character(sample_id)]
  samples[, condition := as.character(condition)]
  samples <- resolve_bw_paths(samples, samples_file)

  if (!all(c(control_name, treatment_name) %in% unique(samples$condition))) {
    stop(
      "sample_info.tsv must contain both control and treatment conditions. Current conditions: ",
      paste(unique(samples$condition), collapse = ", "),
      call. = FALSE
    )
  }

  message("[2/9] Preparing target gene bodies")
  if (!is.null(genes_file) && !is.null(gene_list_file)) {
    stop("Use either --genes OR --gene_list + --gtf, not both.", call. = FALSE)
  }
  if (!is.null(genes_file)) {
    genes <- read_bed_any(genes_file, type = "gene")
  } else {
    if (is.null(gene_list_file) || is.null(gtf_file)) {
      stop("Provide either --genes target_genes.bed OR both --gene_list target_genes.txt and --gtf annotation.gtf.", call. = FALSE)
    }
    gene_obj <- make_genes_from_gtf(gtf_file, gene_list_file, gene_key = gene_key)
    genes <- gene_obj$genes
    generated_bed <- paste0(out_prefix, ".target_genes.generated_from_gene_list.bed")
    generated_tsv <- paste0(out_prefix, ".target_genes.generated_from_gene_list.tsv")
    write_bed6(genes, generated_bed, name_col = "id")
    fwrite(genes, generated_tsv, sep = "\t")
    if (length(gene_obj$unmatched) > 0) {
      fwrite(data.table(input_gene = gene_obj$unmatched), paste0(out_prefix, ".gene_list.unmatched.tsv"), sep = "\t")
      message("  - Unmatched genes written to: ", out_prefix, ".gene_list.unmatched.tsv")
    }
    message("  - Generated target gene BED6: ", generated_bed)
    message("  - Matched genes: ", nrow(genes), " / ", length(gene_obj$input_genes))
  }

  message("[3/9] Harmonizing chromosome style to first bigWig")
  genes <- harmonize_chr_style(genes, samples$bw_file[1])
  cdx2 <- harmonize_chr_style(cdx2, samples$bw_file[1])
  write_bed6(genes, paste0(out_prefix, ".target_genes.used.bed"), name_col = "id")
  fwrite(genes, paste0(out_prefix, ".target_genes.used.tsv"), sep = "\t")

  message("[4/9] Selecting last CDX2 peak along transcription direction")
  selected <- select_last_cdx2(genes, cdx2, tss_exclude = tss_exclude)
  last <- selected$last
  skipped1 <- selected$skipped
  if (nrow(last) == 0) stop("No valid last CDX2 sites found. Check gene list/GTF, CDX2 BED, and chromosome naming.", call. = FALSE)

  message("[5/9] Defining pre-CDX2 and post-CDX2 regions")
  pp <- make_pre_post_regions(
    last,
    tss_exclude = tss_exclude,
    cdx2_buffer = cdx2_buffer,
    min_pre = min_pre,
    min_post = min_post
  )
  regions <- pp$regions
  sites <- pp$sites
  skipped2 <- pp$skipped
  if (nrow(regions) == 0) stop("No valid pre/post regions after filtering. Try smaller --min_pre/--min_post or check CDX2 positions.", call. = FALSE)

  message("[6/9] Writing fixed CDX2 sites and pre/post BED")
  last_bed <- sites[, .(
    chr,
    start0 = cdx2_start0,
    end0 = cdx2_end0,
    name = paste0(gene_id, "|", cdx2_id, "|lastCDX2"),
    score = "0",
    strand
  )]
  fwrite(last_bed, paste0(out_prefix, ".last_cdx2_sites.bed"), sep = "\t", col.names = FALSE)
  write_bed6(regions, paste0(out_prefix, ".pre_post_regions.bed"), name_col = "region_id")

  site_table <- sites[, .(
    gene_id,
    chr,
    strand,
    gene_start0,
    gene_end0,
    tss0,
    tes0,
    cdx2_id,
    cdx2_start0,
    cdx2_end0,
    cdx2_center0,
    pre_start0,
    pre_end0,
    pre_len,
    post_start0,
    post_end0,
    post_len
  )]
  fwrite(site_table, paste0(out_prefix, ".last_cdx2_sites_and_regions.tsv"), sep = "\t")

  skipped <- rbindlist(list(skipped1, skipped2), use.names = TRUE, fill = TRUE)
  if (nrow(skipped) > 0) fwrite(skipped, paste0(out_prefix, ".skipped_genes.tsv"), sep = "\t")

  message("[7/9] Quantifying mean0 signal from each bigWig")
  all_sig <- list()
  for (i in seq_len(nrow(samples))) {
    message("  - ", samples$sample_id[i], " [", samples$condition[i], "]")
    sig <- bw_mean0_regions(samples$bw_file[i], regions)
    tmp <- copy(regions)[, .(
      chr,
      start0,
      end0,
      region_id,
      gene_id,
      region,
      strand,
      region_len = end0 - start0,
      cdx2_id,
      cdx2_start0,
      cdx2_end0,
      cdx2_center0
    )]
    tmp[, `:=`(
      sample_id = samples$sample_id[i],
      condition = samples$condition[i],
      bw_file = samples$bw_file[i],
      mean0 = sig
    )]
    all_sig[[i]] <- tmp
  }
  sig_dt <- rbindlist(all_sig, use.names = TRUE)
  fwrite(sig_dt, paste0(out_prefix, ".per_sample_region_mean0_signal.tsv"), sep = "\t")

  message("[8/9] Calculating per-sample log2(post/pre) ratios")
  ratio_sample <- dcast(
    sig_dt,
    sample_id + condition + bw_file + gene_id + strand + cdx2_id + cdx2_start0 + cdx2_end0 + cdx2_center0 ~ region,
    value.var = c("mean0", "region_len")
  )
  ratio_sample[, log2_post_pre_ratio := log2((mean0_post + pseudo) / (mean0_pre + pseudo))]
  fwrite(ratio_sample, paste0(out_prefix, ".per_sample_log2_post_pre_ratio.tsv"), sep = "\t")

  message("[9/9] Calculating condition means and treatment-control delta")
  ratio_condition <- ratio_sample[, .(
    mean_log2_post_pre_ratio = mean(log2_post_pre_ratio, na.rm = TRUE),
    sd_log2_post_pre_ratio = sd(log2_post_pre_ratio, na.rm = TRUE),
    n_samples = .N,
    mean_pre_signal = mean(mean0_pre, na.rm = TRUE),
    mean_post_signal = mean(mean0_post, na.rm = TRUE)
  ), by = .(gene_id, condition)]

  ratio_delta <- dcast(
    ratio_condition,
    gene_id ~ condition,
    value.var = c(
      "mean_log2_post_pre_ratio",
      "sd_log2_post_pre_ratio",
      "n_samples",
      "mean_pre_signal",
      "mean_post_signal"
    )
  )

  control_col <- paste0("mean_log2_post_pre_ratio_", control_name)
  treatment_col <- paste0("mean_log2_post_pre_ratio_", treatment_name)
  if (!all(c(control_col, treatment_col) %in% colnames(ratio_delta))) {
    stop("Cannot find control/treatment columns after dcast. Check condition names.", call. = FALSE)
  }
  ratio_delta[, delta_log2_post_pre_ratio := get(treatment_col) - get(control_col)]

  ratio_delta <- merge(
    ratio_delta,
    unique(site_table[, .(
      gene_id,
      chr,
      strand,
      gene_start0,
      gene_end0,
      tss0,
      tes0,
      cdx2_id,
      cdx2_start0,
      cdx2_end0,
      cdx2_center0,
      pre_len,
      post_len
    )]),
    by = "gene_id",
    all.x = TRUE
  )

  fwrite(ratio_condition, paste0(out_prefix, ".condition_mean_log2_post_pre_ratio.tsv"), sep = "\t")
  fwrite(ratio_delta, paste0(out_prefix, ".delta_treatment_vs_control.tsv"), sep = "\t")

  stat_file <- paste0(out_prefix, ".delta_stats.txt")
  deltas <- ratio_delta$delta_log2_post_pre_ratio
  deltas <- deltas[is.finite(deltas)]
  sink(stat_file)
  cat("Analysis: last intronic/gene-body CDX2 site-centered Pol II post/pre ratio\n")
  cat("Control condition:", control_name, "\n")
  cat("Treatment condition:", treatment_name, "\n")
  cat("tss_exclude:", tss_exclude, "bp\n")
  cat("cdx2_buffer:", cdx2_buffer, "bp\n")
  cat("min_pre:", min_pre, "bp\n")
  cat("min_post:", min_post, "bp\n")
  cat("pseudocount:", pseudo, "\n")
  cat("Number of genes with valid delta:", length(deltas), "\n")
  if (length(deltas) > 0) {
    cat("Mean delta:", mean(deltas), "\n")
    cat("Median delta:", median(deltas), "\n")
    cat("Fraction delta < 0:", mean(deltas < 0), "\n")
  }
  if (length(deltas) >= 3) {
    cat("\nWilcoxon test, two-sided, delta vs 0:\n")
    print(wilcox.test(deltas, mu = 0, alternative = "two.sided"))
    cat("\nWilcoxon test, alternative = less, delta < 0:\n")
    print(wilcox.test(deltas, mu = 0, alternative = "less"))
    cat("\nWilcoxon test, alternative = greater, delta > 0:\n")
    print(wilcox.test(deltas, mu = 0, alternative = "greater"))
  }
  sink()

  message("Plotting paired violin (control vs treatment)")
  plot_violin(ratio_delta, opt)

  message("Done.")
  message("Main outputs:")
  message("  ", out_prefix, ".target_genes.used.bed")
  message("  ", out_prefix, ".pre_post_regions.bed")
  message("  ", out_prefix, ".per_sample_region_mean0_signal.tsv")
  message("  ", out_prefix, ".per_sample_log2_post_pre_ratio.tsv")
  message("  ", out_prefix, ".delta_treatment_vs_control.tsv")
  message("  ", out_prefix, ".delta_stats.txt")
  message("  Violin plots: ", opt$outdir)
}

main()
