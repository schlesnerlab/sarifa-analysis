source("/cwd/code/km.R")
source("/cwd/code/pairsplot.R")

ggtitle <- function(x) {
  return(ggplot2::ggtitle(""))
}

z_transform_row_wise <- function(x) {
  return(t(scale(t(x))))
  # mns <- colMeans(x, na.rm = na.rm)
  # sds <- apply(x, 2, sd, na.rm = na.rm)
  # x <- sweep(x, 2, mns, "-")
  # x <- sweep(x, 2, sds, "/")
  # return(x)
}

plot_count_density <- function(dds, path) {
  dds_assay <- as.data.frame(assay(dds))
  dds_assay <- dds_assay[sample(nrow(dds_assay), 1000), ]
  df <- data.frame(row_sums = rowSums(dds_assay))
  df["gene_id"] <- c(1:nrow(df))

  plot <- ggplot(df, aes(row_sums, fill = "Dichte")) +
    geom_density() +
    xlim(0, 500000) +
    xlab("Sum of summed counts per gene") +
    ylab("Count size magnitude") +
    ggtitle("Plot showing density of summed gene counts in corresponding count size magnitude")
  ggsave(paste0(path, "count_density.pdf"), plot = plot, height = 7, width = 14)
}

plot_pca <- function(dds, path) {
  taxonomies <- colnames(colData(dds))

  vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  pca <- PCAtools::pca(vsd, metadata = colData(dds), removeVar = 0.95)

  for (i in 2:length(taxonomies)) {
    tax <- taxonomies[i]
    try({
      num_categroies <- length(levels(factor(colData(dds)[, tax])))

      if (num_categroies > 1 && num_categroies < 15) {
        pairsplot(pca,
          components = getComponents(pca, seq_len(4)),
          triangle = FALSE,
          pointSize = 1.35,
          colby = tax,
          legendPosition = "bottom",
          legendLabSize = 10,
          legendTitleSize = 0,
          legendIconSize = 4,
        )
        dir.create(paste0(path, "pca"))
        ggsave(paste0(path, "pca/", tax, ".png"), width = 10, height = 7)
      }
    })
  }
}

plot_pca_mg <- function(dds, path) {
  taxonomies <- colnames(colData(dds))
  vsd <- DESeq2::vst(dds, blind = FALSE)

  for (i in 2:length(taxonomies)) {
    tax <- taxonomies[i]
    try({
      categories <- levels(factor(colData(dds)[, tax]))
      num_categroies <- length(categories)

      if (num_categroies > 1 && num_categroies < 15) {
        for (j in 1:num_categroies) {
          category <- categories[j]

          vsd_ss <- vsd[, which(vsd[[tax]] == category)]
          vsd_ss_se <- SummarizedExperiment::assay(vsd_ss)
          pca <- PCAtools::pca(vsd_ss_se, metadata = colData(vsd_ss), removeVar = 0.95)
          PCAtools::pairsplot(pca, triangle = FALSE, pointSize = 1.0, colby = "sarifa_group", legendPosition = "right", title = paste0("PCA for SARIFA annotated TCGA-STAD set for \n", tax, "=", category, " grouped by NON_SARIFA/SARIFA"))
          dir.create(paste0(path, "pca_mg/"), showWarnings = FALSE)
          ggsave(paste0(path, "pca_mg/", tax, "_", category, ".pdf"), width = 21, height = 14)
        }
      }
    })
  }
}

plot_pca_scree <- function(dds, path) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  pcaobj <- prcomp(t(SummarizedExperiment::assay(vsd)))
  pcaExplorer::pcascree(pcaobj, pc_nr = 20, type = "pev", title = "PCA Scree Plot: principal component\n proportion of SARIFA annotated TCGA-STAD set")
  ggsave(paste0(path, "pca_scree.pdf"))
}

plot_enhanced_volcano <- function(dds, res, path) {
  # Better to use padjust instead of pvalue
  EnhancedVolcano::EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "padj", title = paste0("DESeq2 results"), subtitle = paste0("NON_SARIFA versus SARIFA"))
  ggsave(paste0(path, "ev.pdf"), width = 14, height = 7)
}

plot_violin <- function(dds, path) {
  taxonomies <- colnames(colData(dds))

  for (i in 2:length(taxonomies)) {
    tax <- taxonomies[i]
    try({
      categories <- levels(factor(colData(dds)[, tax]))
      num_categroies <- length(categories)

      if (num_categroies > 5) {
        data <- as.data.frame(colData(dds)[c("sarifa_group", tax)])
        if (is.factor(data[, tax])) {
          data[, tax] <- unfactor(data[, tax])
        }
        data[, tax] <- as.double(data[, tax])
        data <- na.omit(data)
        ggplot(data, aes_string(x = "sarifa_group", y = tax, fill = "sarifa_group")) +
          geom_violin() +
          geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") +
          stat_compare_means(vjust = -2) +
          theme(
            axis.text = element_text(size = rel(1.35)),
            axis.title = element_text(size = rel(1.25)),
            legend.text = element_text(size = rel(1.35)),
            legend.title = element_text(size = rel(1.35))
          ) +
          ggtitle(paste0("Violin plot for numerical/categorical meta data values for\n ", tax))
        dir.create(paste0(path, "violin/"), showWarnings = FALSE)
        ggsave(paste0(path, "violin/", tax, ".png"))
      }
    })
  }
}

# plot_km <- function(dds, res) {
#   # res <- res[order(res$padj),]
#   res_sig <- subset(res, res$padj < 0.05)
#   # goi <- rownames(res[1:100, ])

#   clinical_patient_data <- TCGAbiolinks::GDCquery_clinic("TCGA-STAD", "clinical")

#   km <- TCGAanalyze_SurvivalKM_override(
#     clinical_patient_data,
#     dataGE=assay(dds),
#     Genelist=rownames(res_sig),
#     Survresult=TRUE,
#     p.cut=0.05,
#     ThreshTop=0.67,
#     ThreshDown=0.33
#   )

#   return(km)
# }

plot_km <- function(dds, res, path) {
  res_sig <- subset(res, res$padj < 0.05)
  res_sig <- subset(res_sig, abs(res_sig$log2FoldChange) > 1.0)
  clinical_patient_data <- TCGAbiolinks::GDCquery_clinic("TCGA-STAD", "clinical")

  # clinical_patient_data$bcr_patient_barcode
  colData(dds)$bcr <- substr(colData(dds)$bcr_patient_barcode, 1, 12)
  colData(dds)$bcr_full <- colData(dds)$bcr_patient_barcode
  clin <- merge(colData(dds), clinical_patient_data, by.x = "bcr", by.y = "bcr_patient_barcode", suffixes = c(".dds-colData", ""))

  for (i in 1:nrow(res_sig)) {
    gene_symbol <- rownames(res_sig)[i]
    gene_counts <- counts(dds)[gene_symbol, ]
    gene_expression_level <- sjmisc::dicho(gene_counts, dich.by = "median")

    colName <- paste0(gene_symbol, " Low-High")

    gene_expression_level_df <- data.frame(gene_expression_level)
    colnames(gene_expression_level_df) <- c(colName)
    gene_expression_level_df$bcr_full <- substr(rownames(gene_expression_level_df), 1, 16)

    clin_with_gene_expression_level <- merge(clin, gene_expression_level_df, by = "bcr_full")

    tryCatch({
      TCGAanalyze_survival_override(
        clin_with_gene_expression_level,
        path = path,
        clusterCol = colName,
        genename = colName,
        main = paste0("Kaplan-Meier Overall Survival Curves for gene ", gene_symbol, " expression level"),
        label = c(paste0(gene_symbol, " below median expression for SARIFA annotated TCGA-STAD set"), 
                  paste0(gene_symbol, " above median experssion for SARIFA annotated TCGA-STAD set"))
      )
    })
  }
}

plot_km_sarifa <- function(dds, path) {
  clinical_patient_data <- TCGAbiolinks::GDCquery_clinic("TCGA-STAD", "clinical")

  # clinical_patient_data$bcr_patient_barcode
  colData(dds)$bcr <- substr(colData(dds)$bcr_patient_barcode, 1, 12)
  colData(dds)$bcr_full <- colData(dds)$bcr_patient_barcode
  clin <- merge(colData(dds), clinical_patient_data, by.x = "bcr", 
                by.y = "bcr_patient_barcode", suffixes = c(".dds-colData", ""))

  dir.create(paste0(path, "km/"), showWarnings = FALSE)
  TCGAbiolinks::TCGAanalyze_survival(
    clin,
    path = path,
    clusterCol = "sarifa_group",
    filename = paste0(path, "km/sarifa.pdf"),
    risk.table = FALSE,
    main = paste0("Kaplan-Meier Overall Survival Curves for SARIFA status"),
    label = c("Annotated as NON_SARIFA for TCGA-STAD set", "Annotated as SARIFA for TCGA-STAD set")
  )
}

plot_cluster_heatmap <- function(dds, res, path) {
  # ntd <- DESeq2::normTransform(dds)
  # assay(ntd) <- z_transform_row_wise(assay(ntd))

  # # Most expressed genes
  # select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
  # Most DE genes
  res_sig <- subset(res, res$padj < 0.05)
  res_sig <- subset(res_sig, abs(res_sig$log2FoldChange) > 1.0)
  res_sig <- res_sig[order(res_sig$padj), ]
  # select <- order(res_sig$padj)[1:nrow(res_sig)]

  subset_vector <- c(
    "sarifa_group",
    "paper_Molecular.Subtype",
    "paper_ABSOLUTE.Purity",
    "paper_Estimated.Leukocyte.Percentage",
    "Stromal_score",
    "Immune_score"
  )
  df <- as.data.frame(colData(dds)[, subset_vector])

  df$`paper_Molecular.Subtype` <- unfactor(df$`paper_Molecular.Subtype`)
  df$`paper_Molecular.Subtype`[is.na(df$`paper_Molecular.Subtype`)] <- "UNDEFINED"
  df$`paper_Molecular.Subtype` <- factor(df$`paper_Molecular.Subtype`)

  df$`paper_ABSOLUTE.Purity` <- as.numeric(unfactor(df$`paper_ABSOLUTE.Purity`))
  df$`paper_Estimated.Leukocyte.Percentage` <- as.numeric(unfactor(df$`paper_Estimated.Leukocyte.Percentage`))

  rownames(df) <- colnames(dds)
  rownames(df) <- substr(rownames(df), 1, 15)
  annotation_col_col_names <- c(
    "SARIFA Group",
    "Molecular Subtype",
    "Absolute Purity",
    "Estimated Leukocyte Percentage",
    "Stromal Score",
    "Immune Score"
  )
  colnames(df) <- annotation_col_col_names

  annotation_col <- df
  annotation_colors <- list(
    "SARIFA Group" = c(NON_SARIFA = "#F8766D", SARIFA = "#03BFC4")
  )



  vst <- DESeq2::vst(dds, blind = FALSE)
  # ntd <- DESeq2::normTransform(dds)
  data <- assay(vst)[rownames(res_sig), order(df$`SARIFA Group`)]
  data <- z_transform_row_wise(data)
  colnames(data) <- substr(colnames(data), 1, 15)

  dir.create(paste0(path, "cluster"), showWarnings = FALSE)

  color_steps <- 100
  color <- colorRampPalette(c("navy", "white", "firebrick3"))(color_steps)
  breaks <- quantile(data, probs = seq(0, 1, length.out = color_steps))
  breaks <- breaks[!duplicated(breaks)]

  pheatmap(
    data,
    color = color,
    breaks = breaks,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_method = "complete",
    scale = "row",
    filename = paste0(path, "cluster/heatmap_custer_cols_by_complete.pdf"),
    width = 14,
    height = 9
  )

  pheatmap(
    data,
    color = color,
    breaks = breaks,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_method = "average",
    scale = "row",
    filename = paste0(path, "cluster/heatmap_custer_cols_by_average.pdf"),
    width = 14,
    height = 9
  )

  pheatmap(
    data,
    color = color,
    breaks = breaks,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_method = "ward.D",
    scale = "row",
    filename = paste0(path, "cluster/heatmap_custer_cols_by_ward_D.pdf"),
    width = 14,
    height = 9
  )

  pheatmap(
    data,
    color = color,
    breaks = breaks,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_method = "ward.D2",
    scale = "row",
    filename = paste0(path, "cluster/heatmap_custer_cols_by_ward_D2.pdf"),
    width = 14,
    height = 9
  )
}

plot_enrich_results <- function(enrich_res, path) {
  dir.create(paste0(path, "gsea"), showWarnings = FALSE)

  dotplot <- clusterProfiler::dotplot(enrich_res, showCategory = 10, font.size = 11)
  ggsave(paste0(path, "gsea/dotplot.png"), plot = dotplot, width = 7, height = 9)

  upsetplot <- enrichplot::upsetplot(enrich_res, n = 10)
  ggsave(paste0(path, "gsea/upsetplot.pdf"), plot = upsetplot, width = 14, height = 7)
}

export_raw_deseq2_results <- function(res, path) {
  res$geneSymbol <- rownames(res)
  res <- res[order(res$log2FoldChange), ]

  dir.create(paste0(path, "export"), showWarnings = FALSE)
  writexl::write_xlsx(as.data.frame(res), paste0(path, "export/non_sarifa_vs_sarifa_raw_results.xlsx"))

  res_sig <- subset(res, res$padj < 0.05)
  res_sig <- subset(res_sig, abs(res_sig$log2FoldChange) > 1.0)
  writexl::write_xlsx(as.data.frame(res_sig), paste0(path, "export/non_sarifa_vs_sarifa_sig_results.xlsx"))
}

export_enrich_results <- function(enrich_res, path) {
  enrich_res <- enrich_res[order(enrich_res$enrichmentScore), ]
  writexl::write_xlsx(as.data.frame(enrich_res), paste0(path, "export/non_sarifa_vs_sarifa_enrich_results.xlsx"))
}

export_butchr_input <- function(dds, res, path) {
  res_sig <- subset(res, res$padj < 0.05)
  res_sig <- subset(res_sig, abs(res_sig$log2FoldChange) > 1.0)
  res_sig <- res_sig[order(res_sig$padj), ]

  df <- as.data.frame(colData(dds)[, c("sarifa_group", "paper_Molecular.Subtype")])
  df$`paper_Molecular.Subtype` <- unfactor(df$`paper_Molecular.Subtype`)
  df$`paper_Molecular.Subtype`[is.na(df$`paper_Molecular.Subtype`)] <- "UNDEFINED"
  df$`paper_Molecular.Subtype` <- factor(df$`paper_Molecular.Subtype`)

  rownames(df) <- colnames(dds)
  rownames(df) <- substr(rownames(df), 1, 15)
  colnames(df) <- c("SARIFA Group", "Molecular Subtype")

  stromal_scores <- readr::read_tsv("https://bioinformatics.mdanderson.org/estimate/tables/stomach_adenocarcinoma_RNAseqV2.txt", show_col_type = FALSE)
  stromal_scores <- as.data.frame(stromal_scores)
  rownames(stromal_scores) <- stromal_scores$ID

  df <- merge(df, stromal_scores, by = 0)
  rownames(df) <- df$Row.names

  annotation_col <- df[, c("SARIFA Group", "Molecular Subtype")]
  annotation_col$sampleId <- df$Row.names
  annotation_col$`Stromal Score` <- df$Stromal_score > 0
  annotation_col$`Immune Score` <- df$Immune_score > 0
  # annotation_col$`Purity` <- df$`paper_ABSOLUTE.Purity`

  annotation_col$`Stromal Score` <- factor(annotation_col$`Stromal Score`)
  annotation_col$`Immune Score` <- factor(annotation_col$`Immune Score`)

  # Define explict order
  annotation_col <- annotation_col[, c("sampleId", "SARIFA Group", "Molecular Subtype", "Stromal Score", "Immune Score")]

  # annotation_col$`Purity` <- factor(annotation_col$`paper_ABSOLUTE.Purity`)

  # vst <- DESeq2::vst(dds, blind=FALSE)
  # data <- assay(vst)[rownames(res_sig), order(df$`SARIFA Group`)]

  counts <- assay(dds_reduced)[rownames(res_sig), ]
  sfs <- DESeq2::estimateSizeFactorsForMatrix(counts)
  counts <- t(t(counts) / sfs)
  counts <- apply(counts + 1, 2, log2)


  colnames(counts) <- substr(colnames(counts), 1, 15)
  # Fix order
  annotation_col <- annotation_col[colnames(counts), ]
  # Fix for shiny butchR export format
  rownames(annotation_col) <- 1:nrow(annotation_col)

  dir.create(paste0(path, "export"), showWarnings = FALSE)
  dir.create(paste0(path, "export/butchr"), showWarnings = FALSE)
  saveRDS(data.matrix(counts), paste0(path, "export/butchr/matrix.rds"))
  saveRDS(annotation_col, paste0(path, "export/butchr/annotation.rds"))

}

plot_deseq_results <- function(dds, res, enrich_res, path) {
  dir.create(path, showWarnings = FALSE)

  plot_violin(dds, path)
  return()

  plot_enhanced_volcano(dds, res, path)
  plot_count_density(dds, path)
  plot_pca_scree(dds, path)
  plot_pca(dds, path)
  plot_pca_mg(dds, path)
  plot_km(dds, res, path)
  plot_km_sarifa(dds, path)
  plot_violin(dds, path)
  plot_cluster_heatmap(dds, res, path)

  export_raw_deseq2_results(res, path)
  export_enrich_results(enrich_res, path)
  export_butchr_input(dds, res, path)
}


plot_export_compare_results <- function(res_reduced, res_strict, path) {
  res_merged <- merge(as.data.frame(res_reduced), 
                      as.data.frame(res_strict), by = 0)

  ggscatter(res_merged,
    x = "padj.x",
    y = "padj.y",
    add = "reg.line",
    conf.int = TRUE,
    size = 0.1,
    cor.coef = TRUE,
    cor.method = "pearson",
    xlab = "strict",
    ylab = "reduced"
  ) +
    geom_density_2d_filled(alpha = 0.5) +
    ggtitle("Plot showing correlation plot and pearson correlation coefficient of adjusted p-values,\nprovided by DESeq2 for strict and reduced annotated SARIFA sets")

  ggsave(paste0(path, "padj_correlation.pdf"))


  ggscatter(res_merged,
    x = "pvalue.x",
    y = "pvalue.y",
    add = "reg.line",
    conf.int = TRUE,
    size = 0.1,
    cor.coef = TRUE,
    cor.method = "pearson",
    xlab = "strict",
    ylab = "reduced"
  ) +
    geom_density_2d_filled(alpha = 0.5) +
    ggtitle("Plot showing correlation plot and pearson correlation coefficient of raw p-values,\nprovided by DESeq2 for strict and reduced annotated SARIFA sets")

  ggsave(paste0(path, "pvalue_correlation.pdf"))


  ggscatter(res_merged,
    x = "log2FoldChange.x",
    y = "log2FoldChange.y",
    add = "reg.line",
    conf.int = TRUE,
    size = 0.1,
    cor.coef = TRUE,
    cor.method = "pearson",
    xlab = "strict",
    ylab = "reduced"
  ) +
    geom_density_2d_filled(alpha = 0.5) +
    ggtitle("Plot showing correlation plot and pearson correlation coefficient of raw p-values,\nprovided by DESeq2 for strict and reduced annotated SARIFA sets")

  ggsave(paste0(path, "lfc_correlation.pdf"))

  res_reduced_sig <- subset(res_reduced, res_reduced$padj < 0.05)
  res_reduced_sig <- subset(res_reduced_sig, 
                            abs(res_reduced_sig$log2FoldChange) > 1.0)

  res_strict_sig <- subset(res_strict, res_strict$padj < 0.05)
  res_strict_sig <- subset(res_strict_sig, 
                           abs(res_strict_sig$log2FoldChange) > 1.0)

  ggvenn(
    list(
      reduced = rownames(res_reduced_sig),
      strict = rownames(res_strict_sig)
    ),
    c("reduced", "strict")
  ) +
    ggtitle("Venn diagram for genes classified as significantly\nexpressed by DESeq2 for strict and reduced annotated SARIFA sets")
  ggsave(paste0(path, "venn.pdf"))

  res_merged_sig <- subset(res_merged, padj.y < 0.05 & padj.x < 0.05)
  writexl::write_xlsx(as.data.frame(res_merged_sig), paste0(path, "non_sarifa_vs_sarifa_sig_merged_results.xlsx"))
}
