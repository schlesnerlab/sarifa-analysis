library("readxl")
library("TCGAbiolinks")
library("DESeq2")
library("ggplot2")
library("purrr")
library("cachem")
library("stats")
library("pcaExplorer")
library("PCAtools")
library("SummarizedExperiment")
library("EnsDb.Hsapiens.v86")
library("stringr")
library("EnhancedVolcano")
library("sjmisc")
library("writexl")
library("plyr")
library("survminer")
library("pheatmap")
library("ggvenn")
library("dplyr")
library("clusterProfiler")
library("msigdbr")
library("cachem")
library("cowplot")

options(error = rlang::entrace)

source("/cwd/code/plot.R")

edb <- EnsDb.Hsapiens.v86
cache <- cachem::cache_disk(dir = "/cwd/data/cache", max_size = Inf)

assign_symbols_to_DESeqDataSet <- function(dds) {
  # Use ensb as symbol
  gene_ids <- stringr::str_replace(rownames(dds), "[.][0-9]*$", "")

  gene_names <- mapIds(edb, keys = gene_ids, keytype = "GENEID",
      column = "SYMBOL", multiVals = "first")

  # Polyfill gene names with ids
  gene_names[sapply(gene_names, 
  function(x) length(x) == 0 || is.na(x))] <- gene_ids[sapply(gene_names,
                                      function(x) length(x) == 0 || is.na(x))]

  rownames(dds) <- make.unique(unlist(gene_names))

  return(dds)
}

create_DESeqDataSet <- function(reduced) {
  path_to_xlsx <- "/cwd/data/stad_tcga_pan_can_atlas_2018_clinical_data_bg.xlsx"
  stad_sarifa_cats <- readxl::read_excel(path_to_xlsx)

  stad_sarifa_list <- stad_sarifa_cats$`SARIFA...4`

  if (reduced) {
    stad_sarifa_list <- stad_sarifa_cats$`SARIFAreduziert`
  }

  names(stad_sarifa_list) <- stad_sarifa_cats$`Patient ID`

  stad_sarifa_list <- stad_sarifa_list[stad_sarifa_list != "NA"]
  stad_sarifa_list <- stad_sarifa_list[stad_sarifa_list != "NV"]

  stad_query <- TCGAbiolinks::GDCquery(
    project = "TCGA-STAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )


  TCGAbiolinks::GDCdownload(stad_query, 
                            method = "api", 
                            directory = "/cwd/data/stad")

  stad_se <- TCGAbiolinks::GDCprepare(stad_query, 
                                      directory = "/cwd/data/stad")
  stad_se <- subset(stad_se, 
                    select = colData(stad_se)$patient %in%
                      names(stad_sarifa_list))
  stad_se <- subset(stad_se, 
                    select = colData(stad_se)$sample_type == "Primary Tumor")

  sarifa_group_list <- purrr::map(colData(stad_se)$patient, 
  function(patient_id) {
    if (is.na(stad_sarifa_list[patient_id])) {
      return("UNDEFINED")
    }

    if (stad_sarifa_list[patient_id] == "0") {
      return("NON_SARIFA")
    }

    if (stad_sarifa_list[patient_id] == "1") {
      return("SARIFA")
    }
  }) # nolint

  col_data_sarifa <- data.frame(sarifa_group = factor(as.vector(unlist(sarifa_group_list))), row.names = rownames(colData(stad_se)))

  col_data <- merge(x = col_data_sarifa, y = colData(stad_se), by = 0)

  # Merge Stromalscores
  estimate_link <- "https://bioinformatics.mdanderson.org/estimate/tables/stomach_adenocarcinoma_RNAseqV2.txt"
  stromal_scores <- readr::read_tsv(estimate_link, show_col_type = FALSE)
  col_data_stromal_scores <- as.data.frame(stromal_scores)
  rownames(col_data_stromal_scores) <- col_data_stromal_scores$ID
  rownames(col_data) <- substr(col_data$barcode, 1, 15)
  col_data <- merge(x = col_data_stromal_scores, y = col_data, by = 0)

  rownames(col_data) <- col_data$barcode
  stad_se <- subset(stad_se, select = colData(stad_se)$barcode %in% col_data$barcode)
  # Preserve order of original colData
  col_data <- col_data[rownames(colData(stad_se)), ]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = assay(stad_se), 
                                        colData = col_data, 
                                        design = ~sarifa_group)

  dds <- dds[rowSums(counts(dds)) >= 200, ]

  dds <- DESeq2::DESeq(dds)
  rm(stad_se)
  gc()
  return(dds)
}

get_deseq2_results <- function(dds) {

  contrast <- c("sarifa_group", "SARIFA", "NON_SARIFA")
  res <- results(dds, contrast = contrast)

  return(res)
}

create_gsea <- function(res) {

  gene_list <- (-log10(res$pvalue)) * res$log2FoldChange
  names(gene_list) <- rownames(res)
  gene_list <- sort(gene_list, decreasing = TRUE)
  # Und wir sollten mal das gesamte C2 testen weil nur CP hat jetzt kaum was ergeben
  # TODO Vielleicht auch noch Hallmark "H"
  gene_sets <- msigdbr::msigdbr(species = "human", 
                                category = "C2", subcategory = "CP:REACTOME")
  gene_sets <- as.data.frame(dplyr::distinct(gene_sets, gs_name, ensembl_gene))

  gene_symbols <- mapIds(edb, keys = names(gene_list),
                         keytype = "SYMBOL", column = "GENEID", 
                         multiVals = "first")
  names(gene_list) <- make.unique(unlist(gene_symbols))

  enrich_results <- clusterProfiler::GSEA(gene_list, TERM2GENE = gene_sets)

  return(enrich_results)
}


dds_strict <- cache$get("dds_strict")
if (is.key_missing(dds_strict)) {
  dds_strict <- create_DESeqDataSet(FALSE)
  cache$set("dds_strict", dds_strict)
}
dds_strict <- assign_symbols_to_DESeqDataSet(dds_strict)

res_strict <- cache$get("res_strict")
if (is.key_missing(res_strict)) {
  res_strict <- get_deseq2_results(dds_strict)
  cache$set("res_strict", res_strict)
}
enrich_res_strict <- create_gsea(res_strict)


plot_deseq_results(dds_strict, res_strict, 
                   enrich_res_strict, "/cwd/data/results/strict/")


# Render markdown document with extra plots
rmarkdown::render("/cwd/code/SARIFA_report.Rmd")