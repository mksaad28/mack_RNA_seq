#' DESeq2
#' 6/19/23
#'
#' Written by:
#'
#' Edited by: Benji Bromberg

#' Get Results Lists from DESeq2 Differential Expression Analysis
#'
#' This function runs the DESeq2 pipeline script on a given design argument and
#' contrast argument. It also maps the DESeq2 results by GeneSymbol to results
#' by GeneSymbol and/or GeneID. The function requires the following variables to
#' be loaded: expr_wide_prep (data frame of RNA Seq counts, ready for DESeq2)
#' and meta (metadata, including passage, protocol, and combination).
#'
#' @param design_arg Formula for the model, based on combinations of columns in
#'  the metadata.
#' @param contrast_arg Vector specifying the comparison for log2-fold changes.
#'  Typically, this will contain three elements:
#'   - The condition to test (usually in the design formula)
#'   - The value in that condition for the numerator
#'   - The value in that condition for the denominator
#'
#' @return results A list containing seven data frames:
#'   - df_counts_norm: Data frame with normalized counts
#'   - df_counts_rld: Data frame with regularized log transformed (rld)
#'   - df_res: Data frame with differential expression results, after shrinkage
#'   - df_res_symbol: df_res, but with GeneSymbol is a column not rowname
#'   - df_res_symbol_GeneID: df_res_symbol with GeneID mapping
#'   - df_res_GeneID: df_res_symbol_GeneID without GeneSymbol mapping
#'   - df_res_GeneID_kegg:
#' @export
#'
#' @examples
#' # example code
getResultsLists <- function(
    preProcessOutputList,
    design_arg_list,
    contrast_arg_list,
    GID_mappings,
    OrgDb,
    pvalueCutoffPathway = 0.01,
    logFCCutoffPathway = 2) {
  # Create empty list to fill with results
  results_for_species <- list()

  for (name in names(design_arg_list)) {
    # Run DESeq2 on inputs
    print(paste0("Starting DESeq2 Processing for ", name))

    # Run DESeq2 by gene symbol
    deseq_res <- runDESeq2(
      preProcessOutputList = preProcessOutputList,
      design_arg = design_arg_list[[name]],
      contrast_arg = contrast_arg_list[[name]]
    )

    # Map results by gene symbol to results by gene symbol, GID, and KEGG GID
    deseq_mappings <- addMappings2DESeqRes(deseq_res$df_res, GID_mappings)

    # Output results to list and return results
    deseq_results <- list(
      df_counts_norm = deseq_res$df_counts_norm,
      df_counts_rld = deseq_res$df_counts_rld,
      df_res = deseq_res$df_res,
      df_res_symbol = deseq_mappings$df_res_symbol,
      df_res_symbol_GeneID = deseq_mappings$df_res_symbol_GeneID,
      df_res_GeneID = deseq_mappings$df_res_GeneID,
      df_res_GeneID_kegg = deseq_mappings$df_res_GeneID_kegg
    )

    # Run ORAs on DESeq2 outputs
    print(
      paste0("Starting GO and Kegg Over-representation Analyses for ", name)
    )

    pathway_results <- runPathwayAnalysis(
      df_res_GID_KEGG = deseq_results$df_res_GeneID_kegg,
      pvalueCutoff = pvalueCutoffPathway,
      logFCCutoff = logFCCutoffPathway,
      OrgDb = OrgDb
    )

    new_name <- paste("desqout", name, sep = "_")
    assign(new_name, deseq_results, envir = .GlobalEnv)
    results_for_species[[new_name]] <- get(new_name)
    rm(list = new_name, envir = .GlobalEnv)
    rm(new_name)

    new_name <- paste("pathway", name, sep = "_")
    assign(new_name, pathway_results, envir = .GlobalEnv)
    results_for_species[[new_name]] <- get(new_name)
    rm(pathway_results)
    rm(deseq_results)
    rm(list = new_name, envir = .GlobalEnv)
    rm(new_name)
  }

  return(results_for_species)
}

#' Run DESeq2 for Differential Expression Analysis
#'
#' This function performs differential expression analysis using DESeq2 on a
#' given dataset and specified inputs.
#'
#' @param design_arg Formula for the model, based on combinations of columns in
#'  the metadata.
#' @param contrast_arg Vector specifying the comparison for log2-fold changes.
#'  Typically, this will contain three elements:
#'   - The condition to test (usually in the design formula)
#'   - The value in that condition for the numerator
#'   - The value in that condition for the denominator
#' @include Data_Wrang_and_PreProcess.R
#'
#' @return output_runDESeq2 A list containing three data frames:
#'   - df_counts_norm: Data frame with normalized counts
#'   - df_counts_rld: Data frame with regularized log transformed (rld) counts
#'   - df_res: Data frame with differential expression results, after shrinkage
#'
#' @details This function requires the following variables to be loaded:
#'   - expr_wide_prep: dataframe of RNA Seq counts, ready for DESeq2
#'   - meta: metadata, currently passage, protocol, and the combination
#' @export
#'
#' @examples
#' # example code
runDESeq2 <- function(preProcessOutputList, design_arg, contrast_arg) {
  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(
    countData = preProcessOutputList$expr_wide_prep,
    colData = preProcessOutputList$meta,
    design = design_arg
  )

  # Run DESeq, which runs estimateSizeFactors, estimateDispersions,
  # and nbinomWaldTest
  dds <- DESeq(dds)

  # The rlog transformation computes the log2-transformed normalized count data
  # with a per-sample shift factor, which is chosen to minimize the dependence
  # of the variance on the mean.
  rld <- rlog(dds, blind = TRUE)

  # Data tables for plotting: Normalized and regularized log transformed (rld)
  df_counts_rld <- assay(rld) |> as.data.frame()
  df_counts_norm <- counts(dds, normalized = TRUE) |> as.data.frame()

  # Plot PCA
  plotPCA(rld, intgroup = contrast_arg[1]) + geom_text(aes(label = name))

  # Specify conditions for differential expression then run
  # Contrast matrix is column of meta matrix to compare, then the two conditions
  # to compare. The condition that comes first will be the numerator in log2FC
  res_unshrunken <- results(dds, contrast = contrast_arg)
  summary(res_unshrunken)

  # Shrinkage of log2 fold changes: avoids overestimates of differences between
  # genes with high dispersion
  res <- lfcShrink(dds,
    contrast = contrast_arg,
    res = res_unshrunken,
    type = "ashr"
  )
  df_res <- res |> as.data.frame()

  # Output results to list and return results
  output_runDESeq2 <- list(
    df_counts_norm = df_counts_norm,
    df_counts_rld = df_counts_rld,
    df_res = df_res
  )

  return(output_runDESeq2)
}
