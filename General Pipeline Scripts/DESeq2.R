# Function: run_DESeq2 performs differential expression on our dataset, using 
# DESeq2 and specified inputs

# Pre-requisites
  # expr_wide_prep: dataframe of RNA Seq counts, ready for DESeq2
  # meta: metadata, currently passage, protocol, and the combination

# INPUTS
  # design: formula for model, based on combinations of columns in meta 
    # (e.g. design_arg = ~ passage)
  # contrast: vector to specify the comparison for log2-fold changes. Typically,
    # this will contain three elements: 1) the condition to test (usually in the
    # design formula), 2) the value in that condition for the numerator, 3) the
    # value in that condition for the denominator
    # e.g., contrast_arg = c("passage", "p77", "p27")

# OUTPUTS: List with three data frames
  # df_counts_norm: data frame with normalized counts
  # df_counts_rld: data frame with counts regularized log transformed (rld)
  # df_res: data frame with differential expression results, after shrinkage

# FUNCTION STARTS
run_DESeq2 <- function(design_arg, contrast_arg) {
  # Create DESeq dataset then run DESeq, which runs estimateSizeFactors,
  # estimateDispersions, and nbinomWaldTest
  dds <- DESeqDataSetFromMatrix(countData = expr_wide_prep, 
                                colData = meta, 
                                design = design_arg)
  dds <- DESeq(dds)
  
  # The rlog transformation computes the log2-transformed normalized count data 
  # with a per-sample shift factor, which is chosen to minimize the dependence 
  # of the variance on the mean.
  rld <- rlog(dds, blind=TRUE)
  
  # Data tables for plotting: Normalized and regularized log transformed (rld)
  df_counts_rld <- assay(rld) |> as.data.frame()
  df_counts_norm <- counts(dds, normalized = TRUE) |> as.data.frame()
  
  # Plot PCA
  plotPCA(rld, intgroup=contrast_arg[1]) + geom_text(aes(label=name))
  
  # Specify conditions for differential expression then run
  # Contrast matrix is column of meta matrix to compare, then the two conditions
  # to compare. The condition that comes first will be the numerator in log2FC
  res_unshrunken <- results(dds, contrast=contrast_arg)
  summary(res_unshrunken)
  
  # Shrinkage of log2 fold changes: avoid overestimates of differences between 
  # genes with high dispersion
  res <- lfcShrink(dds, 
                   contrast=contrast_arg, 
                   res=res_unshrunken, 
                   type = 'ashr')
  df_res <- res |> as.data.frame()
  
  
  output_run_DESeq2 <- list(df_counts_norm = df_counts_norm, 
                            df_counts_rld = df_counts_rld, 
                            df_res = df_res)
  
  return(output_run_DESeq2)
}