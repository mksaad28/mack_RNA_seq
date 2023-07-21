createVolcanoPlot <- function(df_res_GeneID, plot_title, receptors_shown) {
  # Load required packages
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  # Define the significance threshold
  sig_threshold <- 0.01

  # Define the logFC threshold
  logFC_threshold <- 1

  # Define differential gene expression data
  DGE <- df_res_GeneID

  # Create a new column to label significant genes
  DGE$significant <- ifelse(DGE$padj <= sig_threshold &
    abs(DGE$log2FoldChange) >= logFC_threshold,
  "Significant", "Not Significant"
  )

  if (receptors_shown == TRUE) {
    # Create a new column to label upregulated receptors
    DGE$receptors <- ifelse(DGE$gene_id %in% testing_receptors$gene_id,
      "Receptors to test", DGE$significant
    )
    DGE$receptor_names <- ifelse(DGE$gene_id %in% testing_receptors$gene_id,
      testing_receptors$names, ""
    )

    # More receptor labeling
    DGE$plot_label <- ifelse(DGE$receptors == "Receptors to test",
      DGE$receptor_names, ""
    )

    # Set up the plot
    ggplot(DGE, aes(
      x = log2FoldChange, y = -log10(padj),
      color = receptors
    )) +
      geom_point(
        size = 1,
        alpha = ifelse(DGE$receptors == "Receptors to test", 1, 0.5)
      ) +
      scale_color_manual(values = c("gray20", "red", "black")) +

      # Add vertical line for log2FoldChange = 1
      geom_vline(
        xintercept = -1, linetype = "dashed", color = "gray40",
        linewidth = 0.5
      ) +

      # Add vertical line for log2FoldChange = -1
      geom_vline(
        xintercept = 1, linetype = "dashed", color = "gray40",
        linewidth = 0.5
      ) +

      # Add horizontal line for padj = 0.01
      geom_hline(
        yintercept = -log10(0.01), linetype = "dashed",
        color = "gray40", linewidth = 0.5
      ) +

      # Add plot labels and title
      labs(
        x = "log2(Fold Change)", y = "-log10(Adjusted p-value)",
        color = "", title = plot_title
      ) +

      # Add receptor labels
      geom_text_repel(aes(label = plot_label),
        size = 4,
        max.overlaps = 1000000000000000000000000000000000,
        show.legend = FALSE,
        force = 20,
        segment.alpha = 0.3
      )
  } else {
    # Set up the plot
    ggplot(DGE, aes(
      x = log2FoldChange, y = -log10(padj),
      color = significant
    )) +
      geom_point(size = 1, alpha = 0.5) +
      scale_color_manual(values = c("gray20", "red", "gray20")) +

      # Add vertical line for log2FoldChange = 1
      geom_vline(
        xintercept = -1, linetype = "dashed", color = "gray40",
        linewidth = 0.5
      ) +

      # Add vertical line for log2FoldChange = -1
      geom_vline(
        xintercept = 1, linetype = "dashed", color = "gray40",
        linewidth = 0.5
      ) +

      # Add horizontal line for padj = 0.01
      geom_hline(
        yintercept = -log10(0.01), linetype = "dashed",
        color = "gray40", linewidth = 0.5
      ) +

      # Add plot labels and title
      labs(
        x = "log2(Fold Change)", y = "-log10(Adjusted p-value)",
        color = "", title = plot_title
      )
  }
}
