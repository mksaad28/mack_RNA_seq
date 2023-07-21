createHeatmap <- function() {
  # Load required packages
  library(pheatmap)
  library(tidyverse)
  library(scico)

  # Subset data
  data <- expr_wide |>
    dplyr::select(-Symbol) |>
    t() |>
    scale() |>
    t() |> # Scale to z-scores
    as.data.frame() |>
    dplyr::select(-contains("m"), -contains("f")) |>
    as.matrix()

  # Beautify column names
  rep_nums <- substr(colnames(data), 6, 6)
  cell_state <- substr(colnames(data), 5, 5)
  passage <- substr(colnames(data), 1, 3)

  # Add in Annotation Information
  annotation_col <- data.frame(CellState = cell_state, Passage = passage)

  ann_colors <- list(
    CellState = c(p = "gray20", s = "darkorange1"),
    Passage = c(p27 = "darkgoldenrod1", p77 = "brown1")
  )

  rownames(annotation_col) <- colnames(data)

  # Create cluster heatmap
  pheatmap(
    mat = data,
    color = scico(200, palette = "vik", direction = 1),
    show_rownames = F,
    show_colnames = T,
    cluster_cols = T,
    cluster_rows = T,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    border_color = F,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cutree_cols = 2,
    angle_col = 315,
    breaks = seq(-4, 4, length.out = 200),
    legend_breaks = c(-4, 0, 4),
    legend_labels = c("Low (-4)", "Medium (0)", "High (4)"),
    treeheight_row = 0
  )
}
