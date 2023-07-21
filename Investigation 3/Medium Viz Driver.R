#' Medium Viz Driver
#' 6/29/23
#'
#' Written by: Benji Bromberg
#'
#' Purpose:
#'  TODO

# Run General Pipeline Scripts ----
# Identify and set working directory
library(here)
script_dir <- here()
setwd(here("General Pipeline Scripts"))

# Load Libraries
source("Libs for Gen Pipeline.R")

# Init HTTPGD Plotting Visualization
hgd()
hgd_browse()

# Pre-Process Data
source("Data_Wrang_and_PreProcess_Sj.R")
source("Data_Wrang_and_PreProcess_Tm.R")

# Load functions and data frames
source("DESeq2.R")
source("GID Mappings.R")
source("Load S japonicus R list.R")
source("Load T maccoyii R list.R")
source("Load Species R list.R")
source("Pathway Enrichment (by GeneID).R")

# Reset working directory
setwd(script_dir)

# Get GID Mappings to KEGG and GO ----
## Set based on species that you are mapping against
GID_mappings_Sj <- dbQueryGIDMappings("org.Sjaponicus.eg.db")
GID_mappings_Tm <- dbQueryGIDMappings("org.Tmaccoyii.eg.db")

# Set design arg and contrast arg lists ----
design_arg_list <- list(
  pro_combo = ~protocol,
  pro_pass_combo = ~ protocol + passage,
  pass_pro_p77 = ~passage_protocol,
  pass_pro_p27 = ~passage_protocol
)

contrast_arg_list <- list(
  pro_combo = c("protocol", "p", "s"),
  pro_pass_combo = c("protocol", "p", "s"),
  pass_pro_p77 = c("passage_protocol", "p77_p", "p77_s"),
  pass_pro_p27 = c("passage_protocol", "p27_p", "p27_s")
)

# Get DESeq and Pathway Expression Analysis Results ----
## Note: Takes a significant amount of time to run

## Get Results for Scomber japonicus
Sj_analysis <- getResultsLists(
  preProcessOutputList = preProcessOutputList_Sj,
  design_arg_list = design_arg_list,
  contrast_arg_list = contrast_arg_list,
  GID_mappings = GID_mappings_Sj,
  OrgDb = org.Sjaponicus.eg.db,
  pvalueCutoffPathway = 0.01,
  logFCCutoffPathway = 2
)

## Get Results for Thunnus maccoyii
Tm_analysis <- getResultsLists(
  preProcessOutputList = preProcessOutputList_Tm,
  design_arg_list = design_arg_list,
  contrast_arg_list = contrast_arg_list,
  GID_mappings = GID_mappings_Tm,
  OrgDb = org.Tmaccoyii.eg.db,
  pvalueCutoffPathway = 0.01,
  logFCCutoffPathway = 2
)

# Species Comparisons (based on total features mapped) ----
## Number Features with Mean Feature Counts over Threshold by Species ----
### Function ----
numFeatOverThreshold <- function(Sj, Tm, subtitle) {
  Sj_mapped_features <- length(Sj$df_res$baseMean)
  Tm_mapped_features <- length(Tm$df_res$baseMean)

  counts <- data.frame(
    Sj_counts = c(
      sum(Sj$df_res$baseMean > 0) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 1) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 10) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 100) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 1000) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 10000) / Sj_mapped_features * 100,
      sum(Sj$df_res$baseMean > 100000) / Sj_mapped_features * 100
    ),
    Tm_counts = c(
      sum(Tm$df_res$baseMean > 0) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 1) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 10) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 100) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 1000) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 10000) / Tm_mapped_features * 100,
      sum(Tm$df_res$baseMean > 100000) / Tm_mapped_features * 100
    )
  )

  counts$labels <- c(
    "> 0", "> 1", "> 10", "> 100", "> 1000", "> 10000", "> 100000"
  )

  # Reorder levels of the variable column
  counts_long <- reshape2::melt(counts, id.vars = "labels")

  counts_comparison <-
    ggplot(data = counts_long, aes(x = variable, y = value, fill = labels)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
      aes(label = paste0(as.character(round(value, 4)), "%")),
      position = position_dodge(width = 0.9),
      vjust = ifelse(counts_long$value < 1, 1.5, -0.5),
      size = 3.5
    ) +
    labs(
      title = paste0(
        "Number Features with Mean Feature Counts over Threshold by Species (",
        subtitle, ")"
      ),
      x = "Threshold",
      y = "(Number of Features / Number of Features Mapped) x 100",
      fill = "Percent of Features with Mean Feature Counts..."
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_x_discrete(
      labels = c(
        "Sj_counts" = "Scomber japonicus",
        "Tm_counts" = "Thunnus maccoyii"
      )
    )

  counts_comparison
}

### Function Calls ----
# All the same!
numFeatOverThreshold(
  Sj_analysis$desqout_pro_combo,
  Tm_analysis$desqout_pro_combo,
  "Pro for P vs S"
)
numFeatOverThreshold(
  Sj_analysis$desqout_pro_pass_combo,
  Tm_analysis$desqout_pro_pass_combo,
  "Pass + Pro for P vs S"
)
numFeatOverThreshold(
  Sj_analysis$desqout_pass_pro_p77,
  Tm_analysis$desqout_pass_pro_p77,
  "Pass_Pro for p77"
)
numFeatOverThreshold(
  Sj_analysis$desqout_pass_pro_p27,
  Tm_analysis$desqout_pass_pro_p27,
  "Pass_Pro for p27"
)

## Number Receptors by Species ----
# Receptor Information
# Link to Receptor list for Tuna:
# https://www.ncbi.nlm.nih.gov/gene/?term=receptor%5BAll+Fields%5D+AND+%22thunnus+maccoyii%22%5Bporgn%5D+AND+alive%5Bprop%5D

# Link to Receptor list for Japonicus:
# https://www.ncbi.nlm.nih.gov/gene/?term=receptor%5BAll+Fields%5D+AND+%22scomber+japonicus%22%5Bporgn%5D+AND+alive%5Bprop%5D

# Search terms for NCBI genes:
# receptor[All Fields] AND "SPECIES"[porgn] AND alive[prop]

### Functions ----
getReceptors <- function(Sj_df_res_GID, Tm_df_res_GID) {
  # Set working directory
  setwd(paste0(script_dir, "/data"))

  # Define a list of receptor genes to compare against
  Tm_recep_list <- read.csv("thunnus_maccoyii_receptors.csv")
  Sj_recep_list <-
    read.table("scomber_japonicus_receptors.txt", sep = "\t", header = TRUE)

  Tm_receptors <- Tm_df_res_GID |> filter(gene_id %in% Tm_recep_list$GeneID)
  Sj_receptors <- Sj_df_res_GID |> filter(gene_id %in% Sj_recep_list$GeneID)

  # Add in metadata to receptor list
  results <- list(
    Tm_receptors = Tm_receptors <-
      transform(Tm_receptors, gene_id = as.numeric(gene_id)) |>
      rename(GeneID = gene_id) |>
      left_join(Tm_recep_list, by = "GeneID"),
    Sj_receptors_meta = Sj_receptors_meta <-
      transform(Sj_receptors, gene_id = as.numeric(gene_id)) |>
      rename(GeneID = gene_id) |>
      left_join(Sj_recep_list, by = "GeneID")
  )

  # Return upregulated receptors list with metadata
  return(results)
}

numMappings <- function(Sj, Tm, subtitle) {
  Sj_mapped_features <- length(Sj$df_res$baseMean)
  Tm_mapped_features <- length(Tm$df_res$baseMean)

  receptors <- getReceptors(Sj$df_res_GeneID, Tm$df_res_GeneID)

  counts <- data.frame(
    Sj_counts = c(
      nrow(subset(Sj$df_res_GeneID, log2FoldChange > 0 & padj <= 0.01)) /
        Sj_mapped_features * 100,
      nrow(subset(Sj$df_res_GeneID, log2FoldChange < 0 & padj <= 0.01)) /
        Sj_mapped_features * 100,
      nrow(subset(Sj$df_res_GeneID, log2FoldChange == 0 & padj <= 0.01)) /
        Sj_mapped_features * 100,
      nrow(subset(receptors$Sj_receptors, log2FoldChange > 0 &
        padj <= 0.01)) / Sj_mapped_features * 100,
      nrow(subset(receptors$Sj_receptors, log2FoldChange < 0 &
        padj <= 0.01)) / Sj_mapped_features * 100,
      nrow(subset(receptors$Sj_receptors, log2FoldChange == 0 &
        padj <= 0.01)) / Sj_mapped_features * 100
    ),
    Tm_counts = c(
      nrow(subset(Tm$df_res_GeneID, log2FoldChange > 0 & padj <= 0.01)) /
        Tm_mapped_features * 100,
      nrow(subset(Tm$df_res_GeneID, log2FoldChange < 0 & padj <= 0.01)) /
        Tm_mapped_features * 100,
      nrow(subset(Tm$df_res_GeneID, log2FoldChange == 0 & padj <= 0.01)) /
        Tm_mapped_features * 100,
      nrow(subset(receptors$Tm_receptors, log2FoldChange > 0 &
        padj <= 0.01)) / Tm_mapped_features * 100,
      nrow(subset(receptors$Tm_receptors, log2FoldChange < 0 &
        padj <= 0.01)) / Tm_mapped_features * 100,
      nrow(subset(receptors$Tm_receptors, log2FoldChange == 0 &
        padj <= 0.01)) / Tm_mapped_features * 100
    )
  )

  counts$labels <- c(
    "Upregulated Genes", "Downregulated Genes", "Genes with no Change",
    "Upregulated Receptor", "Downregulated Receptor", "Receptors with no Change"
  )

  df <- data.frame(
    Labels = rep(counts$labels, 2),
    Counts = c(counts$Sj_counts, counts$Tm_counts),
    Species = rep(c("Scomber japnoicus", "Thunnus maccoyii"),
      each = length(counts$labels)
    )
  )

  ggplot(df, aes(x = Labels, y = Counts, fill = Species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
      aes(label = paste0(as.character(round(Counts, 4)), "%")),
      position = position_dodge(width = 0.9),
      vjust = -0.5,
      size = 3.5
    ) +
    labs(
      title = paste0(
        "Number Features with Mean Feature Counts over Threshold by Species (",
        subtitle, ")"
      ),
      x = "Features with Significance <= 0.01",
      y = "(Number of Features / Number of Features Mapped) x 100",
      fill = "Species"
    ) +
    theme_minimal()
}

### Function Calls ----
# All the same!
numMappings(
  Sj_analysis$desqout_pro_combo,
  Tm_analysis$desqout_pro_combo,
  "Pro for P vs S"
)
numMappings(
  Sj_analysis$desqout_pro_pass_combo,
  Tm_analysis$desqout_pro_pass_combo,
  "Pass + Pro for P vs S"
)
numMappings(
  Sj_analysis$desqout_pass_pro_p77,
  Tm_analysis$desqout_pass_pro_p77,
  "Pass_Pro for p77"
)
numMappings(
  Sj_analysis$desqout_pass_pro_p27,
  Tm_analysis$desqout_pass_pro_p27,
  "Pass_Pro for p27"
)

# Volcano Plots ----
# Load function to build plots
source(paste0(
  script_dir, "/Investigation 3/Investigation 3 Pipeline Scripts",
  "/createVolcanoPlot.R"
))

# Create volcano plots
volcano_combined_with_passage <-
  createVolcanoPlot(
    desqout_pro_pass_combo$df_res_symbol_GeneID,
    "Mack1 DGE During Proliferation (~ Protocol + Passage)",
    FALSE
  )

volcano_p77 <-
  createVolcanoPlot(
    desqout_pass_pro_p77$df_res_GeneID,
    "Mack1 p77 DGE",
    FALSE
  )

volcano_p27 <-
  createVolcanoPlot(
    desqout_pass_pro_p27$df_res_GeneID,
    "Mack1 p27 DGE",
    FALSE
  )

# Display the plots
volcano_combined_with_passage
volcano_p77
volcano_p27

# Save the plots (have to paste in console)
ggsave("volcano_combined_with_passage.png",
  plot = volcano_combined_with_passage, width = 10, height = 6, dpi = 300
)
ggsave("volcano_p77.png",
  plot = volcano_p77, width = 6, height = 6, dpi = 300
)
ggsave("volcano_p27.png",
  plot = volcano_p27, width = 6, height = 6, dpi = 300
)

# Cluster Heatmap ----
# Load function to build plots
source(paste0(
  script_dir, "/Investigation 3/Investigation 3 Pipeline Scripts",
  "/createHeatmap.R"
))

my_heatmap <- createHeatmap()
my_heatmap

# Save the heatmap (have to paste in console)
ggsave("heatmap.png", plot = my_heatmap, width = 10, height = 6, dpi = 300)

# Bar plots ----

# Dot plots ----

#
