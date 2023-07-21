#' Run Pathway Analysis
#' 
#' This function performs pathway analysis on gene expression data, including
#' GO, KEGG, and KEGG Module pathway analysis.
#' 
#' @param df_res_GID_KEGG The input gene expression data frame.
#' @param pvalueCutoff The p-value cutoff for enrichment analysis. Default is
#'  0.01.
#' @param logFCCutoff The log fold change cutoff for defining
#'  differentially-expressed genes. Default is 2.
#' @param OrgDb An AnnotationHub or AnnotationForge OrgDb object containing 
#'  information on genes and pathways.
#' 
#' @return A list containing the results of the pathway analysis.
#' @export
#' 
#' @examples
#' # example code
runPathwayAnalysis <- function(
    df_res_GID_KEGG,
    pvalueCutoff = 0.01,
    logFCCutoff = 2,
    OrgDb) {
  # Load in data frames
  pathway_dfs <- list(
    GO = df_filtered <-
      df_res_GID_KEGG[complete.cases(df_res_GID_KEGG[, 6]), ],
    KOKEGG = df_filtered_KEGG <-
      df_res_GID_KEGG[complete.cases(df_res_GID_KEGG[, 7]), ] |>
      filter(kegg_ko != "KNA")
  )

  # Init Empty list to fill with results
  output_run_pathway <- list(
    EGO = runPathwayGO(pathway_dfs, OrgDb, pvalueCutoff, logFCCutoff),
    KOKEGG = runPathwayKOKEGG(pathway_dfs, OrgDb, pvalueCutoff, logFCCutoff),
    MKEGG = runPathwayMKEGG(pathway_dfs, OrgDb, pvalueCutoff, logFCCutoff)
  )

  return(output_run_pathway)
}

#' Run GO Over-representation Analysis
#' 
#' This function performs GO over-representation analysis on gene
#' expression data.
#' 
#' @param pathway_dfs A list of data frames containing gene expression data for 
#'  different types of enrichment analysis.
#' @param OrgDb An AnnotationHub or AnnotationForge OrgDb object containing 
#'  information on genes and pathways.
#' @param pvalueCutoff The p-value cutoff for enrichment analysis. 
#'  Default is 0.01.
#' @param logFCCutoff The log fold change cutoff for defining 
#'  differentially-expressed genes. Default is 2.
#' 
#' @return A list containing the results of the KEGG pathway analysis.
#' @export
#' 
#' @examples
#' # example code
runPathwayGO <- function(
    pathway_dfs,
    OrgDb,
    pvalueCutoff,
    logFCCutoff) {
  # Note that this is being run
  print("Running GO Over-representation Analysis")
  
  # Create ego List
  egoList <- list()
  
  for (name in names(pathway_dfs)) {
    # Collect GO Results
    if (name == "GO") {
      # Set selection criteria
      pathway_criteria <- list(
        all = pathway_dfs[[name]],
        up = pathway_dfs[[name]] |> filter(log2FoldChange > 0),
        down = pathway_dfs[[name]] |> filter(log2FoldChange < 0)
      )
      
      # Create geneList for universe
      geneListUniverse <- pathway_criteria$all[, 2]
      names(geneListUniverse) <- as.character(pathway_criteria$all[, 6])
      geneListUniverse <- sort(geneListUniverse, decreasing = TRUE)

      # Collect GO Results for each selection criteria
      for (subset_name in names(pathway_criteria)) {
        # Define differentially-expressed genes
        geneList <- pathway_criteria[[subset_name]][, 2]
        names(geneList) <- as.character(pathway_criteria[[subset_name]][, 6])
        geneList <- sort(geneList, decreasing = TRUE)
        de <- names(geneList)[abs(geneList) > logFCCutoff]

        # GO over-representation analyses
        egoALL <- enrichGO(
          gene = de,
          OrgDb = OrgDb,
          ont = "ALL",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          universe = names(geneListUniverse),
          readable = TRUE
        )

        egoCC <- enrichGO(
          gene = de,
          OrgDb = OrgDb,
          ont = "CC",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          universe = names(geneListUniverse),
          readable = TRUE
        )

        egoBP <- enrichGO(
          gene = de,
          OrgDb = OrgDb,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          universe = names(geneListUniverse),
          readable = TRUE
        )

        egoMF <- enrichGO(
          gene = de,
          OrgDb = OrgDb,
          ont = "MF",
          pAdjustMethod = "BH",
          pvalueCutoff = pvalueCutoff,
          universe = names(geneListUniverse),
          readable = TRUE
        )

        # Save results to output
        new_name <- paste("egoALL", name, subset_name, sep = "_")
        assign(new_name, egoALL, envir = .GlobalEnv)
        egoList[[new_name]] <- get(new_name)
        rm(egoALL)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)

        new_name <- paste("egoCC", name, subset_name, sep = "_")
        assign(new_name, egoCC, envir = .GlobalEnv)
        egoList[[new_name]] <- get(new_name)
        rm(egoCC)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)

        new_name <- paste("egoBP", name, subset_name, sep = "_")
        assign(new_name, egoBP, envir = .GlobalEnv)
        egoList[[new_name]] <- get(new_name)
        rm(egoBP)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)

        new_name <- paste("egoMF", name, subset_name, sep = "_")
        assign(new_name, egoMF, envir = .GlobalEnv)
        egoList[[new_name]] <- get(new_name)
        rm(egoMF)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)
      }
      
      # Save universe geneList
      name <- "geneListUniverse"
      egoList[[name]] <- geneListUniverse
    }
  }
  return(egoList)
}

#' Run KEGG Pathway Over-representation Analysis
#' 
#' This function performs KEGG pathway over-representation analysis on gene
#' expression data.
#' 
#' @param pathway_dfs A list of data frames containing gene expression data for 
#'  different types of enrichment analysis.
#' @param OrgDb An AnnotationHub or AnnotationForge OrgDb object containing 
#'  information on genes and pathways.
#' @param pvalueCutoff The p-value cutoff for enrichment analysis. 
#'  Default is 0.01.
#' @param logFCCutoff The log fold change cutoff for defining 
#'  differentially-expressed genes. Default is 2.
#' 
#' @return A list containing the results of the KEGG pathway analysis.
#' @export
#' 
#' @examples
#' # example code
runPathwayKOKEGG <- function(
    pathway_dfs,
    OrgDb,
    pvalueCutoff,
    logFCCutoff) {
  # Note that this is being run
  print("Running KEGG Pathway Over-representation Analysis")
  
  # Create KEGG ko List
  kokeggList <- list()
  
  for (name in names(pathway_dfs)) {
    # Collect KEGG ko Results
    if (name == "KOKEGG") {
      # Set selection criteria
      pathway_criteria <- list(
        all = pathway_dfs[[name]],
        up = pathway_dfs[[name]] |> filter(log2FoldChange > 0),
        down = pathway_dfs[[name]] |> filter(log2FoldChange < 0)
      )
      
      # Create geneList for universe
      geneListUniverse <- pathway_criteria$all[, 2]
      names(geneListUniverse) <- as.character(pathway_criteria$all[, 7])
      geneListUniverse <- sort(geneListUniverse, decreasing = TRUE)

      # Collect KEGG Results for each selection criteria
      for (subset_name in names(pathway_criteria)) {
        # Define differentially-expressed genes
        geneList <- pathway_criteria[[subset_name]][, 2]
        names(geneList) <- as.character(pathway_criteria[[subset_name]][, 7])
        geneList <- sort(geneList, decreasing = TRUE)
        de <- names(geneList)[abs(geneList) > logFCCutoff]

        # KEGG pathway over-representation analysis
        ko <- enrichKEGG(
          gene = de,
          organism = "ko",
          keyType = "kegg",
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = "BH",
          universe = names(geneListUniverse)
        )

        # Save results to output
        new_name <- paste("koKEGG", subset_name, sep = "_")
        assign(new_name, ko, envir = .GlobalEnv)
        kokeggList[[new_name]] <- get(new_name)
        rm(ko)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)
      }
      
      # Save universe geneList
      name <- "geneListUniverse"
      kokeggList[[name]] <- geneListUniverse
    }
  }
  return(kokeggList)
}

#' Run KEGG Module Over-representation Analysis
#' 
#' This function performs KEGG module over-representation analysis on gene
#' expression data.
#' 
#' @param pathway_dfs A list of data frames containing gene expression data for 
#'  different types of enrichment analysis.
#' @param OrgDb An AnnotationHub or AnnotationForge OrgDb object containing 
#'  information on genes and pathways.
#' @param pvalueCutoff The p-value cutoff for enrichment analysis. 
#'  Default is 0.01.
#' @param logFCCutoff The log fold change cutoff for defining 
#'  differentially-expressed genes. Default is 2.
#' 
#' @return A list containing the results of the KEGG pathway analysis.
#' @export
#' 
#' @examples
#' # example code
runPathwayMKEGG <- function(
    pathway_dfs,
    OrgDb,
    pvalueCutoff,
    logFCCutoff) {
  # Note that this is being run
  print("Running KEGG Module Over-representation Analysis")
  
  # Create KEGG Module List
  mkeggList <- list()
  
  for (name in names(pathway_dfs)) {
    # Collect KEGG Module Results
    if (name == "MKEGG") {
      # Set selection criteria
      pathway_criteria <- list(
        all = pathway_dfs[[name]],
        up = pathway_dfs[[name]] |> filter(log2FoldChange > 0),
        down = pathway_dfs[[name]] |> filter(log2FoldChange < 0)
      )
      
      # Create geneList for universe
      geneListUniverse <- pathway_criteria$all[, 2]
      names(geneListUniverse) <- as.character(pathway_criteria$all[, 7])
      geneListUniverse <- sort(geneListUniverse, decreasing = TRUE)

      # Collect KEGG Results for each selection criteria
      for (subset_name in names(pathway_criteria)) {
        # Define differentially-expressed genes
        geneList <- pathway_criteria[[subset_name]][, 2]
        names(geneList) <- as.character(pathway_criteria[[subset_name]][, 7])
        geneList <- sort(geneList, decreasing = TRUE)
        de <- names(geneList)[abs(geneList) > logFCCutoff]

        # KEGG module over-representation analysis
        mKEGG <- enrichMKEGG(
          gene = de,
          organism = "ko",
          keyType = "kegg",
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = "BH",
          universe = names(geneListUniverse)
        )

        # Save results to output
        new_name <- paste("mKEGG", subset_name, sep = "_")
        assign(new_name, mKEGG, envir = .GlobalEnv)
        mkeggList[[new_name]] <- get(new_name)
        rm(mKEGG)
        rm(list = new_name, envir = .GlobalEnv)
        rm(new_name)
      }
      
      # Save universe geneList
      name <- "geneListUniverse"
      mkeggList[[name]] <- geneListUniverse
    }
  }
  return(mkeggList)
}
