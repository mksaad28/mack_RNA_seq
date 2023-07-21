#' GID Mappings
#' 6/20/23
#'
#' Written by: Benji Bromberg
#'
#' Functions Implemented:
#'  dbQueryGIDMappings()
#'  addMappings2DESeqRes()
#'
#' Purpose:
#'  At the beginnning of this pipeline, the featurecounts results are sorted by
#'  Gene Symbol rather than GeneID. To process the DESeq2 results downstream
#'  using the CLusterProfiler package, Gene Symbols must be mapped to their
#'  proper GeneIDs. This is the functionality of addMappings2DESeqRes().
#'  dbQueryGIDMappings() is a helper function.
#'
#'  Note:
#'  Many tRNAs map to many different GeneIDs. Our featurecounts data does not
#'  account for this, and approximately 1,300 entries are added when mapping
#'  occurs. To compensate for this, add_GeneID_to_df() filters out duplicate
#'  entries by Gene Symbol and keeps only one entry per duplicate. This means
#'  that each tRNA is mapped to only one GeneID, but the GeneID it has been
#'  mapped to has been arbitrarily chosen.

#' Add GID and KEGG GID Mappings to DESeq2 Results by Gene Symbol
#'
#' This function takes a data frame containing DESeq2 results, where the row
#' names are Gene Symbols, and adds mappings to Gene IDs (GIDs) and KEGG Gene
#' IDs (KEGG GIDs) using a provided mapping table. Duplicate entries for Gene
#' Symbols are filtered out, keeping only one entry per duplicate. The function
#' returns a list of data frames with the following mappings: Gene Symbol, GID,
#' and KEGG GID.
#'
#' @param df_from_DESeq A data frame containing DESeq2 results, where the row
#'  names are Gene Symbols.
#' @param GID_Mappings A list with two data frames: gene_id_symbol (mapping of
#'  Gene Symbols to GIDs) and gene_id_kegg (mapping of GIDs to KEGG GID).
#'
#' @return A list of data frames with the following mappings: Gene Symbol, GID,
#' and KEGG GID. The list includes the following data frames: df_res_symbol,
#' df_res_symbol_GeneID, df_res_GeneID, and df_res_GeneID_kegg.
#'
#' @examples
#' # example code
#'
#' @seealso
#'   \code{\link{dbQueryGIDMappings}}
#'
#' @importFrom dplyr left_join select
#' @importFrom tibble rownames_to_column
#'
#' @export
addMappings2DESeqRes <- function(df_from_DESeq, GID_Mappings) {
  # Convert Gene Symbol rownames to a column
  df_res_symbol <- df_from_DESeq |>
    rownames_to_column(var = "Symbol")

  # Map Gene Symbols to GIDs using gene_id_symbol as a reference
  df_res_symbol_GeneID <- df_res_symbol |>
    left_join(GID_Mappings$gene_id_symbol, by = "Symbol")

  # For duplicate Gene Symbol entries, keep only one entry per duplicate
  df_res_symbol_GeneID <-
    df_res_symbol_GeneID[!duplicated(df_res_symbol_GeneID$Symbol), ]

  # Create a version of $df_res only mapped to GIDs
  df_res_GeneID <- df_res_symbol_GeneID |> dplyr::select(-Symbol)

  # Map GIDs to KEGG GIDs using gene_id_kegg as a reference
  df_res_GeneID_kegg <- df_res_GeneID |>
    left_join(GID_Mappings$gene_id_kegg, by = "gene_id")

  # Make sure no genes are lost or gained
  print(paste0("The number of rows using gene symbol is ", nrow(df_res_symbol)))
  print(paste0("The number of rows using gene ID is ", nrow(df_res_GeneID)))

  # Output results to list and return results
  output_df <- list(
    df_res_symbol = df_res_symbol,
    df_res_symbol_GeneID = df_res_symbol_GeneID,
    df_res_GeneID = df_res_GeneID,
    df_res_GeneID_kegg =
      df_res_GeneID_kegg
  )
  return(output_df)
}

#' Query Gene Symbol Mappings to GIDs and KEGG GIDs
#'
#' This function queries the Gene ID (GID) mappings from a specified organism
#' database (OrgDb) and returns a list of data frames containing the mappings.
#' The mappings include Gene Symbols to GIDs and GIDs to KEGG Gene IDs
#' (KEGG GID).
#'
#' @param OrgDb The organism database to query. This should be a valid
#' Bioconductor package name (e.g., "org.Hs.eg.db" for Homo sapiens).
#'
#' @return A list containing two data frames:
#'   - gene_id_symbol: A data frame containing the mappings of Gene Symbols to
#'   GIDs.
#'   - gene_id_kegg: A data frame containing the mappings of GIDs to KEGG GIDs.
#'
#' @examples
#' db_mappings <- dbQueryGIDMappings("org.Hs.eg.db")
#'
#' @seealso
#'   \code{\link{addMappings2DESeqRes}}
#'
#' @importFrom BiocManager BiocManager
#' @importFrom dplyr left_join select rename mutate
#' @export
dbQueryGIDMappings <- function(OrgDb) {
  # Check if library is installed
  if (!requireNamespace(OrgDb, quietly = TRUE)) {
    # Install the library
    BiocManager::install(OrgDb)
  }

  # Load library
  lib_call <- paste0("library(", OrgDb, ")")
  eval(parse(text = lib_call))

  # Convert OrgDb into function call
  org_conn_call <- paste0(gsub("\\.db", "", OrgDb), "_dbconn()")

  # Query the OrgDb
  genes <- dbGetQuery(
    eval(parse(text = org_conn_call)),
    "SELECT * FROM genes"
  )
  gene_info <- dbGetQuery(
    eval(parse(text = org_conn_call)),
    "SELECT * FROM gene_info"
  )
  kegg_ko <- dbGetQuery(
    eval(parse(text = org_conn_call)),
    "SELECT * FROM kegg_ko"
  )

  # Join queries
  results <- list(
    gene_id_symbol = gene_id_symbol <- merge(
      genes, gene_info,
      by.x = "_id", by.y = "_id", all.x = TRUE
    )[, c("GID", "SYMBOL")] |>
      rename(Symbol = SYMBOL, gene_id = GID),
    gene_id_kegg = gene_id_kegg <-
      merge(
        genes, kegg_ko,
        by.x = "_id", by.y = "_id", all.x = TRUE
      )[, c("GID", "KEGG_ko")] |>
      rename(kegg_ko = KEGG_ko, gene_id = GID) |>
      mutate(kegg_ko = paste0("K", kegg_ko))
  )

  return(results)
}
