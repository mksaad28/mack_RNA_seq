write_to_db <- function(filename, data_long, int_limma_res, kegg_categories, genes_and_brite, int_enricher_res_tbl, pca_expr) {
  if ("gene_id" %in% colnames(int_limma_res)) {
    int_limma_res <- int_limma_res |>
      mutate(gene_name = gene_id)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = filename)

  copy_to(con,
    data_long,
    "data_long",
    temporary = FALSE,
    overwrite = TRUE,
    indexes = list(
      "sample",
      "gene_id"
    )
  )

  copy_to(con,
    int_limma_res,
    "res",
    temporary = FALSE,
    overwrite = TRUE,
    indexes = list(
      "gene_id"
    )
  )

  copy_to(con,
    kegg_categories,
    "kegg_categories",
    temporary = FALSE,
    overwrite = TRUE,
    indexes = list(
      "category",
      "brite"
    )
  )


  copy_to(con,
    genes_and_brite,
    "genes_and_brite",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "brite"
    )
  )

  copy_to(con,
    int_enricher_res_tbl,
    "gene_ontologies",
    temporary = FALSE,
    overwrite = TRUE
  )

  copy_to(con,
    pca_expr,
    "pca",
    temporary = FALSE,
    overwrite = TRUE
  )


  DBI::dbDisconnect(con)
  # Saving the environment to made intermediate objects available for further
  # analysis by addressing the returned object, i.e. ret$envir$expr_fit for the
  # lmFit.
  return(list(
    filename = filename,
    data_long = data_long,
    int_limma_res = int_limma_res,
    kegg_categories = int_limma_res,
    genes_and_brite = genes_and_brite,
    int_enricher_res_tbl = int_enricher_res_tbl,
    pca_expr = pca_expr
  ))
}
