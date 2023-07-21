#' Load Scomber japonicus Data Tables
#'
#' Load all tables from the Sjaponicus database into a list of data frames.
#' Each data frame is stored as an element in the list with the table name
#' as the variable name.
#'
#' @return A list of data frames containing the Sjaponicus data tables.
#' @export
#' @examples
#' loadSJList()
loadSJList <- function() {
  tables <- dbListTables(org.Sjaponicus.eg_dbconn())

  x <- 1
  sj_list <- list() # create an empty list to store queries

  while (x <= length(tables)) {
    # Query tables[x]
    query_result <- dbGetQuery(
      org.Sjaponicus.eg_dbconn(),
      paste("SELECT * FROM ", noquote(tables[x]))
    )

    # Add result to list with variable name equal to table name
    sj_list[[tables[x]]] <- as.data.frame(query_result)

    x <- x + 1
  }

  return(sj_list)
}
