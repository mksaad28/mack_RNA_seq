#' Load Thunnus maccoyii Data Tables
#'
#' Load all tables from the Tmaccoyii database into a list of data frames.
#' Each data frame is stored as an element in the list with the table name
#' as the variable name.
#'
#' @return A list of data frames containing the Tmaccoyii data tables.
#' @export
#' @examples
#' loadTMList()
loadTMList <- function() {
  tables <- dbListTables(org.Tmaccoyii.eg_dbconn())

  x <- 1
  tm_list <- list() # create an empty list to store queries

  while (x <= length(tables)) {
    # Query tables[x]
    query_result <- dbGetQuery(
      org.Tmaccoyii.eg_dbconn(),
      paste("SELECT * FROM ", noquote(tables[x]))
    )

    # Add result to list with variable name equal to table name
    tm_list[[tables[x]]] <- as.data.frame(query_result)

    x <- x + 1
  }

  return(tm_list)
}
