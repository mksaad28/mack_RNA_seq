#' Load Species Data Tables
#'
#' Load all tables from the the given orgDB database into a list of data frames.
#' Each data frame is stored as an element in the list with the table name
#' as the variable name.
#' @param orgDB A string representation of an orgDB database.
#'
#' @return A list of data frames containing the orgDB data tables.
#' @export
#' @examples
#' loadSpeciesList("org.Hs.eg.db")
loadSpeciesList <- function(orgDB) {
  # Check if library is installed
  if (!requireNamespace(orgDB, quietly = TRUE)) {
    # Install the library
    BiocManager::install(orgDB)
  }
  
  # Load library
  lib_call <- paste0('library(', orgDB, ')')
  eval(parse(text = lib_call))
  
  # Convert orgDB into function call
  func_call <- paste0(gsub("\\.db", "", orgDB), '_dbconn()')
  
  tables <- dbListTables(eval(parse(text = func_call)))
  
  x <- 1
  list <- list() # create an empty list to store queries

  while (x <= length(tables)) {
    # Query tables[x]
    query_result <- dbGetQuery(
      eval(parse(text = func_call)),
      paste("SELECT * FROM ", noquote(tables[x]))
    )

    # Add result to list with variable name equal to table name
    list[[tables[x]]] <- as.data.frame(query_result)

    x <- x + 1
  }

  return(list)
}
