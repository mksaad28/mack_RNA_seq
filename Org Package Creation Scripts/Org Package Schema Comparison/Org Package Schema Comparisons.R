#' Org Package Schema Comparisons
#' 6/28/23
#' 
#' Written by: Benji Bromberg
#' 
#' Purpose:
#'  To compare user-generated databases against the curated databases hosted
#'  on AnnotationHub. See `Schema Comparison Notes` for notes.

# Identify and set working directory
library(here)
script_dir <- here()
setwd(here("General Pipeline Scripts"))

# Load all required packages
source("Libs for Gen Pipeline.R")
source("Load Species R list.R")

# Scomber japonicus (NCBI) ----

# # Load Scomber japonicus NCBI Org Package
# script_dir <- 
#   dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# setwd(paste0(script_dir, "/Japonicus/SJaponicusFromNCBI/"))
# install.packages("org.Sjaponicus.eg.db/", repos = NULL, type = "source")
# 
# library(org.Sjaponicus.eg.db)
# tables = dbListTables(org.Sjaponicus.eg_dbconn())
# 
# x = 1;
# Sj_NCBI_list = list() # create an empty list to store queries
# 
# while (x <= length(tables)) {
#   # Query tables[x]
#   query_result = dbGetQuery(org.Sjaponicus.eg_dbconn(), paste("SELECT * FROM ",
#                                                             noquote(tables[x])))
# 
#   # Add result to list with variable name equal to table name
#   Sj_NCBI_list[[tables[x]]] = as.data.frame(query_result)
# 
#   x = x + 1
# }

# Scomber japonicus (Eggnog) ----

# # Load Scomber japonicus Eggnog Org Package
# script_dir <- 
#   dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# setwd(paste0(script_dir, "/Japonicus/SJaponicusFromEggnog/"))
# install.packages("org.Sjaponicus.eg.db/", repos = NULL, type = "source")

library(org.Sjaponicus.eg.db)
tables = dbListTables(org.Sjaponicus.eg_dbconn())

x = 1;
Sj_Eggnog_list = list() # create an empty list to store queries

while (x <= length(tables)) {
  # Query tables[x]
  query_result = dbGetQuery(org.Sjaponicus.eg_dbconn(), paste("SELECT * FROM ", 
                                                            noquote(tables[x])))
  
  # Add result to list with variable name equal to table name
  Sj_Eggnog_list[[tables[x]]] = as.data.frame(query_result)
  
  x = x + 1
}

# Scomber colias (NCBI) ----

# # Load Scomber colias NCBI Org Package
# script_dir <- 
#   dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# setwd(paste0(script_dir, "/Colias/SColiasFromNCBI/"))
# install.packages("org.Scolias.eg.db/", repos = NULL, type = "source")
# 
# library(org.Scolias.eg.db)
# tables = dbListTables(org.Scolias.eg_dbconn())
# 
# x = 1;
# Sc_NCBI_list = list() # create an empty list to store queries
# 
# while (x <= length(tables)) {
#   # Query tables[x]
#   query_result = dbGetQuery(org.Scolias.eg_dbconn(), paste("SELECT * FROM ", 
#                                                            noquote(tables[x])))
#   
#   # Add result to list with variable name equal to table name
#   Sc_NCBI_list[[tables[x]]] = as.data.frame(query_result)
#   
#   x = x + 1
# }

# Scomber colias (Eggnog) ----

# # Load Scomber colias Eggnog Org Package
# script_dir <- 
#   dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# setwd(paste0(script_dir, "/Colias/SColiasFromEggnog/"))
# install.packages("org.Scolias.eg.db/", repos = NULL, type = "source")

library(org.Scolias.eg.db)
tables = dbListTables(org.Scolias.eg_dbconn())

x = 1;
Sc_Eggnog_list = list() # create an empty list to store queries

while (x <= length(tables)) {
  # Query tables[x]
  query_result = dbGetQuery(org.Scolias.eg_dbconn(), paste("SELECT * FROM ",
                                                           noquote(tables[x])))

  # Add result to list with variable name equal to table name
  Sc_Eggnog_list[[tables[x]]] = as.data.frame(query_result)

  x = x + 1
}

# C elegans
Ce_List <- loadSpeciesList("org.Ce.eg.db")

# M musculus
Mm_list <- loadSpeciesList("org.Mm.eg.db")

# H sapiens
Hs_list <- loadSpeciesList("org.Hs.eg.db")

# S scrofa
Ss_list <- loadSpeciesList("org.Ss.eg.db")

# Chicken
Gg_list <- loadSpeciesList("org.Gg.eg.db")

# Canine
Cf_list <- loadSpeciesList("org.Cf.eg.db")

# Rat
Rn_list <- loadSpeciesList("org.Rn.eg.db")

# Bovine
Bt_list <- loadSpeciesList("org.Bt.eg.db")

# Xenopus
Xl_list <- loadSpeciesList("org.Xl.eg.db")

# Anopheles
Ag_list <- loadSpeciesList("org.Ag.eg.db")

# Fly
Dm_list <- loadSpeciesList("org.Dm.eg.db")

