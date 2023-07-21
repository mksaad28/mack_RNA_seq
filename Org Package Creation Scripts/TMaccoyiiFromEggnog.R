#' TMaccoyiiFromEggnog
#' 7/11/23
#'
#' Written by: Benji Bromberg
#'
#' Purpose:
#'  TODO

# Load all required packages ----
script_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(script_dir)
source("General Pipeline Scripts/Libs for Gen Pipeline.R")
#' Note: Expect an error here as the orgDB is not loaded. This will stop the
#' script. You will have to manually continue it.
library(AnnotationHub) # Bioconductor package

# Load Thunnus maccoyii Org Package (NCBI) ----
#' Collect orgdb from the AnnotationHub since it is available
#' Additionally, create an R list() object of the data in the Org package
ah <- AnnotationHub()

# Search for AnnotationHub for Thunnus maccoyii TaxID
query(ah, pattern = "8240")

res <- ah[["AH112021"]] # Note this query will likely need to be updated
conn <- dbconn(res)

tables <- dbListTables(conn)

x <- 1
tm_list <- list() # create an empty list to store queries

while (x <= length(tables)) {
  # Query tables[x]
  query_result <- dbGetQuery(conn, paste("SELECT * FROM ", noquote(tables[x])))

  # Add result to list with variable name equal to table name
  tm_list[[tables[x]]] <- as.data.frame(query_result)

  x <- x + 1
}

# Clean up unnecessary variables
rm(query_result, x, tables)

# makeOrgPackage Function Prep ----
#' Create data frames for data that will be transferred from NCBI Org package
#' to manually created Org package.
#'
#' All tables must be in terms of "GID" not "_id"

GID2accessions <- merge(tm_list$accessions, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "ACCNUM")]

GID2alias <- merge(tm_list$alias, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "ALIAS")]

GID2chromosomes <- merge(tm_list$chromosome, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "CHR")]

GID2entrez <- merge(tm_list$entrez_genes, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "ENTREZID")]

GID2info <- merge(tm_list$gene_info, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "GENENAME", "SYMBOL")]

GID2refseq <- merge(tm_list$refseq, tm_list$genes,
  by.x = "_id", by.y = "_id",
  all.x = TRUE
)[, c("GID", "REFSEQ")]

# Collect refseq accessions ----
#' refseq accessions -> BatchEntrez (protein) -> fasta results ->
#' Eggnog Mapper -> GO and Kegg Annotations
#'
#' Expecting XM_ files to not be read by BatchEntrez
#'
#' Kept running BatchEntrez over and over until the # of retrived records
#' no longer grew larger (49507 UIDs)
#'
#' Downloaded fasta results multiple times and took largest result as
#' sometimes the files cut off during download

allrefseq <- GID2refseq[, c("REFSEQ")]

allrefseq <- gsub("\\.1$", "", allrefseq)

# Output file path
setwd(paste0(
  script_dir,
  "/Maccoyii/TMaccoyiiFromEggnog/Eggnog Mappings/Eggnog Input/"
))
output_file <- "batchentrez_refseqs_query.txt"

# Write refseq accessions to file
writeLines(allrefseq, output_file)

# eggNOG Mapping ----
#' Collect FASTA sequences into a single .fasta file and process online at:
#' http://eggnog-mapper.embl.de/
#'
#' Used default settings.
#'
#' Link (may be deprecated):
#' http://eggnog-mapper.embl.de/job_status?jobname=MM_1s2xnstt

# Process eggNOG Mapping Data ----
#' eggNOG mapper output collected as .emapper.annotations.xlsx file
eggnog <- readxl::read_xlsx(
  paste0(
    script_dir, "/Maccoyii/TMaccoyiiFromEggno",
    "g/Eggnog Mappings/Eggnog Output/MM_1s2xnstt.emapper.annotations.xlsx"
  ),
  skip = 2
)

# Save GO and Kegg annotation data as tables to add to the Org package
#' Most Org packages only seem to have Kegg Pathway annotation data (unsure
#' why this is exactly)
GID2GO <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "GOs")] |>
  unique() |>
  separate_rows(GOs, sep = ",") |>
  mutate(EVIDENCE = ifelse(!is.na(GOs), "IEA", NA)) |>
  #' GO Evidence Codes set to IEA for all GO annotations since they were
  #' made using eggNOG Mapper
  unique() |>
  rename(GO = GOs) |>
  na.omit() |>
  filter(!grepl("-", GO))

GID2BRITE <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "BRITE")] |>
  unique() |>
  separate_rows(BRITE, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", BRITE))

GID2KEGG_ko <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_ko")] |>
  unique() |>
  separate_rows(KEGG_ko, sep = ",") |>
  mutate(KEGG_ko = gsub("ko:K", "", KEGG_ko)) |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_ko))

GID2KEGG_Pathway <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_Pathway")] |>
  unique() |>
  separate_rows(KEGG_Pathway, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_Pathway))

GID2KEGG_Module <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_Module")] |>
  unique() |>
  separate_rows(KEGG_Module, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_Module))

GID2KEGG_Reaction <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_Reaction")] |>
  unique() |>
  separate_rows(KEGG_Reaction, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_Reaction))

GID2KEGG_rclass <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_rclass")] |>
  unique() |>
  separate_rows(KEGG_rclass, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_rclass))

GID2KEGG_TC <- merge(GID2refseq, eggnog,
  by.x = "REFSEQ", by.y = "query",
  all.x = TRUE
)[, c("GID", "KEGG_TC")] |>
  unique() |>
  separate_rows(KEGG_TC, sep = ",") |>
  unique() |>
  na.omit() |>
  filter(!grepl("-", KEGG_TC))

# Check for duplicates
sum(duplicated(eggnog))
# [1] 13
#' This may be because the 13 "YP_" sequences are actually processed by
#' BatchEntrez.
sum(duplicated(GID2refseq))
# [1] 0
sum(duplicated(GID2GO))
# [1] 0
sum(duplicated(GID2BRITE))
# [1] 0

# Create Org package ----
# Remove NCBI orgDB package
remove.packages(pkgs = "org.Tmaccoyii.eg.db")

# Restart R session
# Note: This will stop the script. You will have to manually continue it.
rstudioapi::restartSession()

# Re-init Packages and directories
script_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(script_dir)
source("General Pipeline Scripts/Libs for Gen Pipeline.R")
#' Note: Expect an error here as the orgDB is not loaded. This will stop the
#' script. You will have to manually continue it.
library(AnnotationHub) # Bioconductor package

# Reset to correct working directory
# note: Setting the correct working directory is critical for a valid install
setwd(paste0(script_dir, "/Maccoyii/TMaccoyiiFromEggnog"))

# Create orgDB
makeOrgPackage(
  accessions = GID2accessions,
  alias = GID2alias,
  brite = GID2BRITE,
  chromosome = GID2chromosomes,
  entrez_genes = GID2entrez,
  gene_info = GID2info,
  go = GID2GO,
  refseq = GID2refseq,
  kegg_ko = GID2KEGG_ko,
  kegg_module = GID2KEGG_Module,
  kegg_pathway = GID2KEGG_Pathway,
  kegg_rclass = GID2KEGG_rclass,
  kegg_reaction = GID2KEGG_Reaction,
  kegg_tc = GID2KEGG_TC,
  version = "0.1",
  maintainer = "Benjamin Bromberg <bromberg.benji@gmail.com>",
  author = "Benjamin Bromberg <bromberg.benji@gmail.com>",
  outputDir = ".",
  tax_id = "8240",
  genus = "Thunnus",
  species = "maccoyii",
  goTable = "go"
)

install.packages("org.Tmaccoyii.eg.db/", repos = NULL, type = "source")

# Validate proper Org package creation ----
library(org.Tmaccoyii.eg.db)

# Collect org.Tmaccoyii.eg.db as an R list() to visualize proper install
tables <- dbListTables(org.Tmaccoyii.eg_dbconn())

x <- 1
tm_list <- list() # create an empty list to store queries

while (x <= length(tables)) {
  # Query tables[x]
  query_result <- dbGetQuery(org.Tmaccoyii.eg_dbconn(), paste(
    "SELECT * FROM ",
    noquote(tables[x])
  ))

  # Add result to list with variable name equal to table name
  tm_list[[tables[x]]] <- as.data.frame(query_result)

  x <- x + 1
}

rm(x, query_result, tables)
