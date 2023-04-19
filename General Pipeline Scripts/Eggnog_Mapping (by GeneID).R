## Eggnog Mapping

## Necessary argument from orthology because neither tuna or mackerel genomes
## are in handy databases like ensembl, gene ontology, or kegg. First we read
## the cds from GCF_910596095.1 and extract the cds accession and the gene id.
## This is necessary because the cds from GCF_910596095.1 were fed into
## eggnogg-mapper. This file was downloaded from NCBI genomes. 

a <- 
  readLines(paste("data/michael_saad/GCF_910596095.1/cds_from_genomic", 
                  ".accessions_only.fna", sep="")) # Read in fasta file

query_to_gene <- tibble()
# Iterates through the fasta file
for (line in a) {                            
  # Processes lines that begin with ">"
  if (!str_starts(line, ">")) next            
  
  # Splits line into vector, based on space
  line_split <- str_split(line, " ")[[1]]     
  accession <- str_remove(line_split[1], ">")  
  # Stores accession as the first element in vector, after removing ">"
  
  if (!str_detect(line_split[2], "gene")) 
    # Ensures there is gene in the second element in vector
    stop("Missing gene symbol in second field") 
  
  gene_symbol <- line_split[2] |>  
    # Pulls out gene_symbol by removing string parts from the 2nd element
    str_remove(fixed("[gene=")) |>
    str_remove(fixed("]"))
  
  if (!str_detect(line_split[3], "db_xref=GeneID")) 
    stop("Missing gene_id in third field")
  GeneID <- line_split[3] |>
    str_remove(fixed("[db_xref=GeneID:")) |>
    str_remove(fixed("]"))
  
  query_to_gene <- bind_rows(query_to_gene, bind_cols(
    query = accession, gene_symbol = gene_symbol,gene_id = GeneID)) 
  # Stores the accession, gene_symbol, and gene_id in query_to_gene
  # Same number of unique values in gene_symbol and gene_id columns
}  


#' Then we read the eggnogg mapper excel file resulting from querying it with
#' GCF_910596095.1/cds_from_genomic.fna and the default settings on
#' http://eggnog-mapper.embl.de
#' 
#' Changed to eggnogg_ID
eggnogg_ID <- 
  readxl::read_xlsx("data/michael_saad/MM_8hli8sc6.emapper.annotations.xlsx", 
                    skip = 2) |>
  left_join(query_to_gene, by = "query") |>
  select(gene_id, everything()) |>
  group_by(gene_id) |>
  slice_head(n = 1)

# Old version:
# eggnogg_symbol <- 
#   readxl::read_xlsx("data/michael_saad/MM_8hli8sc6.emapper.annotations.xlsx", 
#                     skip = 2) |>
#   left_join(query_to_gene, by = "query") |>
#   select(gene_symbol, everything()) |>
#   group_by(gene_symbol) |>
#   slice_head(n = 1)

#' this is a little rough, taking only the first representative per gene when
#' there are multiple transcripts/proteins in the cds dataset. To me (Michael),
#' it seems like the highest scores are always first, so likely the most 
#' accurate mapping?

# Break the eggnogg results into long format
genes_and_brite <- eggnogg_ID |>
  select(gene_id, BRITE) |>
  rename(brite = BRITE) |>
  rowwise() |>
  mutate(brite = str_split(brite, ",")) |>
  unnest(cols = c(brite))

# Delete this or genes_and_brite (whichever isn't used)
term_to_brite <- genes_and_brite[c('brite','gene_id')]

# get gene ontologies out of eggnogg.
term_to_gene <- eggnogg_ID |>
  select(gene_id, GOs) |>
  filter(GOs != "-") |>
  separate_rows(GOs, sep = ",") |>
  rename(GO = GOs) |>
  select(GO, gene_id)

# kegg categories from yanky br08902 file.
kegg_categories <- jsonlite::read_json("data/br08902.json") |>
  unlist() |>
  as_tibble() |>
  separate(value, into = c("brite", NA), remove = FALSE, extra = "drop") |>
  mutate(brite = str_replace(brite, "Protein", NA_character_)) |>
  rename(category = value)