## Eggnog Mapping

## Necessary argument from orthology because neither tuna or mackerel genomes
## are in handy databases like ensembl, gene ontology, or kegg. First we read
## the cds from GCF_910596095.1 and extract the cds accession and the gene id.
## This is necessary because the cds from GCF_910596095.1 were fed into
## eggnogg-mapper. This file was downloaded from NCBI genomes. 

gtf <- read.csv(
  'data/michael_saad/GCF_910596095.1_fThuMac1.1_genomic.gtf',
  sep="\t",
  skip=4,
  header = F
)

# remove last row
gtf <- gtf[-nrow(gtf),]

# isolate GeneIDs
GeneIDs <- stringr::str_match(gtf$V9, "xref GeneID\\s*(.*?)\\s*;")
GeneIDs <- GeneIDs[,2]
GeneIDs <- gsub(":","",GeneIDs)
GeneIDs <- GeneIDs[GeneIDs!=""]
# length = 1523631
# 33656 unique IDs
# Row 1523506 is Na

# isolate gene symbols
gene_sym <- stringr::str_match(gtf$V9, "gene_id\\s*(.*?)\\s*;")
gene_sym <- gene_sym[,2]
gene_sym <- gsub(":","",gene_sym)
#gene_sym <- gene_sym[gene_sym!=""]
# length = 1523631
# 33656 unique IDs
# Row 1523506 is blank ""


Protein_id <- stringr::str_match(gtf$V9, "protein_id\\s*(.*?)\\s*;")
Protein_id <- Protein_id[,2]
# length = 1523632
# 49508 unique IDs


gene_id_symbol <- data.frame(
  gene_id = GeneIDs, Symbol = gene_sym) |>
  group_by(gene_id) |>
  slice_head(n = 1) |>
  drop_na()

gene_symbol_protein <- data.frame(
  Symbol = gene_sym, protein_id = Protein_id) |>
  group_by(protein_id) |>
  slice_head(n = 1) |>
  group_by(Symbol) # unique proteins 49508 (extra proteins for isoforms)


gene_id_symbol_protein <- gene_id_symbol |>
  left_join(gene_symbol_protein, by = "Symbol")


# sum(!is.na(gene_id_symbol_protein$protein_id)) = 24659. So every symbol with 
# a protein ID from the gtf got one




#' Then we read the eggnogg mapper excel file resulting from querying it with
#' GCF_910596095.1/cds_from_genomic.fna and the default settings on
#' http://eggnog-mapper.embl.de
#' 
#' Change to eggnogg_id once the feature_counts is by gene_id
#' 
#' 

eggnogg <- 
  readxl::read_xlsx('data/michael_saad/MM_8hli8sc6.emapper.annotations.xlsx', 
                    skip = 2) 

eggnogg <- eggnogg |> 
  mutate(protein_id = 
           str_match(eggnogg$query, ".*_cds_(.*)(?=_)[^_]*")[,2])
# 48179 unique protein ids

eggnogg_no_protein_id <-  eggnogg[!complete.cases(eggnogg),]
genes_sym_no_protein <- unique(eggnogg_no_protein_id$Preferred_name)
# 11 Preferred names with no protein ID. I don't think we lose much here.

eggnogg_symbol <- eggnogg |>
  left_join(gene_id_symbol_protein, by = "protein_id") |>
  group_by(gene_id) |>
  slice_head(n = 1)

#' this is a little rough, taking only the first representative per gene when
#' there are multiple transcripts/proteins in the cds dataset. To me (Michael),
#' it seems like the highest scores are always first, so likely the most 
#' accurate mapping?

# Break the eggnogg results into long format
term_to_brite <- eggnogg_symbol |>
  select(BRITE, gene_id) |>
  rename(brite = BRITE) |>
  rowwise() |>
  mutate(brite = str_split(brite, ",")) |>
  unnest(cols = c(brite))


# get gene ontologies out of eggnogg.
term_to_gene <- eggnogg_symbol |>
  select(GOs, gene_id) |>
  filter(GOs != "-") |>
  separate_rows(GOs, sep = ",") |>
  rename(GO = GOs) 


# kegg categories from yanky br08902 file.
kegg_categories <- jsonlite::read_json("data/br08902.json") |>
  unlist() |>
  as_tibble() |>
  separate(value, into = c("brite", NA), remove = FALSE, extra = "drop") |>
  mutate(brite = str_replace(brite, "Protein", NA_character_)) |>
  rename(category = value)


