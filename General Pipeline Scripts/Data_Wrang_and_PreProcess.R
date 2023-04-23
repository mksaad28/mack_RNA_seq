## Data wrangling and processing of gene expression data

## Add a header to explain this file


#' # Read in the gene expression data from a TSV file, skip the first row, 
#' select all columns except Chr, Start, End, Strand, and Length, rename the 
#' Geneid column to gene_id, pivot the data from wide to long format, and 
#' remove the "STAR2/output/" and "Aligned.sortedByCoord.out.bam" strings 
#' from the sample column.
expr_long <- read_tsv("data/michael_saad/featurecounts_results.txt", 
                      skip = 1) |>
  select(-Chr, -Start, -End, -Strand, -Length) |>
  rename(gene_symbol = Geneid) |>
  pivot_longer(-gene_symbol, names_to = "sample", values_to = "count") |>
  mutate(sample = str_remove(sample, "STAR2/output/") |> 
           str_remove("Aligned.sortedByCoord.out.bam"))

# Create a sample_lookup table from the sample names in the expr_long and
# breaks in down into passage, protocol, and replicate
sample_lookup <- expr_long |>
  select(sample) |>
  unique() |>
  separate(sample, into = c("passage", "protocol_replicate"), remove = FALSE) |>
  mutate(
    protocol = str_sub(protocol_replicate, 1, 1),
    replicate = str_sub(protocol_replicate, 2, 2)
  ) |>
  select(-protocol_replicate)

#' Group the gene expression data by gene symbol, calculate the proportion of 
#' samples with count 0 for each gene, and create a histogram 
#' of the proportion values.
expr_long |>
  group_by(gene_symbol) |>
  summarize(prop_count0 = mean(count == 0)) |>
  ggplot() +
  aes(x = prop_count0) +
  geom_histogram(binwidth = 0.1)

#' Group the gene expression data by gene symbol, add a column for the 
#' proportion of samples with count 0 for each gene, filter out genes where the 
#' proportion is always 1, remove the additional column, and ungroup the data.
expr_clean <- expr_long |>
  group_by(gene_symbol) |>
  mutate(gene_prop_count0 = mean(count == 0)) |>
  filter(gene_prop_count0 != 1) |>
  select(-gene_prop_count0) |>
  ungroup() ## ^^ Get rid of the genes which are always 0.

# Pivot the gene expression data from long to wide format, with one column 
# for each sample.
expr_wide <- expr_clean |>
  pivot_wider(names_from = sample, values_from = count) |>
  rename(Symbol = gene_symbol)

# Make rowname of expr_wide_prep the gene symbol (needed for DESeq2)
expr_wide_prep <- expr_wide |>
  column_to_rownames(var = "Symbol")

# Create metadata table for protocol, passage, or a combination of the two
meta <- sample_lookup |> 
  column_to_rownames(var = "sample") |>
  select(-replicate) |>
  mutate(passage_protocol = paste(passage,protocol, sep = '_'))