# Be sure to set the right wd
source("api.R")
library(limma)
library(edgeR)
library(tidyverse)
library(clusterProfiler)
library(readxl)

##### ------------ eggnogg ------

## Necessary argument from orthology because neither tuna or mackerel genomes
## are in handy databases like ensembl, gene ontology, or kegg. First we read
## the cds from GCF_910596095.1 and extract the cds accession and the gene id.
## This is necessary because the cds from GCF_910596095.1 were fed into
## eggnogg-mapper. This file was downloaded from NCBI genomes.
a <- readLines("data/michael_saad/GCF_910596095.1/cds_from_genomic.accessions_only.fna")

query_to_gene <- tibble()
for (line in a) {
  if (!str_starts(line, ">")) next
  line_split <- str_split(line, " ")[[1]]
  accession <- str_remove(line_split[1], ">")
  # print(accession)
  # print(line_split[2])
  if (!str_detect(line_split[2], "gene")) stop("Missing gene in second field")
  gene <- line_split[2] |>
    str_remove(fixed("[gene=")) |>
    str_remove(fixed("]"))
  query_to_gene <- bind_rows(query_to_gene, bind_cols(query = accession, gene_id = gene))
  # print("--")
}



## Then we read the eggnogg mapper excel file resulting from querying it with
## GCF_910596095.1/cds_from_genomic.fna and the default settings on
## http://eggnog-mapper.embl.de
eggnogg <- readxl::read_xlsx("data/michael_saad/MM_8hli8sc6.emapper.annotations.xlsx", skip = 2) |>
  left_join(query_to_gene, by = "query") |>
  select(gene_id, everything()) |>
  group_by(gene_id) |>
  slice_head(n = 1)
## this is a little rough, taking only the first representative per gene when
## there are multiple transcripts/proteins in the cds dataset.


## Break the eggnogg results into long format
genes_and_brite <- eggnogg |>
  select(gene_id, BRITE) |>
  rename(brite = BRITE) |>
  rowwise() |>
  mutate(brite = str_split(brite, ",")) |>
  unnest(cols = c(brite))

term_to_brite <- genes_and_brite[c('brite','gene_id')]

## get gene ontologies out of eggnogg.
term_to_gene <- eggnogg |>
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




##### ------------ raw data ------
gene_lookup <- read_tsv("data/michael_saad/gene_result.txt",
  col_types = cols(GeneID = col_character()),
  guess_max = 10^10
) |>
  select(Symbol, description) |>
  unique() |>
  group_by(`Symbol`) |>
  summarize(
    description = paste(description, collapse = ";"),
    .groups = "drop"
  ) |> ## Assure a 1:1 mapping even when there are multiple entries for 1 ensembl ID
  rename(gene_id = Symbol)

gene_lookup |>
  group_by(gene_id) |>
  summarize(n = n()) |>
  filter(n != 1)

expr_long <- read_tsv("data/michael_saad/featurecounts_results.txt", skip = 1) |>
  select(-Chr, -Start, -End, -Strand, -Length) |>
  rename(gene_id = Geneid) |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") |>
  mutate(sample = str_remove(sample, "STAR2/output/") |> str_remove("Aligned.sortedByCoord.out.bam"))

sample_lookup <- expr_long |>
  select(sample) |>
  unique() |>
  separate(sample, into = c("passage", "protocol_replicate"), remove = FALSE) |>
  mutate(
    protocol = str_sub(protocol_replicate, 1, 1),
    replicate = str_sub(protocol_replicate, 2, 2)
  ) |>
  select(-protocol_replicate)

expr_long |>
  group_by(gene_id) |>
  summarize(prop_count0 = mean(count == 0)) |>
  ggplot() +
  aes(x = prop_count0) +
  geom_histogram(binwidth = 0.1)



expr_clean <- expr_long |>
  group_by(gene_id) |>
  mutate(gene_prop_count0 = mean(count == 0)) |>
  filter(gene_prop_count0 != 1) |>
  select(-gene_prop_count0) |>
  ungroup() ## ^^ Get rid of the genes which are always 0.


expr_wide <- expr_clean |>
  pivot_wider(names_from = sample, values_from = count)


expr_wide_cpm <- expr_wide |>
  select(-gene_id) |>
  cpm() |>
  as_tibble() |>
  mutate(gene_id = expr_wide$gene_id) # counts per million

expr_long_cpm <- expr_wide_cpm |>
  pivot_longer(-gene_id,
    names_to = "sample",
    values_to = "cpm"
  )



expr_withmeta <- expr_clean |>
  left_join(sample_lookup, by = "sample")

# Run CalcNormFactors and compare cpm to that without normalization --> conclusion, normalization doesn't change much

# first create DGE object
expr_wide_DGE <- expr_wide |>
    select(-gene_id) |> DGEList()

# Use without interaction design for removing zero or low counts
expr_design <- model.matrix(~ passage + protocol,
                            data = sample_lookup
)

keep <- filterByExpr(expr_wide |> select(-gene_id), expr_design)


expr_wide_DGE <- expr_wide_DGE[keep,keep.lib.sizes=FALSE]


# Then, run normalization
expr_wide_DGE_norm <- calcNormFactors(expr_wide_DGE)

expr_wide_norm_cpm <- expr_wide_DGE_norm$counts |>
    cpm() |>
    as_tibble() |>
    mutate(gene_id = expr_wide[keep,keep.lib.sizes=FALSE]$gene_id) # counts per million

expr_long_norm_cpm <- expr_wide_norm_cpm |>
    pivot_longer(-gene_id,
                 names_to = "sample",
                 values_to = "cpm"
    )

expr_long_comp_norm <- data.frame(
    wo_norm = expr_long_cpm$cpm,
    with_norm = expr_long_norm_cpm$cpm
)


linear_model_comp_norm <- lm(wo_norm ~ with_norm, data = expr_long_comp_norm)
summary(linear_model_comp_norm)

plot( expr_long_comp_norm$wo_norm, expr_long_comp_norm$with_norm )
abline( linear_model_comp_norm)

# Save normalized dataframe
write.csv(expr_wide_norm_cpm, "output/dataframes/counts_norm_cpm.csv")


# PCA ---------------------------------------------------------------------


pca <- prcomp(expr_wide |> select(-gene_id) |> t() |> cpm())

pca_expr <- predict(pca) |>
  as_tibble(rownames = "sample") |>
  select(sample, PC1, PC2) |>
  left_join(sample_lookup)

pca_expr |>
  ggplot() +
  aes(x = PC1, y = PC2, col = passage, shape = protocol) +
  geom_point()



# limma without interaction -------
do_without_interaction <- function(filename) {
  expr_design <- model.matrix(~ passage + protocol,
    data = sample_lookup
  )

  keep <- filterByExpr(expr_wide |> select(-gene_id), expr_design) 
  # filterByExpr: Determine which genes have sufficiently large counts to be retained in a statistical analysis.


  expr_voomed <- voom(
    expr_wide[keep, sample_lookup$sample],
    expr_design,
    plot = TRUE
  ) # Voom: Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observational-level weights. The data are then ready for linear modelling.


  expr_fit <- lmFit(expr_voomed, expr_design)
  # lmFit: Fit linear model for each gene given a series of arrays


  expr_fit_ebayes <- eBayes(expr_fit)
  # eBayes: Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.

  expr_limma_res <- topTable(expr_fit_ebayes,
    coef = "passagep77",
    genelist = expr_wide |> filter(keep) |> select(gene_id),
    number = 20000000,
    confint = TRUE
  ) |>
    as_tibble() |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(hit_num = row_number())
# topTable: Extract a table of the top-ranked genes from a linear model fit.

  # Gene ontologies ---------------

  enricher_obj <- enricher(
    gene = expr_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = expr_limma_res |> pull(gene_id),
    TERM2GENE = term_to_gene
  )
  # enricher: A universal enrichment analyzer
  # genes are signficant genes, universe is all genes. Later used for Fisher's exact test?

  enricher_res <- enricher_obj@result |> as_tibble()

  go_lookup <- go2term(enricher_res$ID) |>
    as_tibble() |>
    left_join(go2ont(enricher_res$ID), by = "go_id") |>
    mutate(Ontology = case_when(
      Ontology == "BP" ~ "Biological Process",
      Ontology == "MF" ~ "Molecular Function",
      Ontology == "CC" ~ "Cellular Component"
    ))

  enricher_res_tbl <- enricher_res |>
    left_join(go_lookup, by = c("ID" = "go_id")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(Term, Ontology, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())

  ## KEGG
  kegg_obj <- enricher(
    gene = expr_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = expr_limma_res |> pull(gene_id),
    TERM2GENE = term_to_brite
  )
  # enricher: A universal enrichment analyzer
  # genes are signficant genes, universe is all genes. Later used for Fisher's exact test?
  
  kegg_res_tbl <- kegg_obj@result |> 
    as_tibble() |>
    inner_join(kegg_categories, by = c("ID" = "brite")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(category, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  write.csv(expr_limma_res, paste("output/dataframes/df",filename,"genes.csv",sep = "_"))
  write.csv(enricher_res_tbl, paste("output/dataframes/df",filename,"pathways.csv",sep = "_"))
  write.csv(kegg_res_tbl, paste("output/dataframes/df",filename,"kegg.csv",sep = "_"))

  
}


without_interaction <- do_without_interaction(filename = "wo_int")



# limma with interactions ----------

#' do an interaction analysis with a specified reference level t
#'
#' @param reference_protocol
#'
#' @return a limma toptable
#' @export
#'
#'
#'
#'
#' @examples
do_interaction_analysis <- function(reference_protocol, filename) {
  int_design <- model.matrix(~ passage * protocol,
    data = sample_lookup |>
      mutate(protocol = relevel(factor(protocol), reference_protocol))
  )

  int_keep <- filterByExpr(expr_wide |> select(-gene_id), int_design)


  int_voomed <- voom(
    expr_wide[int_keep, sample_lookup$sample],
    int_design,
    plot = TRUE
  )

  int_fit <- lmFit(int_voomed, int_design)


  int_fit_ebayes <- eBayes(int_fit)

  int_limma_res <- topTable(int_fit_ebayes,
    coef = "passagep77",
    genelist = expr_wide |> filter(int_keep) |> select(gene_id),
    number = 20000000,
    confint = TRUE
  ) |>
    as_tibble() |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(hit_num = row_number())


  # Gene ontologies ---------------

  int_enricher_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_gene
  )

  int_enricher_res <- int_enricher_obj@result |> as_tibble()

  go_lookup <- go2term(int_enricher_res$ID) |>
    as_tibble() |>
    left_join(go2ont(int_enricher_res$ID), by = "go_id") |>
    mutate(Ontology = case_when(
      Ontology == "BP" ~ "Biological Process",
      Ontology == "MF" ~ "Molecular Function",
      Ontology == "CC" ~ "Cellular Component"
    ))

  int_enricher_res_tbl <- int_enricher_res |>
    left_join(go_lookup, by = c("ID" = "go_id")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(Term, Ontology, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())


  data_long <- expr_long_cpm |>
    left_join(sample_lookup, by = "sample") |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(gene_name = gene_id)

  ## KEGG
  kegg_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_brite
  )
  # enricher: A universal enrichment analyzer
  # genes are signficant genes, universe is all genes. Later used for Fisher's exact test?
  
  kegg_res_tbl <- kegg_obj@result |> 
    as_tibble() |>
    inner_join(kegg_categories, by = c("ID" = "brite")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(category, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  write.csv(int_limma_res, paste("output/dataframes/df",filename,"genes.csv",sep = "_"))
  write.csv(int_enricher_res_tbl, paste("output/dataframes/df",filename,"pathways.csv",sep = "_"))
  write.csv(kegg_res_tbl, paste("output/dataframes/df",filename,"kegg.csv",sep = "_"))
  
}


# Call, run, and write limma with interactions ----------------------
interaction_f <- do_interaction_analysis(
  reference_protocol = "f",
  filename = "f"
)

interaction_m <- do_interaction_analysis(
  reference_protocol = "m",
  filename = "m"
)

interaction_p <- do_interaction_analysis(
  reference_protocol = "p",
  filename = "p"
)

interaction_s <- do_interaction_analysis(
  reference_protocol = "s",
  filename = "s"
)




# ----- D
do_medium_analysis <- function(filename, reference_passage) {
  int_design <- model.matrix(~ passage * protocol,
    data = sample_lookup |>
      mutate(
        protocol = relevel(factor(protocol), "p"),
        passage = relevel(factor(passage), "p77")
      )
  )

  int_keep <- filterByExpr(expr_wide |> select(-gene_id), int_design)


  int_voomed <- voom(
    expr_wide[int_keep, sample_lookup$sample],
    int_design,
    plot = TRUE
  )

  int_fit <- lmFit(int_voomed, int_design)


  int_fit_ebayes <- eBayes(int_fit)

  int_limma_res <- topTable(int_fit_ebayes,
    coef = "protocols",
    genelist = expr_wide |> filter(int_keep) |> select(gene_id),
    number = 20000000,
    confint = TRUE
  ) |>
    as_tibble() |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(hit_num = row_number())


  # Gene ontologies ---------------

  int_enricher_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_gene
  )

  int_enricher_res <- int_enricher_obj@result |> as_tibble()

  go_lookup <- go2term(int_enricher_res$ID) |>
    as_tibble() |>
    left_join(go2ont(int_enricher_res$ID), by = "go_id") |>
    mutate(Ontology = case_when(
      Ontology == "BP" ~ "Biological Process",
      Ontology == "MF" ~ "Molecular Function",
      Ontology == "CC" ~ "Cellular Component"
    ))

  int_enricher_res_tbl <- int_enricher_res |>
    left_join(go_lookup, by = c("ID" = "go_id")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(Term, Ontology, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())


  data_long <- expr_long_cpm |>
    left_join(sample_lookup, by = "sample") |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(gene_name = gene_id)

  
  ## KEGG
  kegg_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_brite
  )
  # enricher: A universal enrichment analyzer
  # genes are signficant genes, universe is all genes. Later used for Fisher's exact test?
  
  kegg_res_tbl <- kegg_obj@result |> 
    as_tibble() |>
    inner_join(kegg_categories, by = c("ID" = "brite")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(category, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  write.csv(int_limma_res, paste("output/dataframes/df",filename,"genes.csv",sep = "_"))
  write.csv(int_enricher_res_tbl, paste("output/dataframes/df",filename,"pathways.csv",sep = "_"))
  write.csv(kegg_res_tbl, paste("output/dataframes/df",filename,"kegg.csv",sep = "_"))

}

medium_p27 <- do_medium_analysis(filename = "media_p27", reference_passage = "p27")
medium_p77 <- do_medium_analysis(filename = "media_p77", reference_passage = "p77")


## Visualizations

