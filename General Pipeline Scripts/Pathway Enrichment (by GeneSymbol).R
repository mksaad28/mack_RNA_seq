# run_pathway function: performs pathway enrichment analysis with enricher --
# using over representation analysis

# INPUTS
# df_DESeq_res: the data frame of the differential expression results
# padj_cutoff: numeric value below which a gene is counted as significant

# OUTPUTS: in list output_run_pathway
# go_res_tbl: table with pathway enrichment based on gene ontology
# kegg_res_tbl: table with pathway enrichment based on KEGG brite

run_pathway <- function(df_DESeq_res, padj_cutoff) {
  
  # enricher: A universal enrichment analyzer
  # genes are significant genes, universe is all genes. 
  # Later used for over representation analysis (ORA)
  enricher_obj <- enricher(
    gene = df_DESeq_res |> filter(padj < padj_cutoff) |> rownames(),
    universe = df_DESeq_res |> rownames(),
    TERM2GENE = term_to_gene #takes gene symbol now; convert to gene_id later
  )
  
  # Sample code to convert to geneID (loses ~350 genes)
  # df_res_geneID <- output_passage_77_27$df_res |> 
  #   rownames() |> 
  #   as.data.frame() |> 
  #   setNames("Symbol") |> 
  #   left_join(gene_id_symbol, by = "Symbol") |> 
  #   as.data.frame()
  
  enricher_res <- enricher_obj@result |> as_tibble()
  
  go_lookup <- go2term(enricher_res$ID) |>
    as_tibble() |>
    left_join(go2ont(enricher_res$ID), by = "go_id") |>
    mutate(Ontology = case_when(
      Ontology == "BP" ~ "Biological Process",
      Ontology == "MF" ~ "Molecular Function",
      Ontology == "CC" ~ "Cellular Component"
    ))
  
  go_res_tbl <- enricher_res |>
    left_join(go_lookup, by = c("ID" = "go_id")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / 
             eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(Term, Ontology, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  
  ## KEGG
  kegg_obj <- enricher(
    gene = df_DESeq_res |> filter(padj < padj_cutoff) |> rownames(),
    universe = df_DESeq_res |> rownames(),
    TERM2GENE = term_to_brite #takes gene symbol now; convert to gene_id later
  )
  
  kegg_res_tbl <- kegg_obj@result |> 
    as_tibble() |>
    inner_join(kegg_categories, by = c("ID" = "brite")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / 
             eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(category, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  output_run_pathway <- 
    list(go_res_tbl = go_res_tbl, kegg_res_tbl = kegg_res_tbl)
  
  return(output_run_pathway)
  
}