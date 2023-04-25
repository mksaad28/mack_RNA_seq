# Convert GeneSymbol to GeneID for ClusterProfiler viz

add_GeneID_to_df <- function(df_from_DESeq) {
  df_res_symbol <- df_from_DESeq |> 
    rownames_to_column(var = 'Symbol')
  
  df_res_symbol_GeneID <- df_res_symbol |> 
    left_join(gene_id_symbol, by = 'Symbol')
  
  df_res_GeneID <- df_res_symbol_GeneID |> select(-Symbol)
  
  df_res_GeneID_for_pathway <- df_res_GeneID |>
    drop_na(padj,gene_id) |>
    column_to_rownames(var = "gene_id")
  
  print(paste0("The number of rows using gene symbol is ",nrow(df_res_symbol)))
  print(paste0("The number of rows using gene ID is ",nrow(df_res_GeneID)))
  
  output_df <- list(df_res_symbol = df_res_symbol, 
                    df_res_symbol_GeneID = df_res_symbol_GeneID, 
                    df_res_GeneID = df_res_GeneID,
                    df_res_GeneID_for_pathway = 
                    df_res_GeneID_for_pathway)
  return(output_df)
}