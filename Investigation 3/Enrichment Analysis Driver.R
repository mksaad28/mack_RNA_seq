# Example: For Humans ----
data(geneList, package = "DOSE")
gene <- names(geneList)[abs(geneList) > 2]
# Entrez gene ID
head(gene)

ego <- enrichGO(
    gene = gene,
    universe = names(geneList),
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
)
head(ego)

# For fish ----

# Filter out cases where gene_id == NA (Maybe try to fix this later)
pathway_dfs <- list(
  GO = df_filtered <- desqout_pro_pass_combo$df_res_GeneID_kegg[complete.cases(
    desqout_pro_pass_combo$df_res_GeneID_kegg[, 6]), ],
  KEGG = df_filtered_KEGG <- 
    desqout_pro_pass_combo$df_res_GeneID_kegg[complete.cases(
      desqout_pro_pass_combo$df_res_GeneID_kegg[, 7]), ] |> 
    filter(kegg_ko != "KNA")
)

my_list2 <- list()

for (name in names(pathway_dfs)) {
if (name == "GO") {
  pathway_criteria <- list(
    all = pathway_dfs[[name]],
    sig = pathway_dfs[[name]] |> 
      filter(abs(log2FoldChange) > 1),
    top500 = pathway_dfs[[name]] |> 
      arrange(padj) |> head(500),
    sig_up = pathway_dfs[[name]] |> 
      filter(log2FoldChange > 0),
    sig_down = pathway_dfs[[name]] |> 
      filter(log2FoldChange < 0)
  )

for (subset_name in names(pathway_criteria)) {
geneList <- pathway_criteria[[subset_name]][, 2]
names(geneList) <- as.character(pathway_criteria[[subset_name]][, 6])

geneList <- sort(geneList, decreasing = TRUE)
de <- names(geneList)[abs(geneList) > 2]

# Create term2gene

# Load sj_list
# sj_list <- loadSJList()
#
# term2gene <- merge(sj_list$entrez_genes, sj_list$go_all,
#                   by.x = "_id", by.y = "_id",
#                   all.x = TRUE)[, c("GOALL", "ENTREZID")]
#
# # Extract a named vector of all GO terms
# goterms <- Term(GOTERM)
#
# term2name <- data.frame("GOID" = names(goterms), "term" = goterms)
#
# ego <- enricher(gene          = de,
#                 pvalueCutoff  = 0.01,
#                 pAdjustMethod = "BH",
#                 universe      = names(geneList),
#                 qvalueCutoff  = 0.05,
#                 TERM2GENE = term2gene,
#                 TERM2NAME = term2name)

egoALL <- enrichGO(
  gene = de,
  OrgDb = org.Sjaponicus.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  universe = names(geneList),
  qvalueCutoff = 0.05,
  readable = TRUE
)

new_name <- paste("egoALL", name, subset_name, sep = "_")
assign(new_name, egoALL, envir = .GlobalEnv)
my_list2[[new_name]] <- get(new_name)
rm(egoALL)
rm(new_name)
}
}
}

egoCC <- enrichGO(
  gene = de,
  OrgDb = org.Sjaponicus.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  universe = names(geneList),
  qvalueCutoff = 0.05,
  readable = TRUE
)

egoBP <- enrichGO(
  gene = de,
  OrgDb = org.Sjaponicus.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  universe = names(geneList),
  qvalueCutoff = 0.05,
  readable = TRUE
)

egoMF <- enrichGO(
  gene = de,
  OrgDb = org.Sjaponicus.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  universe = names(geneList),
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(egoALL)

goplot(egoCC, showCategory = 20)
goplot(egoBP, showCategory = 20)
goplot(egoMF, showCategory = 20)


dotplot(egoALL,
        showCategory = 20,
        font.size = 8,
        title = "Dotplot for ORA"
)

egox <- setReadable(egoALL, org.Sjaponicus.eg.db, keyType = "ENTREZID")
cnetplot(egox,
         showCategory = 1,
         color.params = list(foldChange = geneList, edge = TRUE),
         categorySize = "pvalue",
         circular = TRUE,
         node_label = "category"
)

heatplot(egox,
         foldChange = geneList,
         showCategory = 5
)

edox2 <- pairwise_termsim(edox)
treeplot(edox2,
         cluster.params = list(method = "complete")
)

emapplot(
  edox2,
  layout.params = list(layout = "grid")
)

library(ggupset)
upsetplot(egoALL)

setwd(paste0(script_dir, "/Investigation 3/pathviews"))
pathview(
  gene.data = Sj_analysis$pathway_pro_pass_combo$KOKEGG$geneListUniverse, 
  pathway.id = "ko05200", 
  species = "ko", 
  limit = list(gene = max(abs(Sj_analysis$pathway_pro_pass_combo$KOKEGG$geneListUniverse)), cpd = 1)
)
