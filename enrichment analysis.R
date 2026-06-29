data <- fread('/path_to_data/data.csv')
# GO
go <- enrichGO(gene          = data$ENTREZID,
               OrgDb         = "org.Hs.eg.db",
               ont           = "ALL",#"BP", "MF", "CC"
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)

GO <- go@result 

# KEGG
R.utils::setOption("clusterProfiler.download.method","auto")

kegg <- enrichKEGG(
  gene = data$ENTREZID,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  organism = "hsa")  
KEGG <- kegg@result
