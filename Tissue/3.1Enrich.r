# Analyse
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

diff <- readRDS("./data/2Diff.rds")

# Key mapping
genes <- bitr(rownames(diff),
		   fromType = "ENSEMBL",
           toType = c("ENTREZID", "SYMBOL"),
           OrgDb = "org.Mm.eg.db")

# GO
go <- enrichGO(gene = genes[, "SYMBOL"],
	OrgDb = "org.Mm.eg.db",
	keyType = "SYMBOL",
	ont = "ALL",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH") %>%
	as.data.frame()
write.table(go, "./out/3.1GO.xls", sep = "\t", quote = F)
saveRDS(go, "./data/3.1GO.rds")
 
# KEGG
kegg <- enrichKEGG(gene = genes[, "ENTREZID"],
	organism= "mmu",
	pAdjustMethod = "fdr",
	pvalueCutoff = 0.9,
	qvalueCutoff = 0.9) %>%
	as.data.frame()
write.table(kegg, "./out/3.2KEGG.xls", sep = "\t", quote = F)
saveRDS(kegg, "./data/3.2KEGG.rds")

