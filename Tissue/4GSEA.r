library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(stringr)
library(KEGGREST)
source("./PlotUtils.r")

# Prepare geneset
all.genes <- readRDS('./data/1AllGenes.rds')
gene.map <- bitr(rownames(all.genes),
		   fromType = "ENSEMBL",
           toType = c("ENTREZID"),
           OrgDb = "org.Mm.eg.db")
entrez <- unique(gene.map$ENTREZID)
genelist <- all.genes[gene.map[match(entrez, gene.map$ENTREZID), "ENSEMBL"], "logFC"]
names(genelist) <- entrez
genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[genelist != 0]

id.show <- c("mmu04080", "mmu04151", "mmu04020", "mmu04060", "mmu04022", "mmu04010", "mmu04724")
id.get <- keggGet(id.show)
enrich.db <- lapply(id.get, function(X) {
	gene <- X$GENE[seq(1, length(X$GENE), 2)] %>% as.numeric()
	name <- strsplit(X$NAME, " - ")[[1]][1]
	return(data.frame(gs_name = rep(name, length(gene)), entrez_gene = gene))
})
enrich.db <- do.call(rbind, enrich.db)

# Run GSEA
res <- GSEA(genelist, TERM2GENE = enrich.db)

library(ggplot2)
source("Ridge.r")

id.description <- c("Neuroactive ligand-receptor interaction",
					"Calcium signaling pathway",
					"cGMP-PKG signaling pathway",
					"Glutamatergic synapse")
p <- my.ridgeplot(res, id.description, id.description) +
	xlab("log2 Fold Change") + ylab(NULL) + labs(fill = "FDR") +
	scale_fill_gradient(low = "#FF165D", high = "#FF9A00") +
	theme_classic() +
	theme(axis.line = element_blank(),
		  axis.text.y = element_blank(),
		  panel.border = element_rect(linewidth = 1, fill = NA))
sized.ggsave(p, "./out/4GSEA.pdf", rel.size = c(0.33, 0), ratio = c(1, 1))

# export data
ID.fmt <- function (X) {
	str_sub(X, 6) %>%
		str_replace_all("_", " ") %>%
		str_to_lower() %>%
		str_to_sentence() %>%
		return()
}
export.res <- as.data.frame(res)
export.res$human.friendly <- sapply(export.res$ID, ID.fmt)
write.table(export.res, "./out/4GSEA.xls", quote = FALSE, sep = "\t")

all.genes <- readRDS("./data/1AllGenes.rds")
gene.map <- bitr(rownames(all.genes),
		   fromType = "ENSEMBL",
           toType = c("ENTREZID", "SYMBOL"),
           OrgDb = "org.Mm.eg.db")
all.genes <- all.genes[rownames(all.genes) %in% gene.map$ENSEMBL,]
all.genes$entrez <- gene.map[match(rownames(all.genes), gene.map$ENSEMBL), "ENTREZID"]
lapply(split(enrich.db, enrich.db$gs_name), function (X) {
	title <- X$gs_name[1]
	X <- X[X$entrez_gene %in% gene.map$ENTREZID, ]
	gene <- all.genes[match(X$entrez_gene, all.genes$entrez), c("logFC", "PValue", "FDR")]
	gene$symbol <- gene.map[match(rownames(gene), gene.map$ENSEMBL), "SYMBOL"]
	gene <- gene %>%
		filter(abs(logFC) > 1, PValue < 0.05) %>%
		arrange(desc(logFC))
	gene$symbol <- factor(gene$symbol, levels = gene$symbol)
	p <- ggplot(gene) +
		geom_col(aes(x = logFC, y = symbol, fill = -log10(PValue))) +
		ggtitle(title) +
		scale_fill_gradient(low = "#FF9A00", high = "#FF165D") +
		theme_classic() +
		theme(axis.line = element_blank(),
			  panel.border = element_rect(linewidth = 1, fill = NA))
	sized.ggsave(p, paste0("./out/4Expr-", title, ".pdf"), rel.size = c(0.5, 0), ratio = c(1, 2))
})

