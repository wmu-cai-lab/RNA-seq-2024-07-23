# visualize
library(ggplot2)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
source("./PlotUtils.r")

all.genes <- readRDS("./data/1AllGenes.rds")
gene.map <- bitr(rownames(all.genes),
		   fromType = "ENSEMBL",
           toType = c("ENTREZID", "SYMBOL"),
           OrgDb = "org.Mm.eg.db")

## GO
go.immune <- str_to_title(c("cell junction assembly", "leukocyte migration"))
go.nerve <- c("regulation of membrane potential", "regulation of synapse structure or activity", "regulation of synapse organization", "postsynaptic membrane", "neuron to neuron synapse", "postsynaptic specialization", "asymmetric synapse", "distal axon", "presynaptic membrane", "gated channel activity")
go.all <- readRDS("./data/3.1GO.rds")
go <- go.all %>%
	filter(pvalue < 0.05) %>%
	arrange(desc(Count)) %>%
	head(30) %>%
	arrange(ONTOLOGY)
go$Description <- str_to_title(go$Description)
go$Description <- factor(go$Description, levels = rev(go$Description))
go$LabY <- rep(1, nrow(go))

p <- ggplot(go) +
	geom_col(aes(x = Description, y = Count, fill = ONTOLOGY)) +
	geom_text(aes(x = Description, y = LabY, label = Description), size = 2.5, fontface = "bold", hjust = "left") +
	xlab("Term") +
	ylab("Gene Count") +
	theme_classic() +
	coord_flip() +
	labs(fill = " ") +
	scale_fill_manual(values = c("#FF165D", "#FF9A00", "#3EC1D3")) +
	theme(axis.text.y = element_blank(),
		  axis.text.x = element_text(size = 9, face = "bold"),
		  axis.title = element_text(size = 12, face = "bold"),
		  panel.border = element_rect(linewidth = 1, fill = NA))
sized.ggsave(p, "./out/3.1GO.pdf", rel.size = c(0.33, 0), ratio = c(1.5, 2))

for (idx in 1:length(go.immune)) {
	i <- go.immune[idx]
	genelist <- strsplit(go[match(i, go$Description), ]$geneID, "/")[[1]]
	genelist <- genelist[genelist %in% gene.map$SYMBOL]
	df <- data.frame(symbol = genelist,
					 logfc = all.genes[gene.map[match(genelist, gene.map$SYMBOL), "ENSEMBL"], "logFC"],
					 pval = all.genes[gene.map[match(genelist, gene.map$SYMBOL), "ENSEMBL"], "PValue"]) %>%
		arrange(desc(abs(logfc))) %>%
		head(30) %>%
		arrange(desc(logfc))
	df$symbol <- factor(df$symbol, levels = df$symbol)
	p <- ggplot(df) +
		geom_col(aes(x = symbol, y = logfc, fill = -log10(pval))) +
		coord_flip() +
		ggtitle(str_to_title(i)) +
		theme_bw() +
		xlab(NULL) +
		ylab("Log FC") +
		theme(legend.position = "none", text = element_text(size = 16, family = "sans")) +
		scale_fill_gradient(low = "#FF9A00", high = "#FF165D") +
		theme(plot.title = element_text(size = 16, face = "bold"),
			  axis.text = element_text(size = 9, face = "bold"),
			  axis.title = element_text(size = 12, face = "bold"),
			  panel.border = element_rect(linewidth = 1))
	sized.ggsave(p, paste0("./out/3.", idx + 1, "GO_", i, ".pdf"), rel.size = c(0.33, 0), ratio = c(1, 1.5))
}

## KEGG
calcGeneRatio <- function (x) {
	args <- strsplit(x, "/")[[1]] %>% as.numeric()
	return(args[1] / args[2])
}
stripPathway <- function (X) {
	return(str_trim(strsplit(X, " - ")[[1]][1]))
}

kegg.all <- readRDS("./data/3.2KEGG.rds")
kegg <- kegg.all %>%
	filter(pvalue < 0.05) %>%
	mutate(fGeneRatio = sapply(GeneRatio, calcGeneRatio)) %>%
	arrange(desc(Count)) %>%
	mutate(sPathway = sapply(Description, stripPathway)) %>%
	head(20) %>%
	mutate(ePathway = factor(sPathway, levels = rev(sPathway)))

p <- ggplot(kegg) +
	theme_classic() +
	geom_point(aes(x = fGeneRatio, y = ePathway, size = Count, color = p.adjust)) +
	geom_text(aes(x = fGeneRatio, y = ePathway, label = ePathway), size = 2, fontface = "bold", hjust = ifelse(kegg$fGeneRatio > 0.07, 1, 0), nudge_x = ifelse(kegg$fGeneRatio > 0.07, -0.003, 0.003)) +
	xlab("Gene Ratio") +
	ylab(NULL) +
	labs(color = "FDR") +
	#theme(text = element_text(size = 16, family = "sans")) +
	scale_color_gradient(low = "#FF165D", high = "#FF9A00") +
	scale_y_discrete(expand = c(0, 0.8)) +
	theme(axis.text.y = element_blank(),
		  axis.text.x = element_text(size = 9, face = "bold"),
		  axis.title = element_text(size = 12, face = "bold"),
		  axis.ticks.y = element_blank(),
		  axis.line = element_blank(),
		  panel.border = element_rect(linewidth = 1, fill = NA),
		  legend.title = element_text(size = 9, face = "bold"),
		  legend.text = element_text(size = 9, face = "bold"))
sized.ggsave(p, "./out/7KEGG.pdf", rel.size = c(0.33, 0), ratio = c(1, 1))

