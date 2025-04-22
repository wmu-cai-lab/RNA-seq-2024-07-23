library(ggplot2)
library(ggtree)
library(dplyr)
library(reshape2)
library(patchwork)
library(KEGGREST)
library(org.Mm.eg.db)

fpkm <- readRDS("./data/0FPKM.rds")
source("./PlotUtils.r")

# heatmap
all.genes <- readRDS("./data/1AllGenes.rds")
diff <- all.genes %>% filter(abs(logFC) > 0.5 & PValue < 0.05)
saveRDS(diff, "./data/2Diff.rds")
fpkm <- fpkm[rownames(diff), ]

# custom heatmap
# perform row-wise normalization
fpkm.scaled <- t(apply(fpkm, 1, function (X) {
	return(scale(X))
}))
colnames(fpkm.scaled) <- colnames(fpkm)

# cluster tree panel
cluster.row <- hclust(dist(fpkm.scaled))
cluster.col <- hclust(dist(t(fpkm.scaled)))

# main panel
long.table <- as.data.frame(fpkm.scaled) %>%
	mutate(gene.name = rownames(fpkm.scaled)) %>%
	mutate(gene.name = factor(gene.name, levels = cluster.row$labels[cluster.row$order])) %>%
	reshape2::melt(id.vars = c("gene.name")) %>%
	mutate(variable = factor(variable, levels = cluster.col$labels[cluster.col$order]))
long.table$design <- sapply(long.table$variable, function(X){
	ifelse(X %in% c("DSS1", "DSS2", "DSS3", "DSS4"), "DSS", "CHMs") %>% return()
})

# group annotation panel
design <- c("DSS", "DSS", "DSS", "CMHs", "CMHs", "CMHs")
design <- design[cluster.col$order]
gen.group.panel.df <- function(design) {
	last.design <- NA
	design.width <- 0
	res <- list()
	for(i in 1:length(design)) {
		if (is.na(last.design) | design[i] == last.design) {
			design.width <- design.width + 1
		} else {
			res <- append(res, list(data.frame(x = i - design.width, w = design.width, ann = last.design)))
			design.width <- 1
		}
		last.design <- design[i]
	}
	res <- append(res, list(data.frame(x = length(design) + 1 - design.width, w = design.width, ann = last.design)))
	return(do.call(rbind, res))
}
group.df <- gen.group.panel.df(design)

# patch panels
tree.p <- ggtree(cluster.row, layout = "rectangular", branch.length = "none", size = 0.2) +
	theme(plot.margin = margin())

main.p <- ggplot(long.table) +
	theme_classic() +
	scale_fill_gradient2(low = "#3EC1D3", mid = "#EEEEEE", high = "#FF165D") +
	scale_x_discrete(expand = c(0, 0)) +
	scale_y_discrete(expand = c(0, 0), position = "right") +
	geom_raster(aes(x = variable, y = gene.name, fill = value)) +
	xlab(NULL) + ylab(NULL) +
	theme(#axis.text.y = element_blank(),
		  axis.text.y = element_blank(),
		  axis.text.x = element_blank(),
		  axis.title = element_text(size = 12, face = "bold"),
		  axis.ticks = element_blank(),
		  axis.line = element_blank(),
		  legend.title = element_blank(),
		  panel.border = element_blank(),
		  legend.text = element_text(size = 9, face = "bold"),
		  plot.margin = margin())

group.p <- ggplot(group.df) +
	geom_tile(aes(x = x, y = 1, width = w, height = 1, fill = ann)) +
	geom_text(aes(x = x, y = 1, label = ann)) +
	scale_x_discrete(expand = c(0, 0)) +
	scale_y_continuous(expand = c(0, 0)) +
	scale_fill_manual(values = c("#3EC1D3", "#FF165D")) +
	theme_classic() +
	theme(axis.text = element_blank(),
		  axis.title = element_blank(),
		  axis.ticks = element_blank(),
		  axis.line = element_blank(),
		  legend.position = "none",
		  plot.margin = margin(),
		  panel.border = element_blank())

layout <- c(
	area(3, 1, 40, 2),
	area(1, 3, 2, 20),
	area(3, 3, 40, 20)
)

res.p <- tree.p + group.p + main.p +
	plot_layout(design = layout, guides = 'collect')
sized.ggsave(res.p, "./out/2.1Heatmap.pdf", rel.size = c(0.33, 0), ratio = c(1, 1))

