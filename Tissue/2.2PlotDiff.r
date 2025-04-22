library(ggplot2)
library(dplyr)

fpkm <- readRDS("./data/0FPKM.rds")
all.genes <- readRDS("./data/1AllGenes.rds")
source("./PlotUtils.r")

# PCA
m <- as.matrix(fpkm) %>%
	t() %>% scale(center = TRUE, scale = TRUE)
m <- m[, colSums(is.na(m)) == 0]
pca <- prcomp(m, scale = FALSE, center = FALSE)
pca.dat <- as.data.frame(pca$x[, 1:2])
pca.dat$Sample <- ifelse(rownames(pca.dat) %in% c("DSS1", "DSS2", "DSS3", "DSS4"), "DSS", "CMHs")
pca.dat$Lab <- rownames(pca.dat)
pca.summary <- summary(pca) 
p <- ggplot(pca.dat) +
	geom_label(aes(x = PC1, y = PC2 + 8, label = Sample), size = 3, family = "sans") +
	geom_point(aes(x = PC1, y = PC2, color = Sample, shape = Sample), size = 3) +
	geom_hline(aes(yintercept = 0), linetype = "dashed") +
	geom_vline(aes(xintercept = 0), linetype = "dashed") +
	xlab(paste0("PC1 (", round(pca.summary$importance["Proportion of Variance", "PC1"] * 100, 2), "%)")) +
	ylab(paste0("PC2 (", round(pca.summary$importance["Proportion of Variance", "PC2"] * 100, 2), "%)")) +
	xlim(-200, 100) +
	theme_classic() +
	theme(legend.position = "none") +
	theme(axis.text = element_text(size = 9, face = "bold"),
		  axis.title = element_text(size = 12, face = "bold"),
		  panel.border = element_rect(linewidth = 1, fill = NA))
sized.ggsave(p, "./out/2.2PCA.pdf", rel.size = c(0.25, 0))


# Volcano
diff.genes <- all.genes %>%
	mutate(status =
		ifelse(PValue < 0.05,
			ifelse(logFC < -0.5,
				"Down",
				ifelse(logFC > 0.5,
					"Up",
					"Stable")),
			"Stable")) %>%
	mutate(logpadj = -1 * log10(PValue))

p <- ggplot(diff.genes) +
	geom_point(aes(x = logFC, y = logpadj, color = status)) +
	geom_hline(aes(yintercept = -1 * log10(0.05)), linetype = "dashed") +
	geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
	geom_vline(aes(xintercept = -0.5), linetype = "dashed") +
	scale_color_manual(values = c(Stable = "gray", Up = "#FF165D", Down = "#3EC1D3")) +
	theme_classic() +
	labs(color = " ") +
	xlab("Log2 (Fold Change)") +
	ylab("-Log10 (p-value)") +
	theme(legend.position = "none") +
	theme(axis.text = element_text(size = 9, face = "bold"),
		  axis.title = element_text(size = 12, face = "bold"),
		  panel.border = element_rect(linewidth = 1, fill = NA))
table(diff.genes$status)
sized.ggsave(p, "./out/2.3Volcano.pdf", rel.size = c(0.25, 0))
