library(limma)
library(dplyr)
library(edgeR)

count.data <- readRDS("./data/0Counts.rds") %>%
	as.matrix() %>%
	avereps() %>%
	normalizeBetweenArrays()
count.data <- count.data[rowMeans(count.data) > 0, ]

group <- c(rep("DSS", 3), rep("CMHs", 3))
design <- model.matrix(~ group)

diff <- DGEList(counts = count.data,group = group) %>%
	calcNormFactors() %>%
	estimateCommonDisp() %>%
	estimateTagwiseDisp() %>%
	exactTest(pair = c("DSS", "CMHs")) %>%
	topTags(n = Inf)

diff %>% as.data.frame() %>% saveRDS("./data/1AllGenes.rds")

