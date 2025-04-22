exp <- read.table("Result/4_expression/4.1_expression/Expression_with_annotation.xls", sep = "\t", header = TRUE, quote = "")
rownames(exp) <- exp[, 1]

counts <- exp[, c(2, 4, 6, 8, 10, 12, 14, 16)]
colnames(counts) <- c("DSS1", "DSS2", "DSS3", "DSS4", "CMHs1", "CMHs2", "CMHs3", "CMHs4")
fpkm <- exp[, c(2, 4, 6, 8, 10, 12, 14, 16) + 1]
colnames(fpkm) <- c("DSS1", "DSS2", "DSS3", "DSS4", "CMHs1", "CMHs2", "CMHs3", "CMHs4")

saveRDS(counts[, c("DSS1", "DSS3", "DSS4", "CMHs2", "CMHs3", "CMHs4")], "./data/0Counts.rds")
saveRDS(fpkm[, c("DSS1", "DSS3", "DSS4", "CMHs2", "CMHs3", "CMHs4")], "./data/0FPKM.rds")

