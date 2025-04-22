library(dplyr)
library(patchwork)
library(ggtree)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

# Neuroactive ligand-receptor interaction
interested <- c("1810009J06Rik", "Adcyap1", "Adcyap1r1", "Adm", "Adm2", "Adora1", "Adora2a", "Adora2b", "Adora3", "Adra1a", "Adra1b", "Adra1d", "Adra2a", "Adra2b", "Adra2c", "Adrb1", "Adrb2", "Adrb3", "Agt", "Agtr1a", "Agtr1b", "Agtr2", "Apela", "Apln", "Aplnr", "Avp", "Avpr1a", "Avpr1b", "Avpr2", "Bdkrb1", "Bdkrb2", "Brs3", "C3", "C3ar1", "C5ar1", "Calca", "Calcb", "Calcr", "Calcrl", "Cck", "Cckar", "Cckbr", "Cga", "Chrm1", "Chrm2", "Chrm3", "Chrm4", "Chrm5", "Chrna1", "Chrna10", "Chrna2", "Chrna3", "Chrna4", "Chrna5", "Chrna6", "Chrna7", "Chrna9", "Chrnb1", "Chrnb2", "Chrnb3", "Chrnb4", "Chrnd", "Chrne", "Chrng", "Cnr1", "Cnr2", "Cort", "Crh", "Crhr1", "Crhr2", "Ctsg", "Cysltr1", "Cysltr2", "Drd1", "Drd2", "Drd3", "Drd4", "Drd5", "Edn1", "Edn2", "Edn3", "Ednra", "Ednrb", "F2", "F2r", "F2rl1", "F2rl2", "F2rl3", "Fpr-rs3", "Fpr-rs4", "Fpr-rs6", "Fpr-rs7", "Fpr1", "Fpr2", "Fpr3", "Fshb", "Fshr", "Gabbr1", "Gabbr2", "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabra5", "Gabra6", "Gabrb1", "Gabrb2", "Gabrb3", "Gabrd", "Gabre", "Gabrg1", "Gabrg2", "Gabrg3", "Gabrp", "Gabrq", "Gabrr1", "Gabrr2", "Gabrr3", "Gal", "Galp", "Galr1", "Galr2", "Galr3", "Gcg", "Gcgr", "Gh", "Ghr", "Ghrh", "Ghrhr", "Ghrl", "Ghsr", "Gip", "Gipr", "Glp1r", "Glp2r", "Glra1", "Glra2", "Glra3", "Glrb", "Gm2663", "Gm3867", "Gnrh1", "Gnrhr", "Gpha2", "Gphb5", "Gpr156", "Gpr35", "Gpr50", "Gpr83", "Gria1", "Gria2", "Gria3", "Gria4", "Grid1", "Grid2", "Grik1", "Grik2", "Grik3", "Grik4", "Grik5", "Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b", "Grm1", "Grm2", "Grm3", "Grm4", "Grm5", "Grm6", "Grm7", "Grm8", "Grp", "Grpr", "Gzma", "Hc", "Hcrt", "Hcrtr1", "Hcrtr2", "Hrh1", "Hrh2", "Hrh3", "Hrh4", "Htr1a", "Htr1b", "Htr1d", "Htr1f", "Htr2a", "Htr2b", "Htr2c", "Htr4", "Htr5a", "Htr5b", "Htr6", "Htr7", "Iapp", "Insl3", "Insl5", "Kiss1", "Kiss1r", "Kng1", "Kng2", "Lep", "Lepr", "Lhb", "Lhcgr", "Lpar1", "Lpar2", "Lpar3", "Lpar4", "Lpar6", "Ltb4r1", "Ltb4r2", "Lynx1", "Lypd6", "Lypd6b", "Mas1", "Mc1r", "Mc2r", "Mc3r", "Mc4r", "Mc5r", "Mchr1", "Mtnr1a", "Mtnr1b", "Nmb", "Nmbr", "Nms", "Nmu", "Nmur1", "Nmur2", "Npb", "Npbwr1", "Npff", "Npffr1", "Npffr2", "Nps", "Npsr1", "Npvf", "Npw", "Npy", "Npy1r", "Npy2r", "Npy4r", "Npy5r", "Npy6r", "Nr3c1", "Nts", "Ntsr1", "Ntsr2", "Ogfr", "Ogfrl1", "Oprd1", "Oprk1", "Oprl1", "Oprm1", "Oxt", "Oxtr", "P2rx1", "P2rx2", "P2rx3", "P2rx4", "P2rx5", "P2rx6", "P2rx7", "P2ry1", "P2ry10", "P2ry10b", "P2ry13", "P2ry14", "P2ry2", "P2ry4", "P2ry6", "Paqr6", "Paqr9", "Pate1", "Pate10", "Pate13", "Pate2", "Pate3", "Pate4", "Pate5", "Pate6", "Pate7", "Pdyn", "Penk", "Pgr15l", "Pgrmc1", "Pgrmc2", "Plg", "Pmch", "Pnoc", "Pomc", "Ppy", "Prl", "Prl5a1", "Prl6a1", "Prlh", "Prlhr", "Prlr", "Prss1", "Prss1l", "Prss2", "Prss3", "Prss3b", "Prss3l", "Ptafr", "Ptgdr", "Ptger1", "Ptger2", "Ptger3", "Ptger4", "Ptgfr", "Ptgir", "Pth", "Pth1r", "Pth2", "Pth2r", "Pyy", "Qrfp", "Qrfpr", "Qrfprl", "Rln1", "Rln3", "Rxfp1", "Rxfp2", "Rxfp3", "Rxfp4", "S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sct", "Sctr", "Slurp1", "Slurp2", "Spx", "Sst", "Sstr1", "Sstr2", "Sstr3", "Sstr4", "Sstr5", "Taar1", "Taar2", "Taar3", "Taar4", "Taar5", "Taar6", "Taar7a", "Taar7b", "Taar7d", "Taar7e", "Taar7f", "Taar8a", "Taar8b", "Taar8c", "Taar9", "Tac1", "Tac2", "Tac4", "Tacr1", "Tacr2", "Tacr3", "Tbxa2r", "Thra", "Thrb", "Tmem97", "Trh", "Trhr", "Trhr2", "Trpv1", "Try10", "Try4", "Try5", "Tshb", "Tshr", "Tspo", "Ucn", "Ucn2", "Ucn3", "Uts2", "Uts2b", "Uts2r", "Vgf", "Vip", "Vipr1", "Vipr2")#{{{#}}}

all.genes <- readRDS('./data/1AllGenes.rds')
fpkm <- readRDS("./data/0FPKM.rds")
gene.map <- bitr(rownames(all.genes),
		   fromType = "ENSEMBL",
           toType = c("SYMBOL"),
           OrgDb = "org.Mm.eg.db")

unique.symbol = unique(gene.map$SYMBOL)
all.genes <- all.genes[gene.map[match(unique.symbol, gene.map$SYMBOL), "ENSEMBL"], ] %>%
	mutate(symbol = unique.symbol)
fpkm <- fpkm[gene.map[match(unique.symbol, gene.map$SYMBOL), "ENSEMBL"], ]
rownames(fpkm) <- unique.symbol

interested.genes <- all.genes %>%
	filter(logFC < 0.5, PValue < 0.05) %>%
	arrange(logFC) %>%
	filter(symbol %in% interested)

fpkm.scaled <- fpkm %>%
	apply(1, scale) %>%
	t()
colnames(fpkm.scaled) <- colnames(fpkm)
fpkm.scaled <- as.matrix(fpkm.scaled[interested.genes$symbol, ])

labels_row = rownames(fpkm.scaled)
labels_row[labels_row %in% c("Tacr1", "Tacr2") == FALSE] <- ""
pheatmap(fpkm.scaled, labels_row = labels_row, color = colorRampPalette(colors = c("#3EC1D3","#eee","#FF165D"))(100), cluster_cols = F)

