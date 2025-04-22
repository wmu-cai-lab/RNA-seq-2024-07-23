setwd("D:/老师的生信处理/蔡老师文件夹/11.3派森诺 RAW and EPI/RAW")
library(limma)
exp <- read.table(file = "matrix.txt",header = T, quote="\"")
#exp <-read.csv("去除批次后.csv",row.names = 1)
exp <- as.matrix(exp)
rownames(exp)<-exp[,1]
exp=exp[,2:ncol(exp)]
mode(exp)
dimnames=list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames) ##变成数字形式
exp <- avereps(exp) ##重复的基因名 平均化
exp[exp==0]<-NA
exp=na.omit(exp)
exp<-log2(exp+1)
range(exp)
boxplot(exp)
exp=normalizeBetweenArrays(exp)
###PCA
dat=as.data.frame(t(exp))
#BiocManager::install("factoextra")#按照包
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra)

# pca的统一操作走起
group=c(rep("Normal",4),rep("Tumor",4))
dat.pca <- PCA(dat, graph = FALSE)

pca_plot <- fviz_pca_ind(dat.pca,
                         
                         geom.ind = "point", # show points only (nbut not "text")
                         
                         col.ind = group, # color by groups
                         
                         palette = c("#016DDB", "#E43401"),
                         
                         addEllipses = T, # Concentration ellipses
                         
                         legend.title = "Groups"
                         
)+theme(panel.grid=element_blank())

pca_plot

##热图
cg=names(tail(sort(apply(exp,1,sd)),1000))

n=exp[cg,]

#绘制热图

annotation_col=data.frame(group=group)

rownames(annotation_col)=colnames(n)

library(pheatmap)

pheatmap(n,
         
         show_colnames = F,
         
         show_rownames = F,
         
         annotation_col=annotation_col,
         
         cluster_cols= F,
         
         scale = "row")

dev.off()

#co=cor(exp)；pheatmap(co)这两句代码可做相关性热图

rm(list = ls())

load(file = "step2output.Rdata")

library(dplyr)

exp1 <- as.data.frame(exp)

exp2 <- mutate(exp1,probe_id=rownames(exp))

exp3 <- inner_join(exp2,ids,by="probe_id")

exp4 = subset(exp3,symbol == "CCL13")

write.csv(exp4,file = "CCL13.csv")


##差异分析

design=model.matrix(~group)

fit=lmFit(exp,design)

fit=eBayes(fit)

deg=topTable(fit,coef=2,number = Inf)

#为deg数据框添加几列

#1.加probe_id列，把行名变成一列

library(dplyr)

deg <- mutate(deg,SYMBOL=rownames(deg))

#或者 deg$probe_id=rownames(deg)

head(deg)


#3.加change列,标记上下调基因
setwd("D:/老师的生信处理/同学/GBM/二次分析/DEGs_PPP1R18")
deg<-read.csv("allDIFF.csv",row.names = 1)
logFC_t=0.5

P.Value_t = 0.05

k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t)#可以用table(k1)看看具体有几个
table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t)
table(k2)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))#或者用case_when()
table(change)
deg <- mutate(deg,change)
write.csv(deg, file = "DEGs_change.csv", row.names = FALSE)

mycolor<- c( "#016DDB","gray80","#E43401")
ggplot(
  deg, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
  geom_point(alpha=0.65, size=3) +
  scale_color_manual(values=mycolor)+
  
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=1.0) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=1.0) +
  
  labs(x="log2(Fold Change)",
       y="-log10(FDR)")+
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+theme(panel.grid=element_blank())
dev.off()

#######GSEA富集分析
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library ("org.Mm.eg.db") #人是Hs；小鼠 Mm；大鼠 Rt；
library(patchwork)
library(WGCNA)
library(GSEABase)
library(GSVA)
fix(deg)
gene.df <- bitr(deg$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db) 
Enrich <- gene.df %>%
  inner_join(deg, by = "SYMBOL")
head(Enrich)
Enrichment <- Enrich %>%
  as.data.frame()%>%
  dplyr::select(ENTREZID,logFC)
Enrichment <- Enrichment[!duplicated(Enrichment[,1]),]
geneList <- Enrichment$logFC
names(geneList) <- Enrichment$ENTREZID
geneList <-sort(geneList,decreasing = T)
GSEA1 <- gseKEGG(geneList     = geneList,
                 organism     = 'mmu',
                 nPerm        = 10000,
                 minGSSize    = 1,
                 maxGSSize    = 10000,
                 pvalueCutoff = 1,
                 pAdjustMethod = "none" ) # hsa人 mmu鼠
write.csv(GSEA1,file="GSEA.csv",quote=F,row.names = F) #
###山峦图展示
library(enrichplot)
library(ggplot2)
Category=c("VEGF signaling pathway","B cell receptor signaling pathway","ECM-receptor interaction", "PD-L1 expression and PD-1 checkpoint pathway in cancer",
           "NF-kappa B signaling pathway","HIF-1 signaling pathway","JAK-STAT signaling pathway","Focal adhesion","PI3K-Akt signaling pathway",
           "Hippo signaling pathway"
) ##选择展示部分
#Category=c("Chemical carcinogenesis - reactive oxygen species","Bacterial invasion of epithelial cells",
  #         "PPAR signaling pathway","ECM-receptor interaction","Focal adhesion")
ridgeplot(GSEA1,
          showCategory = Category,##如果写数字
          fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = T
)+
  theme(axis.text.y = element_text(size=8))
## Picking joint bandwidth of 0.212
??ridgeplot()

###GO and KEGG
deg2 = deg %>% 
  filter( abs( logFC ) > 0.5 & adj.P.Val < 0.05 )# &：和；|：或。
  #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选
gene.df <- bitr(deg2$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db) 
gene <- gene.df$ENTREZID
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1,ont="all",readable =T)

write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #
write.csv(go,file="GO.csv",quote=F,row.names = F) #
##可视化
##条形图
pdf(file="GO-barplot.pdf",width = 10,height = 15)
barplot(go, drop = TRUE, showCategory =5,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')+theme(panel.grid=element_blank())
dev.off()

##气泡图
pdf(file="GO-bubble.pdf",width = 10,height = 15)
dotplot(go,showCategory = 5,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')+theme(panel.grid=element_blank())
dev.off()

#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "mmu", pvalueCutoff =1, qvalueCutoff =1, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         
write.csv(kk,file="KEGG.csv",quote=F,row.names = F) #
##可视化
##条形图
KEGGcategory<-c("Cytokine-cytokine receptor interaction","ECM-receptor interaction","Focal adhesion","Neuroactive ligand-receptor interaction","PI3K-Akt signaling pathway",
                "Chemokine signaling pathway","NF-kappa B signaling pathway","MAPK signaling pathway","Cell adhesion molecules","HIF-1 signaling pathway",
                "JAK-STAT signaling pathway","cGMP-PKG signaling pathway","Apoptosis","B cell receptor signaling pathway","Hippo signaling pathway"
)
keggcategory<-c("Chemical carcinogenesis - reactive oxygen species","Bacterial invasion of epithelial cells",
                "PPAR signaling pathway","ECM-receptor interaction","Focal adhesion","PI3K-Akt signaling pathway",
                "cGMP-PKG signaling pathway","TGF-beta signaling pathway","Cell adhesion molecules")

pdf(file="KEGG-barplot.pdf",width = 10,height = 13)
barplot(kk, drop = TRUE, showCategory = KEGGcategory,label_format=1000)
dev.off()

##气泡图
pdf(file="KEGG-bubble.pdf",width = 10,height = 13)
dotplot(kk, showCategory = KEGGcategory,label_format=1000)+theme(panel.grid=element_blank())
dev.off()

####提取KEGG/GSEA通路的DEGs
kure<-kk@result
View(kure)
GS<-kk@geneSets
View(GS)
a<-GS$mmu04668
PI3K_symbol <- bitr(a,
                    fromType = "ENTREZID",
                    toType = c("SYMBOL"),
                    OrgDb = org.Mm.eg.db)
View(PI3K_symbol)
PI3K_symbol<-PI3K_symbol$SYMBOL
PI3K_symbol<-as.data.frame(PI3K_symbol)
PI3K_symbolexp<-deg2[deg2$SYMBOL %in% PI3K_symbol$PI3K_symbol, ]###提取通路的差异分析

write.csv(PI3K_symbolexp, file = "ROS_DEGs.csv", row.names = FALSE) ##导出数据

PI3K_symbolexpress<-rt[rt$SYMBOL %in% PI3K_symbol$PI3K_symbol, ]###表达矩阵
rownames(PI3K_symbolexp)<-NULL
rownames(PI3K_symbolexpress)<-NULL
rownames(PI3K_symbolexpress)=PI3K_symbolexpress[,1]
exp=PI3K_symbolexpress[,2:ncol(PI3K_symbolexpress)]
dimnames=list(rownames(exp),colnames(exp))
PI3K_symbolexpress=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
PI3K_symbolexpress=avereps(PI3K_symbolexpress)
rownames(PI3K_symbolexp)=PI3K_symbolexp[,1]