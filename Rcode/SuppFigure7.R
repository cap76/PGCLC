library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
#library(enrichR)
library(pheatmap)

set.seed(1)
saveext1 = "./Ara_UberAlign18G_wKO6_v6_rerun/"

saveext = "./Ara_UberAlign18G_wKO6_v6_rerun_rr_check/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Marker/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

mammal.combined <- readRDS(file=paste(saveext,"FinalTFAP2AKO_integrated_dataset",".rds",sep=""))

Idents(mammal.combined) <- mammal.combined$Cl5

NewLabs <- paste(mammal.combined$species,mammal.combined$Cl5,sep="_")
Idents(mammal.combined) <- NewLabs
AE <- AverageExpression(mammal.combined)
Ae <- AE$RNA
#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("4) TFAP2A_10"), ident.1 = c("1) EBs_10"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01))
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)

Ae5$ColInd <- 0
Ae5[ which( (Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
genes.to.label5a = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- c("HMGA1","YBX1","CENPX","HMG3","NR6A1","TEAD2","DRAP1","PIN1","TFDP1","THYN1","FOXH1","PRMT3","FOXM1","TGIF1","TFDP2","CREB3",
"SOX2","HESX1","TERF1","POU5F1","GATA2","SP8","MBD2","SOX11","IRF1","MEIS3","KLF4","SOX4","IRX2","PRDM14","ETV1","TET3","YY1","TFAP2C","SP3","ZIC2",
"POU3F1","MBD3","NANOG","ZIC5","SALL4","KMT2B","FOXO1","TEAD1","FOXN3","KDM2A","SALL3","GATAD2A","ZNF292","TET1","PRDM2",
"SON","REST","POU2F1","PBX1","TCF3","ZFX")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_PSC_ESC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Add in MA plot for review
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0
Ae5[ which( (Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'4) TFAP2A_10' + Ae5$'1) EBs_10') +1)

p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
#geom_vline(aes(xintercept = 0)) #+ geom_vline(aes(xintercept = -log(1.5) )) 
geom_hline(aes(yintercept = 0))
# geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_PSC_ESC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_PSC_ESC_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Idents(mammal.combined) <- paste(mammal.combined$species,mammal.combined$Cl5,sep="_")
AE <- AverageExpression(mammal.combined)

#g2 <- c("POU5F1","NANOG","SOX2","SOX3","PRDM14","KLF2","KLF4","TFCP2L1","THY1","DUSP6","ZIC2","ZIC3","ZIC5","CD24","SNAI2","PDGFRA","GATA6","FOXF1",
#"VTCN1","WNT6","ISL1","TFAP2A","GATA3","SOX17","PRDM1","TFAP2C","PDPN","NANOS3")

exp1 <- c("1) EBs_10","4) TFAP2A_10","1) EBs_2","2) Parental_2","4) TFAP2A_2","1) EBs_8","2) Parental_8","4) TFAP2A_8","1) EBs_9","2) Parental_9")
C1 <- cor( log2(AE$integrated[,exp1]+1), log2(AE$integrated[,exp1]+1),  ,method="pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20),border_color = NA,cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/Corr",".pdf",sep="") ,width=5,height=5)


#g1 <- readtable()
g1 <- rownames(GetAssayData(mammal.combined,assay="integrated"))

exp1 <- c("1) EBs_10","4) TFAP2A_10","1) EBs_2","2) Parental_2","4) TFAP2A_2","1) EBs_8","2) Parental_8","4) TFAP2A_8","1) EBs_9","2) Parental_9")

C1 <- cor( log2(AE$RNA[g1,exp1]+1), log2(AE$RNA[g1,exp1]+1),  ,method="pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, border_color = NA,cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/Corr","_RNAmarkers.pdf",sep="") ,width=5,height=5)

#Now generate some scatter plots
Idents(mammal.combined) <- mammal.combined$Cells
mammal.combined <- subset(mammal.combined,idents=c("PGCLC_24h","PGCLC_18h","PGCLC_D4","PGCLC_12h_r2","PGCLC_D4_r2"))
DefaultAssay(mammal.combined) <- "RNA"
p<-FeatureScatter(object = mammal.combined, feature1 = "NANOS3", feature2 = "PDPN")
ggsave(filename=paste(saveext,"SCATTER_lab",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "NANOS3", feature2 = "TFAP2C")
ggsave(filename=paste(saveext,"SCATTER_lab1",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "NANOS3", feature2 = "TFAP2A")
ggsave(filename=paste(saveext,"SCATTER_lab2",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "TFAP2A", feature2 = "TFAP2C")
ggsave(filename=paste(saveext,"SCATTER_lab3",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "TFAP2A", feature2 = "VTCN1")
ggsave(filename=paste(saveext,"SCATTER_lab4",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "NANOS3", feature2 = "VTCN1")
ggsave(filename=paste(saveext,"SCATTER_lab5",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "SOX17", feature2 = "VTCN1")
ggsave(filename=paste(saveext,"SCATTER_lab6",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)
p<-FeatureScatter(object = mammal.combined, feature1 = "TFAP2A", feature2 = "VTCN1")
ggsave(filename=paste(saveext,"SCATTER_lab7",".pdf",sep=""),width = 20, height = 15, limitsize = FALSE, p)


