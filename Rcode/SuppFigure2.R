library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
#library(enrichR)

set.seed(1)

#Ara_UberAlign18G_withMagdaEmb_Tb2
#saveext = "./Ara_UberAlign18G_withMagdaEmb_Tb_noSYS2
saveext = "./Ara_UberAlign18G/"

dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Marker/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

mammal.combined <- readRDS(file = paste(saveext,"mammal.combined.KO5.rds",sep=""))

saveext = "./Ara_UberAlign18G_PreMERR3/"

dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Marker/",sep=""))

#Get Ara's data ready
Idents(mammal.combined) <- mammal.combined$species
mammal.combined <- subset(mammal.combined,idents=c("Human (Ara EBs)"))

Idents(mammal.combined) <- mammal.combined$Cl9

p <- DimPlot(mammal.combined, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"tEB.Cl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)

#Now do scatter plots
TF<-read.table("TF.txt",header = F)
TF <- TF$V1
SIGNAL<-read.table("LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

#Idents(mammal.combined) <- mammal.combined$Cl5
AE <- AverageExpression(mammal.combined)
Ae <- AE$RNA
#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("2"), ident.1 = c("6"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("POU5F1","TFAP2A","NKX1-2","TERF1","NFE2L3","TERF1","NANOG","TFAP2C","SOX17",
"MIXL1","GATA4","MYC","FOXF1","SNAI1","HOXB5","MESP2","HOXB6","HOXB2","TCF4","GATA6","MESP1","SNAI2","HAND1")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_2_6.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'2' + Ae5$'6') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_6.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_6_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("2"), ident.1 = c("0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("POU5F1","MESP1","MILX1","FOXH1","SOX17","HMGA1","GATA6","JUNB","HMGA2","CREB3L1","PITX1","PITX2","HAND1","SNAI2","HOXB6")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_2_0.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'2' + Ae5$'0') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_0.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_0_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("2"), ident.1 = c("13"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("YBX1","SP5","MIXL1","ETV4","FOXH1","EOMES","CDX1","MESP1","SOX17",
"GATA2","TFAP2C","TFAP2A","DLX5","ISL1","MAFB","MEIS2","GATA3","PPARG","HAND1","EPAS1","TSC22D1")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_2_13.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'2' + Ae5$'13') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_13.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_13_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#14,1
#19,1
#3,22
#14,22
#19,22
#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("2"), ident.1 = c("8","5","7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("FOXH1","SP5","ETV4","CENPX","MESP1","POU5F1","EOMES","YBX1","NKX1-2","THYN1","YBX3","DRAP1","NANOG","SOX17",
"TFAP2A","SNAI2","GATA2","SOX7","GATA3","DLX5","TBX3","PHOX2A","HEY1","ISL1","MAFB","EPAS1","HAND1","HOXB6")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_2_857.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'2' + (0.3166870*Ae5$'5'+ 0.4299635*Ae5$'8'+ 0.2533496*Ae5$'7')  ) +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_857.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_857_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(mammal.combined, ident.2 = c("2"), ident.1 = c("17","20","28"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("POU5F1","HMGA1","CXXC5","EOMES","MESP1","FOXA2","HNF1B","CDX2","GATA6","SOX17","TSC22D1",
"HOXB6","HAND1","PITX1")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_2_172028.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'2' + (0.70247934*Ae5$'17'+0.09917355*Ae5$'20'+0.19834711*Ae5$'28')  ) +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_172028.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_2_172028_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(mammal.combined, ident.2 = c("10"), ident.1 = c("13"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("HMGA1","YBX1","CDX1","ETV4","SOX17","EPAS1","PPARG","MEIS2","TSC22D1","HAND1","MSX2","CEBPA","GATA3")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_10_13.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'10' + Ae5$'13') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_10_13.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_10_13_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(mammal.combined, ident.2 = c("10"), ident.1 = c("9"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("HMGA1","SOX3","CDX1","HMGA2","MSX1","MESP1","CEBPZ","SOX17","TFAP2C","DMRTA2","NANOG","GTF3A","PRDM1","TSC22D1","MAFB","SOX15","TET1","POU5F1")
"HMGA1","YBX1","CDX1","ETV4","SOX17","EPAS1","PPARG","MEIS2","TSC22D1","HAND1","MSX2","CEBPA","GATA3")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_10_9.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'10' + Ae5$'9') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_10_9.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_10_9_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(mammal.combined, ident.2 = c("13"), ident.1 = c("9"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("PPARG","EPSA1","CEBPA","DLX5",",MEIS2","GATA3","ISL1","GATA3","REST","DMRTA2","POU5F1","PRDM1","GTF2I","NANOG","TFAP2C","CEBPZ","NKX1-2","SOX17")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_13_9.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'13' + Ae5$'9') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_13_9.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_13_9_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Idents(mammal.combined) <- mammal.combined$Cl9
mammal.combined2 <- subset(mammal.combined,idents=c(9))
mammal.combined2 <- FindNeighbors(mammal.combined2, reduction = "pca", dims = 1:20)

mammal.combined2  <- FindClusters(mammal.combined2, resolution = .5)

p <- DimPlot(mammal.combined2, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"tEBrefcl",".pdf",sep=""),width = 8, height = 8, limitsize = FALSE,p)

AE <- AverageExpression(mammal.combined2)
Ae <- AE$RNA

Cl1 <- FindMarkers(mammal.combined2, ident.2 = c("5"), ident.1 = c("0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("EPAS1","TSC22D1","HAND1","ISL1","SOX4","DLX5","YBX1","MYBL2","CENPA","FOXM1","FOXH1","ETV4","NR2F6")
#"PPARG","EPSA1","CEBPA","DLX5",",MEIS2","GATA3","ISL1","GATA3","REST","DMRTA2","POU5F1","PRDM1","GTF2I","NANOG","TFAP2C","CEBPZ","NKX1-2","SOX17")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_0_5_pgc.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
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
Ae5$AvExp <- log(0.5*(Ae5$'0' + Ae5$'5') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_0_5_pgc.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_0_5_pgc_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

