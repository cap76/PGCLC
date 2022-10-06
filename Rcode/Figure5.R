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

cType <- c(0,1,3,10,5,11,2,8,4,7,9,6,12,13)
BaseCol <- c("#e8a843","#da8f05","#ce8820","#574c9b","#4f99c3","#2891bf","#1381c4","#0e9739","#e83970","#e01551","#c8114a","lightgrey","lightgrey","lightgrey")

colind <- integer( length( levels(mammal.combined$Cl5) )  )
for (i in 1:length( levels(mammal.combined$Cl5) ) ) {
  colind[i] <- which(cType==levels(mammal.combined$Cl5)[i])
}
coluse1 <- BaseCol[colind]

p <- DimPlot(mammal.combined, cols = coluse1, shape.by = "cell.orig", pt.size = 2,reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_split.KO.lab_Cl5",".pdf",sep=""),width = 45, height = 15, limitsize = FALSE, p)


newOrder <- as.character(mammal.combined$Cells)
newOrder[which(newOrder=="PGCLC_12h_r2")] <- "1) 12h"
newOrder[which(newOrder=="PGCLC_18h")] <- "2) 18h"
newOrder[which(newOrder=="PGCLC_32h")] <- "3) 32h"
newOrder[which(newOrder=="PGCLC_40h")] <- "4) 40h"
newOrder[which(newOrder=="PGCLC_48h")] <- "5) 48h"
newOrder[which(newOrder=="PGCLC_D4")] <- "6) 96h"
newOrder[which(newOrder=="18h_cont" )] <- "7) 18h C"
newOrder[which(newOrder=="18h_TFAP2a")] <- "8) 18h T"
newOrder[which(newOrder=="D4_cont" )] <- "9) 96h C"
newOrder[which(newOrder=="D4_TFAP2A")] <- "99) 96h T"

Ano <- mammal.combined$CellLabs

require(tidyverse)
Dat <- data.frame(x=newOrder[which(newOrder!="hESC_W5_E8")],y=Ano[which(newOrder!="hESC_W5_E8")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
ggsave(filename=paste(saveext,"MarmTable","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


NewLabs <- paste(mammal.combined$species,mammal.combined$Cl5,sep="_")
Idents(mammal.combined) <- NewLabs

#Plot expression
DefaultAssay(mammal.combined) <- "RNA"
#
FeaturePlot(mammal.combined, features = "NANOS3", pt.size = 2, cols = c("lightgrey", "#dd1c77"), split.by="species") # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_NANOS3_nosplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "SOX17", pt.size = 2, cols = c("lightgrey", "#dd1c77"), split.by="species") # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_SOX17_nosplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "PDGFRA", pt.size = 2, cols = c("lightgrey", "#dd1c77"), split.by="species") # + xlim(c(-13, 13)) + ylim(c(-13, 13))##
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_PDGFRA_nosplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "VTCN1", pt.size = 2, cols = c("lightgrey", "#dd1c77"), split.by="species") # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_VTCN1_nosplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)

#Now do scatter plots
TF<-read.table("TF.txt",header = F)
TF <- TF$V1

SIGNAL<-read.table("LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

Idents(mammal.combined) <- mammal.combined$Cl5
NewID <- paste(mammal.combined$species, mammal.combined$Cells, mammal.combined$Cl5,sep="_")
NewID[which(mammal.combined$Cells=="hESC_W5_E8")] <- "hESC_W5"
NewID[which(NewID=="4) TFAP2A_D4_TFAP2A_8")] <- "TFAP2Ad"
Idents(mammal.combined) <- NewID

AE <- AverageExpression(mammal.combined)
Ae <- AE$RNA

#Note FindMarkers Returns logFC as natural log
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("4) TFAP2A_D4_TFAP2A_10"), ident.1 = c("2) Parental_D4_cont_9"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
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
genes.to.label5 <- c("TERF1","SOX2","ZIC2","ZIC5","SALL3","TCF7L1","SOX11","HMGA1","TCF7L2","OTX2","GLIS3","MSC","ZHX2","CDX1","TCFP2L1","VENTX","NANOG","BARX1","SOX15","SOX18","IRX6","ARID5B","PRDM1","SOX17")
p1 <- ggplot(Ae5[,], aes(FC1,Pval1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
 geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = -log(1.5) )) + geom_vline(aes(xintercept = log(1.5))) +
 geom_hline(aes(yintercept = -log2(0.01))) 
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta", y = "-log Pval")
ggsave(filename=paste(saveext,"VolanoPlot_PSC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
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
Ae5$AvExp <- log(0.5*(Ae5$'4) TFAP2A_D4_TFAP2A_10' + Ae5$'2) Parental_D4_cont_9') +1)

p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
#geom_vline(aes(xintercept = 0)) #+ geom_vline(aes(xintercept = -log(1.5) )) 
geom_hline(aes(yintercept = 0))
# geom_hline(aes(yintercept = -log2(0.01)))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_PSC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
geom_hline(aes(yintercept = 0))
p1 <- LabelPoints(plot = p1, points = genes.to.label5a, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_PSC_PGCLC_FL.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#AE2 <- AverageExpression(mammal.combined2)
Idents(mammal.combined) <- paste(mammal.combined$species,mammal.combined$Cl5,sep="_")
AE <- AverageExpression(mammal.combined)

g2 <- c("POU5F1","NANOG","SOX2","SOX3","PRDM14","KLF2","KLF4","TFCP2L1","THY1","DUSP6","ZIC2","ZIC3","ZIC5","CD24","SNAI2","PDGFRA","GATA6","FOXF1",
"VTCN1","WNT6","ISL1","TFAP2A","GATA3","SOX17","PRDM1","TFAP2C","PDPN","NANOS3")

exp1 <- c("1) EBs_10","4) TFAP2A_10","1) EBs_2","2) Parental_2","4) TFAP2A_2","1) EBs_8","2) Parental_8","4) TFAP2A_8","1) EBs_9","2) Parental_9")
C1 <- log2(AE$RNA[g2,exp1]+1)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20), border_color = NA, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/Heatmap",".pdf",sep="") ,width=5,height=10)


#
ExpLevel <- log(1)


S17P <- WhichCells(mammal.combined, expression = SOX17 > ExpLevel )

l0 <- colnames(mammal.combined)[which(mammal.combined$Cells=="PGCLC_12h_r2")]
l1 <- colnames(mammal.combined)[which(mammal.combined$Cells=="PGCLC_18h")]
l2 <- colnames(mammal.combined)[which(mammal.combined$Cells=="18h_cont")]
l3 <- colnames(mammal.combined)[which(mammal.combined$Cells=="18h_TFAP2a")]

Counts <- data.frame(x=c(length(intersect(S17P,l0)) / length(l0),
length(intersect(S17P,l1)) / length(l1),
length(intersect(S17P,l2)) / length(l2),
length(intersect(S17P,l3)) / length(l3) ), y = c("12h","18h","18h Parental", "18h KO"), type=c("SOX17+","SOX17+","SOX17+","SOX17+"))

ggplot(Counts, aes(fill=type, x=y, y=x)) + geom_bar(position="stack", stat="identity")
ggsave(filename=paste(saveext,"S17Tab","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

