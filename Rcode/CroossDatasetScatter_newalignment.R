library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/AlignAraDataCross"
fileID <- saveext
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

fileID2 <- "/Volumes/Overflow/Uber18G/"
mammal.combined <- readRDS(paste(fileID2,'mammal.combined.KO5.rds',sep=""))

labs <- as.character( mammal.combined$Cells )
labs[which(labs=="MeLC_R1")] <- "PreME"
labs[which(labs=="MeLC_R2")] <- "PreME"
labs[which(labs=="MeLC_R3")] <- "PreME"
labs[which(labs=="MeLC_R4")] <- "PreME"
labs[which(labs=="MELC_KO")] <- "PreME"
labs[which(labs=="ESC_R1")] <- "ESC"
labs[which(labs=="ESC_R2")] <- "ESC"
labs[which(labs=="ESC_R3")] <- "ESC"
labs[which(labs=="ESC_R4")] <- "ESC"
labs[which(labs=="ESC_KO")] <- "ESC"
labs[which(labs=="PGCLC_D1_R1")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R2")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R3")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R4")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_KO")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D2_R1")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R2")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R3")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R4")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_KO")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D3_R1")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R2")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R3")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R4")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_KO")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D4_R1")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R2")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R3")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R4")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_KO")] <- "PGCLC_96h"
labs[which(labs=="hESC")] <- "ESC"
labs[which(labs=="hPGCLC")] <- "PGCLC"
labs[which(labs=="twAMLC")] <- "AMLC"
labs[which(labs=="MELC1")] <- "MELC1"
labs[which(labs=="MELC2")] <- "MELC2"
labs[which(labs=="PreME_10h")] <- "PreME"
labs[which(labs=="hESC_W5_E8")] <- "ESC"
labs[which(labs=="PGCLC_D4_r2")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_12h_r2")] <- "PGCLC_12h"
labs[which(labs=="TE")] <- "Tb"
labs[which(labs=="PGC_wk6.5_invivo")] <- "PGC"
labs[which(labs=="ME_r2")] <- "ME"

mammal.combined$ID2 <- as.factor(labs)
Idents(mammal.combined) <- as.factor(labs)

#Rejig IDs and order
Species2 <- as.character(mammal.combined$species)
Species2[which(Species2=="Human (Ara EBs)")] <- "1) Human EBs (in vitro)"
Species2[which(Species2=="Human (Amnioid)")] <- "2) Human Amnioids (in vitro)"
Species2[which(Species2=="Human (CS7)")] <- "5) Human CS7 (in vivo)"
Species2[which(Species2=="4) Human (in vitro)")] <- "4) Human (in vitro 1)"
Species2[which(Species2=="human_invivo")] <- "5) Human (in vitro 2)"
Species2[which(Species2=="3) Cynomologous (in vitro)")] <- "6) Cynomolgus (in vitro)"
Species2[which(Species2=="3) Marmoset")] <- "7) Marmoset (in vivo)"
Species2[which(Species2=="10X_1")] <- "4) Human EBs (Clark1)"
Species2[which(Species2=="10X_2")] <- "4) Human EBs (Clark1)"
Species2[which(Species2=="10X_3")] <- "4) Human EBs (Clark2)"
Species2[which(Species2=="10X_4")] <- "4) Human EBs (Clark2)"
Species2[which(Species2=="10X_5")] <- "8) Human EBs (Clark)"
Species2[which(Species2=="Gastruloid2D")] <- "3) Gastruloid"
Species2[which(mammal.combined$Cells=="PGC_wk6.5_inviv")] <- "1) Human EBs (in vivo)"


mammal.combined$species2 <- Species2
ID1 <- paste(mammal.combined$Cl9,mammal.combined$species2,sep="_")

ID2 <- ID1

ID2[which(labs=="PreME")] <- paste(ID2[which(labs=="PreME")],"PreME",sep="_")
ID2[which(labs=="ESC")] <- paste(ID2[which(labs=="ESC")],"ESC",sep="_")
ID2[which(labs=="4i")] <- paste(ID2[which(labs=="4i")],"4i",sep="_")


Idents(mammal.combined) <- ID2


Idents(mammal.combined) <- mammal.combined$Cl9


SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1


#mammal.combined$NewestID <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$NewestID

p <- DimPlot(mammal.combined,  reduction = "umap",  label = TRUE, split.by = "species", repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl11_Am.pdf",sep=""),width = 8*12, height = 8, p, limitsize = FALSE)


p <- FeaturePlot(mammal.combined, feature = "SOX2") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl11_Am.pdf",sep=""),width = 40, height = 30, p, limitsize = FALSE)


#1 is ESC
#3 is PreME
#2 is PS
#10 iis early am
#9,4 PGC
#13,5,7,8 am
#6 early mes
#0 late mes
#17,,20 DE

DefaultAssay(mammal.combined) <- "RNA"
Ae <- AverageExpression(mammal.combined)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)

#saveRDS(mammal.combined2, file = paste(saveext,"mammal_combined_3D.rds",sep=""))

#ESC vs PS 
#ESC PGCLC
#PGC v PS
#Am vs preme
#Mes vs preme
#Endo vs preme

#Ara CS7
#Ara amn
#Ara gast



#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("17_1) Human EBs (in vitro)","20_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("17_HumanCS7","20_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which(( Ae5$`17_1) Human EBs (in vitro)`>0.1 | Ae5$`20_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`17_HumanCS7`>0.1 | Ae5$`20_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_DE_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_DE_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("2_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("2_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`2_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`2_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_PS_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_PS_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("4_1) Human EBs (in vitro)","9_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("4_HumanCS7","9_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`4_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`4_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_PGC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_PGC_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#AdMe vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("0_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`0_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`0_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_AME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_AME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#eMes vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#1 is ESC
#3 is PreME
#2 is PS
#10 iis early am
#9,4 PGC
#13,5,7,8 am
#6 early mes
#0 late mes
#17,,20 DE

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("13_1) Human EBs (in vitro)","5_1) Human EBs (in vitro)","7_1) Human EBs (in vitro)","8_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("13_HumanCS7","5_HumanCS7","7_HumanCS7","8_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`13_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`13_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_LAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_LAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("10_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("10_HumanCS7"), ident.2 = c("1_HumanCS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`10_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`10_HumanCS7`>0.1 | Ae5$`1_HumanCS7`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " CS7 ")
ggsave(filename=paste(saveext,"/Ara_CS7_EAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#######################################
#Now gastruloids

#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("2_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("2_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`2_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`2_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_PS_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_PS_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("4_1) Human EBs (in vitro)","9_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("4_3) Gastruloid","9_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`4_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`4_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_PGC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_PGC_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("17_1) Human EBs (in vitro)","20_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("17_3) Gastruloid","20_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`17_1) Human EBs (in vitro)`>0.1 | Ae5$`20_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`20_3) Gastruloid`>0.1 | Ae5$`17_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_DE_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_DE_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#AdMe vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("0_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`0_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`0_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_AME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_AME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#eMes vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("13_1) Human EBs (in vitro)","5_1) Human EBs (in vitro)","7_1) Human EBs (in vitro)","8_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("13_3) Gastruloid","5_3) Gastruloid","7_3) Gastruloid","8_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`13_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`13_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_LAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_LAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("10_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("10_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`10_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`10_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_EAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#######################################
#Now amnioids

#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("2_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("2_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`2_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`2_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_PS_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_PS_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("4_1) Human EBs (in vitro)","9_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("4_2) Human Amnioids (in vitro)","9_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`4_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`4_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_PGC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_PGC_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#AdMe vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("0_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`0_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`0_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_AME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_AME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#eMes vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_2) Human Amnioids (in vitro)_ESC`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("13_1) Human EBs (in vitro)","5_1) Human EBs (in vitro)","7_1) Human EBs (in vitro)","8_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("13_2) Human Amnioids (in vitro)","5_2) Human Amnioids (in vitro)","7_2) Human Amnioids (in vitro)","8_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`13_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`13_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_LAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_LAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("10_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("10_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`10_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`10_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_EAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("17_1) Human EBs (in vitro)","20_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("17_2) Human Amnioids (in vitro)","20_2) Human Amnioids (in vitro)"), ident.2 = c("1_2) Human Amnioids (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`17_1) Human EBs (in vitro)`>0.1 | Ae5$`20_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`20_2) Human Amnioids (in vitro)`>0.1 | Ae5$`17_2) Human Amnioids (in vitro)`>0.1 | Ae5$`1_2) Human Amnioids (in vitro)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_DE_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Amn ")
ggsave(filename=paste(saveext,"/Ara_Amn_DE_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()







#Now Clark1

#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("2_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("2_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`2_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`2_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)_ESC`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_PS_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_PS_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("4_1) Human EBs (in vitro)","9_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("4_4) Human EBs (Clark1)","9_2) Human Amnioids (in vitro)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`4_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`4_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_PGC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_PGC_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#AdMe vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("0_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`0_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`0_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_AME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_AME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#eMes vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("13_1) Human EBs (in vitro)","5_1) Human EBs (in vitro)","7_1) Human EBs (in vitro)","8_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("13_4) Human EBs (Clark1)","5_4) Human EBs (Clark1)","7_4) Human EBs (Clark1)","8_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`13_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`13_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_LAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_LAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("10_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("10_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`10_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`10_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_EAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("17_1) Human EBs (in vitro)","20_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("17_4) Human EBs (Clark1)","20_4) Human EBs (Clark1)"), ident.2 = c("1_4) Human EBs (Clark1)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`17_1) Human EBs (in vitro)`>0.1 | Ae5$`20_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`20_4) Human EBs (Clark1)`>0.1 | Ae5$`17_4) Human EBs (Clark1)`>0.1 | Ae5$`1_4) Human EBs (Clark1)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_DE_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark1 ")
ggsave(filename=paste(saveext,"/Ara_Clark1_DE_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#Now Clark2

#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("2_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("2_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`2_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`2_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_PS_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_PS_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("4_1) Human EBs (in vitro)","9_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("4_4) Human EBs (Clark2)","9_2) Human Amnioids (in vitro)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`4_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`4_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_PGC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_PGC_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#AdMe vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("0_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`0_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`0_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_AME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_AME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#eMes vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()






#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("6_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("6_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`6_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`6_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EME_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EME_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("13_1) Human EBs (in vitro)","5_1) Human EBs (in vitro)","7_1) Human EBs (in vitro)","8_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("13_4) Human EBs (Clark2)","5_4) Human EBs (Clark2)","7_4) Human EBs (Clark2)","8_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`13_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`13_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_LAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_LAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("10_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("10_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`10_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`10_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EAm_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_EAm_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





#PGCLC vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("17_1) Human EBs (in vitro)","20_1) Human EBs (in vitro)"), ident.2 = c("1_1) Human EBs (in vitro)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("17_4) Human EBs (Clark2)","20_4) Human EBs (Clark2)"), ident.2 = c("1_4) Human EBs (Clark2)_ESC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.5) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$`17_1) Human EBs (in vitro)`>0.1 | Ae5$`20_1) Human EBs (in vitro)`>0.1 | Ae5$`1_1) Human EBs (in vitro)_ESC`>0.1 | Ae5$`20_4) Human EBs (Clark2)`>0.1 | Ae5$`17_4) Human EBs (Clark2)`>0.1 | Ae5$`1_4) Human EBs (Clark2)`>0.1)  ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.5))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.5))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_DE_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.5),linetype="dashed") + geom_hline(yintercept=-log(1.5),linetype="dashed")  + geom_vline(xintercept = log(1.5), linetype="dashed") + geom_vline(xintercept = -log(1.5), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " EBs ", y = " Clark2 ")
ggsave(filename=paste(saveext,"/Ara_Clark2_DE_ESC_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()







#PS vs ESC
M1 <- FindMarkers(mammal.combined, ident.1 = c("0_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1 = c("18_3) Gastruloid"), ident.2 = c("1_3) Gastruloid"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Now do the comparison
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M1)[which(abs(M1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(M1)[which(abs(M1$avg_logFC)>log(1.2))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M1),"X"] <- M1$avg_logFC
Ae5[rownames(M2),"Y"] <- M2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M1)[which(M1$p_val_adj<0.1)],rownames(M2)[which(M2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which(Ae5$`0_3) Gastruloid`>0.1 | Ae5$`18_3) Gastruloid`>0.1 | Ae5$`1_3) Gastruloid`>0.1   ) ,]
genes.to.label1 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.2))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))]), TF)
genes.to.label2 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.2))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))]),SIGNAL1)
genes.to.label3 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.2))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))]),SIGNAL2)
genes.to.label4 = intersect(c(rownames(M1)[which(abs(M1$avg_logFC)>log(1.2))],rownames(M2)[which(abs(M2$avg_logFC)>log(1.2))]),SIGNAL3)
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " Gast ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_0_18_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " Gast ", y = " Gast ")
ggsave(filename=paste(saveext,"/Ara_Gast_0_18_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
