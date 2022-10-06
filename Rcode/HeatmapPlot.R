library(plotly)
library(Seurat)
library(stringr)
library(ggplot2)

saveext <- "/Volumes/Overflow/Uber18G/"
fileID <- "/Volumes/Overflow/Uber18G/"
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

mammal.combined <- readRDS(paste(fileID,'mammal.combined.KO5.rds',sep=""))

#PrincessAra
#LoginForMyCuiteAra

Markers <- c("SOX17",
"TFAP2C",
"PRDM14",
"TFCP2L1",
"KLF4",
"CD38",
"PDPN",
"NANOS3",
"UTF1",
"ARID5B",
"FGFR3",
"DAZL",
"DDX4",
"MAEL",
"NANOG",
"SOX2",
"KLF4",
"POU5F1",
"PRDM14",
"MIXL1",
"MESP1",
"TBXT",
"CDH1",
"FST",
"CDX1",
"CDX2",
"MSX1",
"MSGN1",
"SNAI2",
"ALCAM",
"GATA6",
"PDGFRA",
"MEST",
"HAND1",
"EOMES",
"GSC",
"ZIC2",
"ZIC3",
"ZIC5",
"HAND1",
"GATA3",
"TFAP2A",
"ISL1",
"VTCN1",
"ID3",
"BMP2",
"BMP4",
"BAMBI",
"CER1",
"FOXA2",
"SOX17",
"WNT1",
"WNT2",
"WNT2B",
"WNT3",
"WNT3A",
"WNT4",
"WNT5A",
"WNT5B",
"WNT6",
"WNT7A",
"WNT7B",
"WNT8A",
"WNT8B",
"WNT9A",
"WNT9B",
"WNT10A",
"WNT11",
"WNT16",
"TGFB1",
"TGFB2",
"TGFB3",
"NODAL")

#New labels etc.
CellCycle <- mammal.combined$Phase

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

Idents(mammal.combined) <- mammal.combined$species
Data1 <- subset(mammal.combined, idents="Human (Ara EBs)")

Idents(Data1) <- Data1$Cl9
uID <- as.character(Data1$Cl9)
uID[which(Data1$ID2=="ESC")] <- "ESC"
uID[which(Data1$ID2=="PreME")] <- "PreME"
uID[which(Data1$ID2=="4i")] <- "4i"
uID[which(Data1$ID2=="DE")] <- "DE1"
uID[which(Data1$ID2=="DE" & Data1$Cl9==15)] <- "DE"
uID[which(Data1$ID2=="DE" & Data1$Cl9==16)] <- "DE"
uID[which(Data1$ID2=="ME")] <- "ME"
uID[which(Data1$ID2=="PGC")] <- "PGC"
Idents(Data1) <- uID

#ord <- c("EB","PGC_CS5","PS","Am_CS5","EmDisc_CS5","Endoderm","Mesoderm","Am_CS7","EmDisc_CS7","Am_CS6","PGC_CS6","EmDisc_CS6")
DimPlot(Data1,   pt.size = 4,  reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_1t.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,   pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_2t.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,    pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_3t.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



Markers <- c(
  "SOX2",
  "PRDM14",
  "POU5F1",
  "NANOG",
  "TBXT",
  "CDH1",
  "FST",
  "CDX1",
  "CDX2",
  "MSX1",
  "MSGN1",
  "EOMES",
  "MIXL1",
  "MESP1",
  "FOXF1",
  "SNAI2",
  "HAND1",
  "ALCAM",
  "MEST",
  "PDGFRA",
  "GATA4",
  "GATA6",
  "HNF1B",
  "APOA1",
  "APOB",
  "FOXA1",
  "FOXA2",
  "TTR",
  "SOX17",
  "PRDM1",
  "TFAP2C",
  "PDPN",
  "TFCP2L1",
  "KLF4",
  "CD38",
  "NANOS3",
  "UTF1",
  "ARID5B",
  "DAZL",
  "DDX4",
  "MAEL",
  "TFAP2A",
  "GATA3",
  "ISL1",
  "WNT6",
  "KRT18",
  "KRT17",
  "VTCN1",
  "ITGB6",
  "GABRP"
)

AvExp <- AverageExpression(Data1)

library(pheatmap)
#redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
#mat_breaks <- seq(-1, 1, length.out = 60)
#pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks,color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_4_norm",".pdf",sep=""),width=4,height=6)


redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 10, length.out = 60)
pheatmap(log2(100*AvExp$RNA[Markers,c("ESC","4i","PreME","ME","DE","2","6","0","28","17","20","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_row = c(4,14,21,32,35), gaps_col = c(3,5,8,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_4",".pdf",sep=""),width=4,height=6)
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("ESC","4i","PreME","ME","DE","2","6","0","28","17","20","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks, color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,5,8,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_4_norm",".pdf",sep=""),width=4,height=6)

List2 <- read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/SigMarkers.csv")

#DefaultAssay(Data1) <- "RNA"
#g1 <- intersect(rownames(Data1),List2$V1)
g1 <- as.character(List2$V1)

X <- AvExp$RNA[g1,]          

X <- X[rowSums(X[, -1])>0, ]

redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 10, length.out = 60)
pheatmap(log2(100*X[,c("ESC","4i","PreME","ME","DE","2","6","0","28","17","20","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_col = c(3,5,8,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_4_3",".pdf",sep=""),width=4,height=12)
#mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(X[,c("ESC","4i","PreME","ME","DE","2","6","0","28","17","20","9","4","PGC","10","13","8","5","7")], color =  redblue1(60), border_color = NA, gaps_col = c(3,5,8,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_4_norm_3",".pdf",sep=""),width=4,height=12)

   

#Load in data
#Ara ESC, early PGCLC, late PGCLC, PGC, Amnion
#Human cs7, amnnionn, epi, gPGC
#Cyno, amnion, 
#Amnioid: PGC, ESC, amnion
#Gastruloid: PGCC, ESC, amnion
#Load in the datasets and get the intersecting genes that are in the gene models
raw_countsA1<-read.table("/Users/christopherpenfold/Desktop/Aracely/featurecountsLi.csv",sep=",",header = T, row.names=1)
EBGenes<-read.table("/Users/christopherpenfold/Desktop/Aracely/Pairwise_PGCLC_vs_Amnioids/Plots/GenesMapping.csv",sep=",",header = T)
EBGenes<-EBGenes$Gene
raw_counts4 <- readRDS("/Users/christopherpenfold/Downloads/expression_values_1.rds")


raw_countsA2<-read.table("/Volumes/GoogleDrive/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Cyno/cyInVitData.csv",sep=",",header = T, row.names=1)
#colnames(raw_counts4)[which(colnames(raw_counts4)=="TBXT")] <- "T"
#rownames(raw_counts1)[which(rownames(raw_counts1)=="ENSMFAG00000033399")] = "NANOG"

#ENSCJAG00000040521==GATA1

#Look at these that map in all experiments
common_genes1 <- intersect(intersect(EBGenes,rownames(raw_countsA1)),colnames(raw_counts4))
common_genes2 <- intersect(EBGenes,rownames(raw_countsA2))
#common_genes <- intersect(common_genes1,raw_counts4)

uID <- as.character(mammal.combined$Cells)
uID[which(uID=="Ectoderm")] <- "Human_CS7_Ectoderm"
uID[which(uID=="Epiblast")] <- "Human_CS7_Epiblast"
uID[which(mammal.combined$ID3=="PGC_CS7" & Species2=="CS7")] <- "Human_CS7_PGC"

uID[which(mammal.combined$ID3=="ESC" & Species2=="AraEBsPreME")] <- "EB_ESC"
uID[which(mammal.combined$Cl9=="4" & Species2=="AraEBsPreME")] <- "EB_lPGCLC"
uID[which(mammal.combined$Cl9=="9" & Species2=="AraEBsPreME")] <- "EB_ePGCLC"
uID[which(mammal.combined$Cl9=="13" & Species2=="AraEBsPreME")] <- "EB_AMLC"
uID[which(mammal.combined$Cl9=="8" & Species2=="AraEBsPreME")] <- "EB_AMLC"
uID[which(mammal.combined$Cl9=="5" & Species2=="AraEBsPreME")] <- "EB_AMLC"
uID[which(mammal.combined$Cl9=="7" & Species2=="AraEBsPreME")] <- "EB_AMLC"
uID[which(mammal.combined$Cl9=="4" & mammal.combined$Cells=="PGC_wk6.5_invivo")] <- "EB_PGC"


uID[which(mammal.combined$Cl9=="1" & Species2=="Cyno")] <- "Cyno_Epiblast"
uID[which(mammal.combined$Cl9=="21" & Species2=="Cyno")] <- "Cyno_Epiblast"
uID[which(mammal.combined$Cl9=="3" & Species2=="Cyno")] <- "Cyno_Epiblast"
uID[which(mammal.combined$Cl9=="19" & Species2=="Cyno")] <- "Cyno_Epiblast"
uID[which(mammal.combined$Cl9=="14" & Species2=="Cyno")] <- "Cyno_Epiblast"

uID[which(uID=="Am_CS5" & Species2=="Cyno")] <- "Cyno_Amnion"
uID[which(uID=="Am_CS6" & Species2=="Cyno")] <- "Cyno_Amnion"
uID[which(uID=="Am_CS7" & Species2=="Cyno")] <- "Cyno_Amnion"
uID[which(uID=="PGC_CS6" & Species2=="Cyno")] <- "Cyno_PGC"
uID[which(uID=="PGC" & Species2=="Human1" & mammal.combined$Cl9=="4")] <- "Human_gPGC"

Idents(mammal.combined) <- uID
AvE <-AverageExpression(mammal.combined)

X1<- cor( log2(AvE$RNA[common_genes1,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1), log2(AvE$RNA[common_genes1,c('Human_CS7_Epiblast','Human_CS7_Ectoderm','Human_CS7_PGC','Human_gPGC')]+1))
X2 <- cor( log2(AvE$RNA[common_genes2,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1),log2(AvE$RNA[common_genes2,c('Cyno_Epiblast','Cyno_Amnion','Cyno_PGC')]+1))        
                  
library(pheatmap)
redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 10, length.out = 60)
pheatmap(X1, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC1_RNA",".pdf",sep=""),width=4,height=4)
pheatmap(X2, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC2_RNA",".pdf",sep=""),width=4,height=4)



X1<- cor( log2(AvE$RNA[,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1), log2(AvE$RNA[,c('Human_CS7_Epiblast','Human_CS7_Ectoderm','Human_CS7_PGC','Human_gPGC')]+1))
X2 <- cor( log2(AvE$RNA[,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1),log2(AvE$RNA[,c('Cyno_Epiblast','Cyno_Amnion','Cyno_PGC')]+1))        

library(pheatmap)
redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 10, length.out = 60)
pheatmap(X1, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC1_RNAb",".pdf",sep=""),width=4,height=4)
pheatmap(X2, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC2_RNAb",".pdf",sep=""),width=4,height=4)


X1<- cor( log2(AvE$integrated[,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1), log2(AvE$integrated[,c('Human_CS7_Epiblast','Human_CS7_Ectoderm','Human_CS7_PGC','Human_gPGC')]+1))
X2 <- cor( log2(AvE$integrated[,c('EB_ESC','EB_AMLC','EB_ePGCLC','EB_lPGCLC','EB_PGC')]+1),log2(AvE$integrated[,c('Cyno_Epiblast','Cyno_Amnion','Cyno_PGC')]+1))        

library(pheatmap)
redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 10, length.out = 60)
pheatmap(X1, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC1_Int",".pdf",sep=""),width=4,height=4)
pheatmap(X2, color =  redblue1(60), border_color = NA,  cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste("~/Desktop/CrossCorr_PGC2_Int",".pdf",sep=""),width=4,height=4)


#Idents(mammal.combined) <- mammal.combined$Cl9

DimPlot(mammal.combined,   pt.size = 2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste("~/Desktop/test3.pdf",sep=""),width = 15*10, height = 15, useDingbats=FALSE,limitsize = FALSE)


#Compare:
#EB_PGC, 
#Human_PGC

#uID[which(Species2=="Cyno")] <- paste("Cyno",Idents(mammal.combined)[which(Species2=="Cyno")],sep = "_")


           