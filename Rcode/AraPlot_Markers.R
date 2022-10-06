library(plotly)
library(Seurat)
library(stringr)

saveext <- "/Volumes/Overflow/Uber18G/"
fileID <- "/Volumes/Overflow/Uber18G/"
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

mammal.combined <- readRDS(paste(fileID,'mammal.combined.KO5.rds',sep=""))

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



mammal.combined$ID2 <- as.factor(labs)
Idents(mammal.combined) <- as.factor(labs)
#Load in IDs for gastruloids
saveext2 <- "/Volumes/Overflow/Gast_2D/"
DE <- readRDS(file=paste(saveext2, "DE.rds",sep=""))
PGCLC <- readRDS(file=paste(saveext2, "PGCLC.rds",sep=""))
MELC1 <- readRDS(file=paste(saveext2, "MELC1.rds",sep=""))
MELC2 <- readRDS(file=paste(saveext2, "MELC2.rds",sep=""))
EPI <- readRDS(file=paste(saveext2, "EPI.rds",sep=""))
Ect <- readRDS(file=paste(saveext2, "Ect.rds",sep=""))
Amn <- readRDS(file=paste(saveext2, "Am.rds",sep=""))
DE <- str_replace(DE,'-1_1', '-1_13')
DE <- str_replace(DE,'-1_2', '-1_14')
PGCLC <- str_replace(PGCLC,'-1_1', '-1_13')
PGCLC <- str_replace(PGCLC,'-1_2', '-1_14')
MELC1 <- str_replace(MELC1,'-1_1', '-1_13')
MELC1 <- str_replace(MELC1,'-1_2', '-1_14')
MELC2 <- str_replace(MELC2,'-1_1', '-1_13')
MELC2 <- str_replace(MELC2,'-1_2', '-1_14')
EPI <- str_replace(EPI,'-1_1', '-1_13')
EPI <- str_replace(EPI,'-1_2', '-1_14')
Ect<- str_replace(Ect,'-1_1', '-1_13')
Ect <- str_replace(Ect,'-1_2', '-1_14')
Amn<- str_replace(Amn,'-1_1', '-1_13')
Amn <- str_replace(Amn,'-1_2', '-1_14')

Idents(mammal.combined, cells=DE) <- "DE"
Idents(mammal.combined, cells=PGCLC) <- "PGCLC"
Idents(mammal.combined, cells=MELC1) <- "MELC1A"
Idents(mammal.combined, cells=MELC2) <- "MELC1B"
Idents(mammal.combined, cells=EPI) <- "EPI"
Idents(mammal.combined, cells=Amn) <- "AMLC"
Idents(mammal.combined, cells=Ect) <- "ECLC"

#
DC7 <- readRDS("/Volumes/GoogleDrive/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Human/Human_CS7.rds")
Labs <- readRDS("/Volumes/GoogleDrive/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Human/Human_CS7_annot_umap.rds")
Hum4 <- CreateSeuratObject(counts = as.data.frame(t(DC7)), assay = "RNA",min.cells = 0, min.features = 0)
Idents(Hum4) <- Labs$cluster_id
Hum4 <- subset(Hum4, idents = c("Hemogenic Endothelial Progenitors","Erythrocytes"), invert = TRUE)
Hum4 <- NormalizeData(Hum4, verbose = FALSE)
Hum4$species <- "CS7"
Hum4 <- FindVariableFeatures(Hum4, selection.method = "vst", nfeatures = 3000)
PGC <- intersect(WhichCells(Hum4, expression = NANOS3 > .1),WhichCells(Hum4, expression = SOX17 > .1))
Amn <- intersect(WhichCells(Hum4, expression = VTCN1 > .1),WhichCells(Hum4, expression = TFAP2A > .1))

Idents(mammal.combined,cells=paste(PGC,"_9",sep="") ) <- "PGC_CS7"
Idents(mammal.combined,cells=paste(Amn,"_9",sep="") ) <- "Am_CS7"
mammal.combined$ID3 <- Idents(mammal.combined)
#Rejig IDs and order
Species2 <- as.character(mammal.combined$species)
Species2[which(Species2=="Human (Ara EBs)")] <- "EBPreMe"
Species2[which(Species2=="Human (Amnioid)")] <- "Amnioids"
Species2[which(Species2=="HumanCS7")] <- "HumanCS7"
Species2[which(Species2=="HumanCS5-7")] <- "10) Human (in vitro 1)"
Species2[which(Species2=="human_invivo")] <- "11) Human (in vitro 2)"
Species2[which(Species2=="Cyno (in vitro)")] <- "8) Cynomolgus (in vitro)"
Species2[which(Species2=="3) Marmoset")] <- "9) Marmoset (in vivo)"
Species2[which(Species2=="10X_1")] <- "EBs (Clark1)"
Species2[which(Species2=="10X_2")] <- "EBs (Clark1)"
Species2[which(Species2=="10X_3")] <- "EBs (Clark2)"
Species2[which(Species2=="10X_4")] <- "EBs (Clark2)"
Species2[which(Species2=="10X_5")] <- "EBs (Clark3)"
Species2[which(Species2=="Gastruloid2D")] <- "Gastruloid"

Species2[which(mammal.combined$Cells=="PGCLC_D4_r2")] <- "EB4i"

mammal.combined$species2 <- Species2

Idents(mammal.combined) <- mammal.combined$species
Data1 <- subset(mammal.combined, idents=c("Human (Ara EBs)","Cyno (in vitro)"))
Data2 <- subset(mammal.combined, idents=c("Human (Ara EBs)","HumanCS7"))
Data3 <- subset(mammal.combined, idents=c("Human (Ara EBs)","HumanCS5-7"))
Data4 <- subset(mammal.combined, idents=c("Human (Ara EBs)","Human (Amnioid)"))
Data5 <- subset(mammal.combined, idents=c("Human (Ara EBs)","10) Human (in vitro 1)"))
Data6 <- subset(mammal.combined, idents=c("Human (Ara EBs)","Gastruloid2D"))
#Data1 <- subset(mammal.combined, idents=c("1) Human EBs (in vitro)","Cyno (in vitro)"))
#Data1 <- subset(mammal.combined, idents=c("1) Human EBs (in vitro)","Cyno (in vitro)"))
Data7 <- subset(mammal.combined, idents=c("Human (Ara EBs)","3) Marmoset"))

############## First dooooooo amnioids
Idents(Data4) <- Data4$Cl5
DimPlot(Data4,   pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/cls.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

uID <- as.character(Data4$ID2)
uID[which(Data4$species=="Human (Ara EBs)")] <- "EB"
uID[which(Data4$Cl5=="14")] <- "DELC"
uID[which(Data4$Cl5=="12")] <- "DELC"
Idents(Data4) <- uID

MarmType <- c("ESC","PGCLC","MELC1","MELC2","AMLC","DELC","EB")
MarmC <- c("#F49739","#E166A3","#419CD6","#1B84C7","#00A652","#5D5D5D","#D3D3D3")

D2<- Idents(Data4)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]

saveextb <- '/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/PlotRedo/Amnioid/'
DimPlot(Data4,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_AMLC_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data4,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset6_AMLC_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data4,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_AMLC_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


DefaultAssay(Data4) <- "RNA"


Markers <- c("NANOS3","SOX17","PRDM1","CD38","TFAP2A","TFAP2C","PDPN","PRDM14","ARID5B",
"SNAI2", "EOMES", "FOXF1", "FOXH1", "OTX2", "MIXL1", "MEST", "T", "ALCAM", "HAND1", "GSC", "CDX2", "TBX6", "ZIC2", "ZIC3", "ZIC5", "GATA4","GATA6",
"DAZL", "DDX4", "DMRT1", "DND1", "MAEL", "SYCP3", "UTF1",
"VTCN1", "GABRP", "ISL1", "ITGB6", "KRT7", "KRT17", "GATA3", "GATA2",
"SOX2", "POU5F1", "NANOG", "SOX1", "KLF2", "KLF4", "TFCP2L1", 
"FOXA1", "FOXA2", "PDGFRA", "APOA1", "APOB",
"BMP2", "BMP4", "BMP8", "FGF2", "FGF4", "FGF5", "FGF10", "NODAL", "TGFB1", "TGFB2", "TGFB3",
"WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10A", "WNT11", "WNT16",
"ID1", "ID2", "ID3",
"FST",
"CER1", "NOG", "CHRD", "GREM1", "BAMBI",
"DKK1", "AXIN", "AXIN2", "LEF", "WIF",
"LEFTY1", "LEFTY2", "CER1",
"TGFBR", "BMPRIA", "BMPRIB", "BMPRII", "ACVR2A", "ACVR2B", "ACVR1")

Idents(Data4) <- Data4$species
Data4 <- subset(Data4, idents=c("Human (Amnioid)"))
Markers1 <- intersect(rownames(Data4),Markers)
for (i in 1:length(Markers1)) {
  #possibleError <-  tryCatch(
  FeaturePlot(Data4, features = as.character(Markers1[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))#, error=function(e) e)
  #if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveextb,"/","AmnioidMarkers", "_", as.character(Markers1[i]),".pdf",sep=""),width = 15, height = 15,limitsize = FALSE, useDingbats=FALSE)
  #graphics.off() 
}


#Idents(Data1)
#DimPlot(Data1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


###### Now do marmoset

saveextb <- '/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/PlotRedo/Marmoset/'

uID <- as.character(Data7$ID2)
uID[which(Data7$species=="Human (Ara EBs)")] <- "EB"
Idents(Data7) <- uID

#testMarm <- readRDS('/Users/christopherpenfold/Desktop/Thorsten/FINAL/AllPlatesCCA_4_dataset/mammal.combiined.rds')
#DimPlot(testMarm,  pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) #+ xlim(c(-13, 13)) + ylim(c(-13, 13))
#ggsave(filename=paste(saveextb,"/TestMarm","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
#Idents(testMarm) <- testMarm$Cl1
#DimPlot(testMarm,  pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) #+ xlim(c(-13, 13)) + ylim(c(-13, 13))
#ggsave(filename=paste(saveextb,"/TestMarm","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

#Endo <- intersect(paste(WhichCells(testMarm,idents="2"),"1",sep="_"),colnames(Data7))
#Mes <- intersect( paste(WhichCells(testMarm,idents=c("8","3","6","9") ),"1",sep="_"),colnames(Data7))
#PS1 <- intersect( paste(WhichCells(testMarm,idents="0"),"1",sep="_"),colnames(Data7))
#PS2 <- intersect( paste(WhichCells(testMarm,idents="7"),"1",sep="_"),colnames(Data7))

#Idents(Data7,cells=Endo) <- "EmDisc_Endo_CS7"
#Idents(Data7,cells=Mes) <- "EmDisc_Meso"
#Idents(Data7,cells=PS1) <- "EmDisc_PS1"
#Idents(Data7,cells=PS2) <- "EmDisc_PS2"
#uID <- Idents(Data7)
#uID[which(Data7$ID2=="Am_CS7")] <- "Am_CS7"
#uID[which(Data7$ID2=="Am_CS6")] <- "Am_CS6"
#uID[which(Data7$ID2=="PGC_CS6")] <- "PGC_CS6"
#Idents(Data7) <- uID

#Human 2
MarmType <- c("EmDisc_PS1","EmDisc_PS2","EmDisc_Endo_CS7","EmDisc_Meso","EB","Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","EPI","PGC","PE")
MarmC <- c("black","black","#5D5D5D","#1B84C7","lightgrey","#F27A40","#BF470D","#C81048","#8D340A","#C81048","#C81048","#54DE62","#10C357","#07971F","#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434","#F0621F","#C81048","#5C5C5C")

D2<- Idents(Data7)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data7,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_marm_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data7,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset6_marm_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data7,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_marm_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


DefaultAssay(Data7) <- "RNA"
Idents(Data7) <- Data7$species
Data7 <- subset(Data7, idents=c("3) Marmoset"))
Markers1 <- intersect(rownames(Data7),Markers)
for (i in 71:length(Markers1)) {
  #possibleError <-  tryCatch(
  FeaturePlot(Data7, features = as.character(Markers1[i]), pt.size = 5, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))#, error=function(e) e)
  #if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveextb,"/","MarmosetMarkers", "_", as.character(Markers1[i]),".pdf",sep=""),width = 15, height = 15,limitsize = FALSE, useDingbats=FALSE)
  #graphics.off() 
}

###### Now do cyno...

saveextb <- '/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/PlotRedo/Cyno/'
uID <- as.character(Data1$ID2)
uID[which(Data1$species=="Human (Ara EBs)")] <- "EB"
Idents(Data1) <- uID

#Human 2

MarmType <- c("EB","Epi_CS3","EmDisc_CS5","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","EPI","PGC","PE","PS","Endoderm","Amnion","Mesoderm")
MarmC <- c("lightgrey","#F27A40","#EF5910","#BF470D","#C81048","#8D340A","#C81048","#C81048","#54DE62","#10C357","#07971F","#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434","#F0621F","#C81048","#5C5C5C","#afdfd6","#edbae5","#0E672F","#419CD6")




D2<- Idents(Data1)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_cyno_","_1_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset6_cyno_","_2_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/UMAP_Dataset_cyno_","_3_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


DefaultAssay(Data1) <- "RNA"
Idents(Data1) <- Data1$species
Data1 <- subset(Data1, idents=c("Cyno (in vitro)"))
Markers1 <- intersect(rownames(Data7),Markers)
for (i in 1:length(Markers1)) {
  #possibleError <-  tryCatch(
  FeaturePlot(Data1, features = as.character(Markers1[i]), pt.size = 5, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))#, error=function(e) e)
  #if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveextb,"/","CynoMarkers", "_", as.character(Markers1[i]),".pdf",sep=""),width = 15, height = 15,limitsize = FALSE, useDingbats=FALSE)
  #graphics.off() 
}

###### Now do human 1?
saveextb <- '/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/PlotRedo/CS5-7/'

saveext1 = "~/Desktop/Thorsten/FINAL/Human_amnion_EmDisc3/"
Am<-readRDS(paste(saveext1,"/CynoAm",".rds",sep=""))
Endo<-readRDS(paste(saveext1,"/HumanEndo",".rds",sep=""))
Meso<-readRDS(paste(saveext1,"/HumanMes",".rds",sep=""))
PS<-readRDS(paste(saveext1,"/HumanPS",".rds",sep=""))

#saveRDS(CynoEndo,paste(saveext,"/CynoEndo",".rds",sep=""))
#saveRDS(MarmEndo,paste(saveext,"/MarmEndo",".rds",sep=""))
#saveRDS(CynoMes,paste(saveext,"/CynoMes",".rds",sep=""))
#saveRDS(MarmMes,paste(saveext,"/MarmMes",".rds",sep=""))

#saveRDS(CynoPS,paste(saveext,"/CynoPS",".rds",sep=""))
#saveRDS(MarmPS,paste(saveext,"/MarmPS",".rds",sep=""))

uID <- as.character(Data3$ID2)
uID[which(Data3$species=="Human (Ara EBs)")] <- "EB"
Idents(Data3) <- uID
Idents(object = Data3, cells =paste(Am,"_10",sep="") ) <- "Amnion"
#Idents(object = Data5, cells =paste(Meso,"_10",sep="")) <- "Mesoderm"
#Idents(object = Data5, cells =paste(Endo,"_10",sep="")) <- "Endoderm"
#Idents(object = Data5, cells =paste(PS,"_10",sep="") ) <- "PS"
#Human 2
MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC","PS","Endoderm","Amnion","Mesoderm")
MarmC <- c("lightgrey","#EF5910","#C81048","#8D340A","black","#afdfd6","#edbae5","#0E672F","#419CD6")



D2<- Idents(Data3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/Dataset_goodhuman_","_1_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/Dataset6_goodhuman_","_2_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveextb,"/Dataset_goodhuman_","_3_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



###### Now do human 3
uID <- as.character(Data3$ID2)
uID[which(Data3$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data3) <- uID

#Human 2
MarmType <- c("EB","EPI","EmDisc_CS6","EmDisc_CS7","PGC")
MarmC <- c("lightgrey","#EF5910","#C81048","#8D340A","black")

D2<- Idents(Data3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_badhuman_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_badhuman_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_badhuman_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


######## Now do CS7
uID <- as.character(Data2$ID2)
uID[which(Data2$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data2) <- uID
Idents(Data2,cells=paste(PGC,"_9",sep="") ) <- "PGC"
Idents(Data2,cells=paste(Amn,"_9",sep="") ) <- "Am"

MarmType <- c("ESC","Primitive Streak","Endoderm","PGC","Ectoderm","Am","Emergent Mesoderm","Nascent Mesoderm","Advanced Mesoderm","Axial Mesoderm","YS Mesoderm","Epiblast","EB")
MarmC <- c("#E96F08","#E73D19","#1B1B1B","#E72185","#0C9542","#0E672F","#419CD6","#1B84C7","#0A67A0","#3293D1","#06476E","#107F87","#D3D3D3")

D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_CS7_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_CS7_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_CS7_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


uID <- as.character(Data3$ID2)
uID[which(Data3$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data3) <- uID

#Human 2
MarmType <- c("EB","Epi_CS3","EmDisc_CS5","EmDisc_CS6","VE_CS4","VE_CS5","EPI","PGC","PE")
MarmC <- c("lightgrey","#F16E2F","#DE520F","#C84A0D","#686868","#828282","#F0621F","#C81048","#5C5C5C")
  
D2<- Idents(Data3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_human2_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_human2_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_human2_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



####Now the gastruloid
uID <- as.character(Data6$ID2)
uID[which(Data6$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data6) <- uID

saveext2 <- "/Volumes/Overflow/Gast_2D/"
DE <- readRDS(file=paste(saveext2, "DE.rds",sep=""))
PGCLC <- readRDS(file=paste(saveext2, "PGCLC.rds",sep=""))
MELC1 <- readRDS(file=paste(saveext2, "MELC1.rds",sep=""))
MELC2 <- readRDS(file=paste(saveext2, "MELC2.rds",sep=""))
EPI <- readRDS(file=paste(saveext2, "EPI.rds",sep=""))
Ect <- readRDS(file=paste(saveext2, "Ect.rds",sep=""))
Amn <- readRDS(file=paste(saveext2, "Am.rds",sep=""))


DE <- str_replace(DE,'-1_1', '-1_13')
DE <- str_replace(DE,'-1_2', '-1_14')
PGCLC <- str_replace(PGCLC,'-1_1', '-1_13')
PGCLC <- str_replace(PGCLC,'-1_2', '-1_14')
MELC1 <- str_replace(MELC1,'-1_1', '-1_13')
MELC1 <- str_replace(MELC1,'-1_2', '-1_14')
MELC2 <- str_replace(MELC2,'-1_1', '-1_13')
MELC2 <- str_replace(MELC2,'-1_2', '-1_14')
EPI <- str_replace(EPI,'-1_1', '-1_13')
EPI <- str_replace(EPI,'-1_2', '-1_14')
Ect<- str_replace(Ect,'-1_1', '-1_13')
Ect <- str_replace(Ect,'-1_2', '-1_14')
Amn<- str_replace(Amn,'-1_1', '-1_13')
Amn <- str_replace(Amn,'-1_2', '-1_14')

Idents(Data6, cells=DE) <- "DELC"
Idents(Data6, cells=PGCLC) <- "PGCLC"
Idents(Data6, cells=MELC1) <- "MELC2"
Idents(Data6, cells=MELC2) <- "MELC2"
Idents(Data6, cells=EPI) <- "EPI"
Idents(Data6, cells=Amn) <- "AMLC"
Idents(Data6, cells=Ect) <- "ECLC"
Data6$ID2 <- Idents(Data6)

#Human 2
MarmType <- c("EPI","PGCLC","ECLC","MELC2","AMLC","DELC","EB")
MarmC <- c("#F49739","#E166A3","#419CD6","#1B84C7","#00A652","#5D5D5D","#D3D3D3")

D2<- Idents(Data6)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data6,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_gastruloid_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data6,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_gastruloid_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data6,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_gastruloid_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



D2<- Idents(Data3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_human2_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_human2_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_human2_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

#Clark1
Data1 <- subset(mammal.combined, idents=c("4) Human EBs (Clark1)"))
MarmType <- c("4i","ESC","PreME","ME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_72h","PGCLC_96h","PGC","DE")
MarmC <- c("#25dabf","#20a4c2","#2581b5","#19137a","#f8e8c9","#f5d59e","#f6bb83","#fa8c5a","#ed6648","#d6301f","#B82213","#991306","#4B882C","#692c88")
Idents(Data1) <- Data1$ID2
D2<- Idents(Data1)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1,  cols = coluse1, pt.size = 1, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark1_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,  cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_clark1_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,   cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark1_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

#Clark2
Data2 <- subset(mammal.combined, idents=c("4) Human EBs (Clark2)"))
MarmType <- c("4i","ESC","PreME","ME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_72h","PGCLC_96h","PGC","DE")
MarmC <- c("#25dabf","#20a4c2","#2581b5","#19137a","#f8e8c9","#f5d59e","#f6bb83","#fa8c5a","#ed6648","#d6301f","#B82213","#991306","#4B882C","#692c88")
Idents(Data2) <- Data2$ID2
D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 1, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark2_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,  cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_clark2_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,   cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark2_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



#Clark3
#Data1 <- subset(mammal.combined, idents=c("8) Human EBs (Clark)"))
MarmType <- c("4i","ESC","PreME","ME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_72h","PGCLC_96h","PGC","DE")
MarmC <- c("#25dabf","#20a4c2","#2581b5","#19137a","#f8e8c9","#f5d59e","#f6bb83","#fa8c5a","#ed6648","#d6301f","#B82213","#991306","#4B882C","#692c88")
Idents(Data8c) <- Data8c$ID2
D2<- Idents(Data8c)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data8c,  cols = coluse1, pt.size = 1, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark3_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data8c,  cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_clark3_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data8c,   cols = coluse1, pt.size = 1, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_clark3_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)





Markers <- c(
  "POU5F1",
  "NANOG",
  "SOX2",
  "SOX3",
  "PRDM14",
  "LIN28",
  "ESRRB",
  "KLF2",
  "KLF4",
  "TFCP2L1",
  "CDH1",
  "FST",
  "KDR",
  "EOMES",
  "T",
  "GSC",
  "MIXL1",
  "LHX1",
  "OTX2",
  "LEFTY2",
  "MESP1",
  "FOXF1",
  "SNAI2",
  "ALCAM",
  "HAND1",
  "MEST",
  "PDGFRA",
  "GATA4",
  "GATA6",
  "HNF1B",
  "APOA1",
  "APOB",
  "FOXA1",
  "TTR",
  "FOXA2",
  "HNF1A",
  "HNF4A",
  "SOX17",
  "PRDM1",
  "TFAP2C",
  "PRDM14",
  "TFCP2L1",
  "CD38",
  "PDPN ",
  "NANOS3",
  "UTF1",
  "ARID5B",
  "FGFR3",
  "DND1",
  "DPPA3",
  "DAZL",
  "DDX4",
  "MAEL",
  "SYCP3",
  "RNF17",
  "TDRD5",
  "TDRD9",
  "TDRD10",
  "PIWIL1",
  "PIWIL2",
  "PIWIL3",
  "PRAME",
  "DND1",
  "DMRT1",
  "DNMT3B",
  "ACE",
  "ADAD1",
  "DDX4",
  "M1AP",
  "SOHLH2",
  "SOX30",
  "SPATA22",
  "TDRD5",
  "ZNF597",
  "ZNF99",
  "ZNF208",
  "ZNF354C",
  "ZNF534",
  "ZNF560",
  "ZNF662",
  "ZNF726",
  "ZNF728",
  "ZNF729",
  "TFAP2A",
  "GATA3",
  "VTCN1",
  "WNT6",
  "ISL1",
  "ITGB6",
  "GABRP",
  "MUC16",
  "KRT17",
  "KRT18",
  "HLA-1",
  "HLA-G",
  "CGA",
  "TP63",
  "KRT7",
  "CK7",
  "HAVCR1",
  "SLC28A3",
  "DAB1",
  "ADAP2",
  "XRCC6P2",
  "HNRNPA2B1",
  "PEBP1P2",
  "TFDP1",
  "ALYREF",
  "U2AF1",
  "HNRNPA3P6",
  "SRSF2",
  "XRCC6",
  "TPI1P1",
  "CGA",
  "TMEM258",
  "CGB5",
  "NUCB2",
  "FKBP2",
  "PLPP5",
  "CGB3",
  "SPCS2P4",
  "AC011472.3",
  "CLTB",
  "TIMP2",
  "SKIL",
  "CALU",
  "ESAM",
  "TRIP4",
  "CALR",
  "TMED3",
  "HTRA4",
  "TSPAN9",
  "PAX6",
  "SOX1",
  "TP73L",
  "NES",
  "VIM",
  "SOX2",
  "FGF5",
  "TUJ1",
  "CK18",
  "CK8",
  "FOXJ3",
  "PAX2",
  "OTX2",
  "GBX2",
  "ASCL1",
  "ASCL2")


DefaultAssay(mammal.combined) <- "RNA"

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(
    FeaturePlot(mammal.combined, features = as.character(Markers[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_", as.character(Markers[i]),"_nosplit.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}

Idents(mammal.combined) <- Species2
Data0 <- subset(mammal.combined, idents=c("AraEBs (in vitro)"))
Idents(Data0) <- Data0$ID2
Data0B <- subset(Data0,idents=c("PreME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_96h"))
species3 <- as.character(Idents(Data0B))
species3[which(Data0B$Cells=="PGCLC_D4_r2")] <- "PGCLC_96h_4i"
Data0B$species3 <- factor(species3)
species3 <- factor(Data0B$species3, levels = c("PreME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_96h","PGCLC_96h_4i"))
Data0B$species3 <- species3

DefaultAssay(Data0B) <- "RNA"

for (i in 1:length(markers$V1)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data0B, features = as.character(markers$V1[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/",as.character(markers$V2[i]), "_", as.character(markers$V1[i]),"_Ara.pdf",sep=""),width = 135, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
    possibleError <-  tryCatch(
    FeaturePlot(Data0B, features = as.character(Markers[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
    
    if(inherits(possibleError, "error")) next
    
    ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_Ara.pdf",sep=""),width = 135, height = 15,limitsize = FALSE)
    graphics.off() 
}




Data2 <- subset(mammal.combined, idents=c("EBs (Clark1)"))
Idents(Data2) <- Data2$ID2
Data2B <- subset(Data2,idents=c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
species3 <- factor(Idents(Data2B), levels = c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
Data2B$species3 <- species3
markers <- read.table(file = "/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/Markers.csv",sep=",")
DefaultAssay(Data2B) <- "RNA"


for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data2B, features = as.character(Markers[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_ClarkHC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(markers$V1)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data2B, features = as.character(markers$V1[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/",as.character(markers$V2[i]), "_", as.character(markers$V1[i]),"_ClarkHC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}



Idents(mammal.combined) <- Species2
Data1 <- subset(mammal.combined, idents=c("EBs (Clark2)"))
Idents(Data1) <- Data1$ID2
Data1B <- subset(Data1,idents=c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
species3 <- factor(Idents(Data1B), levels = c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
Data1B$species3 <- species3
markers <- read.table(file = "/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/Markers.csv",sep=",")
DefaultAssay(Data1B) <- "RNA"

for (i in 1:length(markers$V1)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1B, features = as.character(markers$V1[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/",as.character(markers$V2[i]), "_", as.character(markers$V1[i]),"_ClarkLC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1B, features = as.character(Markers[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_ClarkLC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}




Idents(mammal.combined) <- Species2
Data1 <- subset(mammal.combined, idents=c("EBs (Clark3)"))
Idents(Data1) <- Data1$ID2
Data1B <- subset(Data1,idents=c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
species3 <- factor(Idents(Data1B), levels = c("PreME","PGCLC_24h","PGCLC_48h","PGCLC_72h","PGCLC_96h"))
Data1B$species3 <- species3
markers <- read.table(file = "/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/PaperFigures/Markers.csv",sep=",")
DefaultAssay(Data1B) <- "RNA"

for (i in 1:length(markers$V1)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1B, features = as.character(markers$V1[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/",as.character(markers$V2[i]), "_", as.character(markers$V1[i]),"_ClarkLC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1B, features = as.character(Markers[i]), split.by = "species3" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_ClarkLC.pdf",sep=""),width = 75, height = 15,limitsize = FALSE)
  graphics.off() 
}






Species3 <- as.character(Species2)
Species3[1:length(Species3)] <- "Wk4+" 
ind1 <- grep("CS5",as.character(mammal.combined$ID2))
ind2 <- grep("CS6",as.character(mammal.combined$ID2))
ind3 <- grep("CS7",as.character(mammal.combined$ID2))
Species3[ind1] <- "CS5"
Species3[ind2] <- "CS6"
Species3[ind3] <- "CS7"
mammal.combined$Species3 <- Species3

Idents(mammal.combined) <- Species2
Data2 <- subset(mammal.combined, idents="2) Human Amnioids (in vitro)")
Data3 <- subset(mammal.combined, idents="7) Human CS7 (in vivo)")
Data4 <- subset(mammal.combined, idents="10) Human (in vitro 1)")
Data5 <- subset(mammal.combined, idents="11) Human (in vitro 2)")
Data6 <- subset(mammal.combined, idents="8) Cynomolgus (in vitro)")
Data7 <- subset(mammal.combined, idents="9) Marmoset (in vivo)")
Data9 <- subset(mammal.combined, idents="3) Gastruloid")

DefaultAssay(Data2) <- "RNA"
DefaultAssay(Data3) <- "RNA"
DefaultAssay(Data4) <- "RNA"
DefaultAssay(Data5) <- "RNA"
DefaultAssay(Data6) <- "RNA"
DefaultAssay(Data7) <- "RNA"
DefaultAssay(Data9) <- "RNA"






for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data2, features = as.character(Markers[i]) , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_Amnioid.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data3, features = as.character(Markers[i]),  pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_CS7.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data9, features = as.character(Markers[i]),  pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_gastruloid.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}




for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data4, features = as.character(Markers[i]), split.by ="Species3", pt.size = 8, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_goodHumanSplit.pdf",sep=""),width = 60, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data5, features = as.character(Markers[i]),  split.by ="Species3", pt.size = 4, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_badHumanSplit.pdf",sep=""),width = 60, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data6, features = as.character(Markers[i]),  split.by ="Species3", pt.size = 8, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_CynoSplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)
  graphics.off() 
}

for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data7, features = as.character(Markers[i]),  split.by ="Species3", pt.size = 8, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_MarmosetSplit.pdf",sep=""),width = 45, height = 15,limitsize = FALSE)
  graphics.off() 
}





Markers <- c(
  "IHH",
  "SHH",
  "DHH",
  "PTCH1",
  "PTCH2",
  "GLI1",
  "GLI2",
  "GLI3",
  "GPC4",
  "ARL13B",
  "SMO",
  "TULP3",
  "FABP55",
  "YBX1",
  "NR2F2",
  "NRF3",
  "PAGE4",
  "JAM3",
  "UBE2S",
  "NFE2L1"
)
DefaultAssay(Data2) <- "RNA"
DefaultAssay(Data9) <- "RNA"


for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data9, features = as.character(Markers[i]),  pt.size = 4, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_2DGastruloid.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}


for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data2, features = as.character(Markers[i]),  pt.size = 4, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_EB.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}


Data1 <- subset(mammal.combined, idents="1) Human EBs (in vitro)")
DefaultAssay(Data1) <- "RNA"
for (i in 1:length(Markers)) {
  possibleError <-  tryCatch(FeaturePlot(Data1, features = as.character(Markers[i]),  pt.size = 4, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  if(inherits(possibleError, "error")) next
  ggsave(filename=paste(saveext,"/Markers/","Markers", "_", as.character(Markers[i]),"_EBs.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}



Phase <- readRDS("/Volumes/Overflow/Uber18G/DM/Phase.rds")


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#
Data1 <- subset(mammal.combined, idents="PGCLC_12h")
DefaultAssay(Data1) <- "RNA"
Data1 <- CellCycleScoring(Data1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Data1 <- FindVariableFeatures(Data1, selection.method = "vst", nfeatures = 20000)
#Data1 <- ScaleData(Data1, verbose = FALSE)
Data1 <- ScaleData(Data1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Data1))
Data1 <- RunPCA(Data1, npcs = 20, verbose = FALSE)
Data1 <- RunUMAP(Data1, reduction = "pca", dims = 1:20)
Data1 <- RunTSNE(Data1, reduction = "pca", dims = 1:20)
Data1 <- FindNeighbors(Data1, reduction = "pca", dims = 1:20)

Data1 <- FindClusters(Data1, resolution = 0.1)
Data1$Cl1 <- Idents(Data1)
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Cl1_RO.pdf",sep=""),width = 8, height = 8)

Data1 <- FindClusters(Data1, resolution = 0.3)
Data1$Cl3 <- Idents(Data1)
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Cl3_RO.pdf",sep=""),width = 8, height = 8)

Data1 <- FindClusters(Data1, resolution = 0.5)
Data1$Cl5 <- Idents(Data1)
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Cl5_RO.pdf",sep=""),width = 8, height = 8)

Idents(Data1) <- Data1$Cl3
markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl3_RO",".csv",sep=""))
Idents(Data1) <- Data1$Cl5
markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl5_RO",".csv",sep=""))




Idents(Data1) <- Data1$Phase
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Phase_RO.pdf",sep=""),width = 8, height = 8)




FeaturePlot(Data1, features = c("MIXL1"), reduction = "umap") 
ggsave(filename=paste(saveext,"/UMAP_12h_","MIXL1.pdf",sep=""),width = 10, height = 10)


FeaturePlot(Data1, features = c("BMPR1A","BMPR1B","BMPR2","ID3","ID4"), reduction = "umap") 
ggsave(filename=paste(saveext,"/UMAP_12h_","BMP.pdf",sep=""),width = 10, height = 15)



Data1 <- subset(mammal.combined, idents="PGCLC_12h")
Data1 <- FindNeighbors(Data1, reduction = "pca", dims = 1:20)



Data1 <- FindClusters(Data1, resolution = 0.5)
Data1$Cl5 <- Idents(Data1)
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Cl5_orig.pdf",sep=""),width = 8, height = 8)

markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl5_orig",".csv",sep=""))

Data1 <- FindClusters(Data1, resolution = 0.3)
Data1$Cl3 <- Idents(Data1)
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Cl3_orig.pdf",sep=""),width = 8, height = 8)

markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl3_orig",".csv",sep=""))


Idents(Data1) <- Data1$Phase
DimPlot(Data1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/UMAP_12h_","Phase.pdf",sep=""),width = 8, height = 8)



DefaultAssay(Data1) <- "RNA"
Idents(Data1) <- Data1$Cl5
markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl5_origRNA",".csv",sep=""))

Idents(Data1) <- Data1$Cl3
markers <- FindAllMarkers(Data1,only.pos = TRUE,test.use = "MAST")
write.table(as.data.frame(markers),paste(saveext,"/DE_12h_Cl3_origRNA",".csv",sep=""))





Data1 <- subset(mammal.combined, idents=c("AraEBs (in vitro)"))
Data2 <- subset(mammal.combined, idents=c("8) Cynomolgus (in vitro)","AraEBs (in vitro)"))
Data3 <- subset(mammal.combined, idents=c("9) Marmoset (in vivo)","AraEBs (in vitro)"))

uID <- as.character(Data2$ID2)
uID[which(Data2$species2=="AraEBs (in vitro)")] <- "EB"
Idents(Data2) <- uID

#Human 2
MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","PS","Endoderm","Amnion","Mesoderm")
MarmC <- c("lightgrey","#EF5910","#BF470D","#C81048","#262323","#000000","#54DE62","#10C357","#07971F","#afdfd6","#edbae5","#0E672F","#419CD6")
           #,"#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434")


#MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC","PS","Endoderm","Amnion","Mesoderm")
#MarmC <- c("lightgrey","#EF5910","#C81048","#8D340A","black","#afdfd6","#edbae5","#0E672F","#419CD6")


#MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","VE_CS4","VE_CS5","EPI","PGC","PE")
#MarmC <- c("lightgrey","#F16E2F","#DE520F","#C84A0D","#686868","#828282","#F0621F","#C81048","#5C5C5C")

D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data2,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



#Human 2

#MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","PS","Endoderm","Amnion","Mesoderm")
#MarmC <- c("lightgrey","#EF5910","#BF470D","#C81048","#262323","#000000","#54DE62","#10C357","#07971F","#afdfd6","#edbae5","#0E672F","#419CD6")

MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","PS","Endoderm","Amnion","Mesoderm")
MarmC <- c("lightgrey","#EF5910","#BF470D","#C81048","#262323","#000000","#54DE62","#10C357","#07971F","#afdfd6","#edbae5","#0E672F","#419CD6")
#,"#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434")


#MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","VE_CS4","VE_CS5","EPI","PGC","PE")
#MarmC <- c("lightgrey","#F16E2F","#DE520F","#C84A0D","#686868","#828282","#F0621F","#C81048","#5C5C5C")


saveext1 = "~/Desktop/Thorsten/FINAL/Human_amnion_EmDisc3/"
Endo<-readRDS(paste(saveext1,"/MarmEndo",".rds",sep=""))
Meso<-readRDS(paste(saveext1,"/MarmMes",".rds",sep=""))
PS<-readRDS(paste(saveext1,"/MarmPS",".rds",sep=""))

uID <- as.character(Data3$ID2)
uID[which(Data3$species2=="AraEBs (in vitro)")] <- "EB"
Idents(Data3) <- uID
Idents(object = Data3, cells =paste(Meso,"_1",sep="")) <- "Mesoderm"
Idents(object = Data3, cells =paste(Endo,"_1",sep="")) <- "Endoderm"
Idents(object = Data3, cells =paste(PS,"_1",sep="") ) <- "PS"

#uID <- as.character(Data3$ID2)
#uID[which(Data3$species2=="AraEBs (in vitro)")] <- "EB"
#Idents(Data3) <- uID

D2<- Idents(Data3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]

ord <- c("EB","PGC_CS5","PS","Am_CS5","EmDisc_CS5","Endoderm","Mesoderm","Am_CS7","EmDisc_CS7","Am_CS6","PGC_CS6","EmDisc_CS6")

DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_marm_","_1_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_marm_","_2_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_marm_","_3_nL.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


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

ord <- c("EB","PGC_CS5","PS","Am_CS5","EmDisc_CS5","Endoderm","Mesoderm","Am_CS7","EmDisc_CS7","Am_CS6","PGC_CS6","EmDisc_CS6")
DimPlot(Data1,   pt.size = 4, order =  ord, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,   pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,    pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_Cl_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



Markers <- c(
  "SOX2",
  "PRDM14",
  "POU5F1",
  "NANOG",
  "TBXT",
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



#redblue1<-colorRampPalette(c("#007d87","#FFFFFF","#870a00"))
#mat_breaks <- seq(0, 1, length.out = 60)
##pheatmap(log2(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_1",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks,color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_1_norm",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#00FFFF","#FFFFFF","#FF0000"))
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks,color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_2_norm",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#007d87","#FFFFFF","#870a00"))
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks,color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_3_norm",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks,color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_4_norm",".pdf",sep=""),width=4,height=6)


redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 1, length.out = 60)
pheatmap(log2(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_2",".pdf",sep=""),width=4,height=6)
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")],color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_2_norm",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#007d87","#FFFFFF","#870a00"))
mat_breaks <- seq(0, 1, length.out = 60)
pheatmap(log2(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_3",".pdf",sep=""),width=4,height=6)
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks, color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_3_norm",".pdf",sep=""),width=4,height=6)

redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))
mat_breaks <- seq(0, 1, length.out = 60)
pheatmap(log2(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")] +1),breaks = mat_breaks,color =  redblue1(60), gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_ara_fig_4",".pdf",sep=""),width=4,height=6)
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","6","0","28","17","20","DE","9","4","PGC","10","13","8","5","7")], breaks = mat_breaks, color =  redblue1(60), border_color = NA, gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,11,14), cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_4_norm",".pdf",sep=""),width=4,height=6)


#



mat_breaks <- seq(-2, 2, length.out = 20)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","7","1","17","13","DE2","9","6","PGC","11","4","5")],breaks = mat_breaks,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_1a_norm",".pdf",sep=""),width=4,height=6)
mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","7","1","17","13","DE2","9","6","PGC","11","4","5")], breaks = mat_breaks,color =  redblue1(60),gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,10,13), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_1b_normss",".pdf",sep=""),width=4,height=7)

redblue1<-colorRampPalette(c("#1EB5BE","#FFFFFF","#D61C1C"))

mat_breaks <- seq(-1, 1, length.out = 60)
pheatmap(AvExp$RNA[Markers,c("4i","ESC","PreME","ME","2","7","1","17","13","DE2","9","6","PGC","11","4","5")], breaks = mat_breaks,color =  redblue1(60),gaps_row = c(4,14,21,32,35), gaps_col = c(3,7,10,13), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/DimRed/HM_ara_fig_1b_normss",".pdf",sep=""),width=4,height=7)










Idents(mammal.combined) <- mammal.combined$species2
Data2 <- subset(mammal.combined, idents=c("AraEBs (in vitro)","7) Human CS7 (in vivo)"))
#Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 100)
Idents(Data2) <- droplevels(Data2$ID3)

Idents(Data2) <- droplevels(Data2$ID2)


Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 200)
KNN_K200 <- Data2@graphs$integrated_nn
SNN_K200 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 100)
KNN_K100 <- Data2@graphs$integrated_nn
SNN_K100 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 50)
KNN_K50 <- Data2@graphs$integrated_nn
SNN_K50 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 30)

KNN_K30<- Data2@graphs$integrated_nn
SNN_K30 <- Data2@graphs$integrated_snn

saveRDS(KNN_K100,file=paste(saveext,"KNN_K100.rds",sep=""))
saveRDS(SNN_K100,file=paste(saveext,"SNN_K100.rds",sep=""))
saveRDS(KNN_K30,file=paste(saveext,"KNN_K30.rds",sep=""))
saveRDS(SNN_K30,file=paste(saveext,"SNN_K30.rds",sep=""))
saveRDS(KNN_K50,file=paste(saveext,"KNN_K50.rds",sep=""))
saveRDS(SNN_K50,file=paste(saveext,"SNN_K50.rds",sep=""))

Idents(Data2) <- droplevels(Data2$ID3)

uID <- as.character(Data2$ID3)
#uID[which(Data2$species2=="3) Marmoset)")] <- as.character(Data2$ID2[which(Data2$species2=="7) Marmoset (in vivo)")])


#Idents(Data2) <- uID

DimPlot(Data2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE,pt.size = 2) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_split2","_cs7map.pdf",sep=""),width = 20, height = 10)
DimPlot(Data2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE,pt.size = 2) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split2","_cs7map.pdf",sep=""),width = 20, height = 10)

uID <- Idents(Data2)
L1 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L1b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P1 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L2 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L2b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P2 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L3 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L3b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P3 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))

L4 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L4b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P4 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))


Lb <- table(Idents(Data2))
nn <- length(Idents(Data2))

for (i in 1:length(Idents(Data2)) ) {
  L1[i,] <- table(Idents(Data2)[which(SNN_K30[i,]>0)])
  L1b[i,] <- L1[i,] / Lb
  P1[i,] <- phyper(L1[i,]-1, Lb, nn-Lb, 30, log = TRUE)
  
  L2[i,] <- table(Idents(Data2)[which(SNN_K50[i,]>0)])
  L2b[i,] <- 0.5*(L2b[i,] / Lb)
  P2[i,] <- phyper(L2[i,]-1, Lb, nn-Lb, 50, log = TRUE)
  
  L3[i,] <- table(Idents(Data2)[which(SNN_K100[i,]>0)])
  L3b[i,] <- (L3b[i,] / Lb)/3
  P3[i,] <- phyper(L3[i,]-1, Lb, nn-Lb, 100, log = TRUE)
  
  L4[i,] <- table(Idents(Data2)[which(SNN_K200[i,]>0)])
  L4b[i,] <- (L4b[i,] / Lb)/3
  P4[i,] <- phyper(L4[i,]-1, Lb, nn-Lb, 200, log = TRUE)  
  
}

colnames(L1)  <- rownames(Lb)
colnames(L1b) <- rownames(Lb)
colnames(P1)  <- rownames(Lb)
colnames(L2)  <- rownames(Lb)
colnames(L2b) <- rownames(Lb)
colnames(P2)  <- rownames(Lb)
colnames(L3)  <- rownames(Lb)
colnames(L3b) <- rownames(Lb)
colnames(P3)  <- rownames(Lb)
rownames(L1) <- colnames(Data2)
rownames(L2) <- colnames(Data2)
rownames(L3) <- colnames(Data2)
rownames(L1b) <- colnames(Data2)
rownames(L2b) <- colnames(Data2)
rownames(L3b) <- colnames(Data2)
rownames(P1) <- colnames(Data2)
rownames(P2) <- colnames(Data2)
rownames(P3) <- colnames(Data2)

rownames(L4) <- colnames(Data2)
rownames(P4) <- colnames(Data2)

fileID <- saveext

saveRDS(L1,paste(fileID,'L1.rds',sep=""))
saveRDS(L1b,paste(fileID,'L1b.rds',sep=""))
saveRDS(P1,paste(fileID,'P1.rds',sep=""))
saveRDS(L2,paste(fileID,'L2.rds',sep=""))
saveRDS(L2b,paste(fileID,'L2b.rds',sep=""))
saveRDS(P2,paste(fileID,'P2.rds',sep=""))
saveRDS(L3,paste(fileID,'L3.rds',sep=""))
saveRDS(L3b,paste(fileID,'L3b.rds',sep=""))
saveRDS(P3,paste(fileID,'P3.rds',sep=""))
saveRDS(L4,paste(fileID,'L4.rds',sep=""))
saveRDS(P4,paste(fileID,'P4.rds',sep=""))

write.table(L1,file=paste(fileID,'L1.csv',sep=""),sep=",")
write.table(L1b,file=paste(fileID,'L1b.csv',sep=""),sep=",")
write.table(P1,file=paste(fileID,'P1.csv',sep=""),sep=",")
write.table(L2,file=paste(fileID,'L2.csv',sep=""),sep=",")
write.table(L2b,file=paste(fileID,'L2b.csv',sep=""),sep=",")
write.table(P2,file=paste(fileID,'P2.csv',sep=""),sep=",")
write.table(L3,file=paste(fileID,'L3.csv',sep=""),sep=",")
write.table(L3b,file=paste(fileID,'L3b.csv',sep=""),sep=",")
write.table(P3,file=paste(fileID,'P3.csv',sep=""),sep=",")

#L1 <- readRDS(paste(fileID,'L1.rds',sep=""))
#L1b <- readRDS(L1b,paste(fileID,'L1b.rds',sep=""))
P1 <- readRDS(paste(fileID,'P1.rds',sep=""))
#L2 <- readRDS(L2,paste(fileID,'L2.rds',sep=""))
#L2b <- readRDS(L2b,paste(fileID,'L2b.rds',sep=""))
P2 <- readRDS(paste(fileID,'P2.rds',sep=""))
#L3 <- readRDS(L3,paste(fileID,'L3.rds',sep=""))
#L3b <- readRDS(L3b,paste(fileID,'L3b.rds',sep=""))
P3 <- readRDS(paste(fileID,'P3.rds',sep=""))
#L4 <- readRDS(L4,paste(fileID,'L4.rds',sep=""))
P4 <- readRDS(paste(fileID,'P4.rds',sep=""))

#2 advance
#3 axial
#5 ectoderm
#6 emergent
#7 endoder
#8 epiblast
#11 nascent
#12 pgc
#21 prim
#22 YS mes

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- (t((P1[,5])))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AdvancedMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,6]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AxialMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,8]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'EmergentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,9]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Endoderm')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,10]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Epiblast')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,13]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'NascentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,2]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PGCC')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,23]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PrimitiveStreak')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,24]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'YSMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,1]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Amnion')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P1[,7]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Ectoderm')

Idents(Data2) <- Data2$species2
cel <- WhichCells(Data2,idents="AraEBs (in vitro)")
FeaturePlot(Data2, features = "AdvancedMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AdvancedMes_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "AxialMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AxialMes_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "NascentMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/NascentMes_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "EmergentMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/EmergentMes_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "YSMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/YSMes_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Endoderm", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Endoderm_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Epiblast", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Epiblast_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PGC", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PGC_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PrimitiveStreak", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PrimitiveStreak_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Ectoderm", cells = cel, split.by = "species2" , pt.size = 2, cols =   c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Ectoderm_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 

FeaturePlot(Data2, features = "Amnion", cells = cel, split.by = "species2" , pt.size = 2, cols = c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Amnion_P1_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 


cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- (t((P2[,5])))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AdvancedMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,6]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AxialMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,8]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'EmergentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,9]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Endoderm')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,10]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Epiblast')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,13]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'NascentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,2]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PGC')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,23]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PrimitiveStreak')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,24]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'YSMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,1]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Amnion')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P2[,7]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Ectoderm')

Idents(Data2) <- Data2$species2
cel <- WhichCells(Data2,idents="AraEBs (in vitro)")
FeaturePlot(Data2, features = "AdvancedMes", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AdvancedMes_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "AxialMes", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AxialMes_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "NascentMes", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/NascentMes_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "EmergentMes", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/EmergentMes_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "YSMes", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/YSMes_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Endoderm", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Endoderm_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Epiblast", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Epiblast_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PGC", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PGC_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PrimitiveStreak", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PrimitiveStreak_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Ectoderm", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Ectoderm_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Amnion", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))   + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Amnion_P2_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- (t((P3[,5])))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AdvancedMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,6]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AxialMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,8]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'EmergentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,9]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Endoderm')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,10]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Epiblast')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,13]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'NascentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,2]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PGC')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,23]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PrimitiveStreak')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,24]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'YSMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,1]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Amnion')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P3[,7]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Ectoderm')

Idents(Data2) <- Data2$species2
cel <- WhichCells(Data2,idents="AraEBs (in vitro)")
FeaturePlot(Data2, features = "AdvancedMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AdvancedMes_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "AxialMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AxialMes_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "NascentMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/NascentMes_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "EmergentMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/EmergentMes_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "YSMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/YSMes_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Endoderm", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Endoderm_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Epiblast", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Epiblast_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PGC", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PGC_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PrimitiveStreak", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PrimitiveStreak_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Ectoderm", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Ectoderm_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Amnion", cells = cel, split.by = "species2" , pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Amnion_P3_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 








cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- (t((P4[,5])))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AdvancedMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,6]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'AxialMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,8]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'EmergentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,9]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Endoderm')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,10]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Epiblast')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,13]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'NascentMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,2]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PGC')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,23]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'PrimitiveStreak')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,24]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'YSMes')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,1]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Amnion')

cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- t((P4[,7]))
cluster_letters <- as.numeric(cluster_letters)
cluster_letters <- p.adjust(1-exp(cluster_letters), method = "bonferroni")
cluster_letters[which(cluster_letters>0.05)] <- 1
names(cluster_letters)=rownames(Data2@meta.data)
Data2 <- AddMetaData(object = Data2, metadata = cluster_letters, col.name = 'Ectoderm')

Idents(Data2) <- Data2$species2
cel <- WhichCells(Data2,idents="AraEBs (in vitro)")
FeaturePlot(Data2, features = "AdvancedMes", cells = cel, split.by = "species2" , pt.size = 2, cols = c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AdvancedMes_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "AxialMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/AxialMes_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "NascentMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/NascentMes_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "EmergentMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/EmergentMes_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "YSMes", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/YSMes_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Endoderm", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Endoderm_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Epiblast", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Epiblast_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PGC", cells = cel, split.by = "species2" , pt.size = 2, cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PGC_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "PrimitiveStreak", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/PrimitiveStreak_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Ectoderm", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Ectoderm_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 
FeaturePlot(Data2, features = "Amnion", cells = cel, split.by = "species2" , pt.size = 2,  cols =  c("#dd1c77","lightgrey")) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/Markers/Amnion_P4_col1.pdf",sep=""),width = 15, height = 15,limitsize = FALSE) 


BestMod <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
ScoreBestMode <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
subID <- Lb[c(1,2,5,6,7,8,9,10,13,23,24)]
for (i in 1:length(Idents(Data2)) ) {
  BestMod[i] <-  which(P1[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P1[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]
  ScoreBestMode[i] <-  P1[i,c(1,2,5,6,7,8,9,10,13,23,24)][which(P1[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P1[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]]
}
cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "None"
TransferL <-  rownames(subID)[BestMod]
cluster_letters[which(ScoreBestMode>-Inf)] <- TransferL[which(ScoreBestMode>-Inf)]
Data2$TestID <- cluster_letters
MarmType <- c("None","ESC","Primitive Streak","Endoderm","PGC_CS7","Ectoderm","Am_CS7","Emergent Mesoderm","Nascent Mesoderm","Advanced Mesoderm","Axial Mesoderm","YS Mesoderm","Epiblast","EB")
MarmC <- c("lightgrey","#E96F08","#E73D19","#1B1B1B","#E72185","#0C9542","#0E672F","#419CD6","#1B84C7","#0A67A0","#3293D1","#06476E","#107F87","#D3D3D3")
Idents(Data2) <- cluster_letters
D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Transfer_CS7_","P1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

Data2$CS7 <- Idents(Data2)
Idents(Data2) <- Data2$species2
Data3 <- subset(Data2,idents="AraEBs (in vitro)")
Idents(Data3) <- Data3$Cells

Data3a <- Data3
Species4 <- droplevels(Data3a$Cells)
Species4  <- factor(Species4, levels = c("4i","ME_r2","DE","hESC_W5_E8","PreME_10h","PGCLC_12h_r2","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_D4","PGCLC_D4_r2","PGC_wk6.5_invivo"))
Data3a$Species4 <- Species4
#Idents(Data3a) <- Data3a$CS7
Idents(Data3a) <- Data3a$CS7
DimPlot(Data3a,  cols = coluse1, pt.size = 4, reduction = "umap", split.by = "Species4",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Transfer_CS7_","Split.pdf",sep=""),width = 210, height = 15, useDingbats=FALSE, limitsize = FALSE)



BestMod <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
ScoreBestMode <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
subID <- Lb[c(1,2,5,6,7,8,9,10,13,23,24)]
for (i in 1:length(Idents(Data2)) ) {
  BestMod[i] <-  which(P2[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P2[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]
  ScoreBestMode[i] <-  P2[i,c(1,2,5,6,7,8,9,10,13,23,24)][which(P2[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P2[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]]
}
cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "None"
TransferL <-  rownames(subID)[BestMod]
cluster_letters[which(ScoreBestMode>-Inf)] <- TransferL[which(ScoreBestMode>-Inf)]
Data2$TestID <- cluster_letters
MarmType <- c("None","ESC","Primitive Streak","Endoderm","PGC_CS7","Ectoderm","Am_CS7","Emergent Mesoderm","Nascent Mesoderm","Advanced Mesoderm","Axial Mesoderm","YS Mesoderm","Epiblast","EB")
MarmC <- c("lightgrey","#E96F08","#E73D19","#1B1B1B","#E72185","#0C9542","#0E672F","#419CD6","#1B84C7","#0A67A0","#3293D1","#06476E","#107F87","#D3D3D3")
Idents(Data2) <- cluster_letters
D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Transfer_CS7_P2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


BestMod <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
ScoreBestMode <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
subID <- Lb[c(1,2,5,6,7,8,9,10,13,23,24)]
for (i in 1:length(Idents(Data2)) ) {
  BestMod[i] <-  which(P3[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P3[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]
  ScoreBestMode[i] <-  P3[i,c(1,2,5,6,7,8,9,10,13,23,24)][which(P3[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P3[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]]
}
cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "None"
TransferL <-  rownames(subID)[BestMod]
cluster_letters[which(ScoreBestMode>-Inf)] <- TransferL[which(ScoreBestMode>-Inf)]
Data2$TestID <- cluster_letters
MarmType <- c("None","ESC","Primitive Streak","Endoderm","PGC_CS7","Ectoderm","Am_CS7","Emergent Mesoderm","Nascent Mesoderm","Advanced Mesoderm","Axial Mesoderm","YS Mesoderm","Epiblast","EB")
MarmC <- c("lightgrey","#E96F08","#E73D19","#1B1B1B","#E72185","#0C9542","#0E672F","#419CD6","#1B84C7","#0A67A0","#3293D1","#06476E","#107F87","#D3D3D3")
Idents(Data2) <- cluster_letters
D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Transfer_CS7_P3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)


BestMod <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
ScoreBestMode <- matrix(, ncol = 1, nrow = length(Idents(Data2)))
subID <- Lb[c(1,2,5,6,7,8,9,10,13,23,24)]
for (i in 1:length(Idents(Data2)) ) {
  BestMod[i] <-  which(P4[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P4[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]
  ScoreBestMode[i] <-  P4[i,c(1,2,5,6,7,8,9,10,13,23,24)][which(P4[i,c(1,2,5,6,7,8,9,10,13,23,24)]==max(P4[i,c(1,2,5,6,7,8,9,10,13,23,24)]))[1]]
}
cluster_letters <- LETTERS[Data2$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "None"
TransferL <-  rownames(subID)[BestMod]
cluster_letters[which(ScoreBestMode>-Inf)] <- TransferL[which(ScoreBestMode>-Inf)]
Data2$TestID <- cluster_letters
MarmType <- c("None","ESC","Primitive Streak","Endoderm","PGC_CS7","Ectoderm","Am_CS7","Emergent Mesoderm","Nascent Mesoderm","Advanced Mesoderm","Axial Mesoderm","YS Mesoderm","Epiblast","EB")
MarmC <- c("lightgrey","#E96F08","#E73D19","#1B1B1B","#E72185","#0C9542","#0E672F","#419CD6","#1B84C7","#0A67A0","#3293D1","#06476E","#107F87","#D3D3D3")
Idents(Data2) <- cluster_letters
D2<- Idents(Data2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Transfer_CS7_P4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



#clusterdatasets
D1 <- subset(mammal.combined,idents=c("AraEBs (in vitro)","7) Human CS7 (in vivo)","11) Human (in vitro 2)"))
Idents(D1) <- D1$ID3
D1 <- subset(D1,idents=c("EPI","DE"),invert=TRUE)
Idents(D1) <- D1$Cl9
D1 <- subset(D1,idents=c("4","9"))
Idents(D1) <- D1$ID3

DimPlot(D1,  pt.size = 0.5, reduction = "umap", split.by = "species2", label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset9_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

DefaultAssay(D1)<- "integrated"


D1 <- RunPCA(D1, npcs = 30, verbose = FALSE)
D1 <- RunUMAP(D1, reduction = "pca", dims = 1:20)
D1 <- RunTSNE(D1, reduction = "pca", dims = 1:20)
D1 <- FindNeighbors(D1, reduction = "pca", dims = 1:20)

DimPlot(D1,  pt.size = 2, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE) #+ xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/PGC_subPCA.pdf",sep=""),width = 45, height = 15, useDingbats=FALSE)

DimPlot(D1,  pt.size = 2, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/PGC_subUMAP.pdf",sep=""),width = 45, height = 15, useDingbats=FALSE)




#Compare to Tb
Idents(mammal.combined) <- mammal.combined$species2

Tb <- readRDS('/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/Tb/TbC.rds')
Tb <- NormalizeData(Tb, verbose = FALSE)
#Tb <- FindVariableFeatures(Tb, selection.method = "vst", nfeatures = 20000)

Data1 <- subset(mammal.combined, idents="AraEBs (in vitro)")
Idents(Data1) <- Data1$ID3

Tbmerge <- merge(Data1, y = c(Tb))
Tbmerge <- FindVariableFeatures(Tbmerge, selection.method = "vst", nfeatures = 20000)
Tbmerge <- ScaleData(Tbmerge) #, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Data1))
Tbmerge <- RunPCA(Tbmerge, npcs = 20, verbose = FALSE)
Tbmerge <- RunUMAP(Tbmerge, reduction = "pca", dims = 1:20)
Tbmerge <- FindNeighbors(Tbmerge, reduction = "pca", dims = 1:20)

Tbmerge$IDX <- Idents(Tbmerge)
IDDX <- as.character(Tbmerge$IDX)
IDDX[which(Idents(Tbmerge)=="5")] <- "lAMLC"
IDDX[which(Idents(Tbmerge)=="8")] <- "lAMLC"
IDDX[which(Idents(Tbmerge)=="7")] <- "lAMLC"

Idents(Tbmerge) <- IDDX
DimPlot(Tbmerge,  pt.size = 2, reduction = "umap", label = TRUE, repel = TRUE) #+ xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/Tb_EB_merge.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



DefaultAssay(Tbmerge) <- "RNA"
Ae <- AverageExpression(Tbmerge, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)
#


SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)

SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(Tbmerge, ident.2 = "lAMLC", ident.1 = "SCT", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(lAMLC,SCT)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"AMLC_SCT.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(Tbmerge, ident.2 = "lAMLC", ident.1 = "EVT", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(lAMLC,EVT)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"AMLC_EVT.pdf",sep=""),width = 13, height = 13, plot = p1)

Ae5 <- Ae
Ae5$Pval <- 1
Cl1 <- FindMarkers(Tbmerge, ident.2 = "lAMLC", ident.1 = "VCT", verbose = FALSE, test.use = "MAST", only.pos = FALSE)
Ae5[rownames(Cl1),"Pval"] <- 0
Ae5$Pval <- as.factor(Ae5$Pval)
Ae5 <- Ae5[order(Ae5$Pval),]
nn <- dim(Ae5)[2]-2
Ae5[,1:nn] <- log1p(Ae5[,1:nn])
genes.to.label = intersect(rownames(Cl1),TF)
p1 <- ggplot(Ae5, aes(lAMLC,VCT)) + geom_point(aes(color=Pval)) + theme_classic() +  scale_color_manual(values=c('black','lightgrey'))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'red')
ggsave(filename=paste(saveext,"AMLC_VCT.pdf",sep=""),width = 13, height = 13, plot = p1)





FeaturePlot(mammal.combined2, features = "SOX17", pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/SOX17Zoom.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)



FeaturePlot(mammal.combined2, features = "DIO3", pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DIO3Zoom.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)


FeaturePlot(mammal.combined2, features = "HAND1", pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/HAND1Zoom.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)


FeaturePlot(mammal.combined2, features = "KRT19", pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/KRT19Zoom.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)









#Ara'sdata
Genes<-read.table("/Users/christopherpenfold/Desktop/Aracely/Pairwise_PGCLC_vs_Amnioids/Plots/GenesMapping.csv",sep=",",header = T)
metaD<-read.table("/Users/christopherpenfold/Dropbox/PGCs/Processed.dir/ALL_meta.tsv",sep="\t",header = T)
D1<-readRDS('/Users/christopherpenfold/Dropbox/PGCs/Processed.dir/ALL_SCE.RDS')
rownames(D1) <- Genes$Gene
humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
#humaninvit_data$species <- "Human (Ara EBs)"
#humaninvit_data$cond <- "Ape"
labs_2 <- as.character(metaD$CellOrigin)
#Idents(humaninvit_data) <- labs_2 #as.character(metaD$CellOrigin)

#Get variance out
ind1 <- which(labs_2=="hESC_W5_E8")
ind2 <- which(labs_2=="4i")
ind3 <- which(labs_2=="PreME_10h")
ind4 <- which(labs_2=="PGCLC_12h_r2")
ind5 <- which(labs_2=="PGCLC_18h")
ind6 <- which(labs_2=="PGCLC_24h")
ind7 <- which(labs_2=="PGCLC_32h")
ind8 <- which(labs_2=="PGCLC_40h")
ind9 <- which(labs_2=="PGCLC_48h")
ind10 <- which(labs_2=="PGCLC_D4")
ind11 <- which(labs_2=="PGCLC_D4_r2")
ind12 <- which(labs_2=="DE")
ind13 <- which(labs_2=="ME_r2")
ind14 <- which(labs_2=="PGC_wk6.5_invivo")

raw.data <- as.matrix(log(GetAssayData(humaninvit_data, slot = "counts")+1))
dmean <-rowMeans(raw.data[,ind1])

mean1 <- rowMeans(raw.data[,ind1]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean2 <- rowMeans(raw.data[,ind2]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean3 <- rowMeans(raw.data[,ind3]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean4 <- rowMeans(raw.data[,ind4]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean5 <- rowMeans(raw.data[,ind5]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)

mean6 <- rowMeans(raw.data[,ind6]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean7 <- rowMeans(raw.data[,ind7]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean8 <- rowMeans(raw.data[,ind8]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean9 <- rowMeans(raw.data[,ind9]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean10 <- rowMeans(raw.data[,ind10]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean11 <- rowMeans(raw.data[,ind11]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean12 <- rowMeans(raw.data[,ind12]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean13 <- rowMeans(raw.data[,ind13]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)
mean14 <- rowMeans(raw.data[,ind14]) #dmean #rowMean(raw.data[,ind1],na.rm=TRUE)

#dataO$var <- dvar #rowMean(raw.data[,ind1],na.rm=TRUE)


var1 <- apply(raw.data[,ind1],1,var) #vars(raw.data)
var2 <- apply(raw.data[,ind2],1,var) #vars(raw.data)
var3 <- apply(raw.data[,ind3],1,var) #vars(raw.data)
var4 <- apply(raw.data[,ind4],1,var) #vars(raw.data)
var5 <- apply(raw.data[,ind5],1,var) #vars(raw.data)
var6 <- apply(raw.data[,ind6],1,var) #vars(raw.data)
var7 <- apply(raw.data[,ind7],1,var) #vars(raw.data)
var8 <- apply(raw.data[,ind8],1,var) #vars(raw.data)
var9 <- apply(raw.data[,ind9],1,var) #vars(raw.data)
var10 <- apply(raw.data[,ind10],1,var) #vars(raw.data)
var11 <- apply(raw.data[,ind11],1,var) #vars(raw.data)
var12 <- apply(raw.data[,ind12],1,var) #vars(raw.data)
var13 <- apply(raw.data[,ind13],1,var) #vars(raw.data)
var14 <- apply(raw.data[,ind14],1,var) #vars(raw.data)

df <- data.frame(
  Type=c("ESC","4i","PreME","12h","28h","24h","32h","40h","48h","96h","96h4i","DE","ME","PGC"), 
  Var=c(mean(var1),mean(var2),mean(var3),mean(var4),mean(var5),mean(var6),mean(var7),mean(var8),mean(var9),mean(var10),mean(var11),mean(var12),mean(var13),mean(var14)))
                  
Type2 <- factor(df$Type, 
                levels = c("ESC","4i","PreME","12h","28h","24h","32h","40h","48h","96h","96h4i","DE","ME","PGC"))          
df$Type2 <- Type2
p<-ggplot(data=df, aes(x=Type2, y=Var)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
p
ggsave(filename=paste(saveext,"/Variance.pdf",sep=""),p, width = 10, height = 10, useDingbats=FALSE)








Idents(mammal.combined) <- Species2
Data1 <- subset(mammal.combined, idents=c("AraEBs (in vitro)","8) Cynomolgus (in vitro)"))

Idents(Data1) <- Data1$ID2

IDD1 <- as.character(Data1$ID2)
IDD1[which(Data1$species=="Human (Ara EBs)")] <- "EB"
IDD1[which(IDD1=="EmDisc_CS5")] <- "EmDisc_CS5"
IDD1[which(IDD1=="Am_CS5")] <- "Am_CS5"

IDD2 <- as.character(Data1$ID2)
IDD2[which(Data1$species=="Human (Ara EBs)")] <- "EB"
IDD2[which(IDD2=="EmDisc_CS6")] <- "EmDisc_CS6"
IDD2[which(IDD2=="Am_CS6")] <- "Am_CS6"
IDD2[which(IDD2=="PGC_CS6")] <- "PGC_CS6"

IDD3 <- as.character(Data1$ID2)
IDD3[which(Data1$species=="Human (Ara EBs)")] <- "EB"
IDD3[which(IDD1=="EmDisc_CS7")] <- "EmDisc_CS7"
IDD3[which(IDD1=="Am_CS7")] <- "Am_CS7"

Data1_1 <- Data1
Idents(Data1_1) <- IDD1 
Data1_1 <- subset(Data1_1, idents=c("EB","EmDisc_CS5","Am_CS5"))

Data1_2 <- Data1
Idents(Data1_2) <- IDD2
Data1_2 <- subset(Data1_2, idents=c("EB","EmDisc_CS6","Am_CS6","PGC_CS6"))

Data1_3 <- Data1
Idents(Data1_3) <- IDD3 
Data1_3 <- subset(Data1_3, idents=c("EB","EmDisc_CS7","Am_CS7"))




#Human 2
MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","EPI","PGC","PE")
MarmC <- c("lightgrey","#BF470D","#C81048","#8D340A","#C81048","#C81048","#54DE62","#10C357","#07971F","#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434","#F0621F","#C81048","#5C5C5C")

D2<- Idents(Data1_1)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1_1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS5_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_cyno_","CS5_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_1,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS5_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)





#Human 2
MarmType <- c("EB","EmDisc_CS6","PGC_CS6","Am_CS6")
MarmC <- c("lightgrey","#BF470D","#C81048","#8D340A","#54DE62")

D2<- Idents(Data1_2)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1_2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS6_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_2,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_cyno_","CS6_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_2,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS6_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



###### Now do cyno...
uID <- as.character(Data1$ID2)
uID[which(Data1$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data1) <- uID



#Human 2
MarmType <- c("EB","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","EPI","PGC","PE")
MarmC <- c("lightgrey","#BF470D","#C81048","#8D340A","#C81048","#C81048","#54DE62","#10C357","#07971F","#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434","#F0621F","#C81048","#5C5C5C")

D2<- Idents(Data1_3)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1_3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS7_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_3,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_cyno_","CS7_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1_3,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","CS7_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



###### Now do cyno...
uID <- as.character(Data1$ID2)
uID[which(Data1$species2=="1) Human EBs (in vitro)")] <- "EB"
Idents(Data1) <- uID



#Human 2
MarmType <- c("EB","Epi_CS3","EmDisc_CS5","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","EPI","PGC","PE")
MarmC <- c("lightgrey","#F27A40","#EF5910","#BF470D","#C81048","#8D340A","#C81048","#C81048","#54DE62","#10C357","#07971F","#5CBBF5","#3CAEF3","#0F99ED","#4E4E4E","#343434","#F0621F","#C81048","#5C5C5C")

D2<- Idents(Data1)
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(MarmType==levels(D2)[i])
}
coluse1 <- MarmC[colind]
DimPlot(Data1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_1.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,  cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset6_cyno_","_2.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)
DimPlot(Data1,   cols = coluse1, pt.size = 4, reduction = "umap",  label = FALSE, repel = TRUE) + NoLegend() + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset_cyno_","_3.pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)



#Do some DE
DefaultAssay(Data1) <- "RNA"
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
#Mesoderm
Markers1 <- FindMarkers(Data1,ident.1="6",ident.2="2", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers2 <- FindMarkers(Data1,ident.1="0",ident.2="6", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#PS
Markers3 <- FindMarkers(Data1,ident.1="10",ident.2="2", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#Amn
Markers4 <- FindMarkers(Data1,ident.1="13",ident.2="10", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers5 <- FindMarkers(Data1,ident.1=c("5","8","7"),ident.2="13", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#PGCs
Markers6 <- FindMarkers(Data1,ident.1="9",ident.2="10", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers6B <- FindMarkers(Data1,ident.1="9",ident.2="2", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

Markers7 <- FindMarkers(Data1,ident.1="4",ident.2="9", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers8 <- FindMarkers(Data1,ident.1="PGC",ident.2="4", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#DE
Markers9 <- FindMarkers(Data1,ident.1=c("28","17","20"),ident.2="2", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#Pre
Markers10 <- FindMarkers(Data1,ident.1="PreME",ident.2="ESC", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers11 <- FindMarkers(Data1,ident.1="2",ident.2="PreME", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#
Markers12 <- FindMarkers(Data1,ident.1="0",ident.2="PreME", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers13 <- FindMarkers(Data1,ident.1="4",ident.2="PreME", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers14 <- FindMarkers(Data1,ident.1=c("20","17","28"),ident.2="PreME", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers15 <- FindMarkers(Data1,ident.1=c("8","5","7"),ident.2="PreME", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#Amnion vs other terminal lineages
Markers16 <- FindMarkers(Data1,ident.1=c("8","5","7"),ident.2="4", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers17 <- FindMarkers(Data1,ident.1=c("8","5","7"),ident.2="0", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers18 <- FindMarkers(Data1,ident.1=c("8","5","7"),ident.2=c("20","17","28"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#DE vs
Markers19 <- FindMarkers(Data1,ident.1=c("20","17","28"),ident.2="0", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers20 <- FindMarkers(Data1,ident.1=c("20","17","28"),ident.2="4", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers21 <- FindMarkers(Data1,ident.1=c("20","17","28"),ident.2=c("8","5","7"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#Mes vs
Markers22 <- FindMarkers(Data1,ident.1="0",ident.2="4", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers23 <- FindMarkers(Data1,ident.1="0",ident.2=c("8","5","7"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers24 <- FindMarkers(Data1,ident.1="0",ident.2=c("20","17","28"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))

#PGC vs
Markers25 <- FindMarkers(Data1,ident.1="4",ident.2="0", test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers26 <- FindMarkers(Data1,ident.1="4",ident.2=c("8","5","7"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))
Markers27 <- FindMarkers(Data1,ident.1="4",ident.2=c("20","17","28"), test.use = "MAST",only.pos = TRUE, logfc.threshold = log(1.5))


write.table(as.data.frame(Markers6B),paste(saveext,"/Markers6B",".csv",sep=""))

write.table(as.data.frame(Markers1),paste(saveext,"/Markers1",".csv",sep=""))
write.table(as.data.frame(Markers2),paste(saveext,"/Markers2",".csv",sep=""))
write.table(as.data.frame(Markers3),paste(saveext,"/Markers3",".csv",sep=""))
write.table(as.data.frame(Markers4),paste(saveext,"/Markers4",".csv",sep=""))
write.table(as.data.frame(Markers5),paste(saveext,"/Markers5",".csv",sep=""))
write.table(as.data.frame(Markers6),paste(saveext,"/Markers6",".csv",sep=""))
write.table(as.data.frame(Markers7),paste(saveext,"/Markers7",".csv",sep=""))
write.table(as.data.frame(Markers8),paste(saveext,"/Markers8",".csv",sep=""))
write.table(as.data.frame(Markers9),paste(saveext,"/Markers9",".csv",sep=""))
write.table(as.data.frame(Markers10),paste(saveext,"/Markers10",".csv",sep=""))
write.table(as.data.frame(Markers11),paste(saveext,"/Markers11",".csv",sep=""))
write.table(as.data.frame(Markers12),paste(saveext,"/Markers12",".csv",sep=""))
write.table(as.data.frame(Markers13),paste(saveext,"/Markers13",".csv",sep=""))
write.table(as.data.frame(Markers14),paste(saveext,"/Markers14",".csv",sep=""))
write.table(as.data.frame(Markers15),paste(saveext,"/Markers15",".csv",sep=""))

write.table(as.data.frame(Markers16),paste(saveext,"/Markers16",".csv",sep=""))
write.table(as.data.frame(Markers17),paste(saveext,"/Markers17",".csv",sep=""))
write.table(as.data.frame(Markers18),paste(saveext,"/Markers18",".csv",sep=""))
write.table(as.data.frame(Markers19),paste(saveext,"/Markers19",".csv",sep=""))
write.table(as.data.frame(Markers20),paste(saveext,"/Markers20",".csv",sep=""))
write.table(as.data.frame(Markers21),paste(saveext,"/Markers21",".csv",sep=""))
write.table(as.data.frame(Markers22),paste(saveext,"/Markers22",".csv",sep=""))
write.table(as.data.frame(Markers23),paste(saveext,"/Markers23",".csv",sep=""))
write.table(as.data.frame(Markers24),paste(saveext,"/Markers24",".csv",sep=""))
write.table(as.data.frame(Markers25),paste(saveext,"/Markers25",".csv",sep=""))
write.table(as.data.frame(Markers26),paste(saveext,"/Markers26",".csv",sep=""))
write.table(as.data.frame(Markers27),paste(saveext,"/Markers27",".csv",sep=""))



write.table(as.data.frame(rownames(Markers6B)),paste(saveext,"/Markers6Bgenes",".csv",sep=""),row.names= FALSE, quote = FALSE)

write.table(as.data.frame(rownames(Markers1)),paste(saveext,"/Markers1genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers2)),paste(saveext,"/Markers2genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers3)),paste(saveext,"/Markers3genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers4)),paste(saveext,"/Markers4genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers5)),paste(saveext,"/Markers5genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers6)),paste(saveext,"/Markers6genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers7)),paste(saveext,"/Markers7genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers8)),paste(saveext,"/Markers8genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers9)),paste(saveext,"/Markers9genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers10)),paste(saveext,"/Markers10genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers11)),paste(saveext,"/Markers11genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers12)),paste(saveext,"/Markers12genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers13)),paste(saveext,"/Markers13genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers14)),paste(saveext,"/Markers14genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers15)),paste(saveext,"/Markers15genes",".csv",sep=""),row.names= FALSE, quote = FALSE)

write.table(as.data.frame(rownames(Markers16)),paste(saveext,"/Markers16genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers17)),paste(saveext,"/Markers17genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers18)),paste(saveext,"/Markers18genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers19)),paste(saveext,"/Markers19genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers20)),paste(saveext,"/Markers20genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers21)),paste(saveext,"/Markers21genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers22)),paste(saveext,"/Markers22genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers23)),paste(saveext,"/Markers23genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers24)),paste(saveext,"/Markers24genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers25)),paste(saveext,"/Markers25genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers26)),paste(saveext,"/Markers26genes",".csv",sep=""),row.names= FALSE, quote = FALSE)
write.table(as.data.frame(rownames(Markers27)),paste(saveext,"/Marker27genes",".csv",sep=""),row.names= FALSE, quote = FALSE)


DefaultAssay(Data1) <- "RNA"
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

uID[which(Data1$Cl9=="2")] <- "PSLC"
uID[which(Data1$Cl9=="9")] <- "ePGCLC"
uID[which(Data1$Cl9=="4")] <- "lPGCLC"

uID[which(Data1$Cl9=="13")] <- "eAMLC"
uID[which(Data1$Cl9=="8")] <- "lAMLC"
uID[which(Data1$Cl9=="5")] <- "lAMLC"
uID[which(Data1$Cl9=="7")] <- "lAMLC"
uID[which(Data1$Cl9=="6")] <- "eMELC"
uID[which(Data1$Cl9=="0")] <- "lMELC"

uID[which(Data1$Cl9=="28")] <- "eDELC"
uID[which(Data1$Cl9=="17")] <- "eDELC"
uID[which(Data1$Cl9=="20")] <- "eDELC"
Idents(Data1) <- uID


#Idents(Data1) <- Data1$Cl5
DimPlot(Data1,   pt.size = 2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Dataset1_","_TestClCl9_Subtypes.pdf",sep=""),width = 17, height = 15, useDingbats=FALSE)


Av <- AverageExpression(Data1)

C1 <- cor(Av$RNA[,c("4i","ESC","PreME","ME","PSLC","eDELC","DE","eMELC","lMELC","eAMLC","lAMLC","ePGCLC","lPGCLC","PGC")])

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.7, 1, length.out = 60)
pheatmap(C1, breaks = mat_breaks,color =  redblue1(60), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Correlation_motif-clusters",".pdf",sep=""),width=4,height=4)



write.table(as.data.frame(Av$RNA),paste(saveext,"/AverageExpressionForMotifs",".csv",sep=""),row.names= FALSE, quote = FALSE)





Phen<-read.table("/Users/christopherpenfold/Desktop/Phenotypes.csv",header = F,sep=",")
Stage <- toupper(Phen$V2)

Phen <- toupper(Phen$V1)

dir.create(paste(saveext,"/Phenotypes/",sep=""))

DefaultAssay(Data1) <- "RNA"

for (i in 1:length(Phen)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1, features = as.character(Phen[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Phenotypes/","PhenMarkers", "_", as.character(Phen[i]),"_", Stage[i],"nosplit.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}




SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)

SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]

dir.create(paste(saveext,"/Receptor/",sep=""))

for (i in 1:length(Phen)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1, features = as.character(SIGNAL1[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Receptor/","PhenMarkers", "_", as.character(SIGNAL1[i]),"_nosplit.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}



dir.create(paste(saveext,"/Ligand/",sep=""))

for (i in 1:length(Phen)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1, features = as.character(SIGNAL2[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Ligand/","PhenMarkers", "_", as.character(SIGNAL2[i]),"_nosplit.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}


dir.create(paste(saveext,"/ECM/",sep=""))

for (i in 1:length(Phen)) {
  possibleError <-  tryCatch(
    FeaturePlot(Data1, features = as.character(SIGNAL3[i]), pt.size = 2, cols = c("lightgrey", "#dd1c77"))  + xlim(c(-13, 13)) + ylim(c(-13, 13)), error=function(e) e)
  
  if(inherits(possibleError, "error")) next
  
  ggsave(filename=paste(saveext,"/Ligand/","ECM", "_", as.character(SIGNAL3[i]),"_nosplit.pdf",sep=""),width = 15, height = 15,limitsize = FALSE)
  graphics.off() 
}





install.packages("devtools")
devtools::install_github("choisy/cutoff")
library(cutoff)
DW <- read.table("/Users/christopherpenfold/Downloads/expression.csv")
hist(as.numeric( DW$V1 ) )




Phen<-read.table("/Users/christopherpenfold/Desktop/Phenotypes.csv",header = F,sep=",")
Stage <- toupper(Phen$V2)
Phen <- toupper(Phen$V1)
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]


#mammal.combined$Cl9
Idents(mammal.combined) <- mammal.combined$species2
Data0 <- subset(mammal.combined, idents=c("AraEBs (in vitro)"))

#Idents(Data0) <- Data0$ID2
Idents(Data0) <- Data0$ID2
Data0B <- subset(Data0,idents=c("ESC","PreME","DE","ME","PGCLC_12h","PGCLC_18h","PGCLC_24h","PGCLC_32h","PGCLC_40h","PGCLC_48h","PGCLC_96h"))
Idents(Data0B) <- Data0B$Cl9
nID <- as.character(Idents(Data0B))
nID[which(Data0B$ID2=="ESC")] <- "ESC"
nID[which(Data0B$ID2=="PreME")] <- "PreME"
nID[which(Data0B$ID2=="ME")] <- "ME"
nID[which(Data0B$ID2=="DE")] <- "DE1"
nID[which(Data0B$ID2=="DE" & Data0B$Cl9=="15")] <- "DE"
nID[which(Data0B$ID2=="DE" & Data0B$Cl9=="16")] <- "DE"
Idents(Data0B) <- nID

DefaultAssay(Data0B) <- "RNA"

ge1 <- Phen$V1[which(Stage=="E9.5")]
ge2 <- Phen$V1[which(Stage=="E12.5")]
ge3 <- Phen$V1[which(Stage=="E15.5")]
ge4 <- Phen$V1[which(Stage=="E18.5")]
ge5 <- Phen$V1[which(Stage=="LIT")]

ExpR <- AverageExpression(Data0B) 

library(reshape2)

ExpR1 <- ExpR$RNA[unique(ge1),]
ExpR2 <- ExpR$RNA[unique(ge2),]
ExpR3 <- ExpR$RNA[unique(ge3),]
ExpR4 <- ExpR$RNA[unique(ge4),]
ExpR5 <- ExpR$RNA[unique(ge5),]


melt1 <- melt(ExpR1)
melt1$ID <- paste("Cl",melt1$variable,"E095",sep="_")
melt2 <- melt(ExpR2)
melt2$ID <- paste("Cl",melt2$variable,"E125",sep="_")
melt3 <- melt(ExpR3)
melt3$ID <- paste("Cl",melt3$variable,"E155",sep="_")
melt4 <- melt(ExpR4)
melt4$ID <- paste("Cl",melt4$variable,"E185",sep="_")

All <- rbind(melt1,melt2,melt3,melt4)
All$logV <- log2(All$value +1)

All$ID1 <- as.factor(All$ID)

ggplot(All, aes(x=ID, y=logV)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste(saveext,"/PhenotypeExp", ".pdf",sep=""),width = 95, height = 5,limitsize = FALSE)

inds <- which(All$ID %in% c("Cl_0_E095","Cl_0_E125","Cl_0_E155","Cl_0_E185"))
ggplot(All[inds,], aes(x=ID, y=logV)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste(saveext,"/PhenotypeExp_AMes", ".pdf",sep=""),width = 10, height = 5,limitsize = FALSE)



inds <- which(All$ID %in% c("Cl_ESC_E095","Cl_ESC_E125","Cl_ESC_E155","Cl_ESC_E185"))
ggplot(All[inds,], aes(x=ID, y=logV)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste(saveext,"/PhenotypeExp_ESC", ".pdf",sep=""),width = 10, height = 5,limitsize = FALSE)


inds <- which(All$ID %in% c("Cl_PreME_E095","Cl_PreME_E125","Cl_PreME_E155","Cl_PreME_E185"))
ggplot(All[inds,], aes(x=ID, y=logV)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste(saveext,"/PhenotypeExp_PreME", ".pdf",sep=""),width = 10, height = 5,limitsize = FALSE)



inds <- which(All$ID %in% c("Cl_4_E095","Cl_4_E125","Cl_4_E155","Cl_4_E185"))
ggplot(All[inds,], aes(x=ID, y=logV)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste(saveext,"/PhenotypeExp_PGCLC", ".pdf",sep=""),width = 10, height = 5,limitsize = FALSE)


#Comparisoons vs
#Compare starting vs end populations?
Idents(mammal.combined) <- mammal.combined$species2
Data1 <- subset(mammal.combined, idents=c("EBPreMe","Gastruloid","Amnioids","EB4i","HumanCS7") )

uID <- as.character(Data1$Cl9)

#Have to run later parts to generat this
#PreMes <- WhichCells(Dat1, idents = "EBPreMe_PGCLC_12h_r2_8")

uID[which(Data1$Cl9=="2")] <- "PS"
uID[which(Data1$Cl9=="10")] <- "lPS"
#uID[which(Data1$Cl9=="10")] <- "lPS"
#uID[which(Data1$Cl9=="10")] <- "eAMLC"
uID[which(Data1$Cl9=="13")] <- "eAMLC"
uID[which(Data1$Cl9=="5")] <- "lAMLC"
uID[which(Data1$Cl9=="8")] <- "lAMLC"
uID[which(Data1$Cl9=="7")] <- "lAMLC"
uID[which(Data1$Cl9=="9")] <- "ePGCLC"
uID[which(Data1$Cl9=="4")] <- "lPGCLC"
uID[which(Data1$Cl9=="28")] <- "eDELC"
uID[which(Data1$Cl9=="17")] <- "eDELC"
uID[which(Data1$Cl9=="20")] <- "eDELC"
uID[which(Data1$Cl9=="15")] <- "lDELC"
uID[which(Data1$Cl9=="16")] <- "lDELC"
uID[which(Data1$Cl9=="6")] <- "eMELC"
uID[which(Data1$Cl9=="0")] <- "lMELC"
uID[which(Data1$Cl9=="18")] <- "ECLC"
uID[which(Data1$Cl9=="12")] <- "pECLC"
uID[which(Data1$ID2=="PreME")] <- "PreME"
uID[which(Data1$ID2=="ESC")] <- "ESC"
uID[which(Data1$ID2=="4i")] <- "4i"
uID[which(Data1$ID2=="ME")] <- "ME"
uID[which(Data1$ID2=="DE")] <- "DE"
uID[which(Data1$ID2=="PGC")] <- "PGC"

uID[which(Data1$ID3=="EPI")] <- "EPI"

uID[which(Data1$species2=="HumanCS7")] <- as.character(Data1$ID3[which(Data1$species2=="HumanCS7")])


uID2 <- paste(Data1$species2,uID,sep="_")

#Get the overlapping genes?
AmnioidGenes <-rownames(readRDS('/Users/christopherpenfold/Desktop/Aracely/GSE134571_Posterior48h_H9_Amnion_Merged'))
EBGenes<-read.table("/Users/christopherpenfold/Desktop/Aracely/Pairwise_PGCLC_vs_Amnioids/Plots/GenesMapping.csv",sep=",",header = T)
EBGenes<-EBGenes$Gene
GastGenes<-read.table("/Users/christopherpenfold/Desktop/Aracely/GastGenes.csv",sep=" ",header = T)
GastGenes<-GastGenes$x
commongenes <- intersect(intersect(AmnioidGenes,EBGenes),GastGenes)
  
DefaultAssay(Data1) <- "RNA"
Idents(Data1) <- uID2
#Idents(object = Data1, cells = PreMes) <- "PreMes"
Ae <- AverageExpression(Data1, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae) #[commongenes,]
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late MELC (PreME) - PreMe", y = "Late MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lMELC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late MELC (PreME) - PreMe", y = "Late MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lMELC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
#p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
#p1 <- p1 + labs(x = "Late MELC (PreME) - PreMe", y = "Late MELC (Amnioid) - ESC")
#ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lMELC_Lig2.pdf",sep=""),width = 13, height = 13, plot = p1)
#dev.off()
#p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
#p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
#p1 <- p1 + labs(x = "Late MELC (PreME) - PreMe", y = "Late MELC (Amnioid) - ESC")
#ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lMELC_Lig3.pdf",sep=""),width = 13, height = 13, plot = p1)
#dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_eMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_eMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPGCLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPGCLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Gastruloid_EPI", ident.1 = "Gastruloid_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Gastruloid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Gastruloid_EBPreME_lPGCLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Gastruloid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Gastruloid_EBPreME_lPGCLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Gastruloid_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Gastruloid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Gastruloid_EBPreME_lPGCLC_TFv2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (Gastruloid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Gastruloid_EBPreME_lPGCLC_Lig1v2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "HumanCS7_Epiblast", ident.1 = "HumanCS7_PGC_CS7", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "PGC (CS7) - Epiblast")
ggsave(filename=paste(saveext,"Volcano_CS7_EBPreME_lPGCLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "PGC (CS7) - Epiblast")
ggsave(filename=paste(saveext,"Volcano_CS7_EBPreME_lPGCLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_PGC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "HumanCS7_Epiblast", ident.1 = "HumanCS7_PGC_CS7", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "PGC (CS7) - Epiblast")
ggsave(filename=paste(saveext,"Volcano_CS7_EBPreME_PGCW4_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "PGC (CS7) - Epiblast")
ggsave(filename=paste(saveext,"Volcano_CS7_EBPreME_PGCW4_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "EBPreMe_4i", ident.1 = "EB4i_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (PreME) -4i")
ggsave(filename=paste(saveext,"Volcano_EB4i_EBPreME_lPGCLC_TFv2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - PreMe", y = "Late PGCLC (PreME) - 4i")
ggsave(filename=paste(saveext,"Volcano_EB4i_EBPreME_lPGCLC_Lig1v2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_ePGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_ePGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early PGCLC (PreME) - PreMe", y = "Early PGCLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_ePGCLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early PGCLC (PreME) - PreMe", y = "Early PGCLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_ePGCLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_eAMLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_eAMLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early AMLC (PreME) - PreMe", y = "Early AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eAMLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early AMLC (PreME) - PreMe", y = "Early AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eAMLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lAMLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_lAMLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late AMLC (PreME) - PreMe", y = "Late AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lAMLC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late AMLC (PreME) - PreMe", y = "Late AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lAMLC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_lPS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_lPS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PS (PreME) - PreMe", y = "Late AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPS_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PS (PreME) - PreMe", y = "Late PS (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPS_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_PS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_PS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early PS (PreME) - PreMe", y = "Early AMLC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_PS_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early PS (PreME) - PreMe", y = "Early PS (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_PS_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("Amnioids_lPGCLC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Amnioids_EBs_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Amnioids_EBs_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_EBs_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_EBs_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_Amnioids_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_Amnioids_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_Gastruloids_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_Gasrtuloids_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_4i_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCvPGC_Volcano_4i_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("HumanCS7_PGC_CS7"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"CS7PGCvPGC_Volcano_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"CS7PGCvPGC_Volcano_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("EB4i_lPGCLC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_4i_EBs_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_4i_EBs_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("Gastruloid_lPGCLC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Gastruloid_EBs_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Gastruloid_EBs_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("Gastruloid_lPGCLC"), ident.1 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Am_Gastruloid_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_Am_Gastruloid_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_EBs_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_EBs_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_4i_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_4is_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_Amnioid_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_Amnioid_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PreME"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl2 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PreME"), ident.1 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST",logfc.threshold = log(1.2), only.pos=TRUE)
Cl3 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PreME"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl4 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PreME"), ident.1 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST",  logfc.threshold = log(1.2), only.pos=TRUE)
Ints1 <- intersect(intersect(intersect(rownames(Cl4),rownames(Cl3)),rownames(Cl2)),rownames(Cl1))
Cl1 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PreME"), ident.2 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl2 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PreME"), ident.2 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST",logfc.threshold = log(1.2), only.pos=TRUE)
Cl3 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PreME"), ident.2 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl4 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PreME"), ident.2 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST",  logfc.threshold = log(1.2), only.pos=TRUE)
Ints2 <- intersect(intersect(intersect(rownames(Cl4),rownames(Cl3)),rownames(Cl2)),rownames(Cl1))

Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PreME"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(Ints1,Ints2)
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCLC_vs_preME_commonall.pdf",sep=""),width = 40, height = 40, plot = p1)
dev.off()


write.table(as.data.frame(Ints1),paste(saveext,"CommonUp_PGCLCvsPreME.scv",sep=""))
write.table(as.data.frame(Ints2),paste(saveext,"CommonDown_PGCLCvsPreME.scv",sep=""))

Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl2 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST",logfc.threshold = log(1.2), only.pos=TRUE)
Cl3 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl4 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST",  logfc.threshold = log(1.2), only.pos=TRUE)

write.table(as.data.frame(Cl1),paste(saveext,"UpGast.csv",sep=""))
write.table(as.data.frame(Cl2),paste(saveext,"UpAm.csv",sep=""))
write.table(as.data.frame(Cl3),paste(saveext,"UpEB.csv",sep=""))
write.table(as.data.frame(Cl4),paste(saveext,"Up4i.csv",sep=""))


Ints3 <- intersect(intersect(intersect(rownames(Cl4),rownames(Cl3)),rownames(Cl2)),rownames(Cl1))
Cl1 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PGC"), ident.2 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl2 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PGC"), ident.2 = c("Amnioids_lPGCLC"), verbose = FALSE, test.use = "MAST",logfc.threshold = log(1.2), only.pos=TRUE)
Cl3 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PGC"), ident.2 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", logfc.threshold = log(1.2), only.pos=TRUE)
Cl4 <- FindMarkers(Data1, ident.1 = c("EBPreMe_PGC"), ident.2 = c("EB4i_lPGCLC"), verbose = FALSE, test.use = "MAST",  logfc.threshold = log(1.2), only.pos=TRUE)
Ints4 <- intersect(intersect(intersect(rownames(Cl4),rownames(Cl3)),rownames(Cl2)),rownames(Cl1))

#http://www.interactivenn.net
write.table(as.data.frame(Cl1),paste(saveext,"DownGast.csv",sep=""))
write.table(as.data.frame(Cl2),paste(saveext,"DownAm.csv",sep=""))
write.table(as.data.frame(Cl3),paste(saveext,"DownEB.csv",sep=""))
write.table(as.data.frame(Cl4),paste(saveext,"Down4i.csv",sep=""))

write.table(as.data.frame(Ints3),paste(saveext,"CommonUp_PGCvsPGCLC.scv",sep=""))
write.table(as.data.frame(Ints4),paste(saveext,"CommonDown_PGCvsPGCLC.scv",sep=""))

Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGC"), ident.1 = c("EBPreMe_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(Ints1,Ints2)
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGCLC_vs_preME_commonall.pdf",sep=""),width = 40, height = 40, plot = p1)
dev.off()


#Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreME_PreMe"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))


Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))



Cl1 <- FindMarkers(Data1, ident.2 = c("HumanCS7_PGC_CS7"), ident.1 = c("Gastruloid_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_Gastruloid_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta lPGC", y = "-log Pval")
ggsave(filename=paste(saveext,"lPGC_Volcano_CS7_Gastruloid_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





#Cross comparisons
#Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_ESC", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl1 <- FindMarkers(Data1, ident.2 = "Amnioids_lPGCLC", ident.1 = "EBPreMe_lPGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_PreME", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = unique(c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1),
                           intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2),
                           intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3)))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - Late PGCLC Amnioid", y = "PreME vs ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPGCLConly_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Late PGCLC (PreME) - Late PGCLC Amnioid", y = "PreME vs ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_lPGCLConly_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PreME", ident.1 = "EBPreMe_eMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Cl2 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "Amnioids_eMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl3 <- FindMarkers(Data1, ident.2 = "Amnioids_lMELC", ident.1 = "EBPreMe_lMELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Cl4 <- FindMarkers(Data1, ident.2 = "Amnioids_ESC", ident.1 = "EBPreMe_ESC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$Pval2 <- 1
Ae5$FC1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl2),"Pval2"] <- -log2(Cl2$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"FC2"] <- Cl2$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.5) & Ae5$Pval1>-log2(0.01))
pospos2 <- which( abs(Ae5$FC2)>log2(1.5) & Ae5$Pval2>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi[pospos2] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos2],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos2],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos2],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos2],SIGNAL3))

p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_Lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_Lig2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,FC2)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC (PreME) - PreMe", y = "Early MELC (Amnioid) - ESC")
ggsave(filename=paste(saveext,"Volcano_Amnoid_EBPreME_eMELC_Lig3.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Idents(Data1) <- Data1$Cl9
#Idents(Data1) <- Data1$species2
#Data1b <- subset(Data1,idents="Gastruloid")
#Idents(Data1b) <- Data1b$Cl9

uID3 <- paste(Data1$species2,Data1$Cells,uID,sep="_")
Idents(Data1) <- uID3
Ae <- AverageExpression(Data1, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae) #[commongenes,]

#TFs and other factors
#18h : 6 vs 2, 6 vs 10?
#24h : 6 vs 2, 6 vs 10, 6 vs 9

#18hs
Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PGCLC_18h_eMELC", ident.1 = "EBPreMe_PGCLC_18h_PS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_18h_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_18h_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#24
Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_24h_eMELC"), ident.1 = "EBPreMe_PGCLC_24h_PS", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_TF_PS.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_Lig_PS.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#24
Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_24h_eMELC"), ident.1 = "EBPreMe_PGCLC_24h_ePGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_TF_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_Lig_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#24
Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_24h_eMELC"), ident.1 = "EBPreMe_PGCLC_24h_eDELC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - DELC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_TF_DE.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - DELC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_Lig_DE.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#24
Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_24h_eDELC"), ident.1 = "EBPreMe_PGCLC_24h_ePGCLC", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early DELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_TF_DEPGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early DELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_24h_Lig_DEPGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Subset the 12h hour set
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindClusters(mammal.combined, resolution = 2.5)
mammal.combined$Cl25 <- Idents(mammal.combined) # $Cells
Idents(mammal.combined) <- mammal.combined$Cells
Dat1 <- subset(mammal.combined,idents=c("PGCLC_12h_r2"))

uID4 <- paste(Dat1$species2,Dat1$Cells,Dat1$Cl25,sep="_")
Idents(Dat1) <- uID4
Ae <- AverageExpression(Dat1, verbose = FALSE)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae) #[commongenes,]

DefaultAssay(Dat1) <- "RNA"
Cl1 <- FindMarkers(Dat1, ident.2 = "EBPreMe_PGCLC_12h_r2_8", ident.1 = "EBPreMe_PGCLC_12h_r2_44", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_TF_44v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_Lig_44v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PGCLC_12h_r2_8", ident.1 = "EBPreMe_PGCLC_12h_r2_21", verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_TF_21v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_Lig_21v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

Cl1 <- FindMarkers(Data1, ident.2 = "EBPreMe_PGCLC_12h_r2_8", ident.1 = c("EBPreMe_PGCLC_12h_r2_21","EBPreMe_PGCLC_12h_r2_44"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_TF_21_44v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Early MELC - PS", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_12h_Lig_21_44v8.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#32 hours cross talk
Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_32h_eMELC","EBPreMe_PGCLC_32h_lMELC"), ident.1 = c("EBPreMe_PGCLC_32h_ePGCLC","EBPreMe_PGCLC_32h_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "MELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_TF_MELC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "MELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_Lig_MELC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_32h_eDELC"), ident.1 = c("EBPreMe_PGCLC_32h_ePGCLC","EBPreMe_PGCLC_32h_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "DELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_TF_DEPGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "DELC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_Lig_DEPGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_32h_mAMLC","EBPreMe_PGCLC_32h_lAMLC"), ident.1 = c("EBPreMe_PGCLC_32h_ePGCLC","EBPreMe_PGCLC_32h_lPGCLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "AMLC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_TF_AMLC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "AMLC - PGCLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_Lig_AMLC_PGCLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(Data1, ident.2 = c("EBPreMe_PGCLC_32h_eMELC","EBPreMe_PGCLC_32h_lMELC"), ident.1 = c("EBPreMe_PGCLC_32h_mAMLC","EBPreMe_PGCLC_32h_lAMLC"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log2(1.2) & Ae5$Pval1>-log2(0.01))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))

p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "MELC - AMLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_TF_MELC_AMLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5[commongenes,], aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "MELC - AMLC", y = "-log Pval")
ggsave(filename=paste(saveext,"Volcano_32h_Lig_MELC_AMLC.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

