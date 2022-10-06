library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(destiny)

set.seed(1)

saveext = "./Ara_UberAlign18G/"
fileID <- saveext
mammal.combined <- readRDS(paste(fileID,'mammal.combined.KO5.rds',sep=""))

ids <- Idents(mammal.combined)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(mammal.combined) <- "RNA"
mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes)
Idents(mammal.combined) <- ids
DefaultAssay(mammal.combined) <- "integrated"

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

mammal.combined$species2 <- Species2
Idents(mammal.combined) <- mammal.combined$species2

Data2 <- subset(mammal.combined, idents=c("1) Human EBs (in vitro)","HumanCS7","Cyno (in vitro)","7) Marmoset (in vivo)","HumanCS5-7","5) Human (in vitro 2)"))

uID1 <- as.character(Data2$Cl9)
uID2 <- as.character(Data2$Cells)
uID3 <- as.character(Data2$species2)
uID4 <- as.character(Data2$Cl9)

uID4[which(uID2=="Ectoderm")] <- "human_Am_CS7"
uID4[which(uID2=="Epiblast")] <- "human_Epiblast_CS7"
uID4[which(uID2=="Primitive Streak" & (uID1=="4" | uID1=="9" ) )] <- "human_PGC_CS7"
uID4[which(uID2=="PGC" & uID3=="HumanCS5-7")] <- "human_gPGC1"
uID4[which(uID2=="PGC" & uID3=="5) Human (in vitro 2)")] <- "human_gPGC2"

uID4[which(uID2=="Am_CS5" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="Am_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="Am_CS7" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="PGC_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="4" | uID1=="9" ) )] <- "marm_PGC_CS5-7"
uID4[which(uID2=="EmDisc_CS5" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS7" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"

uID4[which(uID2=="Am_CS5" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="Am_CS6" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="Am_CS7" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="PGC_CS6" & uID3=="Cyno (in vitro)" & (uID1=="4" | uID1=="9" ) )] <- "cyno_PGC_CS5-7"
uID4[which(uID2=="EmDisc_CS5" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS6" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS7" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID3=="1) Human EBs (in vitro)" & (uID1=="4" | uID1=="9" ) )] <- "PGCLC"
uID4[which(uID3=="1) Human EBs (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "AMLC"

uID4[which(uID2=="PGC_wk6.5_invivo")] <- "gPGC"
uID4[which(uID2=="hESC_W5_E8")] <- "hESC"
Idents(Data2) <- uID4

list1 <- c("PGCLC","gPGC",
           "marm_PGC_CS5-7","cyno_PGC_CS5-7","human_PGC_CS7","human_gPGC1","human_gPGC2")

Data3 <- subset(Data2,idents=list1)

subs <- Data3 #subset(Data1, idents=c("PreME","ESC","PGCLC","ME"), invert = FALSE)
DN1 <- as.data.frame(subs[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors
rownames(DM) <- colnames(subs)
subs[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(subs))
DimPlot(subs, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC",".pdf",sep=""),width = 12, height = 12)

uID1 <- as.character(Data2$Cl9)
uID2 <- as.character(Data2$Cells)
uID3 <- as.character(Data2$species2)
uID4 <- as.character(Data2$Cl9)
uID4[which(uID2=="Ectoderm")] <- "human_Am_CS7"
uID4[which(uID2=="Epiblast")] <- "human_Epiblast_CS7"
uID4[which(uID2=="Primitive Streak" & (uID1=="4" | uID1=="9" ) )] <- "human_PGC_CS7"
uID4[which(uID2=="PGC" & uID3=="HumanCS5-7" & (uID1=="4" | uID1=="9" ))] <- "human_gPGC1"
uID4[which(uID2=="PGC" & uID3=="5) Human (in vitro 2)" & (uID1=="4" | uID1=="9" ))] <- "human_gPGC2"

uID4[which(uID2=="Am_CS5" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="Am_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="Am_CS7" & uID3=="7) Marmoset (in vivo)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "marm_Am_CS5-7"
uID4[which(uID2=="PGC_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="4" | uID1=="9" ) )] <- "marm_PGC_CS5-7"
uID4[which(uID2=="EmDisc_CS5" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS6" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS7" & uID3=="7) Marmoset (in vivo)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "marm_Epiblast_CS5-7"

uID4[which(uID2=="Am_CS5" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="Am_CS6" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="Am_CS7" & uID3=="Cyno (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "cyno_Am_CS5-7"
uID4[which(uID2=="PGC_CS6" & uID3=="Cyno (in vitro)" & (uID1=="4" | uID1=="9" ) )] <- "cyno_PGC_CS5-7"
uID4[which(uID2=="EmDisc_CS5" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS6" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID2=="EmDisc_CS7" & uID3=="Cyno (in vitro)" & (uID1=="1" | uID1=="21" | uID1=="3" | uID1=="19" | uID1=="14" ) )] <- "cyno_Epiblast_CS5-7"
uID4[which(uID3=="1) Human EBs (in vitro)" & (uID1=="4" | uID1=="9" ) )] <- "PGCLC"
uID4[which(uID3=="1) Human EBs (in vitro)" & (uID1=="13" | uID1=="5" | uID1=="8" | uID1=="7" | uID1=="24" ) )] <- "AMLC"

uID4[which(uID2=="PGC_wk6.5_invivo" & (uID1=="4" | uID1=="9" ))] <- "gPGC"
uID4[which(uID2=="hESC_W5_E8")] <- "hESC"
Idents(Data2) <- uID4

Data3 <- subset(Data2,idents=list1)

subs <- Data3 #subset(Data1, idents=c("PreME","ESC","PGCLC","ME"), invert = FALSE)

DN1 <- as.data.frame(subs[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors
rownames(DM) <- colnames(subs)
subs[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(subs))
DimPlot(subs, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC2",".pdf",sep=""),width = 12, height = 12)

subs2 <- subset(subs,idents=c("PGCLC","gPGC"))

Idents(subs) <- subs$Cells
DimPlot(subs, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC2_lab",".pdf",sep=""),width = 12, height = 12)


DimPlot(subs, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE, dims = c(2,3)) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC2_2_3",".pdf",sep=""),width = 12, height = 12)

DimPlot(subs, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE, dims = c(1,3)) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC2_1_3",".pdf",sep=""),width = 12, height = 12)


Idents(subs2) <- subs2$Phase
DimPlot(subs2, reduction = "dm", pt.size = 2, label.size = 2,  label = TRUE, repel = TRUE, dims = c(1,2)) #+ NoLegend()
ggsave(filename=paste(saveext,"DM_Seurat_PGC_Phase",".pdf",sep=""),width = 12, height = 12)

Idents(Data2) <- Data2$Phase
DimPlot(Data2, reduction = "umap", pt.size = 2, label.size = 2, split.by="species2",  label = TRUE, repel = TRUE, dims = c(1,2)) #+ NoLegend()
ggsave(filename=paste(saveext,"ump_Seurat_PGC_Phase",".pdf",sep=""),width = 12*6, height = 12,limitsize = FALSE)
