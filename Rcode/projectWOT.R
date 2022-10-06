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
Species2[which(mammal.combined$Cells=="4i")] <- "EB4i"
Species2[which(mammal.combined$Cells=="ME_r2")] <- "EB4i"
Species2[which(mammal.combined$Cells=="PGC_wk6.5_invivo")] <- "EB4i"
mammal.combined$species2 <- Species2
Idents(mammal.combined) <- mammal.combined$species
mammal.combined <- subset(mammal.combined,idents = "Human (Ara EBs)")
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID

#saveRDS(mammal.paste(fileID,'subsetforDm.rds',sep=""))

#2,10,9,4,
#13,8,,5,7

Idents(mammal.combined) <- mammal.combined$species2
mammal.combined <- subset(mammal.combined, idents = "EBPreMe")
Idents(mammal.combined) <- mammal.combined$Cells
Idents(mammal.combined) <- mammal.combined$Cl9
D1 <- GetAssayData(mammal.combined, assay = "RNA")



saveext <- "/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOTPlots/"

#Load forward pseudotimes for PGCLC
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "PGCLC"
AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")
D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]
DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_PGCLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

baseID1 <- as.character(mammal.combined$Cl9)
baseID1[which(baseID1==2)] <- "PS1"
baseID1[which(baseID1==10)] <- "PS2"
baseID1[which(baseID1==9)] <- "PGCLC"
baseID1[which(baseID1==4)] <- "PGCLC"
baseID1[which(baseID1==6)] <- "eMELC"
baseID1[which(baseID1==0)] <- "aMELC"
baseID1[which(baseID1==13)] <- "eAMLC"
baseID1[which(baseID1==5)] <- "AMLC"
baseID1[which(baseID1==8)] <- "AMLC"
baseID1[which(baseID1==7)] <- "AMLC"

baseID1[which(baseID1==17)] <- "DELC1"
baseID1[which(baseID1==20)] <- "DELC1"
baseID1[which(baseID1==15)] <- "DELC2"
baseID1[which(baseID1==16)] <- "DELC2"

baseID1[which(baseID1==1)] <- "Basal1"
baseID1[which(baseID1==3)] <- "Basal2"


Idents(mammal.combined) <- baseID1
DimPlot(mammal.combined,   pt.size = 2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))


Idents(mammal.combined) <- mammal.combined$Cl9
DimPlot(mammal.combined,   pt.size = 2,  reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))

#Cl2 = basal, Cl10 = basal2, Cl9,4=PGCLC, Cl 6=MELC1, Cl0-AdMELC,, Cl 13,5,8,7=AMLC

ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl9",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Idents(mammal.combined) <- mammal.combined$Cl9
FM <- FindMarkers(mammal.combined, ident.1 = c("9","4") ,ident.2 = c("13","8","5","7") , test.use = "MAST")

######################
#First compare PGCLC vs Am
#Load PGCLC
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M1 <- merge(x = uID, y = Pt, by = 0, all.x = TRUE)
cellindexes <- which(mammal.combined$Cl9 %in% c(2,10,9,4)) #Subset by cluster also

#Now the amnion
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,10,13,8,5,7))

#
DEG <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/PGCLCdirect_DEgenes.tsv",sep="\t",header=T)
deg <- DEG$external_gene_name[which(abs(DEG$logFC)>log(1.5) )]
deg1 <- intersect(deg,TF)
deg2 <- intersect(deg,SIGNAL$V1)

time <- M1$`Mes1$dpt_pseudotime`[cellindexes]
keepcellindexes <- cellindexes[which(!is.na(time))]
time2 <- M1$`Mes1$dpt_pseudotime`[keepcellindexes]

DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/AMLCdirect_DEgenes.tsv",sep="\t",header=T)
deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)>log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1)) #This is TF
DifG2 <- unique(c(deg2_2,deg2)) #This is other stuff

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]

library(pheatmap)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap1_PGCLC_AMLC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale1_PGCLC_AMLC",".pdf",sep=""),width=5,height=7 )

pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap2_PGCLC_AMLC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale2_PGCLC_AMLC",".pdf",sep=""),width=5,height=15 )

#################################
#Reload the data - do we need to do these? This loads in desired genes
#PGCLC
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M1 <- merge(x = uID, y = Pt, by = 0, all.x = TRUE)
cellindexes <- which(mammal.combined$Cl9 %in% c(2,10,9,4))

#
#Now the amnion
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,10,13,8,5,7))




DEG <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/PGCLCdirect_DEgenes.tsv",sep="\t",header=T)
deg <- DEG$external_gene_name[which(abs(DEG$logFC)>log(1.5))]
deg1 <- intersect(deg,TF)
deg2 <- intersect(deg,SIGNAL$V1)

time <- M1$`Mes1$dpt_pseudotime`[cellindexes]
keepcellindexes <- cellindexes[which(!is.na(time))]
time2 <- M1$`Mes1$dpt_pseudotime`[keepcellindexes]

DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/AMLCdirect_DEgenes.tsv",sep="\t",header=T)
deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)>log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

#Get subset of genes for plotting
DifG1 <- unique(c(deg1_2,deg1,intersect(rownames(FM),TF), deg2_2,deg2,c("T", "GATA2", "MIXL1", "MESP1", "SNAI2", "PDGFRA" ,"TFAP2A", "GATA3", "VTCN1", "SOX17", "TFAP2C", "PRDM1", "ISL1", "HAND1","POU5F1", "SOX2", "NANOG","LHX1","OTX2","LEFTY", "TBXT","MESP1", "PDGFRA", "HAND1","FOXF1", "SNAI2", "PDGFRA",  "GATA6", "BMP4","FOXA2", "FOXA1", "PRDM1", "GATA4", "GATA6", "HNF4A","CER1","SOX3","PDPN","GABAPR","WNT3A","BMP2","SOX17","NANOS3","PRDM1", "CD38", "TFAP2A", "TFAP2C", "PDPN", "PRDM14","ITGB6","KRT7","KRT17")))
DifG1 <- intersect(rownames(mammal.combined@assays$RNA),DifG1)
D1B <- D1[DifG1,keepcellindexes]
D1B2 <- D1[DifG1,keepcellindexes2]
D1B <- D1B[,order(time2)]
D1B2 <- D1B2[,order(time3)]

write.csv(as.data.frame(D1B), file=paste(saveext,"/PGCs_pt.csv",sep=""))
write.csv(as.data.frame(D1B2), file=paste(saveext,"/Am_pt.csv",sep=""))

XX1 <- data.frame(x= D1B[c("SOX17"),] , y= time2, z = seq(1,length(time2) ), lab="SOX17" )
XX2 <- data.frame(x= D1B[c("NANOS3"),] , y= time2, z = seq(1,length(time2) ) , lab="NANOS3")
XX3 <- data.frame(x= D1B[c("PRDM1"),] , y= time2, z = seq(1,length(time2) ) , lab="PRDM1")
XX4 <- data.frame(x= D1B[c("CD38"),] , y= time2, z = seq(1,length(time2) ) , lab="CD38")
XX5 <- data.frame(x= D1B[c("TFAP2A"),] , y= time2, z = seq(1,length(time2) ) , lab="TFAP2A")
XX6 <- data.frame(x= D1B[c("TFAP2C"),] , y= time2, z = seq(1,length(time2) ) , lab="TFAP2C")
XX7 <- data.frame(x= D1B[c("PDPN"),] , y= time2, z = seq(1,length(time2) ) , lab="PDPN")
XX8 <- data.frame(x= D1B[c("PRDM14"),] , y= time2, z = seq(1,length(time2) ) , lab="PRDM14")
XX <- rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_Merge",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_Merge_1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX1 <- data.frame(x= D1B[c("VTCN1"),] , y= time2, z = seq(1,length(time2) ), lab="VTCN1" )
XX2 <- data.frame(x= D1B[c("GABRP"),] , y= time2, z = seq(1,length(time2) ) , lab="GABRP")
XX3 <- data.frame(x= D1B[c("ISL1"),] , y= time2, z = seq(1,length(time2) ) , lab="ISL1")
XX4 <- data.frame(x= D1B[c("ITGB6"),] , y= time2, z = seq(1,length(time2) ) , lab="ITGB6")
XX5 <- data.frame(x= D1B[c("KRT7"),] , y= time2, z = seq(1,length(time2) ) , lab="KRT7")
XX6 <- data.frame(x= D1B[c("KRT17"),] , y= time2, z = seq(1,length(time2) ) , lab="KRT17")
XX7 <- data.frame(x= D1B[c("GATA3"),] , y= time2, z = seq(1,length(time2) ) , lab="GATA3")
XX8 <- data.frame(x= D1B[c("GATA2"),] , y= time2, z = seq(1,length(time2) ) , lab="GATA2")
XX <- rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,2.5)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_MergeAm",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,2.5)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_Merge_Am1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX1 <- data.frame(x= D1B[c("SOX2"),] , y= time2, z = seq(1,length(time2) ), lab="SOX2" )
XX2 <- data.frame(x= D1B[c("POU5F1"),] , y= time2, z = seq(1,length(time2) ) , lab="POU5F1")
XX3 <- data.frame(x= D1B[c("NANOG"),] , y= time2, z = seq(1,length(time2) ) , lab="NANOG")
XX4 <- data.frame(x= D1B[c("KLF4"),] , y= time2, z = seq(1,length(time2) ) , lab="KLF4")
XX <- rbind(XX1,XX2,XX3,XX4)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,4)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_MergePl",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,4)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_Merge_Pl1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX1 <- data.frame(x= D1B2[c("SOX17"),] , y= time3, z = seq(1,length(time3) ), lab="SOX17" )
XX2 <- data.frame(x= D1B2[c("NANOS3"),] , y= time3, z = seq(1,length(time3) ) , lab="NANOS3")
XX3 <- data.frame(x= D1B2[c("PRDM1"),] , y= time3, z = seq(1,length(time3) ) , lab="PRDM1")
XX4 <- data.frame(x= D1B2[c("CD38"),] , y= time3, z = seq(1,length(time3) ) , lab="CD38")
XX5 <- data.frame(x= D1B2[c("TFAP2A"),] , y= time3, z = seq(1,length(time3) ) , lab="TFAP2A")
XX6 <- data.frame(x= D1B2[c("TFAP2C"),] , y= time3, z = seq(1,length(time3) ) , lab="TFAP2C")
XX7 <- data.frame(x= D1B2[c("PDPN"),] , y= time3, z = seq(1,length(time3) ) , lab="PDPN")
XX8 <- data.frame(x= D1B2[c("PRDM14"),] , y= time3, z = seq(1,length(time3) ) , lab="PRDM14")
XX <- rbind(XX1,XX2,XX3,XX4,XX5,XX5,XX7,XX8)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_Merge",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_Merge1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX1 <- data.frame(x= D1B2[c("VTCN1"),] , y= time3, z = seq(1,length(time3) ), lab="VTCN1" )
XX2 <- data.frame(x= D1B2[c("GABRP"),] , y= time3, z = seq(1,length(time3) ) , lab="GABRP")
XX3 <- data.frame(x= D1B2[c("ISL1"),] , y= time3, z = seq(1,length(time3) ) , lab="ISL1")
XX4 <- data.frame(x= D1B2[c("ITGB6"),] , y= time3, z = seq(1,length(time3) ) , lab="ITGB6")
XX5 <- data.frame(x= D1B2[c("KRT7"),] , y= time3, z = seq(1,length(time3) ) , lab="KRT7")
XX6 <- data.frame(x= D1B2[c("KRT17"),] , y= time3, z = seq(1,length(time3) ) , lab="KRT17")
XX7 <- data.frame(x= D1B2[c("GATA3"),] , y= time3, z = seq(1,length(time3) ) , lab="GATA3")
XX8 <- data.frame(x= D1B2[c("GATA2"),] , y= time3, z = seq(1,length(time3) ) , lab="GATA2")
XX <- rbind(XX1,XX2,XX3,XX4,XX5,XX5,XX7,XX8)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,2.5)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_MergeAm",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,2.5)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_MergeAm1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX1 <- data.frame(x= D1B2[c("SOX2"),] , y= time3, z = seq(1,length(time3) ), lab="SOX2" )
XX2 <- data.frame(x= D1B2[c("POU5F1"),] , y= time3, z = seq(1,length(time3) ) , lab="POU5F1")
XX3 <- data.frame(x= D1B2[c("NANOG"),] , y= time3, z = seq(1,length(time3) ) , lab="NANOG")
XX4 <- data.frame(x= D1B2[c("KLF4"),] , y= time3, z = seq(1,length(time3) ) , lab="KLF4")
XX <- rbind(XX1,XX2,XX3,XX4)
XX$ID <- as.factor(XX$lab)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,4)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_MergePl",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
ggplot(XX, aes(x = z, y = x, color = lab, group=lab)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,4)) 
ggsave(filename=paste(saveext,"/DimRed/AMLC_MergePl1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)




XX <- data.frame(x= D1B[c("SOX17","NANOS3","PRDM1", "CD38", "TFAP2A", "TFAP2C", "PDPN", "PRDM14"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.25,2))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SOX17_multiplots",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)



XX <- data.frame(x= D1B[c("SOX17"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.25,2))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SOX17",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("SOX2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.01,0.5)) #ylim(-0.01, 0.5)
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SOX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("TFAP2A"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.2,1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_TFAP2A",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("TFAP2C"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.2,1.8)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_TFAP2C",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B[c("GATA3"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_GATA3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B[c("ISL1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,2)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_ISL1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("PRDM1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,0.6)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_PRDM1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B[c("NANOG"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,2)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_NANOG",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("POU5F1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.3,4)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_POU5F1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("VTCN1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,0.8)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_VTCN1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
#XX <- data.frame(x= D1B[c("FOXA2"),] , y= time2, z = seq(1,length(time2) ) )
#ggplot(XX, aes(x = z, y = x)) + geom_point(shape = 1) + geom_smooth(color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + ylim(-0.2,2)
#ggsave(filename=paste(saveext,"/DimRed/PGCLC_FOXA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("SOX3"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SOX3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("GATA2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_GATA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("GATA4"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.75)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_GATA4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B[c("CER1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_CER1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("BMP2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()  + coord_cartesian(ylim=c(-0,0.25)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_BMP2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("WNT3A"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.25)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_WNT3A",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B2[c("SOX17"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.25,2))#+ coord_cartesian(ylim=c(-0.25,2))
ggsave(filename=paste(saveext,"/DimRed/AMCLC_SOX17",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("SOX2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.01,0.5))
ggsave(filename=paste(saveext,"/DimRed/AMCLC_SOX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("TFAP2A"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.2,1)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_TFAP2A",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("TFAP2C"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.2,1.8)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_TFAP2C",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA3"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_GATA3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("ISL1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_ISL1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("PRDM1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,0.6)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_PRDM1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("NANOG"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,2)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_NANOG",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("POU5F1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()  + coord_cartesian(ylim=c(-0.3,4)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_POU5F1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("VTCN1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.8)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_VTCN1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B2[c("FOXA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.3)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_FOXA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("SOX3"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,0.35))
ggsave(filename=paste(saveext,"/DimRed/AMCLC_SOX3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_GATA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA4"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_GATA4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("CER1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_CER1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("BMP2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.25)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_BMP2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("WNT3A"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,0.25)) 
ggsave(filename=paste(saveext,"/DimRed/AMCLC_WNT3A",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)






XX <- data.frame(x= D1B[c("LHX1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_LHX1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("OTX2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_OTX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("LEFTY2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,1.2)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_LEFTY2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("TBXT"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,2.5)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_TBXT",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("MESP1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.1,3.1))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_MESP1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("PDGFRA"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.5)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_PDGFRA",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("HAND1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,3)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_HAND1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("SNAI2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SNAI2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("GATA6"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_GATA6",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("BMP4"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.4)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_BMP4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)



###############################################################
#Now PGCLC vs MELC
#Now the PGCLC
FM <- FindMarkers(mammal.combined, ident.1 = c("9","4") ,ident.2 = c("6","0") , test.use = "MAST")
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_MELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,6,0))

#MELC
DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/MELC_DEgenes.tsv",sep="\t",header=T)

deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)>log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1)) #This is TF
DifG2 <- unique(c(deg2_2,deg2)) #This is other stuff

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]

#DifG1 <- unique(c(deg1_2,deg1,intersect(rownames(FM),TF), deg2_2,deg2,c("T", "GATA6", "MIXL1", "MESP1", "SNAI2", "PDGFRA" ,"TFAP2A", "GATA3", "VTCN1", "SOX17", "TFAP2C", "PRDM1", "ISL1", "HAND1","POU5F1", "SOX2", "NANOG","LHX1","OTX2","LEFTY", "TBXT","MESP1", "PDGFRA", "HAND1","FOXF1", "SNAI2", "PDGFRA",  "GATA6", "BMP4","FOXA2", "FOXA1", "PRDM1", "GATA4", "GATA6", "HNF4A")))
#DifG1 <- intersect(rownames(mammal.combined@assays$RNA),DifG1)
#D1B2 <- D1[DifG1,keepcellindexes2]
#D1B2 <- D1B2[,order(time3)]

write.csv(as.data.frame(D1B2), file=paste(saveext,"/MELC_pt.csv",sep=""))


#Do heatmaps?
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap1_PGCLC_MELC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale1_PGCLC_MELC",".pdf",sep=""),width=5,height=7 )

pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap2_PGCLC_MELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale2_PGCLC_MELC",".pdf",sep=""),width=5,height=15 )


###############################################################
#Now PGCLC vs MELC - reload for line lpots
#Now the PGCLC
FM <- FindMarkers(mammal.combined, ident.1 = c("9","4") ,ident.2 = c("6","0") , test.use = "MAST")
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_MELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,6,0))


#MELC
DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/MELC_DEgenes.tsv",sep="\t",header=T)


deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)>log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]


#Get subset of genes for plotting
DifG1 <- unique(c(deg1_2,deg1,intersect(rownames(FM),TF), deg2_2,deg2,c("T", "GATA2", "MIXL1", "MESP1", "SNAI2", "PDGFRA" ,"TFAP2A", "GATA3", "VTCN1", "SOX17", "TFAP2C", "PRDM1", "ISL1", "HAND1","POU5F1", "SOX2", "NANOG","LHX1","OTX2","LEFTY", "TBXT","MESP1", "PDGFRA", "HAND1","FOXF1", "SNAI2", "PDGFRA",  "GATA6", "BMP4","FOXA2", "FOXA1", "PRDM1", "GATA4", "GATA6", "HNF4A","CER1","SOX3","PDPN","LEFTY2","GABAPR","WNT3A","BMP2","MEST","GSC","ZIC3","ZIC5")))
DifG1 <- intersect(rownames(mammal.combined@assays$RNA),DifG1)
D1B <- D1[DifG1,keepcellindexes]
D1B2 <- D1[DifG1,keepcellindexes2]
D1B <- D1B[,order(time2)]
D1B2 <- D1B2[,order(time3)]


XX1 <- data.frame(x= D1B2[c("TBXT"),] , y= time3, z = seq(1,length(time3) ), ID="T" )
XX2 <- data.frame(x= D1B2[c("SNAI2"),] , y= time3, z = seq(1,length(time3) ), ID="SNAI2" )
XX3 <- data.frame(x= D1B2[c("PDGFRA"),] , y= time3, z = seq(1,length(time3) ), ID="PDGFRA" )
XX4 <- data.frame(x= D1B2[c("MIXL1"),] , y= time3, z = seq(1,length(time3) ), ID="MIXL1" )
XX5 <- data.frame(x= D1B2[c("HAND1"),] , y= time3, z = seq(1,length(time3) ) , ID="HAND1")
XX6 <- data.frame(x= D1B2[c("MESP1"),] , y= time3, z = seq(1,length(time3) ) , ID="MESP1")
XX7 <- data.frame(x= D1B2[c("GATA4"),] , y= time3, z = seq(1,length(time3) ) , ID="GATA4")
XX7 <- data.frame(x= D1B2[c("GATA6"),] , y= time3, z = seq(1,length(time3) ) , ID="GATA6")
XX8 <- data.frame(x= D1B2[c("BMP4"),] , y= time3, z = seq(1,length(time3) ) , ID="BMP4")
XX9 <- data.frame(x= D1B2[c("LHX1"),] , y= time3, z = seq(1,length(time3) ) , ID="LHX1")
XX10 <- data.frame(x= D1B2[c("EOMES"),] , y= time3, z = seq(1,length(time3) ) , ID="EOMES")
XX11 <- data.frame(x= D1B2[c("FOXF1"),] , y= time3, z = seq(1,length(time3) ) , ID="FOXF1")
XX12 <- data.frame(x= D1B2[c("FOXH1"),] , y= time3, z = seq(1,length(time3) ) , ID="FOXH1")
XX13 <- data.frame(x= D1B2[c("OTX2"),] , y= time3, z = seq(1,length(time3) ) , ID="OTTX2")
XX14 <- data.frame(x= D1B2[c("MEST"),] , y= time3, z = seq(1,length(time3) ) , ID="MEST")
XX15 <- data.frame(x= D1B2[c("ALCAM"),] , y= time3, z = seq(1,length(time3) ) , ID="ALCAM")
XX16 <- data.frame(x= D1B2[c("GSC"),] , y= time3, z = seq(1,length(time3) ) , ID="GSC")
XX17 <- data.frame(x= D1B2[c("CDX2"),] , y= time3, z = seq(1,length(time3) ) , ID="CDX2")
XX18 <- data.frame(x= D1B2[c("ZIC2"),] , y= time3, z = seq(1,length(time3) ) , ID="ZIC2")
XX19 <- data.frame(x= D1B2[c("ZIC3"),] , y= time3, z = seq(1,length(time3) ) , ID="zic3")
XX20 <- data.frame(x= D1B2[c("ZIC5"),] , y= time3, z = seq(1,length(time3) ) , ID="ZIC5")

XX <- rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX10,XX11,XX12,XX13,XX14,XX15,XX16,XX17,XX18,XX19,XX20)
XX$ID <- as.factor(XX$ID)
#ggplot(XX, aes(x=z, y=x)) + geom_point()
#XX <- data.frame(t (D1B2[c("T","SNAI2","PDGFRA","MIXL1","HAND1","MESP1","GATA6","BMP4","LHX1"),] ) )
#XX$time <- time3
#XX$Ind <-  seq(1,length(time3) )
#library(reshape2)
#YY <- melt(XX)

ggplot(XX, aes(x = z, y = x, color = ID, group=ID)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_Merge",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


ggplot(XX, aes(x = z, y = x, color = ID, group=ID)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_Merge1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


#library(scales)
#create exp(x)-1 transformation, the inverse of log(1+p)
#expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
#
#qplot(x, y, data=sites) + stat_smooth(method="loess") +
#  scale_y_continuous(trans=log1p_trans()) +
#  coord_trans(ytrans=expm1_trans())

#ggplot(XX, aes(x = z, y = x, color = ID, group=ID)) + geom_smooth(method = 'loess', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,3)) 
#ggsave(filename=paste(saveext,"/DimRed/MECLC_Merge2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


#geom_smooth(color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 



XX <- data.frame(x= D1B2[c("CER1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35)) 
ggsave(filename=paste(saveext,"/DimRed/MELC_CER1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)


XX <- data.frame(x= D1B2[c("LHX1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,1))
ggsave(filename=paste(saveext,"/DimRed/MECLC_LHX1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("OTX2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_OTX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("LEFTY2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.2)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_LEFTY2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("TBXT"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35)+ coord_cartesian(ylim=c(-0.05,1.2)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_TBXT",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("MESP1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.1,3.1)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_MESP1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

library(scales)
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")

XX <- data.frame(x= D1B2[c("MESP1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.1,2.5)) 


+ ylim(0, 30000) #+ scale_y_continuous(trans=log1p_trans()) + coord_trans(ytrans=expm1_trans())

ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'glm', family='binomial', family='quasibinomial', color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.1,2.5)) + scale_y_log10()+coord_trans(ytrans="pow10") # + scale_y_continuous(trans=log1p_trans()) + coord_trans(ytrans=expm1_trans())



XX <- data.frame(x= D1B2[c("PDGFRA"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.5))
ggsave(filename=paste(saveext,"/DimRed/MECLC_PDGFRA",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("HAND1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,3))
ggsave(filename=paste(saveext,"/DimRed/MECLC_HAND1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("SNAI2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.3))
ggsave(filename=paste(saveext,"/DimRed/MECLC_SNAI2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA6"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1))
ggsave(filename=paste(saveext,"/DimRed/MECLC_GATA6",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("BMP4"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.4))
ggsave(filename=paste(saveext,"/DimRed/MECLC_BMP4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)




XX <- data.frame(x= D1B2[c("SOX3"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35))
ggsave(filename=paste(saveext,"/DimRed/MELC_SOX3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1))
ggsave(filename=paste(saveext,"/DimRed/MELC_GATA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA4"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.75)) 
ggsave(filename=paste(saveext,"/DimRed/MELC_GATA4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("LHX1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,1))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_LHX1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("OTX2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_OTX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("LEFTY2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.2))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_LEFTY2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("TBXT"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,1.2)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_TBXT",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("MESP1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.1,3.1)) 
ggsave(filename=paste(saveext,"/DimRed/PGCLC_MESP1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("PDGFRA"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.5))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_PDGFRA",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("HAND1"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,3))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_HAND1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("SNAI2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.3))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_SNAI2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("GATA6"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_GATA6",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("BMP4"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1.4))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_BMP4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B[c("FOXA2"),] , y= time2, z = seq(1,length(time2) ) )
ggplot(XX, aes(x=z, y=x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.3))
ggsave(filename=paste(saveext,"/DimRed/PGCLC_FOXA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("SOX17"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.25,2))
ggsave(filename=paste(saveext,"/DimRed/MECLC_SOX17",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("SOX2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.01,0.5)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_SOX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("TFAP2A"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.2,1)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_TFAP2A",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("TFAP2C"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.2,1.8)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_TFAP2C",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("GATA3"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_GATA3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("ISL1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,1.3)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_ISL1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("PRDM1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.6)) 
ggsave(filename=paste(saveext,"/DimRed/MECLC_PRDM1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("NANOG"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,2))
ggsave(filename=paste(saveext,"/DimRed/MECLC_NANOG",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("POU5F1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.3,4))
ggsave(filename=paste(saveext,"/DimRed/MECLC_POU5F1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("VTCN1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,0.8)) 
ggsave(filename=paste(saveext,"/DimRed/MELC_VTCN1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("FOXA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.3))
ggsave(filename=paste(saveext,"/DimRed/MELC_FOXA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)



#Now DELC?
FM <- FindMarkers(mammal.combined, ident.1 = c("9","4") ,ident.2 = c("17","20","28") , test.use = "MAST")

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_DELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,20,17,28))

DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/DELC_DEgenes.tsv",sep="\t",header=T)

deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)>log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1)) #This is TF
DifG2 <- unique(c(deg2_2,deg2)) #This is other stuff

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]

#DifG1 <- unique(c(deg1_2,deg1,intersect(rownames(FM),TF), deg2_2,deg2,c("T", "GATA6", "MIXL1", "MESP1", "SNAI2", "PDGFRA" ,"TFAP2A", "GATA3", "VTCN1", "SOX17", "TFAP2C", "PRDM1", "ISL1", "HAND1","POU5F1", "SOX2", "NANOG","LHX1","OTX2","LEFTY", "TBXT","MESP1", "PDGFRA", "HAND1","FOXF1", "SNAI2", "PDGFRA",  "GATA6", "BMP4","FOXA2", "FOXA1", "PRDM1", "GATA4", "GATA6", "HNF4A")))
#DifG1 <- intersect(rownames(mammal.combined@assays$RNA),DifG1)
#D1B2 <- D1[DifG1,keepcellindexes2]
#D1B2 <- D1B2[,order(time3)]

write.csv(as.data.frame(D1B2), file=paste(saveext,"/DELC_pt.csv",sep=""))


#Do heatmaps?
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap1_PGCLC_DLC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale1_PGCLC_DELC",".pdf",sep=""),width=5,height=7 )

pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap2_PGCLC_DELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale2_PGCLC_DELC",".pdf",sep=""),width=5,height=15 )

#Now reload and do line plots
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_DELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,20,17,28))




DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/DELC_DEgenes.tsv",sep="\t",header=T)

deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)> log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1,intersect(rownames(FM),TF), deg2_2,deg2,c("T", "GATA6", "MIXL1", "MESP1", "SNAI2", "PDGFRA" ,"TFAP2A", "GATA3", "VTCN1", "SOX17", "TFAP2C", "PRDM1", "ISL1", "HAND1","POU5F1", "SOX2", "NANOG","LHX1","OTX2","LEFTY", "TBXT","MESP1", "PDGFRA", "HAND1","FOXF1", "SNAI2", "PDGFRA",  "GATA6", "BMP4","FOXA2", "FOXA1", "PRDM1", "GATA4", "GATA6", "HNF4A","FOXA2")))
DifG1 <- intersect(rownames(mammal.combined@assays$RNA),DifG1)
D1B2 <- D1[DifG1,keepcellindexes2]
D1B2 <- D1B2[,order(time3)]



XX1 <- data.frame(x= D1B2[c("EOMES"),] , y= time3, z = seq(1,length(time3) ), ID="EOMES" )
XX2 <- data.frame(x= D1B2[c("PRDM1"),] , y= time3, z = seq(1,length(time3) ), ID="PRDM1" )
XX3 <- data.frame(x= D1B2[c("SOX17"),] , y= time3, z = seq(1,length(time3) ), ID="SOX17" )
XX4 <- data.frame(x= D1B2[c("FOXA1"),] , y= time3, z = seq(1,length(time3) ), ID="FOXA1" )
XX5 <- data.frame(x= D1B2[c("FOXA2"),] , y= time3, z = seq(1,length(time3) ) , ID="FOXA2")
XX <- rbind(XX1,XX2,XX3,XX4,XX5)
XX$ID <- as.factor(XX$ID)
ggplot(XX, aes(x = z, y = x, color = ID, group=ID)) + geom_smooth(method = 'gam', se=TRUE, level = 0.99, span=0.35, fullrange=TRUE) + theme_bw()+ coord_cartesian(ylim=c(-0.6,2)) 
ggsave(filename=paste(saveext,"/DimRed/DELC_Merge",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)




#XX <- data.frame(x= D1B2[c("VTCN1"),] , y= time3, z = seq(1,length(time3) ) )
dd1 <- readRDS("/Users/christopherpenfold/Desktop/Ara/AraDM.rds")



XX <- data.frame(x= D1B2[c("FOXA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.3)) 
ggsave(filename=paste(saveext,"/DimRed/DELC_FOXA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("VTCN1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0.05,0.8)) 
ggsave(filename=paste(saveext,"/DimRed/DELC_VTCN1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("SOX17"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x))+ geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.25,2))
ggsave(filename=paste(saveext,"/DimRed/DELC_SOX17",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("SOX2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35)  + theme_bw()+ coord_cartesian(ylim=c(-0.01,0.5))
ggsave(filename=paste(saveext,"/DimRed/DELC_SOX2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("SOX3"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.05,0.35))
ggsave(filename=paste(saveext,"/DimRed/DELC_SOX3",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA2"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x)) + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw()+ coord_cartesian(ylim=c(-0,1)) 
ggsave(filename=paste(saveext,"/DimRed/DECLC_GATA2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)
XX <- data.frame(x= D1B2[c("GATA4"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x = z, y = x))  + geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0,0.75)) 
ggsave(filename=paste(saveext,"/DimRed/DECLC_GATA4",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)

XX <- data.frame(x= D1B2[c("POU5F1"),] , y= time3, z = seq(1,length(time3) ) )
ggplot(XX, aes(x=z, y=x))+ geom_smooth(method = 'loess',color = "blue", fill = "lightblue", se = TRUE, level = 0.99, span=0.35) + theme_bw() + coord_cartesian(ylim=c(-0.3,4)) 
ggsave(filename=paste(saveext,"/DimRed/DELC_POU5F1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE, limitsize = FALSE)



write.csv(as.data.frame(D1B2), file=paste(saveext,"/DELC_pt.csv",sep=""))


sadsadsadasdasd



#
#2,10,9,4,
#13,8,,5,7
#D1 <- GetAssayData(mammal.combined, assay = "RNA")

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M1 <- merge(x = uID, y = Pt, by = 0, all.x = TRUE)
cellindexes <- which(mammal.combined$Cl9 %in% c(2,10,9,4))

#Now the amnion
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_MELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,6,0))

#
DEG <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/PGCLCdirect_DEgenes.tsv",sep="\t",header=T)
deg <- DEG$external_gene_name[which(abs(DEG$logFC)> log(1.5) )]
deg1 <- intersect(deg,TF)
deg2 <- intersect(deg,SIGNAL$V1)

time <- M1$`Mes1$dpt_pseudotime`[cellindexes]
keepcellindexes <- cellindexes[which(!is.na(time))]
time2 <- M1$`Mes1$dpt_pseudotime`[keepcellindexes]

DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/MELC_DEgenes.tsv",sep="\t",header=T)
deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)> log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1))
DifG2 <- unique(c(deg2_2,deg2))

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]


redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

#annotation_col = data.frame(Stage = factor(ano))
#rownames(annotation_col) <- colnames(AA)

#names(mycolors) <- colnames(X)
#anno_colors <- list(Stage = mycolors)


pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap1_PGCLC_MELC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale1_PGCLC_MELC",".pdf",sep=""),width=5,height=7 )


pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap2_PGCLC_MELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale2_PGCLC_MELC",".pdf",sep=""),width=5,height=15 )


deg1 <- intersect(rownames(FM),TF)
deg2 <- intersect(rownames(FM),SIGNAL$V1)

DifG1 <- unique(c(deg1))
DifG2 <- unique(c(deg2))

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap3_PGCLC_MELC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale3_PGCLC_MELC",".pdf",sep=""),width=5,height=7 )

pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap4_PGCLC_MELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale4_PGCLC_MELC",".pdf",sep=""),width=5,height=15 )







Idents(mammal.combined) <- mammal.combined$Cl9
FM <- FindMarkers(mammal.combined, ident.1 = c("9","4") ,ident.2 = c("17","20") , test.use = "MAST")

#
#2,10,9,4,
#13,8,,5,7
#D1 <- GetAssayData(mammal.combined, assay = "RNA")

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M1 <- merge(x = uID, y = Pt, by = 0, all.x = TRUE)
cellindexes <- which(mammal.combined$Cl9 %in% c(2,10,9,4))

#Now the amnion
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_DELC_pseudotime.tsv",header = T)
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- as.numeric(-100000)
uID <- as.data.frame(uID)
rownames(uID) <- colnames(mammal.combined)
Pt2 <- as.data.frame(Mes1$dpt_pseudotime)
rownames(Pt2) <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
M2 <- merge(x = uID, y = Pt2, by = 0, all.x = TRUE)
cellindexes2 <- which(mammal.combined$Cl9 %in% c(2,17,20))

#
DEG <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/PGCLCdirect_DEgenes.tsv",sep="\t",header=T)
deg <- DEG$external_gene_name[which(abs(DEG$logFC)> log(1.5) )]
deg1 <- intersect(deg,TF)
deg2 <- intersect(deg,SIGNAL$V1)

time <- M1$`Mes1$dpt_pseudotime`[cellindexes]
keepcellindexes <- cellindexes[which(!is.na(time))]
time2 <- M1$`Mes1$dpt_pseudotime`[keepcellindexes]

DEG2 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/DGE/DELC_DEgenes.tsv",sep="\t",header=T)
deg_2 <- DEG2$external_gene_name[which(abs(DEG2$logFC)> log(1.5) )]
deg1_2 <- intersect(deg_2,TF)
deg2_2 <- intersect(deg_2,SIGNAL$V1)

time_2 <- M2$`Mes1$dpt_pseudotime`[cellindexes2]
keepcellindexes2 <- cellindexes2[which(!is.na(time_2))]
time3 <- M2$`Mes1$dpt_pseudotime`[keepcellindexes2]

DifG1 <- unique(c(deg1_2,deg1))
DifG2 <- unique(c(deg2_2,deg2))

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]


redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

#annotation_col = data.frame(Stage = factor(ano))
#rownames(annotation_col) <- colnames(AA)

#names(mycolors) <- colnames(X)
#anno_colors <- list(Stage = mycolors)


pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap1_PGCLC_DELC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale1_PGCLC_DELC",".pdf",sep=""),width=5,height=7 )


pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap2_PGCLC_DELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale2_PGCLC_DELC",".pdf",sep=""),width=5,height=15 )


deg1 <- intersect(rownames(FM),TF)
deg2 <- intersect(rownames(FM),SIGNAL$V1)

DifG1 <- unique(c(deg1))
DifG2 <- unique(c(deg2))

D1B <- D1[DifG1,keepcellindexes]
D1C <- D1[DifG2,keepcellindexes]

D1B2 <- D1[DifG1,keepcellindexes2]
D1C2 <- D1[DifG2,keepcellindexes2]

D1B <- D1B[,order(time2)]
D1C <- D1C[,order(time2)]

D1B2 <- D1B2[,order(-time3)]
D1C2 <- D1C2[,order(-time3)]

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

AA <- cbind(D1B2,D1B)
BB <- cbind(D1C2,D1C)

ano <- c(rep('Am', dim(D1B2)[2]),rep('PGCLC', dim(D1B)[2]))

pheatmap(AA,color =  redblue1(320), border_color = NA, gaps_col=c(dim(D1B2)[2]), cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap3_PGCLC_DELC",".pdf",sep="") ,width=5,height=7)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(AA,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale3_PGCLC_DELC",".pdf",sep=""),width=5,height=7 )

pheatmap(BB,color =  redblue1(320), gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Fullheamtap4_PGCLC_DELC",".pdf",sep="") ,width=5,height=15)
mat_breaks <- seq(-2, 2, length.out = 30)
pheatmap(BB,color =redblue1(30), breaks=mat_breaks, gaps_col=c(dim(D1B2)[2]), border_color = NA, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/DimRed/Fullheamtapscale4_PGCLC_DELC",".pdf",sep=""),width=5,height=15 )








#Rest of the code ...
#Load forward pseudotimes
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "AMLC"
uID <- as.character(Idents(mammal.combined))
inds <- which(uID=="AMLC" & mammal.combined$Cl9 %in% c(2,10,13,8,5,7))
uID[1:length(uID)] <- "Other"
uID[inds] <- "AMLC"
Idents(mammal.combined) <- uID

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_AMLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)


Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "PGCLC"
uID <- as.character(Idents(mammal.combined))
inds <- which(uID=="PGCLC" & mammal.combined$Cl9 %in%  c(2,10,9,4) )
uID[1:length(uID)] <- "Other"
uID[inds] <- "PGCLC"
Idents(mammal.combined) <- uID


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_PGCLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)


Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_MELC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "MELC"
uID <- as.character(Idents(mammal.combined))
inds <- which(uID=="MELC" & mammal.combined$Cl9 %in% c(2,6,0))
uID[1:length(uID)] <- "Other"
uID[inds] <- "MELC"
Idents(mammal.combined) <- uID

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_MELC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)






#Standard plots
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_EarlyMELCTrajectory_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
Idents(mammal.combined,cells = wcells) <- "MELC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_MELC1",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_LateMELCTrajectory_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "MELC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_MELC2",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_AMLCTrajectory_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "AMLC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_AMLC",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_DELCTrajectory_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "DELC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_DELC",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_PGCLCTrajectory_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "PGCLC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_PGCLC",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)









Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "AMLC"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_AMLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)


Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_MELC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "MELC"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_MELC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)




Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_Basal_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "Basal"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Basal_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)




Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_Basal_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "Basal"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Basal_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_DELC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "DELC"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_DELC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)





Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_PGCLC-AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "PGCLC"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_PGCLC-AMLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)





Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_Basal_avecME-sansESC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "MELC1"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_avecME_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)




Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_Basal_sansESC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "Basal"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Basal_ESC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)



 
 






Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Pseudotime/PGCLC_Terminal_AMLC_pseudotime.tsv",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')
uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"
Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells) <- "AMLC"


AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","#5C5C5C")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_AMLC_forward",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)





Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Annotations/WOT_PGCLC_lineages.txt",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')

ind1 <- which(Mes1$Lineage=="AMLC")# & mammal.combined$Cl9 %in% c(2,10,13,8,5,7))
ind2 <- which(Mes1$Lineage=="Basal") 
ind3 <- which(Mes1$Lineage=="Branch")# & mammal.combined$Cl9 %in% c(2) ) 
ind4 <- which(Mes1$Lineage=="DELC")# & mammal.combined$Cl9 %in% c(2,20,17,28)) 
ind4b <- which(Mes1$Lineage=="PGCLC")# & mammal.combined$Cl9 %in% c(2,20,17,28))
#ind5 <- which(Mes1$Lineage=="Early.MELC")
ind6 <- which(Mes1$Lineage=="MELC")# & mammal.combined$Cl9 %in% c(2,6,0))
ind7 <- which(Mes1$Lineage=="PGCLC")# & mammal.combined$Cl9 %in% c(2,10,9,4))

uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"

Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells[ind1]) <- "AMLC"
Idents(mammal.combined,cells = wcells[ind2]) <- "Basal"
Idents(mammal.combined,cells = wcells[ind3]) <- "Branch"
Idents(mammal.combined,cells = wcells[ind4]) <- "DELC"
Idents(mammal.combined,cells = wcells[ind4b]) <- "DELC"
#Idents(mammal.combined,cells = wcells[ind5]) <- "MELC1"
Idents(mammal.combined,cells = wcells[ind6]) <- "MELC"
Idents(mammal.combined,cells = wcells[ind7]) <- "PGCLC"

uIDb <- Idents(mammal.combined)
mammal.combined$Lineage <-uIDb

ind1 <- which(mammal.combined$Lineage=="AMLC" & mammal.combined$Cl9 %in% c(2,10,13,8,5,7))
ind2 <- which(mammal.combined$Lineage=="Basal") 
ind3 <- which(mammal.combined$Lineage=="Branch" & mammal.combined$Cl9 %in% c(2) ) 
ind4 <- which(mammal.combined$Lineage=="DELC" & mammal.combined$Cl9 %in% c(2,20,17,28)) 
ind4b <- which(mammal.combined$Lineage=="PGCLC" & mammal.combined$Cl9 %in% c(2,20,17,28))
#ind5 <- which(Mes1$Lineage=="Early.MELC")
ind6 <- which(mammal.combined$Lineage=="MELC" & mammal.combined$Cl9 %in% c(2,6,0))
ind7 <- which(mammal.combined$Lineage=="PGCLC" & mammal.combined$Cl9 %in% c(2,10,9,4))


uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"

Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = colnames(mammal.combined)[ind1]) <- "AMLC"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind2]) <- "Basal"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind3]) <- "Branch"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind4]) <- "DELC"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind4b]) <- "DELC"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind6]) <- "MELC"
Idents(mammal.combined,cells = colnames(mammal.combined)[ind7]) <- "PGCLC"

Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Annotations/WOT_PGCLC_lineages.txt",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')

ind1 <- which(Mes1$Lineage=="AMLC") #& mammal.combined$Cl9 %in% c(2,10,13,8,5,7))
ind2 <- which(Mes1$Lineage=="Basal") 
ind3 <- which(Mes1$Lineage=="Branch") #& mammal.combined$Cl9 %in% c(2) ) 
ind4 <- which(Mes1$Lineage=="DELC") #& mammal.combined$Cl9 %in% c(2,20,17,28)) 
ind5 <- which(Mes1$Lineage=="Early.MELC")
ind6 <- which(Mes1$Lineage=="MELC") #& mammal.combined$Cl9 %in% c(2,6,0))
ind7 <- which(Mes1$Lineage=="PGCLC") # & mammal.combined$Cl9 %in% c(2,10,9,4))

uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"

Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells[ind1]) <- "AMLC"
Idents(mammal.combined,cells = wcells[ind2]) <- "Basal"
Idents(mammal.combined,cells = wcells[ind3]) <- "Branch"
Idents(mammal.combined,cells = wcells[ind4]) <- "DELC"
Idents(mammal.combined,cells = wcells[ind5]) <- "MELC1"
Idents(mammal.combined,cells = wcells[ind6]) <- "MELC"
Idents(mammal.combined,cells = wcells[ind7]) <- "PGCLC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal","Branch")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","darkgrey","blue")


D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]



DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_WOT_times1",".pdf",sep=""),width = 17, height = 15, useDingbats=FALSE,limitsize = FALSE)

Idents(mammal.combined) <- mammal.combined$Cl9
DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste("~/Desktop/UMAP_WOT_Cl",".pdf",sep=""),width = 17, height = 15, useDingbats=FALSE,limitsize = FALSE)



uID2 <- paste(Idents(mammal.combined),mammal.combined$Cells,sep="_")

Idents(mammal.combined) <- uID2
AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal","Branch")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","darkgrey","blue")

AmType <- c("DELC_PGCLC_32h", "DELC_PGCLC_40h" , "DELC_PGCLC_48h",
  "Other_DE","Other_PGCLC_12h_r2","Other_PGCLC_18h",   
            "Other_PGCLC_24h","Other_PGCLC_32h","Other_PGCLC_40h","Other_PGCLC_48h","Other_PGCLC_D4","Other_PreME_10h", 
"AMLC_PGCLC_24h",
"AMLC_PGCLC_32h",
"AMLC_PGCLC_40h",
"AMLC_PGCLC_48h",
"AMLC_PGCLC_D4",
"Basal_DE",
"Basal_hESC_W5_E8",
"Basal_PGCLC_12h_r2",
"Basal_PGCLC_18h",
"Basal_PreME_10h",
"Branch_PGCLC_12h_r2",
"Branch_PGCLC_18h",
"Branch_PGCLC_24h",
"Branch_PGCLC_32h",
"Branch_PreME_10h",
"DELC_DE",
"DELC_PGCLC_D4",
"MELC_PGCLC_12h_r2",
"MELC_PGCLC_18h",
"MELC_PGCLC_24h",
"MELC_PGCLC_32h",
"MELC_PGCLC_40h",
"MELC_PGCLC_48h",
"MELC_PGCLC_D4",
"PGCLC_DE",
"PGCLC_PGCLC_12h_r2",
"PGCLC_PGCLC_32h",
"PGCLC_PGCLC_40h",
"PGCLC_PGCLC_48h",
"PGCLC_PGCLC_D4")

AmCol <- c("#737373",
  "#525252",
  "#252525",
  "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey",
           "#66c2a4",
"#41ae76",
"#238b45",
"#006d2c",
"#00441b",
"#828282",
"#662506",
"#d4b9da",
"#c994c7",
"#df65b0",
"#d4b9da",
"#c994c7",
"#df65b0",
"#e7298a",
"#cc4c02",
"#686868",
"#000000",
"#d0d1e6",
"#a6bddb",
"#74a9cf",
"#3690c0",
"#0570b0",
"#045a8d",
"#023858",
"#1A1A1A",
"#d4b9da",
"#e7298a",
"#ce1256",
"#980043",
"#67001f")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_WOT_times",".pdf",sep=""),width = 17, height = 15, useDingbats=FALSE,limitsize = FALSE)



DimPlot(mammal.combined,   pt.size = 2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_WOT_times",".pdf",sep=""),width = 17, height = 15, useDingbats=FALSE,limitsize = FALSE)


DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste("~/Desktop/UMAP_WOT",".pdf",sep=""),width = 17, height = 15, useDingbats=FALSE,limitsize = FALSE)




Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/ReverseTrajectories/PGCLC_Terminal_trajectoryLineages.txt",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')

ind1 <- which(Mes1$AMLC==TRUE)
ind2 <- which(Mes1$Not.Lineagee==TRUE) 
#ind3 <- which(Mes1$Lineage=="Branch") 
ind4 <- which(Mes1$DELC==TRUE) 
ind5 <- which(Mes1$Early.MELC==TRUE)
ind6 <- which(Mes1$Late.MELC==TRUE)
ind7 <- which(Mes1$PGCLC==TRUE)

uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"

Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells[ind1]) <- "AMLC"
#Idents(mammal.combined,cells = wcells[ind2]) <- "Basal"
#Idents(mammal.combined,cells = wcells[ind3]) <- "Branch"
Idents(mammal.combined,cells = wcells[ind4]) <- "DELC"
Idents(mammal.combined,cells = wcells[ind5]) <- "MELC1"
Idents(mammal.combined,cells = wcells[ind6]) <- "MELC"
Idents(mammal.combined,cells = wcells[ind7]) <- "PGCLC"

AmType <- c("Other","ESC","PGCLC","AMLC","MELC1","MELC","DELC","Basal","Branch")
AmCol <- c("lightgrey","#ED5810","#ED1E5D","#07971F","#0D8BD7","#0C81C7","#5C5C5C","darkgrey","blue")

D2<- Idents(mammal.combined) #$ID2
D2 <- droplevels(D2)
colind <- integer( length( levels(D2) )  )
for (i in 1:length( levels(D2) ) ) {
  colind[i] <- which(AmType==levels(D2)[i])
}
coluse1 <- AmCol[colind]

#DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_WOT",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)

DimPlot(mammal.combined,   pt.size = 2, cols = coluse1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + xlim(c(-13, 13)) + ylim(c(-13, 13))
ggsave(filename=paste(saveext,"/UMAP_WOT_rev",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE,limitsize = FALSE)


#Need to load WOT here first
Mes1 <-read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/WOT\ final\ calculations/Annotations/WOT_PGCLC_lineages.txt",header = T)
wcells <- str_replace(as.character(Mes1$Sample),'-1', '-1_7')

ind1 <- which(Mes1$Lineage=="AMLC") #& mammal.combined$Cl9 %in% c(2,10,13,8,5,7))
ind2 <- which(Mes1$Lineage=="Basal") 
ind3 <- which(Mes1$Lineage=="Branch") #& mammal.combined$Cl9 %in% c(2) ) 
ind4 <- which(Mes1$Lineage=="DELC") #& mammal.combined$Cl9 %in% c(2,20,17,28)) 
ind5 <- which(Mes1$Lineage=="Early.MELC")
ind6 <- which(Mes1$Lineage=="MELC") #& mammal.combined$Cl9 %in% c(2,6,0))
ind7 <- which(Mes1$Lineage=="PGCLC") # & mammal.combined$Cl9 %in% c(2,10,9,4))

uID <- as.character(mammal.combined$species)
uID[1:length(uID)] <- "Other"

Idents(mammal.combined) <- uID
Idents(mammal.combined,cells = wcells[ind1]) <- "AMLC"
Idents(mammal.combined,cells = wcells[ind2]) <- "Basal"
Idents(mammal.combined,cells = wcells[ind3]) <- "Branch"
Idents(mammal.combined,cells = wcells[ind4]) <- "DELC"
Idents(mammal.combined,cells = wcells[ind5]) <- "MELC1"
Idents(mammal.combined,cells = wcells[ind6]) <- "MELC"
Idents(mammal.combined,cells = wcells[ind7]) <- "PGCLC"
Av1 <- AverageExpression(mammal.combined)

baseID1 <- as.character(mammal.combined$Cl9)
baseID1[which(baseID1==2)] <- "PS1"
baseID1[which(baseID1==10)] <- "PS2"
baseID1[which(baseID1==9)] <- "PGCLC"
baseID1[which(baseID1==4)] <- "PGCLC"
baseID1[which(baseID1==6)] <- "eMELC"
baseID1[which(baseID1==0)] <- "aMELC"
baseID1[which(baseID1==13)] <- "eAMLC"
baseID1[which(baseID1==5)] <- "AMLC"
baseID1[which(baseID1==8)] <- "AMLC"
baseID1[which(baseID1==7)] <- "AMLC"
baseID1[which(baseID1==17)] <- "DELC1"
baseID1[which(baseID1==20)] <- "DELC1"
baseID1[which(baseID1==15)] <- "DELC2"
baseID1[which(baseID1==16)] <- "DELC2"
baseID1[which(baseID1==1)] <- "Basal1"
baseID1[which(baseID1==3)] <- "Basal2"
Idents(mammal.combined) <- baseID1
Av2 <- AverageExpression(mammal.combined)


#intgenes <- rownames(avexp$RNA)
a <- log2(Av1$integrated[,c("Basal","Branch","MELC","DELC","AMLC","PGCLC")])
b <- log2(Av2$integrated[,c("Basal1","Basal2","PS1","PS2","eMELC","aMELC","DELC1","DELC2","eAMLC","AMLC","PGCLC")])

C1 <- cor(a,b, method = "pearson",use="complete.obs")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 2, length.out = 20)
pheatmap((C1),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/WOT-vs-anno.pdf",,width=16,height=10)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, scale="row", filename = "~/Desktop/HM_ara_logged_scaled.pdf",width=10,height=16)


baseID1 <- as.character(mammal.combined$Cl9)
baseID1[which(baseID1==2)] <- "PS"
baseID1[which(baseID1==10)] <- "PS"
baseID1[which(baseID1==9)] <- "PGCLC"
baseID1[which(baseID1==4)] <- "PGCLC"
baseID1[which(baseID1==6)] <- "MELC"
baseID1[which(baseID1==0)] <- "MELC"
baseID1[which(baseID1==13)] <- "AMLC"
baseID1[which(baseID1==5)] <- "AMLC"
baseID1[which(baseID1==8)] <- "AMLC"
baseID1[which(baseID1==7)] <- "AMLC"
baseID1[which(baseID1==17)] <- "DELC"
baseID1[which(baseID1==20)] <- "DELC"
baseID1[which(baseID1==15)] <- "DELC"
baseID1[which(baseID1==16)] <- "DELC"
baseID1[which(baseID1==1)] <- "Basal"
baseID1[which(baseID1==3)] <- "Basal"
Idents(mammal.combined) <- baseID1
Idents(mammal.combined) <- baseID1
Av3 <- AverageExpression(mammal.combined)

a <- log2(Av1$integrated[,c("Basal","Branch","MELC","DELC","AMLC","PGCLC")])
b <- log2(Av3$integrated[,c("Basal","PS","MELC","DELC","AMLC","PGCLC")])

C1 <- cor(a,b, method = "pearson",use="complete.obs")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 1, length.out = 20)
pheatmap((C1),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/WOT-vs-anno2.pdf",,width=12,height=10)


