library(ggplot2)
library(gg3D)
library(plotly)
library(lme4)
library(geosphere)
library(reshape2)
library(heatmap3)
library(ComplexHeatmap)
library(colorRamp2)
library(dplyr)
library(car)
library(eks)
library(colorspace)
library(sf)
library(ggspatial)
library(maps)
library(cowplot)

###Heterozygosity vs the distance######
het <- read.table("./maizeEU.all.filteredQual_exp_het", header=TRUE, sep=",")
head(het)

#give a general description for the two groups
het$East <- het$Longitude > 2.35
het$North <- het$Latitude > 45.43

#distance from specific places. 
#spanish etry point : Tuy
latS <- het[30,2]
lonS <- het[30,3]
tmp <- NULL
for (i in 1: nrow(het)){
  tmp <- c(tmp, distm(c(het[i,3], het[i,2]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
tmp
het$DistSpain <- tmp

#german entry point: Barisis
latG <- het[2,2]
lonG <- het[2,3]
tmp <- NULL
for (i in 1: nrow(het)){
  tmp <- c(tmp, distm(c(het[i,3], het[i,2]), c(lonG, latG), fun = distHaversine) / 1000)  #transform them to km
}
tmp
het$DistGerman <- tmp

cor.test(het$DistSpain, het$mean.expected.heterozygosity, method="spearman")
cor.test(het$DistGerman, het$mean.expected.heterozygosity, method="spearman")

cor.test(het[het$East == FALSE,5], het[het$East == FALSE,8], method="spearman") #spain, within spain
cor.test(het[het$East == TRUE,5], het[het$East == TRUE,8], method="spearman") #spain, within CEu

cor.test(het[het$East == FALSE,5], het[het$East == FALSE,9], method="spearman") #german, within spain
cor.test(het[het$East == TRUE,5], het[het$East == TRUE,9], method="spearman") #german, within CEu


tiff("heter_entryPointSpain.tiff", units="in", width=5, height=5, res=300)
ggplot(het, aes(x=DistSpain, y=mean.expected.heterozygosity, col=East)) + geom_point(size=3, show.legend = F) + theme_bw() + geom_smooth(method = "lm",  fill="NA", show.legend = F) + ylab("Expected Heterozygosity") + xlab("Distance from entry point in Spain (km)") #theme(text=element_text(size=15))
dev.off()

tiff("heter_entryPointNorth.tiff", units="in", width=5, height=5, res=300)
ggplot(het, aes(x=DistGerman, y=mean.expected.heterozygosity, col=East)) + geom_point(size=3, show.legend = F) + theme_bw() + geom_smooth(method = "lm", fill="NA", show.legend = F) + ylab("Expected Heterozygosity") +  xlab("Distance from entry point in Central Europe (km)") # theme(text=element_text(size=15))
dev.off()

###Plot PCA####
pcaA <- read.csv("PCA/PCA/maizeEU.all.filteredQual_pca_without_DH.csv")
head(pcaA)
metadata <- read.csv("./master_thesis_scripts copy/metadata_for_PCA.csv")
metadata2 <- metadata 
metadata2 <- metadata2[-c(14,20,29),]


tiff("pca_maizeEU.noDH.tiff", units="in", width=5, height=5, res=300)
ggplot(pcaA, aes(pc1, pc2, col=PopID, shape=Type)) + geom_point(alpha=0.5, size=3, show.legend=FALSE) + theme_bw() +  scale_color_manual(values = metadata2$hex_col) + xlab("PC1\n4.3%") + ylab("PC2\n2.8%") + theme(text=element_text(size=13))
dev.off()

tiff("pca_maizeEU.noDH.withLegend.tiff", units="in", width=20, height=10, res=300)
ggplot(pcaA, aes(pc1, pc2, col=PopID, shape=Type)) + geom_point(alpha=0.5, size=3) + theme_bw() +  theme(legend.position = "bottom", legend.text = element_text(size=20),) + guides(color=guide_legend(nrow=2,byrow=TRUE), shape=guide_legend(nrow=2,byrow=T)) +
  scale_color_manual(values = metadata2$hex_col) + xlab("PC1\n4.3%") + ylab("PC2\n2.8%") + theme(text=element_text(size=15))
dev.off()

pcaA$Type <- as.factor(pcaA$Type)
pcaA$Type <- factor(pcaA$Type, c("landrace", "elite"))

###Plot map populations#####
library(rnaturalearth)
library(sf)
library(ggplot2)
library(ggrepel)
library(ggalt)

world <- ne_countries(scale = 'medium', returnclass = 'sf')
theme_set(theme_bw())

tiff("map_maizeEU.tiff",units="in", width=10, height=10, res=300)
ggplot(world) + 
  geom_sf(color="ivory2", fill="grey80")  +
  coord_sf(xlim=c(-15,25), ylim=c(35,57)) +
  theme(legend.position = "none") + 
  geom_point(data=metadata2, aes(x=Longitude, y=Latitude, col=Abb.), size=2, alpha=0.8) + 
  geom_text_repel(data=metadata2, aes(x=Longitude, y=Latitude,label=Abb.,  col=Abb.),max.overlaps=25) +
  scale_color_manual(name=NULL, values=metadata2$hex_col) +
  theme(text=element_text(size=20)) +
  geom_vline(xintercept=2.3514, col="darkblue")
dev.off()  

##
#####selscan result#####
#
setwd("./vcf")
filelist <- list.files(pattern="ihs.out")
selchrPE <- NULL
selchrKE <- NULL
for (file in filelist){
  if (length(strsplit(file, "[.]")[[1]]) > 7){
    tmp <- read.table(file, header=F)
    if (grepl("KE", file, fixed = TRUE)){tmp$Chr <- strsplit(file, "[.]")[[1]][4]; selchrKE <- rbind(selchrKE, tmp)}
    else if (grepl("PE", file, fixed = TRUE)){tmp$Chr <- strsplit(file, "[.]")[[1]][4]; selchrPE <- rbind(selchrPE, tmp)}
  } 
}
colnames(selchrPE) <- c("PhysicalPos", "geneticPos", "FreqAlt", "ihh1", "ihh0", "unstandiHS", "iHS", "Selection", "Chr")
summary(selchrPE)
levels(as.factor(selchrPE$Chr))
selchrPE$PosPlot <- 1:nrow(selchrPE)

colnames(selchrKE) <- c("PhysicalPos", "geneticPos", "FreqAlt", "ihh1", "ihh0", "unstandiHS", "iHS", "Selection", "Chr")
summary(selchrKE)
selchrKE$Chr <- as.factor(selchrKE$Chr)
selchrKE$PosPlot <- 1:nrow(selchrKE)

hist(abs(selchrPE$unstandiHS))

#empirical cutoff for both 
selchrPE$pval <- rank(-abs(selchrPE$iHS))/length(selchrPE$iHS)
hist(selchrPE$pval)
nrow(selchrPE[selchrPE$pval < 0.01,])
nrow(selchrPE[abs(selchrPE$iHS) > 2,])
nrow(selchrPE[selchrPE$pval < 0.01 & abs(selchrPE$iHS) > 2,])
selchrPE$significant <- selchrPE$PhysicalPos %in% selchrPE[selchrPE$pval < 0.01 & abs(selchrPE$iHS) > 2,1]

tmp <- selchrPE[selchrPE$pval < 0.01 & abs(selchrPE$iHS) > 2,]
summary(abs(tmp$iHS))

selchrKE$pval <- rank(-abs(selchrKE$iHS))/length(selchrKE$iHS)
hist(selchrKE$pval)
nrow(selchrKE[selchrKE$pval < 0.01,])
nrow(selchrKE[abs(selchrKE$iHS) > 2,])
nrow(selchrKE[selchrKE$pval < 0.01 & abs(selchrKE$iHS) > 2,])
selchrKE$significant <- selchrKE$PhysicalPos %in% selchrKE[selchrKE$pval < 0.01 & abs(selchrKE$iHS) > 2,1]

selchrKE$pval.adj <- p.adjust(selchrKE$pval, method = "bonferroni", n = length(selchrKE$pval))
nrow(selchrKE[selchrKE$pval.adj < 0.01,])

tmp <- selchrKE[selchrKE$pval < 0.01 & abs(selchrKE$iHS) > 2,]
summary(abs(tmp$iHS))

selchrPE$Chr = factor(selchrPE$Chr, levels=c("chr1", "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", "chr10" ))
levels(selchrPE$Chr) <- c("1", "2",  "3", "4",  "5",  "6",  "7",  "8",  "9", "10" )
tiff("maizeEU.DH.PE.iHS.tiff", units="in", width=10, height=5, res=600)
ggplot(selchrPE, aes(x=geneticPos/1000000, y=abs(iHS), col=significant)) + geom_point(alpha=0.5, size=3, show.legend = FALSE) + theme_bw() + geom_hline(yintercept = 3.298) + facet_grid( ~ Chr, scales="free_x", space = "free_x") + xlab("Position (Mbp)") + scale_color_manual(values = c("lightgrey", "navyblue")) + ylab("iHS") + theme(text=element_text(size=15), axis.text.x=element_text(size=11, angle=45))
dev.off()

selchrKE$Chr = factor(selchrKE$Chr, levels=c("chr1", "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", "chr10" ))
levels(selchrKE$Chr) <- c("1", "2",  "3", "4",  "5",  "6",  "7",  "8",  "9", "10" )
tiff("maizeEU.DH.KE.iHS.tiff", units="in", width=10, height=5, res=600)
ggplot(selchrKE %>% arrange(significant), aes(x=geneticPos/1000000, y=abs(iHS), col=significant)) + geom_point(alpha=0.5, size=3, show.legend = FALSE) + theme_bw() + geom_hline(yintercept = 3.747)  + facet_grid( ~ Chr, scales="free_x", space = "free_x") + xlab("Position (Mbp)") + scale_color_manual(values = c("lightgrey", "navyblue")) + ylab("iHS") + theme(text=element_text(size=15), axis.text.x=element_text(size=11, angle=45))
dev.off()

###
####Selscan Genes GOs Plot#####
setwd("/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/")
genesKE <- read.table("iHS.KE.closest.gene")
head(genesKE)
nrow(genesKE)

#check the genes that are closest to the SNP
genesKE <- genesKE[genesKE$V6 == "gene",]
nrow(genesKE)
tmp <- genesKE[!duplicated(genesKE[c(1,2,3)]),]     
nrow(tmp)
genesKE <- tmp

tmp <- NULL
for (i in genesKE$V12){
  tmp2 <- strsplit(i, ";")[[1]][1]
  tmp2 <- strsplit(tmp2, ":")[[1]][2]
  tmp2 <- strsplit(tmp2, "_")[[1]][1]
  tmp <- c(tmp, tmp2)
}
tmp <- unique(tmp) 
length(tmp)

write.table(tmp, file="genes.KE.selected.wind", col.names = F, row.names = F)

##Load the results from AgriGo GO enrichment
setwd("../")
gos.pe <- read.table("gos.PE.sign", header=T, sep="\t")
head(gos.pe)

ggdata <- gos.pe
ggdata$GO.Term <- factor(ggdata$GO.Term, levels = rev(ggdata$GO.Term))  # Fix order

# Plotting using ggplot2
gg1 <- ggplot(ggdata,
              aes(x = Description, y = -log10(p.value), size = -log10(p.value), fill = -log10(p.value))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5, 12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'PE DH',
    subtitle = 'Terms Ordered',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
  ) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),colour = c("black", "black", "black"),size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12), axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'), axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),legend.key.size = unit(1, "cm"),legend.text = element_text(size = 16),
    title = element_text(size = 12)) + coord_flip()

# Display and save the plot
print(gg1)
ggsave("GO_results_pe.png", plot = gg1, width = 11, height = 10, dpi = 300)

#PE
genesPE <- read.table("iHS.PE.closest.gene")
head(genesPE)
nrow(genesPE)

genesPE <- genesPE[genesPE$V6 == "gene",]
nrow(genesPE)
tmp <- genesPE[!duplicated(genesPE[c(1,2,3)]),]     
nrow(tmp)
genesPE <- tmp

tmp <- NULL
for (i in genesPE$V12){
  tmp2 <- strsplit(i, ";")[[1]][1]
  tmp2 <- strsplit(tmp2, ":")[[1]][2]
  tmp2 <- strsplit(tmp2, "_")[[1]][1]
  tmp <- c(tmp, tmp2)
}
tmp <- unique(tmp) 
length(tmp)
genesPE <- tmp

write.table(genesPE, file="genes.PE.selected.wind", col.names = F, row.names = F)

#plot for PE
gos.kl <- read.table("gos.KE.sign", header=T, sep="\t")
head(gos.kl)

ggdata <- gos.kl
ggdata$GO.Term <- factor(ggdata$GO.Term, levels = rev(ggdata$GO.Term))  # Fix order

# Plotting using ggplot2
gg1 <- ggplot(ggdata,
              aes(x = Description, y = -log10(p.value), size = -log10(p.value), fill = -log10(p.value))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5, 12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'KL DH',
    subtitle = 'Terms Ordered',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
  ) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),colour = c("black", "black", "black"),size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12), axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'), axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),legend.key.size = unit(1, "cm"),legend.text = element_text(size = 16),
    title = element_text(size = 12)) + coord_flip()

# Display and save the plot
print(gg1)
ggsave("GO_results_kl.png", plot = gg1, width = 11, height = 10, dpi = 300)

###
######Admixture#####
setwd("/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/")
metadata <- read.csv("./master_thesis_scripts copy/metadata_for_PCA.csv")
metadata2 <- metadata 
metadata2 <- metadata2[-c(14,20,29),]

k35 <- read.table("./admixture.maizeEU/admixture.maizeEU/admixture_all_k35.csv", sep=";", header = TRUE)
k35 <- melt(k35)

#merge with the metadata
k35 <- merge(k35, metadata, by.x="Site", by.y="Abb.", all.x=TRUE)

#reoorder factors based on Fst
poporder <- c("MD", "RT", "RM", "PF", "SO", "SF", "RO", "WA", "GB", "CO", "MB", "OM", "PE", "KL", "KR", "PL", "ND", "TR", "CA", "LC", "BU","BA", "MO", "AN", "SA", "VI","TU","ML","RD","GA","LD","LL","OE","KN","GL")

k35$Site <- as.factor(k35$Site)
k35$Site <- factor(k35$Site, levels = poporder)

tiff("maizeEU.k35.tiff", units="in", width=12, height=5, res=600)
ggplot() + geom_col(data= k35, width=1, mapping = aes(y=value, x=reorder(Individual ,Site), fill=variable), show.legend = FALSE) + theme_bw() + xlab("") + ylab("") + theme(text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +  scale_fill_manual(values = metadata2$hex_col) + geom_vline(xintercept = clusters, col="yellow", linetype="dashed")
dev.off()

####
######Fst#####
fst <- read.table("fsts_maizeEU.txt")
head(fst)

#fix the mistake in the printing
tmp <- NULL
for (i in fst$V2){
  tmp <- rbind(tmp, strsplit(i, "[\\]")[[1]])
}
tmp
fst <- cbind(fst, tmp)
colnames(fst) <- c("Pop1", "ce", "SE", "Pop2", "Fst")
fst$Fst <- as.numeric(as.character(fst$Fst))
#replace the negative values with 0s
fst[fst < 0] <- 0
summary(fst)

#fix the names of the populations
tmp <- NULL
for (i in fst$Pop1){
  tmp <- rbind(tmp, strsplit(i, "_")[[1]][1])
}
fst$Pop1 <- tmp

tmp <- NULL
for (i in fst$Pop2){
  tmp <- rbind(tmp, strsplit(i, "_")[[1]][1])
}
fst$Pop2 <- tmp

hist(fst$Fst)
ggplot(fst, aes(x=Pop1,y=Pop2,fill=Fst)) + geom_tile()
ggplot(fst, aes(x=Pop1, y=Fst)) + geom_point()

length(levels(as.factor(fst$Pop1)))
#create the matrix for the heatmap
fst2 <- fst[,c("Pop1", "Pop2", "Fst")]
fst2 <-acast(fst2, Pop1 ~ Pop2)

#fill up the reverse
for (i in 1:ncol(fst2)){
  fst2[c((i+1):ncol(fst2)),i] <- fst2[i,c((i+1):ncol(fst2))] 
}
#fst2[upper.tri(fst2)]=NA
heatmap3(fst2, symm=T, cexRow = 0.4, cexCol= 0.4, na.rm=F)
col2 = colorRamp2(c(0,0.25,0.6), c( "darkblue", "white","darkred"))

pdf("maizeEU_Fst.pdf")
Heatmap(fst2, rect_gp = gpar(type = "none"), column_dend_side = "bottom",col=col2, border=F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        })
dev.off()

tiff("maizeEU.Fst.tiff", units="in", width=7, height=7, res=600)
Heatmap(fst2, rect_gp = gpar(type = "none"), column_dend_side = "bottom",col=col2, border=F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        })
dev.off()

#compare with distance/cluster
metadata <- read.csv("./master_thesis_scripts copy/metadata_for_PCA.csv")
head(metadata)

fst3 <- merge(fst, metadata, by.x="Pop1", by.y="Abb.")
nrow(fst)
fst3 <- merge(fst3, metadata, by.x="Pop2", by.y="Abb.")
nrow(fst3)
colnames(fst3)
fst3 <- fst3[,c("Pop2", "Pop1", "Fst", "Latitude.x", "Longitude.x", "Latitude.y", "Longitude.y")]
head(fst3)

#estimate distance
tmp <- NULL
for (i in 1:nrow(fst3)){
  tmp <- c(tmp, distm(c(fst3[i,5], fst3[i,4]), c(fst3[i,7], fst3[i,6]), fun = distHaversine) / 1000)  #transform them to km
}
fst3$Distance <- tmp
ggplot(fst3, aes(x=Distance, y=Fst)) + geom_point()

#group by cluster
fst3$East.x <- fst3$Longitude.x > 2.35
fst3$East.y <- fst3$Longitude.y > 2.35
fst3$Cluster <- fst3$East.y==fst3$East.x
#check if the pair is in the same cluster or not
tiff("maizeEU.Fst.Distance.tiff", units="in", width=7, height=7, res=600)
ggplot(fst3, aes(x=Distance, y=Fst, col=Cluster)) + geom_point(size=3, show.legend = F) + theme_bw() + geom_smooth(method = "lm",  fill="NA", show.legend = F) + ylab("Fst") + xlab("Distance between the populations (km)") + theme(text=element_text(size=20))
dev.off()

tiff("maizeEU.Fst.Density.tiff", units="in", width=7, height=7, res=600)
ggplot(fst3, aes(x=Fst, fill=Cluster)) + geom_density(alpha=0.7, show.legend = F) + 
  xlab("Fst") + theme_bw() + theme(text=element_text(size=20))
dev.off()

cor.test(fst3$Distance, fst3$Fst, method="spearman")
cor.test(fst3[fst3$East.x == FALSE & fst3$East.y == FALSE,3], fst3[fst3$East.x == FALSE & fst3$East.y == FALSE,8], method="spearman") #west, within west
cor.test(fst3[fst3$East.x == TRUE & fst3$East.y == TRUE,3], fst3[fst3$East.x == TRUE & fst3$East.y == TRUE,8], method="spearman") #east within east
cor.test(fst3[fst3$Cluster == TRUE,3], fst3[fst3$Cluster == TRUE,8], method="spearman") #within cluster
cor.test(fst3[fst3$Cluster == FALSE,3], fst3[fst3$Cluster == F,8], method="spearman") #between cluster

summary(fst3[fst3$East.x == FALSE & fst3$East.y == FALSE,3]) #within west
summary(fst3[fst3$East.x == TRUE & fst3$East.y == TRUE,3]) #within east
summary(fst3[fst3$East.x == TRUE,3]) #within cluster
summary(fst3[fst3$Cluster == FALSE,3]) #outside cluster
ks.test(fst3[fst3$Cluster == FALSE,3], fst3[fst3$East.x == TRUE,3])

#highest differentiation of populations
tapply(fst3$Fst, fst3$Pop2, median)
tapply(fst3$Fst, fst3$Pop1, mean)
summary(fst3[fst3$Pop2 == "PL" | fst3$Pop1 == "PL",3])

##Slatkin's Distance
fst3$SlatkinDistance <- (fst3$Fst / (1 - fst3$Fst))
hist(fst3$SlatkinDistance)
ggplot(fst3, aes(x=Fst, y=SlatkinDistance)) + geom_point() + geom_smooth(method="lm")

tiff("maizeEU.slatkinsFst.tiff", units="in", width=7, height=7, res=600)
ggplot(fst3, aes(x=Distance, y=SlatkinDistance)) + geom_point() + geom_smooth(method="lm") + theme(text=element_text(size=15)) + ylab("Slatkin's Linearized Fst") + xlab("Distance between the populations (km)") 
dev.off()

cor.test(fst3$Distance, fst3$SlatkinDistance, method="spearman")
cor.test(fst3$Fst, fst3$SlatkinDistance, method="spearman")
cor.test(fst3$Distance, fst3$Fst, method="spearman")
##
###GERP based genetic load per individual####
chr1 <- read.table("./gerp_scores_1", header=TRUE)
chr2 <- read.table("./gerp_scores_2", header=TRUE)
chr3 <- read.table("./gerp_scores_3", header=TRUE)
chr4 <- read.table("./gerp_scores_4", header=TRUE)
chr5 <- read.table("./gerp_scores_5", header=TRUE)
chr6 <- read.table("./gerp_scores_6", header=TRUE)
chr7 <- read.table("./gerp_scores_7", header=TRUE)
chr8 <- read.table("./gerp_scores_8", header=TRUE)
chr9 <- read.table("./gerp_scores_9", header=TRUE)
chr10 <- read.table("./gerp_scores_10", header=TRUE)

chrs <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

##need to make the negative values to 0
chrs[chrs < 0] <- 0
tmp <- NULL
for (i in 2:ncol(chrs)){
  tmp <- c(tmp, sum(is.na(chrs[,i])))
}
summary(tmp)
nrow(chrs)
tmp <- (tmp / nrow(chrs)) * 100 
summary(tmp)
quantile(tmp, c(0.90))
tmp <- as.data.frame(tmp)
tmp$indvs <- colnames(chrs)[2:ncol(chrs)]
colnames(chrs)[2]
toremove <- tmp[tmp$tmp >= 1,]
nrow(toremove)

gerp <- NULL
for (i in 2:ncol(chrs)){
  gerp <- c(gerp, sum(chrs[,i], na.rm=T)) #when imputed data are used and there is no missing information
 }
gerp <- as.data.frame(gerp)
gerp$indvs <- colnames(chrs)[2:ncol(chrs)]
summary(gerp)

pops <- rep("elite", 155)
pops <- c(pops, substr(gerp[156:nrow(gerp),2], start = 1, stop = 2))
gerp$pops <- pops

load <- gerp
colnames(load)[1] <- "load"
ggplot(load, aes(x=pops, y=load)) + geom_boxplot()

#read in the metadata
metadata <- read.csv("../../maizeEU/master_thesis_scripts copy/metadata_for_PCA.csv")
head(metadata)
load2 <- merge(load, metadata, by.x="pops", by.y="Abb.", all.x=TRUE)
load2 <- load2[!(load2$indvs %in% toremove$indvs),]

load2$East <- load2$Longitude > 2.35

load <- as.data.frame(load)
medias <- tapply(load2$load, load2$pops, median)
medias <- as.data.frame(medias)

load2$Elites <- load2$pops == "elite"
summary(load2)

ks.test(load2[load2$Elites == FALSE,2], load2[load2$Elites == TRUE,2])

#distance from specific places. 
#spanish entry point : Tuy
latS <- metadata[41,3]
lonS <- metadata[41,4]
tmp <- NULL
for (i in 1: nrow(load2)){
  tmp <- c(tmp, distm(c(load2[i,6], load2[i,5]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
load2$DistSpain <- tmp

cor.test(load2$load, load2$DistSpain, method="spearman") ##significant!!
#german entry point: Barisis
latS <- metadata[3,3]
lonS <- metadata[3,4]
tmp <- NULL
for (i in 1: nrow(load2)){
  tmp <- c(tmp, distm(c(load2[i,6], load2[i,5]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
load2$DistGermany <- tmp
cor.test(load2$load, load2$DistGermany, method="spearman") #significant?
cor.test(load2$load, load2$DistSpain, method="spearman")
cor.test(load2$load, load2$DistGermany, method="spearman")

tapply(load2$load, load2$East, median)

##add the elites info
elites <- read.table("../../maizeEU/elites")
head(elites)
colnames(elites) <- c("pop", "elite")

tmp <- load2[load2$Elites == TRUE,]
nrow(tmp)
tmp1 <- NULL
for (i in tmp$indvs){
  tmp1 <- c(tmp1, strsplit(i, "_")[[1]][1])
}
tmp1
tmp$type <- tmp1
tmp2 <- load2[load2$Elites == FALSE,]
tmp2$type <- tmp2$pops
load2 <- rbind(tmp2, tmp)
elites2 <-read.table("../../maizeEU/elites2")
head(elites2)
colnames(elites2) <- c("pop", "elite")

load2 <- merge(load2, elites2, by.x="type", by.y="pop", all.x=TRUE)
head(load2)
load2[1089,19] <- "dent" #fix an error in the dataset
load2[1190,19] <- "dent"
ggplot(load2, aes(x=elite, y=load, fill=East)) + geom_boxplot() + aes(x=reorder(elite, Longitude)) + ylab("sum(gpns)")

#AL, CA, PL, ND, TR as dent like landraces
tmp <- c("dent landrace", rep("landrace", 3), "dent landrace", "landrace", "dent", "landrace", "flint", rep("landrace", 13), "dent landrace", rep("landrace", 4), "dent landrace", rep("landrace", 8), "dent landrace", "landrace",  "unknown", "landrace", "landrace")
load2$elite2 <- as.factor(load2$elite)
levels(load2$elite2) <- tmp

load2 <- load2[!(load2$elite2 == "unknown"),]
ggplot(load2, aes(y=load,x=elite2, fill=elite2)) + geom_boxplot()
tapply(load2$load, load2$elite2, mean)

ks.test(load2[load2$elite2 == "landrace",3], load2[load2$elite2 == "flint",3])
ks.test(load2[load2$elite2 == "landrace",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "flint",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "dent landrace",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "flint",3], load2[load2$elite2 == "dent landrace",3])

jpeg("loadGerpElitesLandraces.jpg")
ggplot(load2, aes(y=load,x=elite2, fill=elite2)) + geom_boxplot(show.legend = F) + theme_bw() +  theme(text=element_text(size=20)) + xlab("") + ylab("Genetic Load") 
dev.off()

jpeg("maizeEU.load.gerp.pops.jpeg", units="in", width=12, height=7, res=600)
ggplot(load2, aes(x=pops, y=load, fill=East)) + geom_boxplot(show.legend = F) + aes(x=reorder(pops, DistSpain)) + ylab("Genetic Load") + xlab("")  + theme(text=element_text(size=20)) + theme_bw()
dev.off() 

##
###highly deleterious load######
ggCom <- NULL
for (i in 1:10){
  gerp10 <- read.table(paste0("../gerps/Zea_mays.v4.chr", i, ".gerpcolRatesPlusOGs.txt",sep=""), header=FALSE)
  gpn10 <- read.table(paste0("./gerp_scores_", i, "all", sep=""), header=FALSE)
  gerp10 <- gerp10[!(gerp10$V6 == "-"),]
  gerp10 <- gerp10[gerp10$V4 >= 5, c(1,2)]
  tmp <- gpn10[gpn10$V1 %in% as.character(gerp10$V2),]
  ggCom <- rbind(ggCom, tmp)
}
nrow(ggCom)
head(ggCom)

gerp <- NULL
for (i in 2:ncol(ggCom)){
  gerp <- c(gerp, sum(ggCom[,i]>0, na.rm=T)) #count the number of those above 0
 }
gerp <- as.data.frame(gerp)
gerp$indvs <- colnames(chrs)[2:ncol(chrs)]
summary(gerp)

pops <- rep("elite", 155)
pops <- c(pops, substr(gerp[156:nrow(gerp),2], start = 1, stop = 2))
gerp$pops <- pops

load <- gerp
colnames(load)[1] <- "load"
ggplot(load, aes(x=pops, y=load)) + geom_boxplot()

#read in the metadata
metadata <- read.csv("../../maizeEU/master_thesis_scripts copy/metadata_for_PCA.csv")
head(metadata)
load2 <- merge(load, metadata, by.x="pops", by.y="Abb.", all.x=TRUE)
load2 <- load2[!(load2$indvs %in% toremove$indvs),]

load2$East <- load2$Longitude > 2.35
ggplot(load2, aes(x=pops, y=load, fill=East)) + geom_boxplot() + aes(x=reorder(pops, Longitude)) + ylab("sum(GREP)")

load <- as.data.frame(load)
medias <- tapply(load2$load, load2$pops, median)
medias <- as.data.frame(medias)

load2$Elites <- load2$pops == "elite"
summary(load2)

ks.test(load2[load2$Elites == FALSE,2], load2[load2$Elites == TRUE,2])

#distance from specific places. 
#spanish entry point : Tuy
latS <- metadata[41,3]
lonS <- metadata[41,4]
tmp <- NULL
for (i in 1: nrow(load2)){
  tmp <- c(tmp, distm(c(load2[i,6], load2[i,5]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
load2$DistSpain <- tmp

cor.test(load2$load, load2$DistSpain, method="spearman")

#german entry point: Barisis
latS <- metadata[3,3]
lonS <- metadata[3,4]
tmp <- NULL
for (i in 1: nrow(load2)){
  tmp <- c(tmp, distm(c(load2[i,6], load2[i,5]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
load2$DistGermany <- tmp
cor.test(load2$load, load2$DistGermany, method="spearman") #significant?

cor.test(load2$load, load2$DistSpain, method="spearman")
cor.test(load2$load, load2$DistGermany, method="spearman")

tapply(load2$load, load2$East, median)

##add the elites info
elites <- read.table("../../maizeEU/elites")
head(elites)
colnames(elites) <- c("pop", "elite")

tmp <- load2[load2$Elites == TRUE,]
nrow(tmp)
tmp1 <- NULL
for (i in tmp$indvs){
  tmp1 <- c(tmp1, strsplit(i, "_")[[1]][1])
}
tmp1
tmp$type <- tmp1
tmp2 <- load2[load2$Elites == FALSE,]
tmp2$type <- tmp2$pops
load2 <- rbind(tmp2, tmp)
elites2 <-read.table("../../maizeEU/elites2")
head(elites2)
colnames(elites2) <- c("pop", "elite")

load2 <- merge(load2, elites2, by.x="type", by.y="pop", all.x=TRUE)
head(load2)

#AL, CA, PL, ND, TR as dent like landraces
tmp <- c("dent landrace", rep("landrace", 3), "dent landrace", "landrace", "dent", "landrace", "flint", rep("landrace", 13), "dent landrace", rep("landrace", 4), "dent landrace", rep("landrace", 8), "dent landrace", "landrace",  "unknown", "landrace", "landrace")
load2$elite2 <- as.factor(load2$elite)
levels(load2$elite2) <- tmp

load2 <- load2[!(load2$elite2 == "unknown"),]

ks.test(load2[load2$elite2 == "landrace",3], load2[load2$elite2 == "flint",3])
ks.test(load2[load2$elite2 == "landrace",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "flint",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "dent landrace",3], load2[load2$elite2 == "dent",3])
ks.test(load2[load2$elite2 == "flint",3], load2[load2$elite2 == "dent landrace",3])

load2 <- load2[!is.na(load2$elite2),]

jpeg("loadGerpElitesLandracesHigh.jpg")
ggplot(load2, aes(y=load,x=elite2, fill=elite2)) + geom_boxplot(show.legend = F) + theme_bw() +  theme(text=element_text(size=20)) + xlab("") + ylab("Highly Deleterious Mutations") 
dev.off()

jpeg("maizeEU.load.gerp.high.pops.jpeg", units="in", width=12, height=7, res=600)
ggplot(load2, aes(x=pops, y=load, fill=East)) + geom_boxplot(show.legend = F) + aes(x=reorder(pops, DistSpain)) + ylab("Highly Deleterious Mutations") + xlab("")  + theme(text=element_text(size=20)) + theme_bw()
dev.off() 

save.image("gerpbasedLoad.RData")

###
###Fixed & segregating load####
setwd("/home/mtakou/Dropbox/margarita/posdocAGStetter/gpn/Scores/")
load("gerpbasedLoad.RData")
head(load2)

write.table(load2, file="load.txt", col.names = T, row.names = F)
write.table(table(load2$pops), file="indiv.txt")
#fixed
ggfix <- NULL
fixedall <- NULL
for (u in c("dentlandrace", "landrace" ,"dent", "flint")){
  ggfix <- NULL
  for (i in 1:10){
    print(u)
    gerp10 <- read.table(paste0("../fixed_load_", i, "_", u, sep=""), header=TRUE)
    gpn10 <- read.table(paste0("./gerp_scores_", i, "all", sep=""), header=FALSE)
    tmp <- gpn10[gpn10$V1 %in% as.character(gerp10$POS),]
    ggfix <- rbind(ggfix, tmp)
  }
  if (u == "dentlandrace"){indc <- load2[load2$elite2 == "dent landrace", 4]}
  else{indc <- load2[load2$elite2 == u, 4]}
  #print(indc)
  ggfix[ggfix < 0] <- 0
  fx <- NULL
  print(fx)
  for (i in 2:ncol(ggfix)){
    fx <- c(fx, sum(as.numeric(ggfix[,i]), na.rm=TRUE)) 
   }
  fx <- as.data.frame(fx)
  fx$indvs <- colnames(chrs)[2:ncol(chrs)]
  fx <- fx[fx$indvs %in% indc,]
  fx$elite2 <- u
  summary(fx)
  fixedall <- rbind(fixedall, fx)
}
summary(fixedall)
nrow(fixedall)
head(fixedall)
fixedall$elite2 <- as.factor(fixedall$elite2)
fixedall$elite2 <- factor(fixedall$elite2, levels = c("dentlandrace", "landrace" ,"dent", "flint"))

fixedallM <- tapply(fixedall$fx, fixedall$elite2, mean)
fixedallM <- as.data.frame(fixedallM)
fixedallM$elites <- as.factor(rownames(fixedallM))
colnames(fixedallM)[1] <- "load"
fixedallM$elites <- factor(fixedallM$elites, levels = c("dentlandrace", "landrace" ,"dent", "flint"))

jpeg("loadGerpElitesLandracesFixed.jpg")
ggplot(fixedallM, aes(y=load,x=elites, fill=elites)) + geom_col(show.legend = F) + theme_bw() +  theme(text=element_text(size=20)) + xlab("") + ylab("Fixed Genetic Load") 
dev.off()

ks.test(fixedall[fixedall$elite2 == "landrace",1], fixedall[fixedall$elite2 == "flint",1])
ks.test(fixedall[fixedall$elite2 == "landrace",1], fixedall[fixedall$elite2 == "dent",1])
ks.test(fixedall[fixedall$elite2 == "flint",1], fixedall[fixedall$elite2 == "dent",1])
ks.test(fixedall[fixedall$elite2 == "dentlandrace",1], fixedall[fixedall$elite2 == "dent",1])
ks.test(fixedall[fixedall$elite2 == "flint",1], fixedall[fixedall$elite2 == "dentlandrace",1])

#segregating
ggfix <- NULL
segall <- NULL
for (u in c("dentlandrace", "landrace" ,"dent", "flint")){
  ggfix <- NULL
  for (i in 1:10){
    gerp10 <- read.table(paste0("../segregating_load_", i, "_", u, sep=""), header=TRUE)
    gpn10 <- read.table(paste0("./gerp_scores_", i, "all", sep=""), header=FALSE)
    tmp <- gpn10[gpn10$V1 %in% as.character(gerp10$POS),]
    ggfix <- rbind(ggfix, tmp)
  }
  if (u == "dentlandrace"){indc <- load2[load2$elite2 == "dent landrace", 4]}
  else{indc <- load2[load2$elite2 == u, 4]}
  fx <- NULL
  for (i in 2:ncol(ggfix)){
    fx <- c(fx, sum(as.numeric(ggfix[,i]>0), na.rm=TRUE)) 
  }
  fx <- as.data.frame(fx)
  fx$indvs <- colnames(chrs)[2:ncol(chrs)]
  fx <- fx[fx$indvs %in% indc,]
  fx$elite2 <- u
  summary(fx)
  segall <- rbind(segall, fx)
}

nrow(segall)
head(segall)
segall$elite2 <- as.factor(segall$elite2)
segall$elite2 <- factor(segall$elite2, levels = c("dentlandrace", "landrace" ,"dent", "flint"))

ggplot(segall, aes(x=elite2, y=fx, fill=elite2)) + geom_boxplot()

jpeg("loadGerpElitesLandracesSeg.jpg")
ggplot(segall, aes(y=fx,x=elite2, fill=elite2)) + geom_boxplot(show.legend = F) + theme_bw() +  theme(text=element_text(size=20)) + xlab("") + ylab("Segregating Genetic Load") 
dev.off()

ks.test(segall[segall$elite2 == "landrace",1], segall[segall$elite2 == "flint",1])
ks.test(segall[segall$elite2 == "landrace",1], segall[segall$elite2 == "dent",1])
ks.test(segall[segall$elite2 == "flint",1], segall[segall$elite2 == "dent",1])
ks.test(segall[segall$elite2 == "dentlandrace",1], segall[segall$elite2 == "dent",1])
ks.test(segall[segall$elite2 == "flint",1], segall[segall$elite2 == "dentlandrace",1])
##
#######Fis######
fis <- read.table("./maizeEU.noDH.filtQual.het", header=T)
head(fis)
hist(fis$F)

metadata <- read.csv("./master_thesis_scripts copy/metadata_for_PCA.csv")
head(metadata)

pops <- rep("elite", 155)
pops <- c(pops, substr(fis[156:nrow(fis),1], start = 1, stop = 2))
fis$pops <- pops

ggplot(fis, aes(x=pops, y=F)) + geom_boxplot()

fis <- merge(fis, metadata, by.x="pops", by.y="Abb.", all.x=TRUE)

fis$East <- fis$Longitude > 2.35
ggplot(fis, aes(x=pops, y=F, fill=East)) + geom_boxplot() + aes(x=reorder(pops, Longitude)) 

fis <- as.data.frame(fis)

fis$Elites <- fis$pops == "elite"
summary(fis)
ggplot(fis, aes(x=Elites, y=F, fill=Elites)) + geom_boxplot()
ks.test(fis[fis$Elites == FALSE,6], fis[fis$Elites == TRUE,6])

#compare with load
load.fis <- merge(load2, fis, by="pops")

ggplot(load.fis, aes(x=F, y=load)) + geom_point() + geom_smooth(method="lm")
cor.test(load.fis$load, load.fis$F, method="spearman")

ggplot(load.fis, aes(x=elite, y=F, fill=East.x)) + geom_boxplot() + aes(x=reorder(elite, Longitude.x)) 

#compare with Het
het.fis <- merge(fis, het, by="Full.name", all.x=TRUE)
head(het.fis)

ggplot(het.fis, aes(x=mean.expected.heterozygosity, y=F)) + geom_point()

cor.test(het.fis$mean.expected.heterozygosity, het.fis$F, method="spearman")

tiff("maizeEU.fis.density.tiff", units="in", width=5, height=5, res=600)
ggplot(het.fis, aes(x=F)) + geom_density() + theme_bw() + xlab("Inbreeding Coefficient") + theme(text=element_text(size=15))
dev.off()

tiff("maizeEU.het.density.tiff", units="in", width=5, height=5, res=600)
ggplot(het.fis, aes(x=mean.expected.heterozygosity)) + geom_density() + theme_bw() + xlab("Mean Expected Heterozygosity") + theme(text=element_text(size=15))
dev.off()

tiff("maizeEU.fis.pops.tiff", units="in", width=12, height=5, res=600)
ggplot(het.fis, aes(x=pops,y=F)) + geom_boxplot() + theme_bw() + ylab("Inbreeding Coefficient") + xlab("Population") + theme(text=element_text(size=15)) + aes(x=reorder(pops, Longitude.x)) 
dev.off()

