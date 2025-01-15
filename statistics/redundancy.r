library(pegas)
library(ggplot2)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(gridExtra)

setwd("./rda")
load("rda.maizeEU.RData")

##per population samples
filelist <- list.files("./", pattern="frq")
filelist

#create the allele frequencies table
allel.freqs <- NULL
pops <- NULL
for (i in 1:length(filelist)){
  tmp <- read.table(filelist[i], header=T, fill=T, row.names=NULL)
  allel.freqs <- cbind(allel.freqs, tmp[,6])
  pops <- c(pops, strsplit(filelist[i],"[.]")[[1]][1])
}
allel.freqs <- as.data.frame(allel.freqs)
colnames(allel.freqs) <- pops
allel.freqs$Chrom <- tmp$row.names
allel.freqs$Pos <- tmp$CHROM
allel.freqs <- allel.freqs[!(allel.freqs$Chrom == "0"),]
allel.freqs2 <- allel.freqs[,1:38]
summary(allel.freqs2)
allel.freqs2 <- na.omit(allel.freqs2)
allel.freqs2 <- t(allel.freqs2)
allel.freqs2 <- as.data.frame(allel.freqs2)

#load metadata
setwd("./wc2.1_30s_bio/")
ras <- stack(list.files("./"))
names(ras) <- unlist(unlist(strsplit(list.files("./"), split=".tif")))

##metadata for the populations
setwd("../")
coords <- read.table("metadata_for_PCA.csv", header=T, sep=",", fill=T)
coords <- coords[coords$Abb. %in% rownames(allel.freqs2),c(1,3,4)]
colnames(coords)[1] <-"Population" 
head(coords)

#extract the environmental data for each value
env <- data.frame(extract(ras, coords[,3:2]))
head(env)
row.names(env) <- c(coords$Population)

save.image("rda.maizeEU.RData")

#load information about degenerate sites and keep the neutral (4-fold)
neutrals <- NULL
for (i in 1:10){
  tmp <- read.table(paste('./degenerate', i, '/degeneracy-all-sites.bed', sep=""), header=FALSE, fill=TRUE)
  tmp <- tmp[tmp$V5 == "4",]
  tmp1 <- allel.freqs[allel.freqs$Chrom == as.character(i),]
  tmp1 <- tmp1[tmp$V3 %in% tmp1$Pos,]
  neutrals <- rbind(neutrals, tmp1)
}
head(neutrals)
rm(tmp)
rm(tmp1)
save.image("rda.maizeEU.RData")

neutrals2 <- neutrals[,1:38]
summary(neutrals2)
neutrals2 <- na.omit(neutrals2)
neutrals2 <- t(neutrals2)
neutrals2 <- as.data.frame(neutrals2)

#Infer population structure
pca <- rda(neutrals2, scale=T)
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

#based on screeplot keep the first 2
PCs <- scores(pca, choices=c(1:2), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(allel.freqs2), PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2")

#gather all the variables together
variables <- data.frame(coords, PopStruct[,-1], env)

## europe range shapefile
range <- read_sf("./Europe/europe_10km.shp") 

##variable selection: forward model building
## Null model
RDA0 <- rda(allel.freqs2 ~ 1,  variables) 

## Full model
RDAfull <- rda(allel.freqs2 ~ wc2.1_30s_bio_1 + wc2.1_30s_bio_2 + wc2.1_30s_bio_3 +wc2.1_30s_bio_4 + wc2.1_30s_bio_5 +
                 wc2.1_30s_bio_6 + wc2.1_30s_bio_7 + wc2.1_30s_bio_8 +  wc2.1_30s_bio_9 + wc2.1_30s_bio_10 + wc2.1_30s_bio_11 +
                 wc2.1_30s_bio_12 + wc2.1_30s_bio_13 + wc2.1_30s_bio_14 + wc2.1_30s_bio_15 + wc2.1_30s_bio_16 + wc2.1_30s_bio_17 + 
                 wc2.1_30s_bio_18 + wc2.1_30s_bio_19, variables)


## Stepwise selection procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
mod$anova

mod2 <- ordiR2step(RDA0, RDAfull, Pin = 0.05, R2permutations = 1000, R2scope = T)
mod2$anova #keep from here the values

##variance partitioning
pRDAfull <- rda(allel.freqs2 ~ Latitude + Longitude + PC1 + PC2 + wc2.1_30s_bio_9 + wc2.1_30s_bio_15 + wc2.1_30s_bio_19 , variables)
RsquareAdj(pRDAfull)
anova(pRDAfull) #significant p=0.001

pRDAclim <- rda(allel.freqs2 ~  wc2.1_30s_bio_9 + wc2.1_30s_bio_15 + wc2.1_30s_bio_19  + Condition(Latitude + Longitude + PC1 + PC2), variables)
RsquareAdj(pRDAclim)
anova(pRDAclim) #p=0.001

pRDAlat <- rda(allel.freqs2 ~ Latitude + Longitude + Condition(PC1 + PC2 + wc2.1_30s_bio_9 + wc2.1_30s_bio_15 + wc2.1_30s_bio_19 ), variables)
RsquareAdj(pRDAlat)
anova(pRDAlat) #significant p=0.011

pRDApcs <- rda(allel.freqs2 ~  PC1 + PC2 + Condition(Latitude + Longitude + wc2.1_30s_bio_9 + wc2.1_30s_bio_15 + wc2.1_30s_bio_19), variables)
RsquareAdj(pRDApcs)
anova(pRDApcs) #significant p=0.001

save.image("rda.maizeEU.RData")

##correlogram
library(GGally)
toPlot <- variables[,c(2:5,12, 16,24,25,26)]
colnames(toPlot)[5:9] <- c("bio15", "bio19", "bio9", "SouthRoute", "NorthRoute")
ggsave(ggpairs(toPlot), file="corr.RDA.jpeg", width=10, height=10, dpi=300)

#genome scan using pRDA
screeplot(pRDAclim)

####Correlate with distance from entry point the significant rdas####
#or based on clusters?
head(variables)

#add distances from too hypothetical points
#spanish etry point : Tuy
latS <- variables[variables$Population == "TU",2]
lonS <- variables[variables$Population == "TU",3]
tmp <- NULL
for (i in 1: nrow(variables)){
  tmp <- c(tmp, distm(c(variables[i,3], variables[i,2]), c(lonS, latS), fun = distHaversine) / 1000)  #transform them to km
}
tmp
variables$DistSpain <- tmp

#german entry point: Barisis
latG <- variables[variables$Population == "BA",2]
lonG <- variables[variables$Population == "BA",3]
tmp <- NULL
for (i in 1: nrow(variables)){
  tmp <- c(tmp, distm(c(variables[i,3], variables[i,2]), c(lonG, latG), fun = distHaversine) / 1000)  #transform them to km
}
tmp
variables$DistGerman <- tmp

cor.test(variables$DistGerman, variables$wc2.1_30s_bio_9, method="spearman")
cor.test(variables$DistGerman, variables$wc2.1_30s_bio_15, method="spearman")
cor.test(variables$DistGerman, variables$wc2.1_30s_bio_19, method="spearman")


cor.test(variables$DistSpain, variables$wc2.1_30s_bio_9, method="spearman")
cor.test(variables$DistSpain, variables$wc2.1_30s_bio_15, method="spearman")
cor.test(variables$DistSpain, variables$wc2.1_30s_bio_19, method="spearman")

variables$East <- variables$Longitude > 3.4
variables[variables$Abb. == "ND",27] <- FALSE
variables[variables$Abb. == "SO",27] <- FALSE
variables[variables$Abb. == "BA",27] <- FALSE

for (i in 6:24){
  print(cor.test(variables[,i], variables$DistSpain, method="spearman"))
}

for (i in 6:24){
  print(cor.test(variables[,i], variables$DistGerman, method="spearman"))
}

clm.ks <- NULL
for (i in 6:24){
  clm.ks <- rbind(clm.ks ,c(ks.test(variables[variables$East == TRUE, i], variables[variables$East == FALSE, i])$p.value, colnames(variables)[i]))
}
head(clm.ks)

variables$East2 <- variables$East
variables$East2 <- as.factor(variables$East2)
levels(variables$East2) <- c("West", "East")

p1 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_9, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Mean temperature of driest quarter")+ theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p2 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_15, fill=East2)) + geom_boxplot() + theme(text=element_text(size=15)) + theme_bw() +  xlab("") +
  ylab("Precipitation Seasonality") + theme(axis.text.x = element_blank()) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p3 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_19, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Precipitation of Coldest Quarter")+ theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p4 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_1, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Annual Mean Temperature") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p5 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_2, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Mean Diurnal Range")+ theme(axis.text.x = element_blank()) + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p4 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_3, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Isothermality") + theme(axis.text.x = element_blank())+ theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p6 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_4, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Temperature Seasonality") + theme(legend.position="none") + theme(text=element_text(size=20)) + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p7 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_5, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Max Temperature of Warmst Month") + theme(axis.text.x = element_blank())+ theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p8 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_6, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Min Temperature of Coldest Month") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p9 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_7, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Temperature Annual Range")+ theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p10 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_8, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Mean Temperature of Wettest Quarter") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p11 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_10, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Mean Temperature of Warmest Quarter") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p12 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_11, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Mean Temperature of Coldest Quarter") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p13 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_12, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Annual Precipitation") + theme(legend.position="none") + theme(text=element_text(size=20)) + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p14 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_13, fill=East2)) + geom_boxplot() + theme_bw() +  xlab("") +
  ylab("Precipitation of Wettest Month") + theme(axis.text.x = element_blank()) + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))

p15 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_14, fill=East2)) + geom_boxplot() + theme_bw() + xlab("Cluster") +
  ylab("Precipitation of Driest Month")  + theme(text=element_text(size=15))  + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p16 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_16, fill=East2)) + geom_boxplot() + theme_bw() + 
  ylab("Precipitation of Wettest Quarter") + xlab("Cluster") + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p17 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_17, fill=East2)) + geom_boxplot() + theme_bw() + 
  ylab("Precipitation of Driest Quarter") + xlab("Cluster")  + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))
p18 <- ggplot(variables, aes(x=East2, y=wc2.1_30s_bio_18, fill=East2)) + geom_boxplot() + theme_bw() + 
  ylab("Precipitation of the Warmest Quarter") + xlab("Cluster") + theme(text=element_text(size=15)) + theme(legend.position="none") + stat_compare_means(aes(label=..p.signif..), hide.ns = FALSE,position = position_nudge(x=0.5))


ggsave(grid.arrange(p1, p2, p3, p4,  p8, p9, p10, p11, ncol=4), file="biov1.png",dpi=300, width=10, height=10)
ggsave(grid.arrange(p5, p7, p12, p14, p15, p16, p17, p18, ncol=4), file="biov3.png",dpi=300, width=10, height=10)
ggsave(grid.arrange(p13, p6,ncol=2), file="Bio.main.png",dpi=300, width=10, height=7)


####Environmental Distance between the populations#####
head(variables)
variables2 <- variables
variables2[,6:24] <- scale(variables2[,6:24], center = TRUE, scale = TRUE)

envirDist <- dist(variables2[,6:24], method = "euclidean", diag = TRUE, upper = TRUE)
envirDist <- as.matrix(envirDist)
colnames(envirDist) <- variables2$Population
rownames(envirDist) <- variables2$Population
envirDist <- as.data.frame(envirDist)
envirDist$Pop2 <- variables2$Population
envirDist <- melt(envirDist, c("Pop2"))
colnames(envirDist) <- c("Pop2", "Pop1", "EnvDist")

envirDist <- merge(envirDist, fst3, by=c("Pop1", "Pop2"))
head(envirDist)

ggplot(envirDist, aes(x=EnvDist, y=Fst)) + geom_point() + geom_smooth(method="lm")
ggplot(envirDist, aes(x=EnvDist, y=SlatkinDistance)) + geom_point() + geom_smooth(method="lm")

envirDist$East <- envirDist$Longitude.x > 3.34
envirDist$East2 <- envirDist$Longitude.y > 3.34
envirDist[envirDist$Pop1 == "BA", 11] <- FALSE
envirDist[envirDist$Pop1 == "SO", 11] <- FALSE
envirDist[envirDist$Pop1 == "ND", 11] <- FALSE

envirDist[envirDist$Pop2 == "BA", 12] <- FALSE
envirDist[envirDist$Pop2 == "SO", 12] <- FALSE
envirDist[envirDist$Pop2 == "ND", 12] <- FALSE
envirDist$Cluster <- envirDist$East==envirDist$East2

tiff("maizeEU.envir.Fst.tiff", units="in", width=7, height=7, res=600)
ggplot(envirDist, aes(x=EnvDist, y=Fst, col=Cluster)) + theme_bw() + geom_point() + geom_smooth(method="lm") + theme(text=element_text(size=20), legend.position = "null") + ylab("Fst") + xlab("Environmental Distance between the populations") 
dev.off()

cor.test(envirDist$Fst, envirDist$EnvDist, method="spearman")
cor.test(envirDist$SlatkinDistance, envirDist$EnvDist, method="spearman")

cor.test(envirDist[envirDist$Cluster == TRUE, 4], envirDist[envirDist$Cluster == TRUE, 3], method="spearman")
cor.test(envirDist[envirDist$Cluster == FALSE, 4], envirDist[envirDist$Cluster == FALSE, 3], method="spearman")

cor.test(envirDist[envirDist$East == TRUE, 4], envirDist[envirDist$East == TRUE, 9], method="spearman")
cor.test(envirDist[envirDist$East == FALSE, 4], envirDist[envirDist$East == FALSE, 9], method="spearman")

##Correlate with environment load
load("gerpbasedAllLoad.RData")

head(load2)
loadc <- merge(load2, variables, by.x="pops", by.y="Population", all.x=TRUE)
summary(loadc)

cor.test(loadc$wc2.1_30s_bio_9, loadc$load, method="spearman")
cor.test(loadc$wc2.1_30s_bio_15, loadc$load, method="spearman")
cor.test(loadc$wc2.1_30s_bio_19, loadc$load, method="spearman")
cor.test(loadc$load, loadc$Longitude.x, method="spearman")

cor.test(loadc$wc2.1_30s_bio_15, loadc$Longitude.x, method="spearman")

cor.test(loadc$DistSpain.x, loadc$load, method="spearman")
cor.test(loadc$DistGermany, loadc$load, method="spearman")