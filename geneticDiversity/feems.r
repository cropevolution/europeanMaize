setwd("/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU")

coords <- read.table("./FEEMS-maizeEU.all.filteredQual/maizeEU.all.filteredQual-FEEMS/coords.raw")
head(coords)
nrow(coords)

samples <- read.table("./FEEMS-maizeEU.all.filteredQual/maizeEU.all.filteredQual-FEEMS/Maize-0.8.pruned.fam")
head(samples)
nrow(samples)
tmp <- samples$V2
samples <- NULL
for (i in tmp){samples <- c(samples, strsplit(i, '[.]')[[1]][1])}
samples <- as.data.frame(samples)

names <- read.table("./FEEMS-maizeEU.all.filteredQual/maizeEU.all.filteredQual-FEEMS/coords.names", sep=";")
head(names)
nrow(names)
names <- names$V1

coords <- cbind(names, coords)
head(coords)
head(samples)
ncoords <- merge(coords, samples, by.x="names", by.y="samples", all.y=TRUE)
nrow(ncoords)
head(ncoords)
ncoords <- ncoords[,2:3]
nrow(ncoords)

tmp <- seq(0.000001, to=0.0009, length.out=nrow(ncoords))
ncoords$V1 <- ncoords$V1 + tmp
ncoords$V2 <- ncoords$V2 + tmp

write.table(ncoords, file="european.fixed.coord", sep="\t", row.names = F, col.names = F)
