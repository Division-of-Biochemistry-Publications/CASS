# https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html
# packages ####
library(tidyverse)
library(circlize)
library(viridis)

# load data ####

dataPath = "../RNASEQ/"
outputPath = paste0(dataPath,"10_CRN/")
figurePath = paste0(outputPath,"FIGURES/CHORDDIAGRAM/")
dir.create(figurePath, recursive = TRUE)

wald = read.csv(paste0(dataPath,"08_DESEQ2/WALD_TISSUE.csv")) %>% rename("Geneid" = X) # SR versus Stem, independent of cluster

clusters = read.csv(paste0(dataPath,"10_CRN/CLUSTERS.csv")) %>% rename("Geneid" = X) %>%
  filter(Cluster != -1)
clusterOrder = c("P", "P/C", "P/X", "C/X", "X")
clusters = clusters[order(match(clusters$clusterNameShort,clusterOrder)),] %>%
  mutate(clusterNameShort = paste(Tissue, clusterNameShort, sep = "\n"))
clusterOrder = c(paste("Stem",clusterOrder, sep = "\n"), paste("SR",clusterOrder, sep = "\n"))




# preparations ####

wald = wald %>% filter((padj < 0.001) & (abs(log2FoldChange) > log2(1.5))) %>%
  split((.$log2FoldChange) > 0) %>% lapply(.,"[[","Geneid")
names(wald) = c("SR↓", "SR↑")
clusters = clusters %>% split(.$clusterNameShort) %>% lapply(.,"[[","Geneid")
clusters = clusters[clusterOrder]
degs = c(clusters,wald)
data = data.frame(matrix(NA, ncol=length(degs), nrow=length(degs)),row.names = names(degs))

names(data) = names(degs)
combs = as.data.frame(combinat::combn2(names(data))) %>% mutate(tmp = paste0(V1,V2))
for(i in combs$tmp) {
  print(i)
  x = combs[combs$tmp == i,"V1"]
  y = combs[combs$tmp == i,"V2"]
  data[x,y] = length(intersect(degs[[x]],degs[[y]]))
  data[y,x] = length(intersect(degs[[x]],degs[[y]]))
  data[x,x] = length(intersect(degs[[x]],degs[[x]]))
}
group = c(rep("SR",5),rep("Stem",5),rep("Wald",2))
names(group) = names(data)


baseGridColor = plasma(n = length(unique(group)),direction = -1)
names(baseGridColor) = unique(group)
gridColor = c(rep(baseGridColor["SR"],5),rep(baseGridColor["Wald"],2),rep(baseGridColor["Stem"],5))
names(gridColor) = NULL


rowColor = c(rep(baseGridColor["SR"],5),rep(baseGridColor["Stem"],5),rep(baseGridColor["Wald"],2))
names(rowColor) = NULL
rowColor[1:11] = gsub("FF", "1A", rowColor[1:11]) # changes transparency for the first 10 colors 1 -> 0.2
circos.clear()
png(paste0(figurePath,"WALD_HIGHLIGHT_UP.png"),units = "cm", width = 6, height = 6, res = 600)
par(cex = 0.5)
circos.clear()
circos.par$track.margin = c(0.0000001, 0.0375)
chordDiagram(as.matrix(data), symmetric = TRUE, self.link = 0, group = group, grid.col = gridColor,
             row.col = rowColor, big.gap = 5, small.gap = 2.5)
dev.off()
circos.clear()


rowColor = c(rep(baseGridColor["SR"],5),rep(baseGridColor["Stem"],5),rep(baseGridColor["Wald"],2))
names(rowColor) = NULL
rowColor[-11] = gsub("FF", "1A", rowColor[-11]) # changes transparency for the first 10 colors 1 -> 0.2
png(paste0(figurePath,"WALD_HIGHLIGHT_DOWN.png"),units = "cm", width = 6, height = 6, res = 600)
par(cex = 0.5)
circos.clear()
circos.par$track.margin = c(0.0000001, 0.0375)
chordDiagram(as.matrix(data),symmetric = TRUE,self.link = 0,group = group,grid.col = gridColor,
             row.col = rowColor, big.gap = 5, small.gap = 2.5)
dev.off()
circos.clear()


rowColor = c(rep(baseGridColor["SR"],5),rep(baseGridColor["Stem"],5),rep(baseGridColor["Wald"],2))
names(rowColor) = NULL
rowColor[c(1:5,11,12)] = gsub("FF", "1A", rowColor[c(1:5,11,12)]) # changes transparency for the first 10 colors 1 -> 0.2
png(paste0(figurePath,"WALD_HIGHLIGHT_SR.png"),units = "cm", width = 6, height = 6, res = 600)
circos.clear()
circos.par$track.margin = c(0.0000001, 0.0375)
par(cex = 0.5)
chordDiagram(as.matrix(data),symmetric = TRUE,self.link = 0,group = group,grid.col = gridColor,
             row.col = rowColor, big.gap = 5, small.gap = 2.5)
dev.off()

