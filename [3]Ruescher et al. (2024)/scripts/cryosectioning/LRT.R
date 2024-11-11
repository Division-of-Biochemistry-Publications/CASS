# load packages ####
library(tidyverse)
library(DESeq2)
library(ggvenn)

# variables ####
print("creating variables")

dataPath = "../RNASEQ/"
outputPath = paste0(dataPath,"08_DESEQ2/")
figurePath = paste0(outputPath,"FIGURES")
dir.create(figurePath, recursive = TRUE)

# load data ####
metaData <- read.csv(paste0(dataPath,"00_META_DATA/META_DATA.csv")) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Distance = factor(Distance, levels = unique(Distance)))
metaDataClusters <- read.csv(paste0(dataPath,"09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv")) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Distance = factor(Distance, levels = unique(Distance)))
featureCounts <- read.delim(paste0(dataPath,"07_FEATURECOUNTS/FEATURECOUNTS_gene_UNIQUE.txt"), 
                            row.names=1, comment.char="#")

# change colnames in featureCounts ####
print("change colnames in featureCounts")
for (i in metaData$Sample) {
  colnames(featureCounts)[grepl(i, colnames(featureCounts))] = i
}

# LRT against cluster ####
lrt_res = list()
for(i in unique(metaDataClusters$Tissue)) {
  print(i)
  samples = metaDataClusters[metaDataClusters$Tissue == i,"Sample"]
  count_matrix = as.matrix(featureCounts[,samples])
  dds = DESeqDataSetFromMatrix(countData = count_matrix, 
                               colData = metaDataClusters[metaDataClusters$Sample %in% colnames(count_matrix),],
                               design = ~ Plant + Cluster)
  deseq = DESeq(dds, test = "LRT", reduced = ~ Plant)
  lrt_res[[i]] = data.frame(Tissue = i,as.data.frame(results(deseq, pAdjustMethod = "fdr", 
                                                             cooksCutoff=FALSE))) %>%
    rownames_to_column("Geneid")
} # i
lrt_res = bind_rows(lrt_res)
write.csv(lrt_res, file = paste0(outputPath,"LRT.csv"), row.names = FALSE)

# venn diagram 

degs = lrt_res[lrt_res$padj <0.001,] %>%
  split(.,f=.$Tissue)
degs = lapply(degs,"[[",1)

ggvenn(
  degs,
  digits = 0,
  fill_color = viridis::plasma(n = 2),
  stroke_size = 0,
  set_name_size = 2,
  text_size = 1
)
ggsave(
  filename = "LRT_VENN.png",
  path = figurePath,
  units = "cm",
  width = 8,
  height = 8,
  dpi = 600,
  device = "png"
)

# Wald test controling for cluster ####
print("Wald test")
count_matrix = as.matrix(featureCounts[,metaDataClusters$Sample])
dds = DESeqDataSetFromMatrix(countData = count_matrix, 
                             colData = metaDataClusters,
                             design = ~ Plant + Cluster + Tissue) 
deseq = DESeq(dds, test = "Wald")
res = as.data.frame(results(deseq,contrast = c("Tissue","SR","Stem"),pAdjustMethod = "fdr",
              cooksCutoff = FALSE))
write.csv(res,file = paste0(outputPath,"WALD_TISSUE.csv"))

# venn

waldDegs = res %>% rownames_to_column("Geneid") %>% filter((padj < 0.001) & (abs(log2FoldChange) > log2(1.5))) %>%
  split((.$log2FoldChange) > 0) %>% lapply(.,"[[","Geneid")
names(waldDegs) = c("SR↓", "SR↑")

combinedDegs = c(degs,waldDegs)

ggvenn(
  combinedDegs,
  digits = 0,
  fill_color = viridis::plasma(n = 4),
  stroke_size = 0,
  set_name_size = 3,
  text_size = 2
)
ggsave(
  filename = "VENN.png",
  path = figurePath,
  units = "cm",
  width = 8,
  height = 8,
  dpi = 600,
  device = "png"
)
