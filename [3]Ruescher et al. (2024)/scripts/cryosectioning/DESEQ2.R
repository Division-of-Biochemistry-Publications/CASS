# load packages ####
library(tidyverse)
library(DESeq2)
library(limma)

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
featureCounts <- read.delim(paste0(dataPath,"07_FEATURECOUNTS/FEATURECOUNTS_gene_UNIQUE.txt"), 
                            row.names=1, comment.char="#")

# change colnames in featureCounts ####
print("change colnames in featureCounts")
for (i in metaData$Sample) {
  print(i)
  colnames(featureCounts)[grepl(i, colnames(featureCounts))] = i
}

# count normalization, transformation and batch correction ####

print("normalization")
count_matrix = as.matrix(featureCounts[,metaData$Sample])
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = metaData,
                             design = ~ Distance + Tissue)
norm = as.data.frame(counts(estimateSizeFactors(dds),normalized=TRUE))
write.csv(norm, file = paste0(outputPath,"NORMALIZED.csv"))

print("VST")
drl = vst(dds)
logvalues = as.data.frame(assay(drl))
write.csv(logvalues, file = paste0(outputPath,"VST.csv"))

print("Batch correction")

corrected = list()

for(tissue in unique(metaData$Tissue)) {
  print(tissue)
  metaDataTissue = metaData[metaData$Tissue == tissue,]
  samples = metaDataTissue$Sample
  corrected[[tissue]] = as.data.frame(removeBatchEffect(
    logvalues[,samples],
    batch = metaDataTissue$Plant,
    design = model.matrix(data=metaDataTissue,~Distance)
  ))
}

corrected = bind_cols(corrected)
corrected = corrected[,metaData$Sample]
write.csv(corrected, file = paste0(outputPath,"VST_LIMMA.csv"))
