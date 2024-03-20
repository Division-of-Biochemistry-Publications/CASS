# load packages ####

library(tidyverse)
library(DESeq2)

# variables ####

inputPath = "../data/RNASEQ/"
outputPath = "../data/RNASEQ/08_DESEQ2/"

dir.create(outputPath, recursive = TRUE)

# prepare meta data ####

attributes = read.csv(paste0(inputPath,"00_META_DATA/BIOSAMPLE_ATTRIBUTES.csv")) %>%
  filter(Tissue %in% c("SoL","Stem_LP","Stem_LC","SR","FR")) %>%
  select(Sample,Tissue,Stage,Replicate) %>%
  mutate(Tissue = factor(Tissue,levels=unique(Tissue)))
levels(attributes$Tissue) = c("SoL","Bark","Wood","SR","FR")

meta_data = attributes[grepl("T",attributes$Sample),]
write.csv(meta_data,paste0(inputPath,"00_META_DATA/META_DATA.csv"),row.names = FALSE)

# load data ####
#load meta data
meta_data <- read.csv(paste0(inputPath,"00_META_DATA/META_DATA.csv")) %>% mutate(
  Tissue = factor(Tissue, levels = unique(Tissue)),
  Stage = factor(Stage, levels = unique(Stage))
  )
# load count table
featureCounts <- read.delim(paste0(inputPath,"07_featureCounts/featureCounts_gene_UNIQUE.txt"), row.names=1, comment.char="#")

# change colnames in featureCounts
for (sample in attributes$Sample) {
  print(sample)
  colnames(featureCounts)[grepl(sample, colnames(featureCounts))] = sample
}

# count normalization and transformation ####

count_matrix = as.matrix(featureCounts[,colnames(featureCounts) %in% meta_data$Sample])
dds = DESeqDataSetFromMatrix(countData = count_matrix,
  colData = meta_data[meta_data$Sample %in% colnames(count_matrix),],
  design = ~ Stage + Tissue)
logvalues = as.data.frame(assay(vst(dds)))
write.csv(logvalues,
  file = paste0(outputPath,"VST.csv"))
norm = as.data.frame(counts(estimateSizeFactors(dds),normalized=TRUE))
write.csv(norm,
  file = paste0(outputPath,"NORMALIZED.csv"))

# LRT for each tissue ####

deseq_res = list()
for(tissue in levels(meta_data$Tissue)) {
  print(tissue)
  
  samples = meta_data[meta_data$Tissue == tissue,]$Sample
  
  dds = DESeqDataSetFromMatrix(
    countData = count_matrix[,samples],
    colData = meta_data[meta_data$Sample %in% samples,],
    design = ~ Stage
  )
  
  deseq_res[[tissue]] = data.frame(
    Tissue = tissue,
    as.data.frame(
      results(DESeq(dds,reduced = ~1, test = "LRT"),pAdjustMethod = "fdr")
    )
  ) %>%
    rownames_to_column(var = "Geneid")
  
} # tissue

deseq_res = bind_rows(deseq_res)

write.csv(
  deseq_res,
  file = paste0(outputPath,"LRT.csv"),
  row.names = FALSE
)



# second SR dataset ####
print("Analysis of the SR reproduction samples")
# variables ####
print("creating variables")

dataPath = "../data/RNASEQ_BULKING/"
outputPath = paste0(dataPath,"08_DESEQ2/")
figurePath = paste0(outputPath,"FIGURES")
dir.create(figurePath, recursive = TRUE)
dir.create(paste0(dataPath,"00_META_DATA"), recursive = TRUE)

# prepare meta data ####

meta_data_bulking = attributes[attributes$Tissue %in% c("SR","FR"),] %>%
  mutate(Experiment = str_sub(Sample,1,1))
write.csv(meta_data_bulking,paste0(dataPath,"00_META_DATA/META_DATA.csv"),row.names = FALSE)

# load data ####
metaData <- meta_data_bulking %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)),
         Tissue = factor(Tissue, levels = unique(Tissue)),
         Stage = factor(Stage, levels = unique(Stage)),
         Experiment = factor(Experiment, levels = unique(Experiment))
  )

# count normalization, transformation and batch correction ####

print("normalization")
count_matrix = as.matrix(featureCounts[,levels(metaData$Sample)])
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = metaData,
                             design = ~ Experiment + Stage + Tissue)
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
    batch = metaDataTissue$Experiment,
    design = model.matrix(data=metaDataTissue,~Stage)
  ))
}

corrected = bind_cols(corrected)
corrected = corrected[,levels(metaData$Sample)]
write.csv(corrected, file = paste0(outputPath,"VST_LIMMA.csv"))

# LRT per Tissue ####

print("LRT for each Tissue")

deseqRes = list()

for(tissue in unique(metaData$Tissue)) {
  print(tissue)
  metaDataTissue = metaData[metaData$Tissue == tissue,]
  samples = metaDataTissue$Sample
  
  dds = DESeqDataSetFromMatrix(countData = count_matrix[,samples], colData = metaDataTissue,
                               design = ~ Experiment + Stage)
  
  deseq = DESeq(dds,test="LRT",reduced = ~ Experiment)
  
  deseqRes[[tissue]] = data.frame(
    Tissue = tissue,
    as.data.frame(results(deseq))
  ) %>%
    rownames_to_column("Geneid")
}

deseqRes = bind_rows(deseqRes)
write.csv(deseqRes, file = paste0(outputPath,"LRT.csv"))
