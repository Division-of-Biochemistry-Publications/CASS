# load packages ####

library(tidyverse)
library(UpSetR)
library(clusterProfiler)

# variables ####

inputPath = "../RNASEQ/"
outputPath = "../RNASEQ/12_PROMOTERS/"

# load data ####
#load meta data
print("load data")
meta_data <- read.csv(paste0(inputPath,"09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv")) %>% mutate(
  Tissue = factor(Tissue, levels = unique(Tissue)),
  Distance = factor(Distance, levels = unique(Distance))
)

clusters <- read.csv(paste0(inputPath,"10_CRN/CLUSTERS.csv")) %>% rename("Geneid" = X) %>%
  filter(Cluster != -1) %>%
  mutate(
    clusterNameShort = paste(Tissue,clusterNameShort,sep=" ")
  )

intersections = read.csv(paste0(inputPath,"10_CRN/WALD_CLUSTER_INTERSECTIONS.csv")) %>% rename("Geneid" = X)

patterns = read.csv(paste0(inputPath,"11_GOI/patternsOfInterest.csv")) %>%
  filter(Pattern != "")

bestPeaks = read.table(paste0(outputPath,"FIMO/best_site.narrowPeak")) %>%
  select(1,4,8,9) %>%
  mutate(V1 = str_sub(V1,1,-4))
names(bestPeaks) = c("Geneid","TF","p","q")
bestPeaks = bestPeaks %>% mutate(TF = paste0(TF,".v8.1"))

synonyms = read.delim(paste0(outputPath,"Mesculenta_671_v8.1.synonym.txt"), na.strings=c("","NA"), header=FALSE) %>%
  mutate(V1 = sub("\\.[^.]*$", "", V1)) %>%
  unique()
names(synonyms) = c("locusName","old1","old2")

TFs = read.table(paste0(outputPath,"Mes_TF_list.txt"),header=TRUE) %>%
  select(-TF_ID) %>%
  filter(!duplicated(Gene_ID)) %>%
  rename("TF" = Gene_ID)

oldLoci = TFs[TFs$TF %in% c(synonyms$V2,synonyms$V3),"TF"]
newLoci = synonyms[synonyms$V2 %in% oldLoci | synonyms$V3 %in% oldLoci,]
newLoci = newLoci[order(match(newLoci$V2,oldLoci)),]
TFs[TFs$TF %in% oldLoci,"TF"] = newLoci$V1
TFs = TFs %>% mutate(TF = paste0(TF,".v8.1"))

TFs = TFs[order(TFs$Family, TFs$TF),]
TFs$Family = gsub("/","_",TFs$Family)
write.csv(TFs,paste0(outputPath,"TF_LIST_CLEAN.csv"), row.names = FALSE)

TFs = TFs %>%
  mutate(name = paste0(TF," (",Family,")")) %>%
  unique()

tfBindingSites = inner_join(bestPeaks,TFs)

# TF family enrichment ####

enrichmentRes = list()

for(pattern in unique(patterns$Pattern)) {
  
  print(pattern)
  
  genes = patterns[patterns$Pattern == pattern,"Geneid"]
  
  enrichment = enricher(
    gene = genes,
    TERM2GENE = TFs[,c("Family","TF")],
    TERM2NAME = TFs[,c("Family","Family")],
    minGSSize = 0,
    maxGSSize = Inf,
    pAdjustMethod = "fdr"
  ) 
  
  enrichmentRes[[pattern]] = data.frame(
    Pattern = pattern,
    enrichment@result,
    row.names = NULL
  ) %>% filter(p.adjust < 0.05)
  
  if(dim(enrichmentRes[[pattern]])[1] != 0) {
    
    
    dotplot(
      enrichment,
      font.size = 8,
      title = pattern
    )
    ggsave(
      path = paste0(outputPath,"FIGURES/ENRICHPLOTS/PATTERNS/TFs"),
      filename = gsub("/", "", paste0(pattern,".png")),
      units = "cm",
      height = 14,
      width = 12,
      dpi = 600,
      device = "png"
    )
    
  }
  
}

enrichmentRes = bind_rows(enrichmentRes)

write.csv(
  enrichmentRes,
  file = paste0(outputPath,"PATTERN_TF_ENRICHMENT.csv"),
  row.names = FALSE
)

# pattern TFBS enrichment ####

enrichmentRes = list()

for(pattern in unique(patterns$Pattern)) {
  
  print(pattern)
  
  genes = patterns[patterns$Pattern == pattern,"Geneid"]
  
  enrichment = enricher(
    gene = genes,
    TERM2GENE = tfBindingSites[,c("TF","Geneid")],
    TERM2NAME = tfBindingSites[,c("TF","name")],
    minGSSize = 0,
    maxGSSize = Inf,
    pAdjustMethod = "fdr"
  )
  
  enrichmentRes[[pattern]] = data.frame(
    Pattern = pattern,
    enrichment@result,
    row.names = NULL
  )
  
  dotplot(
    enrichment,
    font.size = 8,
    title = pattern
  )
  ggsave(
    path = paste0(outputPath,"FIGURES/ENRICHPLOTS/PATTERNS/TFBS"),
    filename = gsub("/", "", paste0(pattern,".png")),
    units = "cm",
    height = 14,
    width = 12,
    dpi = 600,
    device = "png"
  )
  
}

enrichmentRes = bind_rows(enrichmentRes) %>% filter(p.adjust < 0.05)

write.csv(
  enrichmentRes,
  file = paste0(outputPath,"PATTERN_TFBS_ENRICHMENT.csv"),
  row.names = FALSE
)

plotTFsData = TFs

plotEnrichmentData = enrichmentRes %>% 
  select(Pattern, ID, Description, p.adjust) 
names(plotEnrichmentData) = c("Pattern","TF","Name","pAdj")

enrichedTFs = sort(unique(plotEnrichmentData$TF))
profiles = unique(patterns$Pattern) 

plotData = list()


for(pattern in profiles) {
  
  print(pattern)
  
  tmp = data.frame(matrix(nrow=length(enrichedTFs),ncol=3),
                   row.names = enrichedTFs)
  names(tmp) = c("Family","Pattern","value")

  for(tf in unique(plotEnrichmentData$TF)) {
  
    print(tf)
    
    tmp[tf,"Family"] = unique(plotTFsData[plotTFsData$TF == tf,"Family"])[1]
    tmp[tf,"Pattern"] = pattern
    
    enriched = tf %in% plotEnrichmentData[plotEnrichmentData$Pattern == pattern,"TF"]
    expressed = tf %in% patterns[patterns$Pattern == pattern,"Geneid"]
    
    value = sum(enriched,expressed, all(enriched,expressed) )
    
    if( value == 1 ) {
      
      value = which.max(c(enriched,expressed))
      
    } # if
    
    tmp[tf,"value"] = value
    
  } # tf
  
  tmp = tmp %>% rownames_to_column("TF")
  plotData[[pattern]] = tmp
  
}

plotData = bind_rows(plotData)

translation = data.frame(
  value = c(0,1,2,3),
  Name = c("None","Enriched","Expressed","Enriched + Expressed")
)

plotData = inner_join(plotData,translation) %>%
  mutate(
    Name = factor(Name, levels = c("None","Expressed","Enriched","Enriched + Expressed")),
    Pattern = factor(Pattern, levels = c("MYB46","PXL", "MYB7","KNOX1","WOX14"))
  )

ggplot(plotData,aes(x=Pattern,y=TF)) +
  geom_tile(aes(fill=Name),color="black", linewidth=.1) +
  facet_grid( rows = vars(Family), scales = "free", space = "free", switch="y") +
  scale_fill_manual(values = c("grey90","blue","gold2","red")) +
  scale_y_discrete(limits=rev) +
  labs(y = "Transcription factor", x = "Expression profile") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust =1 , vjust =1),
    plot.title = element_text(size = 10, face = "bold", hjust = .5, vjust = 0.5),
    legend.position = "bottom",
    legend.justification = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 8, angle=0, hjust = 1, vjust = .5, face = "bold"),
    strip.placement = "outside",
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
ggsave(
  path = paste0(outputPath,"FIGURES/"),
  filename = "Enriched_TFS.png",
  units = "cm",
  height = 24,
  width = 10,
  dpi = 600,
  device = "png"
)
