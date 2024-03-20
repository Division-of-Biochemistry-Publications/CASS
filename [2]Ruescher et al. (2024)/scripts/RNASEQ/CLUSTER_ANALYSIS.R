# load packages ####

library(tidyverse)
library(UpSetR)
library(clusterProfiler)

# variables ####

inputPath = "../data/RNASEQ/"
outputPath = "../data/RNASEQ/10_CRN/"

# load data ####
#load meta data
meta_data <- read.csv(paste0(inputPath,"00_META_DATA/META_DATA.csv")) %>% mutate(
  Tissue = factor(Tissue, levels = unique(Tissue)),
  Stage = factor(Stage, levels = unique(Stage))
)

GO <- read.csv(paste0(inputPath,"00_META_DATA/Mesculenta_BEST_HIT_ARAPORT11_GO.csv"))
KO = read.csv(paste0(inputPath,"00_META_DATA/Mesculenta_KEGG_PATHWAY.csv"))
clusters <- read.csv(paste0(inputPath,"10_CRN/CLUSTERS.csv")) %>% rename("Geneid" = X) %>%
  mutate(
    Cluster = paste(Tissue,Cluster,sep="_")
  )

# generate upset plot ####
print("Cluster Upset Plot")
deg = split(
  clusters,
  clusters$Cluster
)
# only Geneids
deg = lapply(deg, "[[", "Geneid")

# open png device
png(
  filename = paste0(outputPath,"FIGURES/UPSET.png"),
  res = 600,
  units = "cm",
  width = 8,
  height = 11
)
# upset plot
upset(
  fromList(deg),
  order.by = "freq",
  nsets = length(deg),
  keep.order = TRUE,
  nintersects = 25,
  show.numbers = FALSE,
  text.scale = c(.9, 0.75, 0.65, 0.75, .9, 0.5),
  sets.x.label = "Cluster size",
  mainbar.y.label = "Cluster intersection size",
  mb.ratio = c(0.6,0.4),
  point.size = 1,
  line.size = 0.5
)
# close ong device
dev.off()

# cluster GO enrichment ####
print("GO Enrichment")
# create empty list
enrichementRes = list()

for(cluster in unique(clusters$Cluster)) {
  
  print(cluster)
  
  genes = clusters[clusters$Cluster == cluster ,"Geneid"]
  
  for(ontology in c("BP","CC","MF")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichement = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichementRes[[cluster]][[ontology]] = data.frame(
        Cluster = cluster,
        Ontology = ontology,
        enrichement@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichement))[1] != 0) {
        dotplot(
          enrichement,
          font.size = 8,
          title = cluster
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS/GO"),
          filename = paste0(cluster,"_",ontology,".png"),
          units = "cm",
          height = 12,
          width = 12,
          dpi = 600,
          device = "png"
        )
      
    }
  }
    
  } # ontology
  
  enrichementRes[[cluster]] = bind_rows(enrichementRes[[cluster]])
  
} # cluster

enrichementRes = bind_rows(enrichementRes)

write.csv(
  enrichementRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"CLUSTER_GO_ENRICHMENT.csv"),
  row.names = FALSE
)



# KO enrichment ####
print("KEGG Pathway Enrichment")
enrichementRes = list()

for(i in unique(clusters$Cluster)) {
  
  print(i)
  
  genes = clusters[clusters$Cluster == i ,"Geneid"]
  
  if(sum(genes %in% KO$Geneid) != 0) {
    
    enrichement = enricher(
      gene = genes,
      TERM2GENE = unique(KO[,c("Pathway","Geneid")]),
      TERM2NAME = unique(KO[,c("Pathway","Name")]),
      minGSSize = 0,
      maxGSSize = Inf,
      pAdjustMethod = "fdr"
    )
    
    enrichementRes[[i]] = data.frame(
      Cluster = i,
      enrichement@result,
      row.names = NULL
    )
    
    if(dim(as.data.frame(enrichement))[1] != 0) {
      
      dotplot(
        enrichement,
        font.size = 8,
        title = i
      )
      ggsave(
        path = paste0(outputPath,"FIGURES/ENRICHPLOTS/KO"),
        filename = paste0(i,".png"),
        units = "cm",
        height = 12,
        width = 12,
        dpi = 600,
        device = "png"
      )
    } # if 2
    
  } #if 1
  
} # i

enrichementRes = bind_rows(enrichementRes)

write.csv(
  enrichementRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"CLUSTER_KO_ENRICHMENT.csv"),
  row.names = FALSE
)

