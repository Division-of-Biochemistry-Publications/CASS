# load packages ####

library(tidyverse)
library(UpSetR)
library(clusterProfiler)
library(ggvenn)

# variables ####

inputPath = "../RNASEQ/"
outputPath = "../RNASEQ/10_CRN/"

# load data ####
#load meta data
print("load data")
meta_data <- read.csv(paste0(inputPath,"09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv")) %>% mutate(
  Tissue = factor(Tissue, levels = unique(Tissue)),
  Distance = factor(Distance, levels = unique(Distance))
)

GO <- read.csv(paste0(inputPath,"00_META_DATA/Mesculenta_BEST_HIT_ARAPORT11_GO.csv"))


clusters <- read.csv(paste0(inputPath,"10_CRN/CLUSTERS.csv")) %>% rename("Geneid" = X) %>%
  filter(Cluster != -1) %>%
  mutate(
    clusterNameShort = paste(Tissue,clusterNameShort,sep=" ")
  )

intersections = read.csv(paste0(inputPath,"10_CRN/WALD_CLUSTER_INTERSECTIONS.csv")) %>% rename("Geneid" = X)

# generate upset plot ####
print("Cluster Upset Plot")
deg = split(
  clusters,
  clusters$clusterNameShort
)
# only Geneids
deg = lapply(deg, "[[", "Geneid")

# open png device
png(
  filename = paste0(outputPath,"FIGURES/UPSET.png"),
  res = 600,
  units = "cm",
  width = 12,
  height = 16
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
  mb.ratio = c(0.7,0.3),
  point.size = 1,
  line.size = 0.5
)
# close ong device
dev.off()

# cluster GO enrichment ####
print("clusterNameShort GO Enrichment")
# create empty list
enrichmentRes = list()

for(cluster in unique(clusters$clusterNameShort)) {
  
  print(cluster)
  
  genes = clusters[clusters$clusterNameShort == cluster ,"Geneid"]
  
  for(ontology in c("BP")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichment = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichmentRes[[cluster]][[ontology]] = data.frame(
        clusterNameShort = cluster,
        Ontology = ontology,
        enrichment@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichment))[1] != 0) {
        dotplot(
          enrichment,
          font.size = 8,
          title = cluster
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS/"),
          filename = gsub("/", "", paste0(cluster,"_",ontology,".png")),
          units = "cm",
          height = 14,
          width = 12,
          dpi = 600,
          device = "png"
        )
        
      }
    }
    
  } # ontology
  
  enrichmentRes[[cluster]] = bind_rows(enrichmentRes[[cluster]])
  
} # cluster

enrichmentRes = bind_rows(enrichmentRes)

write.csv(
  enrichmentRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"CLUSTER_GO_ENRICHMENT.csv"),
  row.names = FALSE
)

# intersections
print("Intersection GO Enrichment")

enrichmentRes = list()

for(cluster in unique(intersections$clusterIntersection)) {
  
  print(cluster)
  
  genes = intersections[intersections$clusterIntersection == cluster ,"Geneid"]
  
  for(ontology in c("BP")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichment = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichmentRes[[cluster]][[ontology]] = data.frame(
        clusterIntersection = cluster,
        Ontology = ontology,
        enrichment@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichment))[1] != 0) {
        dotplot(
          enrichment,
          font.size = 8,
          title = cluster
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS/INTERSECTIONS"),
          filename = gsub("/", "", paste0(cluster,"_",ontology,".png")),
          units = "cm",
          height = 14,
          width = 12,
          dpi = 600,
          device = "png"
        )
        
      }
    }
    
  } # ontology
  
  enrichmentRes[[cluster]] = bind_rows(enrichmentRes[[cluster]])
  
} # cluster

enrichmentRes = bind_rows(enrichmentRes)

write.csv(
  enrichmentRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"INTERSECTION_GO_ENRICHMENT.csv"),
  row.names = FALSE
)

# wald intersections
print("Intersection GO Enrichment")

enrichmentRes = list()

for(cluster in unique(intersections$waldClusterIntersection)) {
  
  print(cluster)
  
  genes = intersections[intersections$waldClusterIntersection == cluster ,"Geneid"]
  
  for(ontology in c("BP")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichment = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichmentRes[[cluster]][[ontology]] = data.frame(
        waldClusterIntersection = cluster,
        Ontology = ontology,
        enrichment@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichment))[1] != 0) {
        dotplot(
          enrichment,
          font.size = 8,
          title = cluster
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS/WALD_INTERSECTIONS"),
          filename = gsub("/", "", paste0(cluster,"_",ontology,".png")),
          units = "cm",
          height = 14,
          width = 12,
          dpi = 600,
          device = "png"
        )
        
      }
    }
    
  } # ontology
  
  enrichmentRes[[cluster]] = bind_rows(enrichmentRes[[cluster]])
  
} # cluster

enrichmentRes = bind_rows(enrichmentRes)

write.csv(
  enrichmentRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"WALD_INTERSECTION_GO_ENRICHMENT.csv"),
  row.names = FALSE
)


# WALD GO enrichment ####

print("WALD GO Enrichment")

outputPath = "../RNASEQ/08_DESEQ2/"

wald = read.csv(paste0(inputPath,"/08_DESEQ2/WALD_TISSUE.csv")) %>% rename("Geneid" = X) %>%
  na.omit()

wald = wald %>% filter((padj < 0.001) & (abs(log2FoldChange) > log2(1.5))) %>%
  split((.$log2FoldChange) > 0) %>% lapply(.,"[[","Geneid")
names(wald) = c("SR ↓", "SR ↑")

enrichmentRes = list()

for (type in names(wald)) {
  
  genes = wald[[type]]
  
  for(ontology in c("BP")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichment = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichmentRes[[type]][[ontology]] = data.frame(
        Type = type,
        Ontology = ontology,
        enrichment@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichment))[1] != 0) {
        dotplot(
          enrichment,
          font.size = 8,
          title = type
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS/"),
          filename = paste0(type,"_",ontology,".png"),
          units = "cm",
          height = 14,
          width = 12,
          dpi = 600,
          device = "png"
        )
        
      } # plot if
      
    } # if
  
  } # ontology
  
  enrichmentRes[[type]] = bind_rows(enrichmentRes[[type]])
    
} # type

enrichmentRes = bind_rows(enrichmentRes)


write.csv(
  enrichmentRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"WALD_GO_ENRICHMENT.csv"),
  row.names = FALSE
)


# patterns enrichment ####

print("Pattern of Interest Enrichment")

outputPath = "../RNASEQ/11_GOI/"

patternsOfInterest = read.csv(paste0(outputPath,"patternsOfInterest.csv")) %>%
  filter(Pattern != "")

enrichmentRes = list()

for(cluster in unique(patternsOfInterest$Pattern)) {
  
  print(cluster)
  
  genes = patternsOfInterest[patternsOfInterest$Pattern == cluster ,"Geneid"]
  
  for(ontology in c("BP")) {
    
    if(sum(genes %in% GO[GO$Ontology == ontology,]$Geneid) != 0) {
      
      print(ontology)
      
      enrichment = enricher(
        gene = genes,
        TERM2GENE = unique(GO[GO$Ontology == ontology,c("GO","Geneid")]),
        TERM2NAME = unique(GO[GO$Ontology == ontology,c("GO","Term")]),
        minGSSize = 0,
        maxGSSize = Inf,
        pAdjustMethod = "fdr"
      )
      
      enrichmentRes[[cluster]][[ontology]] = data.frame(
        Pattern = cluster,
        Ontology = ontology,
        enrichment@result,
        row.names = NULL
      )
      
      if(ontology == "BP" & dim(as.data.frame(enrichment))[1] != 0) {
        dotplot(
          enrichment,
          font.size = 8,
          title = cluster
        )
        ggsave(
          path = paste0(outputPath,"FIGURES/ENRICHPLOTS"),
          filename = gsub("/", "", paste0(cluster,"_",ontology,".png")),
          units = "cm",
          height = 14,
          width = 12,
          dpi = 600,
          device = "png"
        )
        
      }
    }
    
  } # ontology
  
  enrichmentRes[[cluster]] = bind_rows(enrichmentRes[[cluster]])
  
}
  
enrichmentRes = bind_rows(enrichmentRes)


write.csv(
  enrichmentRes %>% filter(p.adjust < 0.05),
  file = paste0(outputPath,"PATTERN_ENRICHMENT.csv"),
  row.names = FALSE
)