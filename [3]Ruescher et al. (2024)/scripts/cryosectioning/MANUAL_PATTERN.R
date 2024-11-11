library(tidyverse)

wald = read.csv("../RNASEQ/08_DESEQ2/WALD_TISSUE.csv")

patterns = read.csv(paste0("../RNASEQ/11_GOI/patternsOfInterest.csv")) %>%
  filter(!Pattern %in% c("","MYB7"))

a = patterns %>% split(.$Pattern)

PATTERN_ENRICHMENT <- read.csv("../RNASEQ/11_GOI/PATTERN_ENRICHMENT.csv")

GO = read.csv("../RNASEQ/00_META_DATA/Mesculenta_BEST_HIT_ARAPORT11_GO.csv")

extractGO = function(term,pattern) {
  
  '
  term: str; GOID (GO:XXXXXXX)
  pattern: str; pattern name
  '
  
  tmp = GO[GO$GO == term,c("Geneid","AT_locusName","AT_Alias")] %>% unique()
  tmp = tmp[tmp$Geneid %in% patterns[patterns$Pattern == pattern,"Geneid"],]
  return(tmp)
  
}

waterDepr = extractGO("GO:0009414","KNOX1")
coldResp = extractGO("GO:0009409","KNOX1")
saltResp = extractGO("GO:0009651","KNOX1")
heatresp = extractGO("GO:0009266","KNOX1")
postEmbryoDev = extractGO("GO:0090698","KNOX1")

sum(waterDepr$Geneid %in% coldResp$Geneid)
sum(saltResp$Geneid %in% waterDepr$Geneid)


cellDiviosn = extractGO("GO:0051301","WOX14") 
geneExpression = extractGO("GO:0010468","WOX14")
mitosis = extractGO("GO:0000278","WOX14")
intracellularProtein = extractGO("GO:0006886","WOX14")

glucose = extractGO("GO:0009749","PXL")
a = extractGO("GO:0032544","PXL")
