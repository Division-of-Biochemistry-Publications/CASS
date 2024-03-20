# set WD ####

setwd("../data/PHYLOGENIES")

# load packages ####

library(tidyverse)

# execute ####

files = list.files(
  path = "./01_CDS/03_TRIMAL/MAFFT",
  pattern = ".fa"
)

dir.create(paste0("./01_CDS/01_SEQUENCES/FILTERED/"))

manual_removal = c(
  "Manes.18G106451.1" # SWEET
  ,"Potri.003G136800.1" # TPT
)


b = list()
c = list()

for (i in files) {
  
  CDS = Biostrings::readDNAStringSet(
    paste0("./01_CDS/01_SEQUENCES/",i)
  )
  
  trimal <- Biostrings::readDNAStringSet(
    paste0("./01_CDS/03_TRIMAL/MAFFT/",i)
  )
  a = as.data.frame(trimal)
  
  for(j in rownames(a)[!grepl("AT",rownames(a))]) {
    
    b[[i]][[j]] = data.frame(
      file = i,
      name = j,
      pct_gap = (nchar(a[j,"x"])-nchar(gsub("[^a-zA-Z]", "", a[j,"x"])))/nchar(a[j,"x"])
    )
    
  }
  
  
  b[[i]] = bind_rows(b[[i]])
  c[[i]] = c(
    manual_removal,
    b[[i]][b[[i]]$pct_gap > 0.50,"name"]
  )
  
  if(i != "VATPase.fa") {
    
    CDS[c[[i]]] = NULL
    
  }
  
  Biostrings::writeXStringSet(
    CDS,
    filepath =  paste0("./01_CDS/01_SEQUENCES/FILTERED/",i)
  )
  
}

c = unique(unlist(c))

blast_names <- read.csv("./05_NOMENCLATURE/BLAST_NAMES.csv")

kept = blast_names %>%
  filter(
    !transcriptName %in% c
  )

write.csv(
  kept,
  file = "./05_NOMENCLATURE/KEPT.csv",
  row.names = FALSE
)

removed = blast_names %>%
  filter(
    transcriptName %in% c
  )

write.csv(
  removed,
  file = "./05_NOMENCLATURE/REMOVED.csv",
  row.names = FALSE
)

