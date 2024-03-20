# set WD ####

setwd("../data/PHYLOGENIES")


library(tidyverse)
#library(Biostrings)
#library(seqRFLP)

# create directories

dir.create(
  "./01_CDS/01_SEQUENCES",
  recursive = TRUE
)

dir.create(
  "./02_PEPTIDE/01_SEQUENCES",
  recursive = TRUE
)

dir.create(
  "./05_NOMENCLATURE",
  recursive = TRUE
)


AT_PHYLOGENY_REFERENCE <- read.csv("./00_META_DATA/AT_PHYLOGENY_REFERENCE.csv")
family = read.csv("./00_META_DATA/FAMILY_PFAM.csv") %>%
  select(family,PFAM)

# prepare list of organisms ####

organisms <- c(
  "Athaliana"
  #,"Brapa"
  ,"Mesculenta"
  #,"Rcommunis"
  #,"Mtruncatula"
  #,"Gmax"
  ,"Ptrichocarpa"
  #,"Pdeltoides"
  #,"Stuberosum"
  #,"Slycopersicum"
  #,"Bvulgaris"
  #,"Osativa"
  #,"Zmays"
  #,"Ppatens"
  #,"Mpolymorpha"
)

genes_list = list()
PEPTIDE_sequences <- list()
CDS_sequences <- list()

for(j in organisms) {
  
  
  print(j)
  
  # PFAM
  # The PFAM files contains the information about the protein domains per full-length peptide of each gene and the domain's position in the peptide
  
  PFAM <- read.table(
    paste0("./00_META_DATA/FUNCTIONAL_ANNOTATION/PFAM/",j,"_PFAM.txt"), 
    quote="\""
  )
  
  PFAM <- separate( # HMMSCAN adds useless information to the PFAM ID after a "."; remove it via separate
    PFAM,
    col = "V2",
    into = c("PFAM", "atachement"),
    sep = "[.]"
  )
  
  PFAM <- data.frame(
    peptideName = PFAM$V4,
    transcriptName = gsub(x = PFAM$V4, pattern = "[.]p", replacement = ""),
    PFAM = PFAM$PFAM,
    start = PFAM$V18, # start of the PFAM domain
    end = PFAM$V19 # end of the PFAM domain
  )
  
  PFAM[,"length"] = PFAM$end-PFAM$start
  PFAM[,"CDS_start"] = (PFAM$start*3)-2 # CDS is 3 times longer, than AA sequence; 1*3=3, but the first NT is 1; hence (1*3)-2=1
  PFAM[,"CDS_end"] = (PFAM$end*3) # CDS is 3 times longer
  PFAM[,"CDS_length"] = PFAM$CDS_end-PFAM$CDS_start
  
  PFAM$transcriptName <- gsub(x = PFAM$transcriptName, pattern = "_P0", replacement = "_T0")
  PFAM$transcriptName <- gsub(x = PFAM$transcriptName, pattern = "XP", replacement = "XM")
  
  
  # BLAST
  
  if(j != "Athaliana"){
    
    
    
    BLASTP <- read.table( # load the blastp versus athaliana file
      paste0("./00_META_DATA/FUNCTIONAL_ANNOTATION/BLAST/",j,"_BLASTP.txt"),
      header = TRUE,
      quote="\""
    )
    
    BLASTP <- BLASTP %>%
      group_by(qseqid) %>% 
      slice(which.min(evalue)) %>% #select only rows containing the lowest e-value per peptideName
      ungroup()
    
    BLASTP <- BLASTP[BLASTP$evalue <= 10**-3 , ] # set minimum e-value
    
    BLASTP <- data.frame( # make new data frame wit all necessary information
      peptideName = BLASTP$qseqid,
      transcriptName = gsub(x = BLASTP$qseqid, pattern = "[.]p", replacement = ""),
      AT_locusName = str_split_fixed(BLASTP$sseqid, pattern = "[.]", n=2)[,1],
      AT_transcriptName = BLASTP$sseqid,
      AT_peptideName = BLASTP$sseqid
    )
    
    BLASTP$transcriptName <- gsub(x = BLASTP$transcriptName, pattern = "_P0", replacement = "_T0") # different annotation in some genomes
    BLASTP$transcriptName <- gsub(x = BLASTP$transcriptName, pattern = "XP", replacement = "XM") # different annotation in some genomes
    
  } else{
    
    print("Athaliana")
    
  } # else
  
  # Sequences
  # Both the fasta files containingpeptides and CDS sequences of the genes will be loaded
  # a data frame will be produced containing the ID, the length and the sequence
  
  PEPTIDE <- Biostrings::readAAStringSet(
    paste0("./00_META_DATA/PEPTIDE/",j,"_PEPTIDE.fa")
  )
  
  PEPTIDE <- data.frame(
    peptideName = str_split_fixed(rownames(as.data.frame(PEPTIDE)), n=Inf, pattern = " ")[,1], # extrats only the first part of the name; the actual ID
    length = PEPTIDE@ranges@width, # contains the information about the length of each sequence
    sequence = as.data.frame(PEPTIDE)$x # as.data frame produces a df of rownames=ID and x=Sequence; only take the sequence (x)
  )
  
  
  CDS <- Biostrings::readDNAStringSet(
    paste0("./00_META_DATA/CDS/",j,"_CDS.fa")
  )
  
  CDS <- data.frame(
    transcriptName = str_split_fixed(rownames(as.data.frame(CDS)), n=Inf, pattern = " ")[,1], # extrats only the first part of the name; the actual ID
    length = CDS@ranges@width, # contains the information about the length of each sequence
    sequence = as.data.frame(CDS)$x # as.data frame produces a df of rownames=ID and x=Sequence; only take the sequence (x)
  )
  
  
  
  for(i in unique(family$family)) {
    
    print(i)
    
    PFAM_ID <- family[family$family == i, "PFAM"]
    
    AT_genes <- AT_PHYLOGENY_REFERENCE[AT_PHYLOGENY_REFERENCE$family == i,"locusName"] # get list of AT_genes for filtering BLASTP
    
    # get gene names per organisms
    
    if(j == "Athaliana") {
      
      
      if(PFAM_ID != "None") {
        
        genes <- AT_PHYLOGENY_REFERENCE[AT_PHYLOGENY_REFERENCE$peptideName %in% PFAM[PFAM$PFAM %in% PFAM_ID,"peptideName"]
                                        & AT_PHYLOGENY_REFERENCE$locusName %in% AT_genes
                                        , ]
        
        remove = AT_PHYLOGENY_REFERENCE[AT_PHYLOGENY_REFERENCE$family == i & !AT_PHYLOGENY_REFERENCE$peptideName %in% genes, ]$peptideName
        
        AT_genes <- AT_genes[!AT_genes %in% remove]
        
      } else {
        
        genes <- AT_PHYLOGENY_REFERENCE[AT_PHYLOGENY_REFERENCE$family == i,]
        
      }
      
    } else {
      
      
      if(PFAM_ID != "None") {
        
        genes <- BLASTP[
          (BLASTP$AT_locusName %in% AT_genes) & # filter peptides that blast to the AT_genes
            (BLASTP$peptideName %in% PFAM[PFAM$PFAM %in% PFAM_ID,"peptideName"] ), # filter peptides with the correct PFAM domain
        ] 
        
      } else {
        
        genes <- BLASTP[
          (BLASTP$AT_locusName %in% AT_genes), # filter peptides that blast to the AT_genes
        ]
        
      }
      
      if(length(genes$peptideName) != 0) {
        genes_list[[j]][[i]] = data.frame(
          Organism = j,
          Family = i,
          genes
        )
      }
    }
    

    
    # Sequences
    
    
    FAMILY_PEPTIDE <- PEPTIDE[PEPTIDE$peptideName %in% genes$peptideName,]
    FAMILY_CDS <- CDS[CDS$transcriptName %in% genes$transcriptName,]
    
    for(k in FAMILY_PEPTIDE$peptideName) {
      
      print(k)
      
      tmp = Biostrings::AAString(
        FAMILY_PEPTIDE[FAMILY_PEPTIDE$peptideName == k,"sequence"],
        start = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$peptideName == k),"start"][1],
        nchar = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$peptideName == k),"length"][1]
      )
      
      FAMILY_PEPTIDE[FAMILY_PEPTIDE$peptideName == k,"PFAM_length"] = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$peptideName == k),"length"][1]
      FAMILY_PEPTIDE[FAMILY_PEPTIDE$peptideName == k,"PFAM_sequence"] = as.character(tmp)
      
      
    } # end of k
    
    
    for(k in FAMILY_CDS$transcriptName) {
      
      print(k)
      
      tmp = Biostrings::DNAString(
        FAMILY_CDS[FAMILY_CDS$transcriptName == k,"sequence"],
        start = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$transcriptName == k),"CDS_start"][1],
        nchar = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$transcriptName == k),"CDS_length"][1]
      )
      
      FAMILY_CDS[FAMILY_CDS$transcriptName == k,"PFAM_length"] = PFAM[(PFAM$PFAM == PFAM_ID) & (PFAM$transcriptName == k),"CDS_length"][1]
      FAMILY_CDS[FAMILY_CDS$transcriptName == k,"PFAM_sequence"] = as.character(tmp)
      
      
    } # end of k
    
    PEPTIDE_sequences[[i]][[j]] <- FAMILY_PEPTIDE
    CDS_sequences[[i]][[j]] <- FAMILY_CDS
    
  } # end of i
  
  genes_list[[j]] = bind_rows(genes_list[[j]])
  
} # end of j

genes_list = bind_rows(genes_list)

genes_list = inner_join(
  genes_list,
  AT_PHYLOGENY_REFERENCE %>% rename(
    AT_locusName = locusName, AT_transcriptName = transcriptName,AT_peptideName = peptideName
    ) %>%
    select(
      AT_locusName,name
    )
) %>% filter(
  Organism == "Mesculenta"
)

genes_list$name = gsub(
  pattern = "At",replacement = "Me",genes_list$name
)

genes_list$name = gsub(
  pattern = "ANAC",replacement = "MeNAC",genes_list$name
)

genes_list$name = gsub(
  pattern = "ATAF",replacement = "MeTAF",genes_list$name
)


genes_list = split(genes_list,genes_list$name)

for(i in names(genes_list)) {
  
  print(i)
  
  tmp = genes_list[[i]]
  
  tmp = tmp[order(tmp$transcriptName),]
  
  for(j in 1:(length(tmp$transcriptName)))
    
    genes_list[[i]]$name[j] = paste0(tmp$name[j],letters[j])
  
}

genes_list = bind_rows(genes_list)

genes_list = genes_list[order(genes_list$name),]
genes_list = genes_list[order(genes_list$Family),]



write.csv(
  genes_list,
  file = "./05_NOMENCLATURE/BLAST_NAMES.csv",
  row.names = FALSE
)







AT_PHYLOGENY_REFERENCE[,"name"] <- paste(
  AT_PHYLOGENY_REFERENCE$peptideName,
  AT_PHYLOGENY_REFERENCE$name,
  sep = "_"
)

for (i in names(PEPTIDE_sequences)) {
  
  tmp <- bind_rows(PEPTIDE_sequences[[i]])
  
  tmp <- left_join(
    tmp,
    AT_PHYLOGENY_REFERENCE[,c("peptideName", "name")]
  )
  
  full_sequence <- tmp[,
    c( "peptideName", "sequence", "name" )
  ]
  
  full_sequence[
    !is.na(full_sequence$name),"peptideName"
  ] <- full_sequence[
    !is.na(full_sequence$name),"name"
  ]
  
  df = full_sequence[,c("peptideName", "sequence")]
  names(df) = c("seq.name", "seq.text")
  
  phylotools::dat2fasta(
    dat = df,
    outfile = paste0("./02_PEPTIDE/01_SEQUENCES/",i,".fa")
  )
  
}




for (i in names(CDS_sequences)) {
  
  tmp <- bind_rows(CDS_sequences[[i]])
  
  tmp <- left_join(
    tmp,
    AT_PHYLOGENY_REFERENCE[,c("transcriptName", "name")]
  )
  
  full_sequence <- tmp[,
    c( "transcriptName", "sequence", "name" )
  ]
  
  full_sequence[
    !is.na(full_sequence$name),"transcriptName"
  ] <- full_sequence[
    !is.na(full_sequence$name),"name"
  ]
  
  df = full_sequence[,c("transcriptName", "sequence")]
  names(df) = c("seq.name", "seq.text")
  
  phylotools::dat2fasta(
    dat = df,
    outfile = paste0("./01_CDS/01_SEQUENCES/",i,".fa")
  )

  
}