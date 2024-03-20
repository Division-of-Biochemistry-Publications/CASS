# load packages ####

library(tidyverse)

# variables ####

data_path = "../data/METABOLITES/METABOLITE_PROFILES/"
figure_path = "../data/METABOLITES/METABOLITE_PROFILES/FIGURES/"

# create directories ####


dir.create(
  path = paste0(data_path),
  recursive = TRUE
)

# load data ####

meta_data <- read.csv(paste0(data_path,"00_META_DATA/META_DATA.csv"))

intensities <- read.csv(paste0(data_path,"01_RAW_DATA/INTENSITIES.csv"), row.names=1)

sugar_starch = read.csv("../data/METABOLITES/SUGAR_STARCH/NSC.csv") %>%
  filter(
    Sample %in% meta_data$Sample
  )


# normalize to internal standard ####

tmp = as.data.frame(t(intensities))

normalized_intensities = as.data.frame(t(tmp/tmp$Ribitol)) %>%
  rownames_to_column("Metabolite")

rm(tmp)

normalized_intensities = normalized_intensities[normalized_intensities$Metabolite != "Ribitol",]

# calculate calibration curve from standard ####
# use normalized intensities

standards = normalized_intensities[,grepl(colnames(normalized_intensities), pattern = "QC_")]   %>%
  mutate(Metabolite = normalized_intensities$Metabolite) %>%
  gather( # make long table
    key = "amount",
    value = "value",
    -Metabolite
  ) %>%
  mutate( 
    amount = as.numeric(gsub("QC_","",amount))
  )


model_data = list()


for (i in unique(standards$Metabolite)) {
  
  print(i)
  
  input = standards[
    standards$Metabolite == i,
  ]
  
  model_coefficients = lm(data = input,formula = value ~  amount)$coefficients
  
  model_data[[i]] = data.frame(
    Metabolite = i,
    intercept = model_coefficients[1],
    slope = model_coefficients[2],
    row.names = NULL
  )
  
  
}


model_data = bind_rows(model_data)

normalized_intensities = left_join(
  model_data, normalized_intensities
)


metabolite_profiles = data.frame( # results in µg
  Metabolite = normalized_intensities$Metabolite,
  ((
    normalized_intensities[,meta_data$Sample]-normalized_intensities$intercept)/ # x = (y-intercept/slope)
     normalized_intensities$slope
  )/1000 # from ng to µg
) %>%
  gather(
    key = "Sample",
    value = "value",
    -Metabolite
  ) %>%
  inner_join(meta_data[c("Sample", "Weight")]) %>%
  mutate(
    value = value/(Weight/1000)
  ) %>%
  spread(
    key = "Metabolite",
    value = "value"
  )  %>%
  select(-Weight)


metabolite_profiles[metabolite_profiles <= 0] = NA # replace negative values with NA

metabolite_profiles = metabolite_profiles[, colSums(is.na(metabolite_profiles)) != nrow(metabolite_profiles)]

# combine with sugar + starch data ####

# sugars overshoots in MS data; replace with our measurement; add starch


metabolite_profiles = inner_join(
  metabolite_profiles[,!colnames(metabolite_profiles) %in% c("Glucose", "Fructose", "Sucrose")],
  sugar_starch[
    ,c("Sample", "Glc", "Frc", "Suc", "Starch")
  ]
) %>%
  column_to_rownames("Sample") 
metabolite_profiles$Pyruvate = NULL # too many missing values
write.csv(
  as.data.frame(t(metabolite_profiles)),
  file = paste0(data_path,"METABOLITES.csv")
)

# ttest 

results = list()
data = as.data.frame((metabolite_profiles))
for(i in unique(meta_data$Tissue)){
  print(i)
  samples = meta_data[meta_data$Tissue == i,"Sample"]
  for(j in names(data)) {
    if(i != "Stem_LC" | j != "Trehalose") {
      print(j)
      y=data[samples,j]
      x=meta_data[meta_data$Tissue == i, "Stage"]
      tmp=t.test(formula=y~x,var.equal=FALSE)
      results[[i]][[j]] = data.frame(Tissue = i,Metabolite = j,
                                     t = tmp$statistic, p = tmp$p.value,row.names = NULL) 
    }
  }#j
  results[[i]] = bind_rows(results[[i]])
  results[[i]]$p_adj = p.adjust(results[[i]]$p,"fdr")
}#i
results = bind_rows(results)
write.csv(
  results,
  file = paste0(data_path,"METABOLITES_T_TEST.csv"),
  row.names = FALSE
)

