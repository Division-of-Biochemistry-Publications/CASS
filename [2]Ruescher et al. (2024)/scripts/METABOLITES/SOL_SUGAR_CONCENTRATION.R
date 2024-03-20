# load packages ####
library(tidyverse)
library(emmeans)
library(car)
# variables ####
data_path = "../data/METABOLITES/SUGAR_STARCH/"
figure_path = "../data/METABOLITES/SUGAR_STARCH/FIGURES/"
# load data ####
data <- read.csv(paste0(data_path,"/LEAF_FW_SUGAR.csv")) %>%
  filter(
    !is.na(Weight)
  )
data$waterVol = ((data$Weight * (1-data$DM))) * (1/1000) # g = ml
# sugar ####
mMolSugarCalc =
  function(
    data,col1,col2,waterVol="waterVol",weightCol="Weight",VolCol="Volume",sugarType=c("Sucrose","Hexose","Glucose","Fructose"),
    VolWellSample=5,VolBuffer=200,r=0.3195
  ) {
    VolWell=VolBuffer+VolWellSample
    d=(VolWell/1000)*(1/(pi*(r**2)))
    extCoef=6220
    intermediateResult=(data[col2]-data[col1]) * (1/(extCoef*d)) * (data[VolCol]/1e6) * (VolWell/VolWellSample) * 
      (1/data[waterVol])
    if(sugarType == "Sucrose"){result = intermediateResult * (1/2)} else{
      result = intermediateResult}
    return(result)
  }
data["Glucose"] = mMolSugarCalc(data = data, col1 = "Blank", col2 = "HXK",sugarType = "Hexose")*1000
data["Fructose"] = mMolSugarCalc(data = data, col1 = "HXK", col2 = "PGI",sugarType = "Hexose")*1000
data["Sucrose"] = mMolSugarCalc(data = data, col1 = "PGI", col2 = "INV",sugarType = "Sucrose")*1000
# statistics ####
dir.create(
  path = paste0(figure_path,"/00_DIAGNOSTICS"),
  recursive = TRUE
)
input = data %>%
  gather(
    key = "Metabolite",
    value = "value",
    Glucose,Fructose,Sucrose
  ) %>%
  mutate(
    Stage = factor(Stage,levels=unique(Stage)),
    Metabolite = factor(Metabolite,levels=unique(Metabolite))
  )
anova_res = list()
tukey_res = list()
for(j in unique(input$Metabolite)) {
  print(j)
  model = glm(
    data = input[input$Metabolite == j,],
    formula = value ~ Stage,
    family= Gamma(link='log')
  )
  png(
    filename = paste0(figure_path,"/00_DIAGNOSTICS/SoL_FW_",j,".png"),
    units = "cm",
    width = 20,
    height = 20,
    res = 300
  )
  par(mfrow=c(2,2))
  plot(model)
  dev.off()
  tmp = as.data.frame(Anova(model,test="F", type = "II"))
  colnames(tmp) = c("Df","SS","F","p")
  tukey_res[[j]] = data.frame(
    Metabolite = j,
    as.data.frame(emmeans(model,pairwise ~ Stage)$contrasts)
  )
  anova_res[[j]] = data.frame(
    Metabolite = j,
    tmp %>% rownames_to_column("Index")
  ) 
} # j
anova_res = bind_rows( anova_res)
tukey_res = bind_rows( tukey_res)
anova_res$p_adj = NA
anova_res[anova_res$Index == "Stage","p_adj"] = p.adjust(anova_res[anova_res$Index == "Stage","p"],method = "fdr")
tukey_res = inner_join(
  anova_res %>% filter(Index == "Stage") %>% select(Metabolite,p_adj),
  tukey_res
) %>%
  rename(
    ANOVA_p = p_adj
  )
tukey_res[tukey_res$ANOVA_p > 0.05,]$p.value = 1
# plots
sig_code = list()
for(j in unique(tukey_res$Metabolite)) {
  tmp = tukey_res[tukey_res$Metabolite == j,]
  sig_code[[j]] = data.frame(
    Metabolite = j,
    rcompanion::cldList(
      comparison = tmp$contrast,
      p.value = tmp$p.value
    )
  ) %>%
    rename(
      Stage = Group
    )
} # j
sig_code = bind_rows(sig_code)
plot_df = data %>%
  gather(
    key = "Metabolite",
    value = "value",
    Glucose,Fructose,Sucrose
  ) 
plot_df = plot_df %>%
  group_by(Metabolite) %>%
  slice(
    which.max(value)
  ) %>%
  mutate(
    sig_dist = value * .1
  ) %>%
  select(
    Metabolite,sig_dist
  ) %>%
  right_join(
    plot_df
  ) %>%
  group_by(
    Metabolite,Stage
  ) %>%
  slice(
    which.max(value)
  ) %>%
  mutate(
    sig_height = value + sig_dist
  ) %>%
  right_join(
    plot_df
  ) %>%
  inner_join(sig_code) %>%
  mutate(
    Stage = factor(Stage,levels=unique(data$Stage)),
    Metabolite = factor(Metabolite,levels=c("Glucose","Fructose","Sucrose"))
  )
ggplot(
  plot_df,
  aes(
    x = Stage,
    y = value,
    fill = Tissue
  )
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.25,
    size = 0.2
  ) +
  geom_boxplot(
    width = 0.75,
    size = .25,
    outlier.shape = NA
  ) +
  geom_jitter(
    shape = 21,
    fill = "white",
    color = "black",
    size = 1,
    width = 0.25
  ) +
  geom_text(
    aes(
      y = sig_height,
      label = Letter
    ),
    size = 2.5
  ) +
  facet_grid(
    Metabolite ~ Tissue,
    scales = "free"
  ) +
  scale_fill_viridis_d(
    option = "plasma"
  ) +
  scale_y_continuous(
    limits = c(0,NA),
    expand = expansion(
      mult = c(0,0.15),
      add = 0
    )
  ) +
  labs(
    y = expression(bold(~"["*mg~g^-1~DW*"]"))
  ) +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 10, face = "bold", hjust = .5, vjust = 0.5),
    legend.position = "none",
    legend.justification = "top",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text.y = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.text.x = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave(
  path = figure_path,
  filename = "LEAF_SUGAR_CONCENTRATION.png",
  units = "cm",
  width = 18,
  height = 15,
  dpi = 600
)



# pblication plot ####



plot_df = data %>%
  filter(
    Tissue %in% c("SoL","Stem_LP","Stem_LC","SR","FR")
  ) %>%
  gather(
    key = "Metabolite",
    value = "value",
    Glucose,Fructose,Sucrose
  ) %>%
  mutate(
    Stage = factor(Stage,levels=unique(data$Stage)),
    Metabolite = factor(Metabolite,levels=c("Glucose","Fructose","Sucrose")),
    Tissue = factor(Tissue,levels=c("SoL","Stem_LP","Stem_LC","SR","FR"))
  )

maxValues = plot_df %>% group_by(Tissue,Metabolite,Stage) %>%  slice(which.max(value)) %>% mutate(sigDist = value + max(plot_df$value)*.1) %>%
  select(Stage,Metabolite,sigDist)
sigCode = inner_join(sig_code,maxValues) %>%
  mutate(
    Stage = factor(Stage,levels=unique(data$Stage)),
    Metabolite = factor(Metabolite,levels=c("Glucose","Fructose","Sucrose")),
    Tissue = factor(Tissue,levels=c("SoL","Stem_LP","Stem_LC","SR","FR"))
  )
levels(plot_df$Metabolite) = c("Glucose","Fructose","Sucrose","Starch") 



ggplot(
  plot_df,
  aes(
    x = Stage,
    y = value,
    fill = Tissue
  )
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.25,
    size = 0.2
  ) +
  geom_boxplot(
    width = 0.9,
    size = .25,
    outlier.shape = NA
  ) +
  geom_jitter(
    shape = 21,
    fill = "white",
    color = "black",
    size = 1,
    width = 0.25
  ) +
  geom_text(
    data = sigCode,
    aes(
      y = sigDist,
      label = Letter
    ),
    size = 2.5
  ) +
  facet_grid(
    cols=vars(Metabolite),
    scales = "free"
  ) +
  scale_fill_viridis_d(
    option = "plasma"
  ) +
  scale_y_continuous(
    limits = c(0,NA),
    expand = expansion(
      mult = c(0.05,0.15),
      add = 0
    )
  ) +
  labs(
    y = expression(bold(~"["*mmol~l^-1*"]"))
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 10, face = "bold", hjust = .5, vjust = 0.5),
    legend.position = "none",
    legend.justification = "top",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text.y = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.text.x = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.background = element_rect(fill = NA, color = NA),
    panel.border = element_blank()
  )
ggsave(
  path = "../../PUBLICATIONS/Transportome_2023/Draft3/SUPPLEMENT/FIGURES",
  filename = "S4_SoL_FW_SUGAR.tiff",
  units = "cm",
  width = 14,
  height = 8,
  dpi = 600
)
ggsave(
  path = "../../PUBLICATIONS/Transportome_2023/Draft3/SUPPLEMENT/FIGURES",
  filename = "S4_SoL_FW_SUGAR.png",
  units = "cm",
  width = 14,
  height = 8,
  dpi = 600
)

