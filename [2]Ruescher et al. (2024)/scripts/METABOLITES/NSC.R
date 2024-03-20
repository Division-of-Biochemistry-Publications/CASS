# load packages ####

library(tidyverse)
library(emmeans)
library(car)

# variables ####

data_path = "../data/METABOLITES/SUGAR_STARCH/"
figure_path = "../data/METABOLITES/SUGAR_STARCH/FIGURES/"


p_thrsh = 0.05

M_hexose = 180.156 # g/mol

M_sucrose = 342.30 # g/mol

ext = 6220 # l/(mol*cm)

r = 0.639/2 # appr., since well is slightly conical; cm

# load data ####

data <- read.csv(paste0(data_path,"/NSC_RAW.csv")) %>%
  filter(
    !is.na(Weight)
  )

# sugar ####

sugar = data[data$Type == "sugar",]

sugar[,"d"] = ((sugar$V_buffer + sugar$V_measured)/(pi*(r**2))) #[(cm³/cm²)] = [cm]

sugar$Glc = (sugar$hxk - sugar$blank) * # no unit
  (1/(ext * sugar$d)) * # [mol/l]
  (sugar$V_total/1000) * # [(mol/l) * l] = [mol]
  ( (sugar$V_buffer + sugar$V_measured) / sugar$V_measured) * # [mol * (ml/ml)] = [mol]
  sugar$Digestion *
  (M_hexose*1000) * # [mol * (mg/mol)] = [mg] in total extract
  (1/sugar$Weight)  # [mg/g] 


sugar$Frc = (sugar$pgi - sugar$hxk) * # no unit
  (1/(ext * sugar$d)) * # [mol/l]
  (sugar$V_total/1000) * # [(mol/l) * l] = [mol]
  ( (sugar$V_buffer + sugar$V_measured) / sugar$V_measured) * # [mol * (ml/ml)] = [mol]
  sugar$Digestion *
  (M_hexose*1000) * # [mol * (mg/mol)] = [mg] in total extract
  (1/sugar$Weight)  # [mg/g] 


sugar$Suc = (0.5*(sugar$inv - sugar$pgi)) * # no unit
  (1/(ext * sugar$d)) * # [mol/l]
  (sugar$V_total/1000) * # [(mol/l) * l] = [mol]
  ( (sugar$V_buffer + sugar$V_measured) / sugar$V_measured) * # [mol * (ml/ml)] = [mol]
  sugar$Digestion *
  (M_hexose*1000) * # [mol * (mg/mol)] = [mg] in total extract
  (1/sugar$Weight)  # [mg/g] 


sugar = sugar %>%
  select(
    Sample,Stage,Tissue,Glc,Frc,Suc
  ) %>%
  group_by(
    Sample,Stage,Tissue
  ) %>%
  summarise_if(
    is.numeric,
    .funs = list(mean)
  )




# starch ####

starch = data[data$Type == "starch",]

starch[,"d"] = ((starch$V_buffer + starch$V_measured)/(pi*(r**2))) #[(cm³/cm²)] = [cm]

starch$Starch = (starch$hxk - starch$blank) * # no unit
  (1/(ext * starch$d)) * # [mol/l]
  (starch$V_total/1000) * # [(mol/l) * l] = [mol]
  ( (starch$V_buffer + starch$V_measured) / starch$V_measured) * # [mol * (ml/ml)] = [mol]
  starch$Digestion *
  (M_hexose*1000) * # [mol * (mg/mol)] = [mg] in total extract
  (1/starch$Weight)  # [mg/g] 

starch = starch  %>%
  select(
    Sample,Stage,Tissue,Starch
  ) %>%
  group_by(
    Sample,Stage,Tissue
  ) %>%
  summarise_if(
    is.numeric,
    .funs = list(mean)
  )



# combine ####

data = inner_join(sugar,starch)

data[data<0] = NA

write.csv(data,paste0(data_path,"/NSC.csv"),row.names = FALSE)

# statistics ####

dir.create(
  path = paste0(figure_path,"/00_DIAGNOSTICS"),
  recursive = TRUE
)

input = data %>%
  gather(
    key = "Metabolite",
    value = "value",
    Glc,Frc,Suc,Starch
  ) %>%
  mutate(
    Stage = factor(Stage,levels=unique(Stage)),
    Metabolite = factor(Metabolite,levels=unique(Metabolite))
  )

anova_res = list()
tukey_res = list()

for(i in unique(input$Tissue)) {
  
  print(i)
  
  for(j in unique(input$Metabolite)) {
    
    print(j)
    
    
    model = glm(
      data = input[input$Tissue == i & input$Metabolite == j,],
      formula = value ~ Stage,
      family= Gamma(link='log')
    )
    
    png(
      filename = paste0(figure_path,"/00_DIAGNOSTICS/",i,"_",j,".png"),
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
    
    tukey_res[[i]][[j]] = data.frame(
      Tissue = i,
      Metabolite = j,
      as.data.frame(emmeans(model,pairwise ~ Stage)$contrasts)
    )
    
    anova_res[[i]][[j]] = data.frame(
      Tissue = i,
      Metabolite = j,
      tmp %>% rownames_to_column("Index")
    ) 
    
  } # j
  
  anova_res[[i]] = bind_rows( anova_res[[i]])
  tukey_res[[i]] = bind_rows( tukey_res[[i]])
  
} # i 

anova_res = bind_rows(anova_res)
anova_res$p_adj = NA
anova_res[anova_res$Index == "Stage","p_adj"] = p.adjust(anova_res[anova_res$Index == "Stage","p"],method = "fdr")


tukey_res = bind_rows(tukey_res)

tukey_res = inner_join(
  anova_res %>% filter(Index == "Stage") %>% select(Tissue,Metabolite,p_adj),
  tukey_res
) %>%
  rename(
    ANOVA_p = p_adj
  )
tukey_res[tukey_res$ANOVA_p > p_thrsh,]$p.value = 1



# plots

sig_code = list()

for(i in unique(tukey_res$Tissue)) {
  
  print(i)
  
  for(j in unique(tukey_res$Metabolite)) {
    
    tmp = tukey_res[tukey_res$Tissue == i & tukey_res$Metabolite == j,]
    
    sig_code[[i]][[j]] = data.frame(
      Tissue = i,
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
  
  sig_code[[i]] = bind_rows(sig_code[[i]])
  
} # i

sig_code = bind_rows(sig_code)



plot_df = data %>%
  gather(
    key = "Metabolite",
    value = "value",
    Glc,Frc,Suc,Starch
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
    Tissue,Metabolite,Stage
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
    Metabolite = factor(Metabolite,levels=c("Glc","Frc","Suc","Starch")),
    Tissue = factor(Tissue,levels=unique(data$Tissue))
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
      mult = c(0.05,0.15),
      add = 0
    )
  ) +
  labs(
    y = expression(bold(~"["*mg~g^-1~DW*"]"))
  ) +
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
    panel.background = element_rect(fill = "transparent", color = "black"),
    #axis.line = element_line(color = 'black', size = 0.25),
    strip.text.y = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.text.x = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave(
  path = figure_path,
  filename = "NSC.png",
  units = "cm",
  width = 14,
  height = 14,
  dpi = 600
)
