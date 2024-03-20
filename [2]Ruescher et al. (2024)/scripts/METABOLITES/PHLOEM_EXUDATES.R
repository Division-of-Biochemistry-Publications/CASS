# set WD ####
setwd("/Users/david/iCloud/PhD/RESULTS/Transport/")
# load packages ####
library(tidyverse)
library(emmeans)
# load data ####
data <- read.csv("./METABOLITES/03_PHLOEM_EXUDATES/PHLOEM_EXUDATES.csv") %>%
  mutate(Plant=factor(Plant),Leaf=factor(Leaf)) %>% gather("Metabolite","value",-Plant,-Leaf) %>%
  group_by(Plant,Metabolite) %>% summarise_at("value",mean,na.rm=TRUE) %>%
  mutate(Metabolite=factor(Metabolite,levels=c("Sucrose","Glucose","Fructose","Raffinose"))) %>%
  ungroup()
# statistics ####
model=glm(value~Metabolite,Gamma(link="log"),data)
par(mfrow=c(2,2))
plot(model)
tukey=as.data.frame(emmeans(model,pairwise~Metabolite)$contrasts)
sigCode=rcompanion::cldList(p.value = tukey$p.value, comparison = tukey$contrast,threshold = 0.001) %>%
  rename(Metabolite=Group)
# plot ####
sigDist=max(data$value)*0.1
plotData = data %>% group_by(Metabolite) %>% slice(which.max(value)) %>% 
  mutate(sigHeight=value+sigDist) %>% select(Metabolite,sigHeight) %>% inner_join(sigCode)
ggplot(data,aes(x=Metabolite,y=value,fill=Metabolite)) +
  stat_boxplot(geom = "errorbar", width = 0.25, size = 0.2) +
  geom_boxplot(width = 0.75, size = .25, outlier.shape = NA) +
  geom_jitter(shape = 21, fill = "white", color = "black", size = 1, width = 0.25) +
  geom_text(data=plotData,aes(y = sigHeight, label = Letter), size = 2.5) +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.15), add = 0)) +
  labs(x="Metabolite",y = expression(bold(~"["*mg~L^-1*"]"))) +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 10, face = "bold", hjust = .5, vjust = 0.5),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text.y = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.text.x = element_text(size = 8, face = "bold", hjust = .5, vjust = 0.5),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave("S2_PHLOEM_EXUDATES.tiff", path="../../PUBLICATIONS/Transportome_2023/SUPPLEMENT/FIGURES/",
       units = "cm", width = 8, height = 8, dpi = 600)
ggsave("S2_PHLOEM_EXUDATES.png", path="../../PUBLICATIONS/Transportome_2023/SUPPLEMENT/JPEGS/",
       units = "cm", width = 8, height = 8, dpi = 600)
