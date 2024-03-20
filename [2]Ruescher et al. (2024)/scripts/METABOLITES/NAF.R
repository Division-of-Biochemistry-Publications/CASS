#load packages ####
library(tidyverse)
library(emmeans)
#load data ####
data = read.csv("./METABOLITES/04_NAF/01_RAW_DATA/NAF.csv") %>%
  mutate(Relative = Relative * 100)
#statistics ####
figure_path = "./00_FIGURES/METABOLITES/04_NAF/"
dir.create(paste0(figure_path,"00_DIAGNOSTICS"),recursive=TRUE)
tukey = list()
sigCode = list()
for(i in unique(data$Metabolite)) {
  input=data[data$Metabolite == i,]
  model = lm(Relative ~ Compartment,input)
  print(summary(model))
  png(paste0(figure_path,"00_DIAGNOSTICS/",i,".png"),units = "cm",width=16,height=16,
      res = 600)
  par(mfrow=c(2,2))
  plot(model)
  dev.off()
  tmp = data.frame(Metabolite = i,emmeans(model,pairwise ~ Compartment)$contrast)
  tukey[[i]] = tmp
  sigCode[[i]] = data.frame(Metabolite= i, rcompanion::cldList(p.value ~ contrast, tmp),threshold = 0.05) %>% 
    rename("Compartment" = "Group")
}
tukey = bind_rows(tukey)
sigCode = bind_rows(sigCode)
#plot####
sigDist=max(data$Relative)*0.1
plotData = data %>% group_by(Compartment) %>% slice(which.max(Relative)) %>% 
  mutate(sigHeight=Relative+sigDist) %>% select(Compartment,sigHeight) %>% inner_join(sigCode)
ggplot(data,aes(x=Compartment,y=Relative,fill=Compartment)) +
  stat_boxplot(geom = "errorbar", width = 0.25, size = 0.2) +
  geom_boxplot(width = 0.75, size = .25, outlier.shape = NA) +
  geom_jitter(shape = 21, fill = "white", color = "black", size = 1, width = 0.25) +
  geom_text(data=plotData,aes(y = sigHeight, label = Letter), size = 2.5) +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.15), add = 0)) +
  labs(x="Compartment",y = expression(bold(~"[%]"))) +
  geom_hline(yintercept = 0) +
  facet_grid(rows = "Metabolite") +
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
ggsave("08_NAF.tiff", path="../../PUBLICATIONS/Transportome_2023/Draft3/FIGURES/",
       units = "cm", width = 8, height = 16, dpi = 600)
ggsave("08_NAF.png", path="../../PUBLICATIONS/Transportome_2023/Draft3/FIGURES/JPEGS/",
       units = "cm", width = 8, height = 16, dpi = 600)  
