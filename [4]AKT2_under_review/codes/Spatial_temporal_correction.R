
################################################################################
############    Spatial and temporal correction of field data   ################
################################################################################
# "A non-rectifying potassium channel increases cassava drought tolerance and storage root yield"
# by Mercedes Thieme


library(tidyverse)
library(data.table)
library(SpATS)    # spatial correction
library(metan)    # correlations
library(lmerTest) # linear mixed models (lme4) + p-values

# custom formula to calculate corrected values & identify outlier
spats_blue <- function(data,
                       genotype = "line", 
                       response = "rfw", 
                       row = "row", 
                       col = "col"){
  
  df <- as.data.frame(data)
  
  df$row_factor <- as.factor(df$row)
  df$col_factor <- as.factor(df$col)
  
  # prepare nested basis 
  # reduce computational cost by using nested bases with only the main effects (Lee et al. 2013) 
  # number of segments should be a divisor of original number of segments
  ratio = c(0.66,0.66) # use only 2/3rds of the original number of segments
  nx <- ceiling(nlevels(df$row_factor) * ratio[1]) - ceiling(nlevels(df$row_factor) * ratio[1])%%2 # original number of rows (64) vs reduced number (2/3) of rows (42)
  ny <- ceiling(nlevels(df$col_factor) * ratio[2]) - ceiling(nlevels(df$col_factor) * ratio[2])%%2
  
  
  m_blue <- SpATS(response = response,  # phenotype of interest
                  genotype = genotype, # variety/line name
                  spatial = ~ PSANOVA(row, col, nseg = c(nx, ny), nest.div = c(2, 2)), # nseg = number of knots for smoothing
                  fixed = ~ 1, 
                  random = ~ row_factor + col_factor, # or c(row_factor, col_factor), 
                  genotype.as.random = FALSE,
                  data = df) 
  
  
  # add BLUEs and SEs
  pr_blue <- predict(m_blue, which= "line")
  blue <- pr_blue %>% dplyr::select(line, BLUE = "predicted.values", BLUE_SE = "standard.errors")
  #m_blue$fitted
  
  df_out <- left_join(df, blue) %>% as.data.frame()
  #df_out <- left_join(df[[t]], blue) %>% as.data.frame()
  
  # add spatial components (residuals, spatial trend & spatially corrected trait)
  df_out$residual <- m_blue$residuals
  df_out$spatial_trend <- df_out[[response]] - df_out$residual - df_out$BLUE
  df_out[[paste0("trait_spatially_corrected")]] <- df_out$BLUE + df_out$residual
  
  # outlier identification
  # residuals that are ±3 standard deviations away form residual mean
  VarE <- m_blue$psi[1]
  u <- 3 * sqrt(VarE)
  
  temp_outlier <- filter(df_out, residual >= abs(u)) 
  
  #col_name_outlier <- paste0("outlier_", response) # create 1 outlier col per trait
  #df_out[[col_name_outlier]] <- ifelse(df_out$plot %in% temp_outlier$plot, yes = "y", no = "n")
  
  df_out$outlier <- ifelse(df_out$plot %in% temp_outlier$plot, yes = "y", no = "n")
  
  output <- tibble(m_blue = list(m_blue), results_table = list(df_out))
  
  return(output)
  
  
}


################################################################################
### SPATIAL CORRECTION 
################################################################################
# per year and trait

setwd("YOUR/PATH")
df <- fread("./input_AKT2.txt") 

toi <- c("sfw", "rfw", "tfw",  "hi", "rs_ratio", "trdm", "dmc") # list traits
years <- unique(df$year) # list years
res <- data.frame() # empty results table
save_location <- "./graphs"

# loop through traits and years
for (t in toi){
  for(y in years){
    print(paste("Processing year:", y, "trait:", t)) 
    
    df_sub <- df %>% filter(year == y)
    
    m <- spats_blue(data = df_sub, 
                    genotype = "line", 
                    response = t,  
                    row = "row", 
                    col = "col")
    
    # append results
    temp_res <- as.data.frame(m$results_table) # results per year
    temp_res$trait <- t

    res <- rbind(res, temp_res) # append results table (for all years)
    
    # OPTIONAL: save spatial distribution figures (per trait & year)
    # png(paste0(save_location, "/", y, "_", t, "_", "spatial.png"), width = 15, height = 10, units = "cm", res = 300)
    # plot(m$m_blue[[1]])
    # dev.off()
    
  }
}

fwrite(res, "results_spats.csv")


###
# plot outliers

out_count_year <- res %>% 
  filter(outlier == "y") %>% 
  group_by(year) %>% 
  summarise(n = n()) 

ggplot(out_count_year, aes(x = factor(year), y = n)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Year", y = "Count", title = "Excluded outliers per year") 


###
# plot BLUEs

res <- res %>% 
  filter(outlier == "n") %>% #exclude outliers
  mutate_at(c('line', 'year', 'accession', 'stock', 'trait', 'line_num'), as.factor) # convert to factors

traits <- c("rfw", "sfw", "tfw", "hi", "rs_ratio", "dmc", "trdm")
cbf_colors <- c("#D81B60", "#1E88E5", "#004D40")

for(t in traits){
  ggplot(filter(res, trait == t), aes(x = line, y = trait_spatially_corrected, color = factor(year))) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    #geom_jitter() +
    #geom_dotplot(binaxis = 'y', stackdir = 'center') +
    scale_color_manual(values = cbf_colors) +
    theme_minimal() +
    labs(x = "Line", y = t, color = "Year") +
    theme(axis.text=element_text(size = 14),
          axis.title=element_text(size = 16),
          legend.text=element_text(size = 14),
          legend.title=element_text(size = 16), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          axis.line = element_line(size = 1),
          legend.position = "bottom")
  ggsave(paste0("./graphs/boxplot_traits_years_", t, ".png"), width = 15, height = 15, unit = "cm", dpi = 300)
}
  
  
################################################################################
### CORRELATIONS 
################################################################################
# per year 
 
res_wide <- res %>% 
  select('line', 'rep', 'trait', 'year', 'trait_spatially_corrected' ) %>% 
  pivot_wider(names_from = trait, 
              values_from = c(trait_spatially_corrected))

for(y in years){
  res_wide_year <- res_wide %>% 
    filter(year %in% y) %>% 
    select(-c(rep, year))

  c <- corr_coef(res_wide_year) # metan package
  plot(c,
       reorder = F,
       col.low = "#2166AC",
       col.mid = "white",
       col.high = "#B2182B", 
       caption = F,  # signif. legend
       legend.position = "bottom") 
  ggsave(paste0("./graphs/", y, "_correlation.png"), width = 12, height = 12, units = "cm", dpi = 300)
}


################################################################################
### TEMPORAL CORRECTION 
################################################################################
# per trait
# reference has to be alphanumeric first item -> create vector control dummy variable

# dummy variable
# average all VCs, but keep replicates (average replicate 1 of all VC lines in 2022, average replicate 2 of all VC lines in 2022, ...)
dummy <- res %>% 
  filter(str_detect(line, "^VC")) %>%     # filter all VCs
  group_by(trait, year, rep) %>%          # mean of all rep 1 in year 22 for trait xyz
  summarise(trait_spatially_corrected = mean(trait_spatially_corrected)) %>% 
  mutate(line = "0000_VC_mean_rep", 
         line_num = "0000")          # add name

res_no_VC <- res %>% 
  filter(!str_detect(line, "^VC")) # exclude single VCs
  
res_dummy <- rbind(dummy, res_no_VC) %>%  # add VC_means to rest of data
  mutate_at(c('trait', 'year','rep', 'line', 'line_num'), as.factor) # make factor (important for model)


# run model
var_out <- data.frame()
res_out <- data.frame()
traits <- c("rfw", "sfw", "tfw", "hi", "rs_ratio", "dmc", "trdm")

for (t in traits){
  df <- res_dummy %>% 
    filter(trait == t)
  
  m_blue <- lmer(trait_spatially_corrected ~ line +(1|year) +(1|rep), data = df) # genotype as fixed to get BLUEs
  m_blup <- lmer(trait_spatially_corrected ~ 1 + (1|line) +(1|year) +(1|rep), data = df) # genotype as random (BLUPs) to get variances 
  
  sum_m_blue <- summary(m_blue)
  
  # Variances
  var_comp <- as.data.frame(VarCorr(m_blup))
  var_rep <- var_comp[1, "vcov"]  # Variance of V1 (rep)
  var_line <- var_comp[2, "vcov"]  # Variance of V2 (line)
  var_year <- var_comp[3, "vcov"]  # Variance of V3 (year)
  var_residual <- var_comp[4, "vcov"]  # Variance of V4 (residual)
  total_var <- var_rep + var_line + var_year + var_residual
  var_line_percent <- (var_line / total_var) * 100 # % reps explain
  var_rep_percent <- (var_rep / total_var) * 100 # % reps explain
  var_year_percent <- (var_year / total_var) * 100 # % year explains
  var_residual_percent <- (var_residual / total_var) * 100 # % reps explain
  
  var <- cbind(t, round(var_line_percent, 1), round(var_year_percent, 1), round(var_rep_percent, 1), round(var_residual_percent, 1))
    var_out <- rbind(var_out, var)
  
  # Plot BLUEs
  res_blue <- as.data.frame(sum_m_blue$coefficients) %>% 
    rownames_to_column("Effect") %>% 
    rename(Std.Error = "Std. Error", t_value = "t value", p_value =  "Pr(>|t|)" ) %>%   # rename cols
    mutate(value = Estimate + Estimate[1]) %>%                                          # sum estimates with ref estimate
    mutate(value = replace(value, value == Estimate[1]*2, Estimate[1])) %>%             # all are summed (also the ref itself) -> change back to original value
    mutate(line = str_remove(Effect, "line")) %>%                                       # remove prefix line
    mutate(line = replace(line, line == "(Intercept)", "0000")) %>%                     # rename Intercept to VC_mean_rep
    mutate(signif = ifelse(p_value < 0.001, "***", 
                           ifelse(p_value < 0.01, "**", 
                                  ifelse(p_value < 0.05, "*", 
                                         no = "")))) %>%
    mutate(trait = t)
  
  res_out <- rbind(res_out, res_blue)
  
  # define position of highlighted control sample
  min_vc <- filter(res_blue, line == "0000")$value - filter(res_blue, line == "0000")$Std.Error
  max_vc <- filter(res_blue, line == "0000")$value + filter(res_blue, line == "0000")$Std.Error
  
  # define position for significance labels
  max_y <- max(res_blue$value + res_blue$Std.Error)
  range_y <- max(res_blue$value + res_blue$Std.Error) - min(res_blue$value + res_blue$Std.Error)
  label_pos <- max_y + (0.01 * range_y)
  
  p <- ggplot(res_blue, aes(x = line, y = value)) + 
    geom_point(size = 1.5)+
    geom_errorbar(width = .4, linewidth = 0.9, aes(ymin = value - Std.Error, ymax = value + Std.Error)) +
    geom_hline(yintercept = filter(res_blue, line == "0000")$value, linetype = "dashed") + # mean of VC-mean as horizontal line
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = min_vc, ymax = max_vc), fill = "cadetblue", alpha = 0.01) + # color rectangle of VC_mean variance
    geom_text(aes(label = signif), y = label_pos, size = 9) + # add label at max value
    theme_minimal() +
    labs(x = "Line", y = t) +
    theme(axis.text=element_text(size = 14),
          axis.title=element_text(size = 16),
          legend.text=element_text(size = 14),
          legend.title=element_text(size = 16), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          axis.line = element_line(size = 1)) #, 
  
   ggsave(paste0("./graphs/plot_corrected_VC_mean_",t, ".png"), width = 15, height = 15, units = "cm", dpi = 300) 
 
}

colnames(var_out) <- c("trait", "var_geno", "var_year", "var_rep", "var_res")
fwrite(var_out, "./temporal_model_plot_corrected_variances.csv")

fwrite(res_out, "./temporal_model_results.csv")


















