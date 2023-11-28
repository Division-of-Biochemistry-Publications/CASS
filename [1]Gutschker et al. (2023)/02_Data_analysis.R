########################################################################################################################
# LIBRARIES
########################################################################################################################
library("dplyr")
library("tibble")
library("DESeq2")
library("clusterProfiler")
library("ggplot2")
library("tidyverse")
library("ppcor")
library("ggpubr")
library("cowplot")
source("http://www.sthda.com/upload/rquery_cormat.r")
library("factoextra")
library("EnhancedVolcano")
library("reshape2")
library("corrplot")
library("GenomicFeatures")
library("regioneR")
library("yarrr")
library("karyoploteR")
library("lme4")
library("car")
library("gridExtra")
library("patchwork")
library("ggcorrplot")

########################################################################################################################
# SET DIRECTORIES
########################################################################################################################
dir.create("./")
WD <- "./"
dir.create(paste0(WD, "XX_FIGURES/"))
figure_path <- paste0(WD, "XX_FIGURES/")
dir.create(paste0(WD, "01_RESULTS/"))
results_path <- paste0(WD, "01_RESULTS/")

########################################################################################################################
# FUNCTIONS AND VARIABLES
########################################################################################################################
# DIFFERENTIAL EXPRESSION ANALYSES
Differential_expression <- function(results_path, output_suffix, dds_dataset, contrast, 
                                    pAdjustMethod, padj_filter, L2FC_filter) {
  message(sprintf("%s comparison is being processed...", output_suffix))
  res <- results(dds_dataset, 
                 contrast = contrast, 
                 pAdjustMethod = pAdjustMethod)
  write.table(res, 
              sprintf("%s%s_DESeq2_results_unfiltered.txt", results_path, output_suffix), 
              row.names = TRUE, 
              sep='\t')
  res <- as.data.frame(res)
  res_up <- dplyr::filter(res, padj <= padj_filter & log2FoldChange >= L2FC_filter)
  write.table(res_up, 
              sprintf("%s%s_DESeq2_results_filtered_UPREGULATED.txt", results_path, output_suffix), 
              row.names = TRUE, 
              sep='\t')
  res_down <- dplyr::filter(res, padj <= padj_filter & log2FoldChange <= -L2FC_filter)
  write.table(res_down, 
              sprintf("%s%s_DESeq2_results_filtered_DOWNREGULATED.txt", results_path, output_suffix), 
              row.names = TRUE, 
              sep='\t')
}

# COUNT NORMALIZATION SEPARATElY FOR WHITE AND YELLOW
Count_cormalization_by_type <- function(data_type, meta_data, count_data, results_path) {
  meta_data_filtered <- meta_data %>% filter(type == data_type)
  count_data_filtered <- count_data[, names(count_data) %in% meta_data_filtered$Sample]
  dds <- DESeqDataSetFromMatrix(countData = count_data_filtered, 
                                colData = meta_data_filtered, 
                                design = ~ 1)
  dds <- DESeq(dds)
  normcounts <- counts(dds, normalized = TRUE)
  write.table(normcounts, 
              sprintf("%s%s_Normcounts.txt", results_path, data_type), 
              row.names = TRUE, 
              sep='\t')
  write.table(as.data.frame(assay(vst(dds))), 
              sprintf("%s%s_VSTcounts.txt", results_path, data_type), 
              row.names = TRUE, 
              sep='\t')
  normcounts_depth_filtered <- normcounts[rowMeans(normcounts) >= 10,]
  write.table(normcounts_depth_filtered, 
              sprintf("%s%s_Normcounts_depth_10.txt", results_path, data_type), 
              row.names = TRUE, 
              sep='\t')
}

# ENRICHMENT
Functional_enrichment <- function(input_path, annotation_file, TERM2GENE, TERM2NAME, results_path, figure_path,
                                  First_component, Second_component, comparison, method) {
  message(sprintf("%s_%s - %s is being processed...", First_component, Second_component, comparison))
  Input <- read.table(input_path, row.names = 1, sep='\t')
  Input$Geneid <- substr(row.names(Input), 1, nchar(row.names(Input))-5)
  Input_annot <- merge(Input, annotation_file, by = "Geneid")
  filename_annot <- sprintf("%s%s_vs_%s_%s_annotated.txt", results_path, First_component, Second_component, comparison)
  write.table(Input_annot, filename_annot, sep="\t", quote = FALSE, row.names = FALSE)
  Input_annot <- dplyr::select(Input_annot, -Symbol, -Full_name) %>% dplyr::distinct()
  Input_annot_enrichment <- enricher(gene = Input_annot$Geneid, 
                                     TERM2GENE = TERM2GENE,
                                     TERM2NAME = TERM2NAME, 
                                     pAdjustMethod = "bonferroni",
                                     pvalueCutoff = 1, 
                                     maxGSSize = 1000, 
                                     minGSSize = 0)
  x <- dotplot(Input_annot_enrichment, showCategory = 15)
  filename_enrich <- sprintf("%s%s_%s_%s_enrichment_%s.txt", 
                             results_path, First_component, Second_component, comparison, method)
  filename_figure <- sprintf("%s%s_%s_%s_enrichment_%s.png",
                             figure_path, First_component, Second_component, comparison, method)
  if (length(Input_annot_enrichment$ID) > 0) {
    write.table(Input_annot_enrichment, 
                file = filename_enrich, 
                row.names = FALSE, 
                sep='\t')
    ggsave(x, file = filename_figure, width = 8, height = 10, dpi = 300)
  } else {
    print(sprintf("NO ENRICHED TERM FOUND FOR: %s_%s - %s", First_component, Second_component, comparison))
    write.table(Input_annot_enrichment, file = filename_enrich, row.names = FALSE, sep='\t')
  }
}

# PARTIAL CORRELATION
Partial_correlation <- function(count_data_file_name, results_path, meta_data) {
  normcounts_depth_filtered <- read.csv(paste0(results_path, count_data_file_name), sep='\t')
  normcounts_filtered <- as.data.frame(t(normcounts_depth_filtered))
  corlist <- list()
  for(i in names(normcounts_filtered)) {
    df <- normcounts_filtered[i]
    df$Sample <- row.names(df)
    df <- drop_na(merge(df, meta_data, by="Sample"))
    res <- pcor.test(df[i], df$STARCH, df$YELLOWNESS, method="spearman")
    corlist[[i]] <- res
  }
  cor_data <- do.call(rbind, corlist)
  cor_data_filtered <- cor_data %>% filter(p.value <= 0.05)
  return(cor_data_filtered)
}


# FUNCTIONAL ENRICHMENT CORRELATION DATA
Functional_enrichment_cor_data <- function(color, cor_direction, method, annotation, TERM2GENE_KEGG, 
                                           TERM2NAME_KEGG, TERM2GENE_GO, TERM2NAME_BP) {
  data_name <- get_data_by_color(color)
  cor_data <- read.delim(paste0(results_path, data_name, ".txt"), header = TRUE)
  filtered_data <- filter_data_by_cor_direction(cor_data, cor_direction)
  message(sprintf("%s_%s is being processed...", cor_direction, data_name))
  filtered_data$Geneid <- substr(row.names(filtered_data), 1, nchar(row.names(filtered_data)) - 5)
  Input_annot <- merge(filtered_data, annotation, by = "Geneid")
  filename_annot <- sprintf("%s%s_%s_annotated.txt", results_path, cor_direction, data_name)
  write.table(Input_annot, 
              filename_annot, 
              sep="\t", 
              quote = FALSE, 
              row.names = FALSE)
  Input_annot <- dplyr::select(Input_annot, -Symbol, -Full_name) %>% dplyr::distinct()
  TERM2GENE <- if (method == "KEGG") TERM2GENE_KEGG else TERM2GENE_GO
  TERM2NAME <- if (method == "KEGG") TERM2NAME_KEGG else TERM2NAME_BP
  Input_annot_enrichment <- enricher(
    gene = Input_annot$Geneid,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 1,
    maxGSSize = 1000,
    minGSSize = 0
  )
  x <- dotplot(Input_annot_enrichment, showCategory = 15)
  filename_enrich <- sprintf("%s%s_%s_enrichment_%s.txt", results_path, cor_direction, data_name, method)
  filename_figure <- sprintf("%s%s_%s_enrichment_%s.png", figure_path, cor_direction, data_name, method)
  if (length(Input_annot_enrichment$ID) > 0) {
    write.table(Input_annot_enrichment, file = filename_enrich, row.names = FALSE, sep='\t')
    ggsave(x, file = filename_figure, width = 8, height = 10, dpi = 600)
  } else {
    print(sprintf("NO ENRICHED TERM FOUND FOR: %s_%s - %s ENRICHMENT", cor_direction, data_name, method))
    write.table(Input_annot_enrichment, file = filename_enrich, row.names = FALSE, sep='\t')
  }
}

# VARIABLES
padj_filter <- 0.01
set.seed(42)
########################################################################################################################
# LOAD AND PREPARE DATA
########################################################################################################################
metabolite_data <- read.delim(sprintf("%s00_DATA/PYT52_metabolites_2021.txt", WD), sep="\t", header = TRUE, row.names=1) %>%
  tibble::rownames_to_column("plot_number")
new_sugar_measures <- read.delim(sprintf("%s00_DATA/PYT52_SugarMeasurements.txt", WD), sep="\t", header = TRUE) %>%
  rename(plot_number = X) %>%
  group_by(plot_number) %>%
  summarise_all(mean)
final_metabolite_data <- merge(metabolite_data, new_sugar_measures, by = "plot_number") %>%
  column_to_rownames("plot_number")
metabolite_data_scaled <- final_metabolite_data %>%
  scale() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot_number")
meta_data <- read.delim(sprintf("%s00_DATA/PYT52_phenotype_data_2021.txt", WD), sep="\t", header=TRUE)
accession_naming <- read.delim(sprintf("%s00_DATA/CS_GT_Plot.txt", WD), sep="\t", header=TRUE)
meta_data_renamed <- merge(meta_data, accession_naming, by=c("accession_name", "plot_number"))
metabolite_data_scaled <- merge(meta_data_renamed, metabolite_data_scaled, by="plot_number") %>%
  mutate(type = ifelse(YELLOWNESS < 25, "WHITE_GT", "YELLOW_GT")) %>%
  dplyr::select(-c(STARCH, DM, YELLOWNESS, RFW, accession_name)) %>%
  column_to_rownames("Sample")
count_data <- read.delim(sprintf("%s00_DATA/PYT52_FEATURECOUNTS_gene_UNIQUE.txt", WD), row.names=1) %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  dplyr::select(matches("CS")) %>%
  dplyr::select(order(names(.)))
meta_data <- inner_join(meta_data, accession_naming, by=c("plot_number", "accession_name")) %>%
  mutate(type = case_when(
    YELLOWNESS < 25  ~ "WHITE_GT",
    YELLOWNESS >= 25 ~ "YELLOW_GT",
    TRUE             ~ NA_character_)) %>%
  filter(!is.na(type))
count_data <- count_data %>%
  dplyr::select(meta_data$Sample)
meta_data_WHITE <- meta_data %>% filter(type == "WHITE_GT")
meta_data_YELLOW <- meta_data %>% filter(type == "YELLOW_GT")

CASSAVA_KO_PATHWAY <- read.csv(paste0(WD, "00_DATA/Mesculenta_KEGG_PATHWAY.csv"), sep=",", stringsAsFactors = FALSE)
TERM2GENE_KEGG <- CASSAVA_KO_PATHWAY %>% dplyr::select(Pathway, locusName) %>% unique() %>% na.omit()
TERM2NAME_KEGG <- CASSAVA_KO_PATHWAY %>% dplyr::select(Pathway, Name) %>% unique() %>% na.omit()
annotation <- read.csv(paste0(WD, "00_DATA/Gene2AT_hits.txt"), sep = "\t")
names(annotation)[names(annotation) == "locusName"] <- "Geneid"
GO_annotation <- read.csv(paste0(WD, "00_DATA/Mesculenta_GO.csv"), sep=",") %>%
  rename(GO_ID = GO) %>%
  dplyr::select(locusName, AT_locusName, AT_Alias, GO_ID, Ontology, Term)
list_GO_data <- split(GO_annotation, GO_annotation$Ontology)
GO_BP <- list_GO_data$BP
GO_MF <- list_GO_data$MF
GO_CC <- list_GO_data$CC
TERM2GENE_GO <- GO_annotation %>% dplyr::select(GO_ID, locusName)
TERM2NAME_BP <- GO_BP %>% dplyr::select(GO_ID, Term)
TERM2NAME_MF <- GO_MF %>% dplyr::select(GO_ID, Term)
TERM2NAME_CC <- GO_CC %>% dplyr::select(GO_ID, Term)

regulations <- c("UPREGULATED", "DOWNREGULATED")
databases <- list(
  KEGG = list(TERM2GENE = TERM2GENE_KEGG, TERM2NAME = TERM2NAME_KEGG),
  GO = list(TERM2GENE = TERM2GENE_GO, TERM2NAME = TERM2NAME_BP)
)

########################################################################################################################
# DIFFERENTIAL EXPRESSION ANALYSES
########################################################################################################################
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = ~ type)
dds <- DESeq(dds)
Differential_expression(results_path, 
                        "WHITE_vs_YELLOW", 
                        dds, c("type", "WHITE_GT", "YELLOW_GT"),
                        "fdr", 
                        0.01, 
                        0)

########################################################################################################################
# COUNT NORMALIZATION SEPARATED BY WHITE AND YELLOW GTs AND DEPTH FILTER
########################################################################################################################
Count_cormalization_by_type("WHITE_GT", 
                            meta_data, 
                            count_data, 
                            results_path)
Count_cormalization_by_type("YELLOW_GT", 
                            meta_data, 
                            count_data, 
                            results_path)

########################################################################################################################
# DEPTH FILTER EXPRESSION CORRELATED TO DM WHITE and YELLOW SEPERATELY
########################################################################################################################
cor_data_white_filtered <- Spearman_correlation("WHITE_GT_Normcounts_depth_10.txt", results_path, meta_data_WHITE)
write.table(as.data.frame(cor_data_white_filtered), paste0(results_path, "WHITE_GT_cor_expression_STARCH_without_carotenoids_effect.txt"), sep="\t", quote = FALSE)
cor_data_yellow_filtered <- Spearman_correlation("YELLOW_GT_Normcounts_depth_10.txt", results_path, meta_data_YELLOW)
write.table(cor_data_yellow_filtered, paste0(results_path, "YELLOW_GT_cor_expression_STARCH_without_carotenoids_effect.txt"), sep="\t", quote = FALSE)

########################################################################################################################
# DEPTH FILTER EXPRESSION CORRELATED TO DM WITH CAROTENOIDS AS COVARIATE WHITE and YELLOW SEPERATELY
########################################################################################################################
cor_data_white_filtered <- Partial_correlation("WHITE_GT_Normcounts_depth_10.txt", 
                                               results_path, 
                                               meta_data_WHITE)
write.table(cor_data_white_filtered, 
            paste0(results_path, "WHITE_GT_cor_expression_STARCH.txt"), 
            sep="\t", 
            quote = FALSE)
cor_data_yellow_filtered <- Partial_correlation("YELLOW_GT_Normcounts_depth_10.txt", 
                                                results_path, 
                                                meta_data_YELLOW)
write.table(cor_data_yellow_filtered, 
            paste0(results_path, "YELLOW_GT_cor_expression_STARCH.txt"), 
            sep="\t", 
            quote = FALSE)

########################################################################################################################
# ENRICHMENT WHITE vs YELLOW
########################################################################################################################
# KEGG and GO
for (reg in regulations) {
  for (db_name in names(databases)) {
    db <- databases[[db_name]]
    file_path <- paste0(results_path, 
                        "WHITE_vs_YELLOW_DESeq2_results_filtered_", 
                        reg, 
                        ".txt")
    Functional_enrichment(file_path, 
                          annotation, 
                          db$TERM2GENE, 
                          db$TERM2NAME, 
                          results_path, 
                          figure_path,
                          "WHITE", 
                          "YELLOW", 
                          reg, 
                          db_name)
  }
}

########################################################################################################################
# ENRICHMENT CORRELATED TO DM WITH CAROTENOIDS AS COVARIATE
########################################################################################################################
get_data_by_color <- function(color) {
  if (color == "WHITE") {
    return("WHITE_GT_cor_expression_STARCH")
  } else {
    return("YELLOW_GT_cor_expression_STARCH")
  }
}
# Filter data based on correlation direction
filter_data_by_cor_direction <- function(data, direction) {
  if (direction == "POS") {
    return(data %>% filter("estimate.rho" >= 0.3))
  } else {
    return(data %>% filter("estimate.rho" <= -0.3))
  }
}
# Loop through all combinations of parameters
for (color in c("WHITE", "YELLOW")) {
  for (cor_direction in c("POS", "NEG")) {
    for (method in c("KEGG", "GO")) {
      Functional_enrichment_cor_data(color, 
                                     cor_direction, 
                                     method,
                                     annotation, 
                                     TERM2GENE_KEGG, 
                                     TERM2NAME_KEGG, 
                                     TERM2GENE_GO, 
                                     TERM2NAME_BP)
    }
  }
}
########################################################################################################################
# CELL_WALL_TO_STARCH_RATIO
########################################################################################################################
meta_data <- read.delim(paste0(WD, "00_DATA/PYT52_phenotype_data_2021.txt"))
meta_data$type[meta_data$YELLOWNESS < 25] <- "WHITE_GT"
meta_data$type[meta_data$YELLOWNESS >= 25] <- "YELLOW_GT"
meta_data <- meta_data[!is.na(meta_data$type),]
meta_data$starch <- (meta_data$FYLD/100)*meta_data$STARCH
meta_data$cell_wall <- meta_data$DYLD-meta_data$starch
meta_data$ratio <- meta_data$cell_wall/meta_data$DYLD
meta_data$ratio <- (meta_data$DYLD-((meta_data$FYLD/100)*meta_data$STARCH))/meta_data$DYLD
BLUE <- lmer(ratio ~ (1|accession_name) + plot_number, data = meta_data)
meta_data <- meta_data %>%
  na.omit() %>%
  mutate(prediction = predict(BLUE))
meta_data2 <- meta_data %>%
  group_by(accession_name, type) %>%
  summarise(mean_prediction = mean(prediction)) %>%
  ungroup()
x <- filter(meta_data2, type == "WHITE_GT")$mean_prediction
y <- filter(meta_data2, type == "YELLOW_GT")$mean_prediction
shapiro.test(x)
ggqqplot(x)
shapiro.test(y)
ggqqplot(y)
ks.test(x,y)
leveneTest(mean_prediction ~ type, data = meta_data2)
t.test(mean_prediction ~ type, data = meta_data2, var.equal = TRUE)

########################################################################################################################
# SAVE SESSION INFO
########################################################################################################################
writeLines(capture.output(sessionInfo()), paste0(WD, "R_sessionInfo.txt"))

