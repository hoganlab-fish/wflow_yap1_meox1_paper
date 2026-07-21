#!/usr/bin/Rscript
# regenerate fig 7F split by genotype EMBO LEC VEC score plots split by genotype
# ../../Saki_data/Revision_analysis/queue_plots/
#   embo_markers_degbarplot_LECs.pdf
#   embo_markers_degbarplot_VECs.pdf
library(patchwork)
library(qs2)
library(Seurat)
library(SeuratObject)
library(tidyverse)

infile_path <- "../paper/D10051_meox1_dataset_Level_03_annotated.qs2"
data <- qs_read(infile_path)

data$Genotype <- gsub("meox1_mutant", "meox1 mutant", data$Genotype)
data$Genotype <- gsub("wildtype", "WT", data$Genotype)
data$Genotype <- factor(
    data$Genotype, 
    levels = c("WT", "meox1 mutant")
)

all_degs <- "../../Saki_data/Revision_analysis/queue_plots/meox1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv"
all_degs <- read_csv(all_degs, show_col_types = F)

lec_vec_labels <- "../../Saki_data/Revision_analysis/data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforlabels_only.csv"
lec_vec_labels <- read_csv(lec_vec_labels)

# calculate the scores
embo_data <- "../../Saki_data/paper/EMBO_JOURNAL_Level_02_seuratObject_LM.RDS"
Level_02_seuratObject <- readRDS(embo_data)

genes_labels <- "data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforlabels_only.csv"
genes_scoring <- "data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforscoring_only.csv"
genes_everything <- "data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf.csv"

# get markers based on the analysis performed for https://link.springer.com/article/10.15252/embj.2022112590
Level_02_seuratObject <- SetIdent(Level_02_seuratObject, value = 'L2_Phenotype_Stage')

# ---- write a function so we can speed things up a bit ----
run_LEC_VEC_markers <- function(seurat = Level_02_seuratObject,
                                timepoint){
  lec_idents <- c(sprintf("LEC_Stage_%idpf", timepoint), sprintf("LEC_prox1a_low_Stage_%idpf", timepoint))
  vec_idents <- c(sprintf("VEC_Stage_%idpf", timepoint))

  tmp_markers <- FindMarkers(object = seurat,
                             ident.1 = lec_idents,
                             ident.2 = vec_idents,
                             assay = "RNA",
                             logfc.threshold = 0.25, #threshold from EMBO paper
                             min.pct = 0.1) #threshold from EMBO paper
  #add expression ratio
  tmp_markers <- tmp_markers %>%
    mutate(., ExpnRatio = ifelse(avg_log2FC < 0, -(1/(pct.1/pct.2)), pct.1/pct.2),
              ExpnRatio = case_when(is.infinite(ExpnRatio) & pct.1 == 0 ~ 10000*pct.2, #replace infinite values with something as they're actually useful
                                    is.infinite(ExpnRatio) & pct.2 == 0 ~ 10000*pct.1,
                                    TRUE ~ ExpnRatio),
              direction = ifelse(avg_log2FC < 0, "VEC", "LEC"),
              timepoint = paste(timepoint, "dpf")) %>%
    rownames_to_column("Gene")

  #filter based on thresholds in EMBO paper
  genes <- tmp_markers %>%
    filter(., abs(ExpnRatio) > 1.5) #both up and down

  return(genes)

}

# ---- run for 3,4,5 dpf ----
list_all_timepoints <- lapply(c(3,4,5), function(TP) run_LEC_VEC_markers(timepoint = TP)) %>%
  bind_rows()


# ---- Get genes for LEC/VEC across all timepoints - to give genes labels ----
# note that this will ignore fold changes etc, this is for marking the genes as LEC or VEC markers

lec_vec <- list_all_timepoints %>%
  group_by(Gene) %>%
  mutate(., always_same_direction = ifelse(all(avg_log2FC > 0)|all(avg_log2FC <0), "yes", "no")) %>% #make sure the directionality is always lec or always vec
  filter(., always_same_direction == "yes") %>%
  select(., Gene, direction) %>%
  distinct() %>%
  filter(., !(grepl("si:|CABZ|zgc:|im:" ,Gene))) #EMBO paper also cleans up more and removes unannotated genes


lec_vec$direction %>% table()
# LEC  VEC
# 2583  321

# ---- Get genes for LEC/VEC across all timepoints - for LEC/VEC scoring ----
lec_vec_scoring_genes <- list_all_timepoints %>%
  filter(., Gene %in% lec_vec$Gene & p_val_adj < 0.05) %>% # use genes from above + add significance threshold
  group_by(Gene) %>%
  reframe(., sum_log2fc = sum(avg_log2FC),
                direction = direction) %>%  #this is just for ranking them
  distinct() %>%
  group_by(direction) %>%
  top_n(., n = 25, wt = abs(sum_log2fc))


# ---- Export for plotting later ----
#only genes
write_csv(file = genes_labels, x = lec_vec, col_names = T)

#for scoring
write_csv(file = genes_scoring, x = lec_vec_scoring_genes, col_names = T)

#everything
write_csv(file = genes_everything, x = list_all_timepoints, col_names = T)


theme_for_deg_plot <- plot_sizing <- theme_bw() + 
    theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

wrap_elements(plot_list[[1]]) + wrap_elements(plot_list[[2]]) + plot_layout(guides = "collect", axes="collect", ncol=1)


# make the actual barplots
plot_list <- lapply(unique(all_degs$group), function(GROUP){
    # MAIN PLOT
    plot <- all_degs %>% 
    filter(., !(grepl("si:|CABZ|zgc:|im:" ,Gene))) %>%  #like EMBO, remove unannotated genes
    filter(., group == GROUP & p_val < 0.05) %>% # raw p-value to reduce density of plotting
    mutate(., plot_order_variable = ifelse(avg_log2FC < 0, "down", "up")) %>% 
    left_join(., lec_vec_labels, by="Gene") %>% 
    mutate(., direction = factor(ifelse(is.na(direction), "none", paste0(direction, " marker")), levels = c("LEC marker", "VEC marker", "none"))) %>% 
    arrange(., plot_order_variable, -avg_log2FC) %>%    # primary key, then secondary key
    mutate(., Gene = factor(Gene, levels = unique(Gene))) %>% 
    ggplot(., aes(x = Gene, y = avg_log2FC, fill= direction)) +
    geom_bar(stat = "identity") +
    theme_for_deg_plot+
    scale_fill_manual(values = colours)

    #PLIE CHARTS
    in_pie <- all_degs %>% 
    filter(., !(grepl("si:|CABZ|zgc:|im:" ,Gene))) %>%  #like EMBO, remove unannotated genes
    filter(., group == GROUP & p_val < 0.05) %>% # raw p-value to reduce density of plotting
    mutate(., plot_order_variable = ifelse(avg_log2FC < 0, "down-\nregulated", "up-\nregulated")) %>% 
    left_join(., lec_vec_labels, by="Gene") %>% 
    mutate(., direction = factor(ifelse(is.na(direction), "none", paste0(direction, " marker")), levels = c("LEC marker", "VEC marker", "none"))) %>% 
    arrange(., plot_order_variable, -avg_log2FC) %>%    # primary key, then secondary key
    mutate(., Gene = factor(Gene, levels = unique(Gene))) %>% 
    group_by(direction, plot_order_variable) %>% 
    summarise(., count = n())
    pie   <- ggplot(in_pie, aes(x="", 
                                y=count, fill=direction)) + geom_bar(width = 1, stat = "identity")+ 
    coord_polar("y", start=0)+
    scale_fill_manual(values = colours) + 
    facet_wrap(~plot_order_variable, scale = "free") +
    theme_void() +
    theme(strip.text = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.text = element_text(size=14))

    #COMBINE AND OUT
    plot_combined <- plot + pie + plot_annotation(title=paste0("DEG results for ", GROUP), 
                                                theme = theme(plot.title = element_text(size = 18, face = "bold"))) + plot_layout(guides = "collect", widths = c(2,1))
    return(plot_combined)
})

# wrap_elements(plot_list[[1]]) + wrap_elements(plot_list[[2]]) + plot_layout(guides = "collect", axes="collect", ncol=1)

ggsave(
    plot = plot_list[[1]],
    filename = "queue_plots/embo_scores/embo_markers_degbarplot_VECs.pdf",
    device = "pdf",
    height = 9,
    width = 16
    )
ggsave(
    plot = plot_list[[1]],
    filename = "queue_plots/embo_scores/embo_markers_degbarplot_VECs.png",
    device = "png",
    height = 9,
    width = 16
    )
ggsave(
    plot = plot_list[[2]],
    filename = "queue_plots/embo_scores/embo_markers_degbarplot_LECs.pdf",
    device = "pdf",
    height = 9,
    width = 16
    )
ggsave(
    plot = plot_list[[2]],
    filename = "queue_plots/embo_scores/embo_markers_degbarplot_LECs.png",
    device = "png",
    height = 9,
    width = 16
    )


