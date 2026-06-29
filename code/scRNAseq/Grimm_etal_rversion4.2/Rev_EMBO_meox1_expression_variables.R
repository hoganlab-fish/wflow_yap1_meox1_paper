#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited May 2023

################################################################################
# LOAD FUNCTIONS
################################################################################
# ---- Load 00_seurat_function ----
# source seurat functions
source("/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_scripts/00_seurat_functions_R4.2.0Core.R")

library(scales)		# for rescaling color in bar plots
library(tidyr)
library(dplyr)
################################################################################
# IN-SCRIPT FUNCTIONS
################################################################################
plot_umap <- function(seurat.object,
                      group,
                      save.dir,
                      cols = NULL,
                      project.name,
                      size_h = 6,
                      size_w = 7){
  #with legend
  umap <- DimPlot(seurat.object,
                  group.by = group,
                  cols = cols,
                  shuffle = T,
                  label = F)
  umapFileName <- sprintf('%s/%s_umap_%s.pdf', save.dir, project.name, group)
  ggsave(umapFileName,
         umap,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  #without legend
  umap.nolegend.label <- DimPlot(seurat.object,
                                 group.by = group,
                                 cols = cols,
                                 shuffle = T,
                                 label = F) + NoLegend() + NoAxes() + ggtitle('')
  #without legend
  umap.nolegend <- DimPlot(seurat.object,
                           group.by = group,
                           cols = cols,
                           shuffle = T,
                           label = F) + NoLegend() + NoAxes() + ggtitle('')
  umapFileName <- sprintf('%s/%s_umap_%s_nolegend.pdf', save.dir, project.name, group)
  ggsave(umapFileName,
         umap.nolegend.label,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  return(umap)
}

make_barplot <- function(DEGList,
                         logfc_t = 0.5,
                         padj_t = 0.05){

  cts <- DEGList$celltype_cluster %>% unique()
  up_ns <- down_ns <- matrix(nrow = length(cts), ncol = 1, dimnames = list( cts %>% sort(),"count")) %>% as.data.frame()

  up_ns_t <- DEGList %>%
    filter(., AvLogFC > logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  down_ns_t <- DEGList %>%
    filter(., AvLogFC < -logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  up_ns[up_ns_t$celltype_cluster, 1] <- up_ns_t$n
  down_ns[down_ns_t$celltype_cluster, 1] <- down_ns_t$n
  all_ns <- cbind(up_ns, down_ns)
  colnames(all_ns) <- c('up', 'down')
  all_ns <- all_ns %>%
    mutate(., celltype_cluster = rownames(all_ns)) %>%
    gather('direction', 'count', -celltype_cluster)
  all_ns[is.na(all_ns)] <- 0
  #workaround to get order
  all_ns_tmp <- all_ns %>%
    group_by(., celltype_cluster) %>%
    summarise(across(count, sum)) %>%
    arrange(. , -count)
  all_ns <- all_ns %>% mutate(celltype_cluster=factor(celltype_cluster, levels=all_ns_tmp$celltype_cluster))

  #make plot
  barplot <- ggplot(all_ns, aes(fill=direction, y=count, x=celltype_cluster)) +
    geom_bar(position="stack", stat="identity")+
    xlab('') + ylab('') +
    ggtitle('Stacked barplot') + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +theme(axis.text=element_text(size=18),
                                                                                                   plot.title = element_text(size = 20))
  return(barplot)

}

make_barplot_flipped <- function(DEGList,
                         logfc_t = 0.5,
                         padj_t = 0.05){

  cts <- DEGList$celltype_cluster %>% unique()
  up_ns <- down_ns <- matrix(nrow = length(cts), ncol = 1, dimnames = list( cts %>% sort(),"count")) %>% as.data.frame()

  up_ns_t <- DEGList %>%
    filter(., AvLogFC > logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  down_ns_t <- DEGList %>%
    filter(., AvLogFC < -logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  up_ns[up_ns_t$celltype_cluster, 1] <- up_ns_t$n
  down_ns[down_ns_t$celltype_cluster, 1] <- down_ns_t$n
  all_ns <- cbind(up_ns, down_ns)
  colnames(all_ns) <- c('up', 'down')
  all_ns <- all_ns %>%
    mutate(., celltype_cluster = rownames(all_ns)) %>%
    gather('direction', 'count', -celltype_cluster)
  all_ns[is.na(all_ns)] <- 0
  #workaround to get order
  all_ns_tmp <- all_ns %>%
    group_by(., celltype_cluster) %>%
    summarise(across(count, sum)) %>%
    arrange(. , -count)
  all_ns <- all_ns %>% mutate(celltype_cluster=factor(celltype_cluster, levels=all_ns_tmp$celltype_cluster))

  #make plot
  barplot <- ggplot(all_ns, aes(fill=direction, y=celltype_cluster, x=count)) +
    geom_bar(position="stack", stat="identity")+
    xlab('') + ylab('') +
    ggtitle('Stacked barplot') + theme_bw() + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) +theme(axis.text=element_text(size=18),
                                                                                                   plot.title = element_text(size = 20))
  return(barplot)

}

make_barplot_flipped_splitupdown <- function(DEGList,
                         logfc_t = 0.5,
                         padj_t = 0.05,
                         order){

  cts <- DEGList$celltype_cluster %>% unique()
  up_ns <- down_ns <- matrix(nrow = length(cts), ncol = 1, dimnames = list( cts %>% sort(),"count")) %>% as.data.frame()

  up_ns_t <- DEGList %>%
    filter(., AvLogFC > logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  down_ns_t <- DEGList %>%
    filter(., AvLogFC < -logfc_t & AdjPVal < padj_t) %>%
    group_by(., celltype_cluster) %>%
    count(celltype_cluster)
  up_ns[up_ns_t$celltype_cluster, 1] <- up_ns_t$n
  down_ns[down_ns_t$celltype_cluster, 1] <- down_ns_t$n
  all_ns <- cbind(up_ns, down_ns)
  colnames(all_ns) <- c('up', 'down')
  all_ns <- all_ns %>%
    mutate(., celltype_cluster = rownames(all_ns)) %>%
    gather('direction', 'count', -celltype_cluster)
  all_ns[is.na(all_ns)] <- 0

  #order based on input
  if (order == "both"){
  all_ns_tmp <- all_ns %>%
    group_by(., celltype_cluster) %>%
    summarise(across(count, sum)) %>%
    arrange(. , -count)
  all_ns <- all_ns %>% mutate(celltype_cluster=factor(celltype_cluster, levels=all_ns_tmp$celltype_cluster))
  } else if (order == "up"){
    up_ns_t <- up_ns_t %>%
      arrange(. , -n)
    all_ns <- all_ns %>% mutate(celltype_cluster=factor(celltype_cluster, levels=up_ns_t$celltype_cluster))
  } else if (order == "down"){
    down_ns_t <- down_ns_t %>%
      arrange(. , -n)
    all_ns <- all_ns %>% mutate(celltype_cluster=factor(celltype_cluster, levels=down_ns_t$celltype_cluster))
  }


  #make plot
  barplot <- ggplot(all_ns, aes(fill = direction, y=celltype_cluster, x=count)) +
    geom_bar(position="stack", stat="identity")+
    xlab('') + ylab('') +
    ggtitle(sprintf('Number of DEG: |log2FC| > %0.2f', logfc_t)) + theme_bw() + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) +theme(axis.text=element_text(size=18),
                                                                                                   plot.title = element_text(size = 20)) + facet_grid(~direction)
  return(barplot)


}

plot_umap_feature <- function(seurat.object,
                              feature,
                              save.dir,
                              cols = NULL,
                              project.name,
                              size_h = 6,
                              size_w = 7){
  #with legend
  umap <- FeaturePlot(seurat.object,
                      features = feature,
                      cols = cols)
  umapFileName <- sprintf('%s/%s_umap_%s.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  #without legend
  umap.nolegend <- FeaturePlot(seurat.object,
                               features = feature,
                               cols = cols) + NoLegend() + NoAxes() + ggtitle('')
  umapFileName <- sprintf('%s/%s_umap_%s_nolegend.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap.nolegend,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  return(umap)
}
plot_umap_enrichment <- function(seurat.object,
                                 feature,
                                 save.dir,
                                 cols = NULL,
                                 project.name,
                                 size_h = 6,
                                 size_w = 7){
  #with legend
  umap <- FeaturePlot(seurat.object,
                      features = feature,
                      cols = cols,
                      min.cutoff = "q1",
                      max.cutoff =  "q99", order = T)
  umapFileName <- sprintf('%s/%s_umap_%s.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  #without legend
  umap.nolegend <- FeaturePlot(seurat.object,
                               features = feature,
                               min.cutoff = "q1",
                               max.cutoff =  "q99",
                               cols = cols, order = T) + NoLegend() + NoAxes() + ggtitle('')
  umapFileName <- sprintf('%s/%s_umap_%s_nolegend.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap.nolegend,
         device="pdf",
         width=size_w,
         height=size_h,
         units="cm",
         scale=1)

  return(umap)
}

plot_umap_enrichment_split_genotype <- function(seurat.object,
                                                feature,
                                                save.dir,
                                                cols = NULL,
                                                project.name,
                                                size_h = 12,
                                                size_w =15,
                                                point_size = 1){
  n <- seurat.object$Sample_Name %>% unique() %>% length()
  #with legend
  umap <- FeaturePlot(seurat.object,
                      features = feature,
                      cols = cols,
                      min.cutoff = "q1",
                      max.cutoff =  "q99",
                      split.by = "Sample_Name",
                      order = T,
                      pt.size = point_size)
  umapFileName <- sprintf('%s/%s_umap_%s_splitgenotype.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap,
         device="pdf",
         width=size_w * n,
         height=size_h,
         units="cm",
         scale=1)

  #without legend
  umap.nolegend <- FeaturePlot(seurat.object,
                               features = feature,
                               cols = cols,
                               min.cutoff = "q1",
                               max.cutoff =  "q99",
                               split.by = "Sample_Name",
                               order = T,
                               pt.size = point_size) + NoLegend() + NoAxes() + ggtitle('')
  umapFileName <- sprintf('%s/%s_umap_%s_splitgenotype_nolegend.pdf', save.dir, project.name, feature)
  ggsave(umapFileName,
         umap.nolegend,
         device="pdf",
         width=size_w * n,
         height=size_h,
         units="cm",
         scale=1)

  return(umap)
}
################################################################################
# SET VARIABLES
################################################################################
# ---- General  ----
# Set save directory
save.dir <- '/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/Figures_for_Revision/'
project.name.level.03_embo <- 'EMBO_Level_03'

# ---- Data  ----
Level_03_seuratObject <- readRDS("/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/R_Data_files/Level_03_seuratObject_LEC_VEC.RDS")
EMBO_Level_03_seuratObject <- readRDS('/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/R_Data_files/EMBO_JOURNAL_Level_03_seuratObject_LM.RDS')

# ---- Variables ----

embo_phenotype_stage_order <- c(
  "LEC_Stage_40hpf",
  'LEC_Stage_3dpf',
  'LEC_Stage_4dpf',
  'LEC_Stage_5dpf',
  'LEC_prox1a_low_Stage_40hpf',
  'LEC_prox1a_low_Stage_3dpf',
  'LEC_prox1a_low_Stage_4dpf',
  'LEC_prox1a_low_Stage_5dpf',
  'muLEC_Stage_40hpf',
  'muLEC_Stage_3dpf',
  'muLEC_Stage_4dpf',
  'muLEC_Stage_5dpf',
  'VEC_Stage_40hpf',
  'VEC_Stage_3dpf',
  'VEC_Stage_4dpf',
  'VEC_Stage_5dpf',
  'VEC_preLEC_Stage_40hpf',
  'VEC_preLEC_Stage_3dpf',
  'VEC_preLEC_Stage_4dpf',
  'VEC_preLEC_Stage_5dpf'
)

################################################################################
# PLOTTING SETTINGS
################################################################################
# ---- Plot settings: colours  ----
genotype_colours <- c(
  "#dbe2c6", #wildtype
  "#657c95" #mutant
)

names(genotype_colours) <- c("WT", "yap1_mutant")

expression_colours <- c("#d9d9d9", "#40004b")

cols_cellcycle <- c( '#ffff66', '#cc85ff','#a0d0e0')
names(cols_cellcycle) <- c("G1/G0", "S", "G2M")

