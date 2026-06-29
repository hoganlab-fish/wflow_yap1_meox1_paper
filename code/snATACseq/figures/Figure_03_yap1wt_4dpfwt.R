#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited November 2023

#   Depends on R/4.2.0.Core + Renv

################################################################################
# SET VARIABLES
################################################################################
#Load variables
source('./Figure_02_yap1wt_4dpfwt_variables.R')

################################################################################
# SUPP  YAP LEVEL_01: DOTPLOT PHENOTYPE
################################################################################
plot <- make_dotplot(
  project=project.yap1.level01,
  groupBy=c("ClustersID"),
  useGroups=ClustersID_ordered_Level_01_dotplot,
  groupName=cClustersID_ordered_Level_01_dotplot,
  markerGenes=dotplot_marker_genes_L1,
  exportFileName="Yap1_Level_01_ClustersID_markergenes",
  exportFilePath = figure.dir.out,
  accessibilityColour = c("#d9d9d9", "#7a0177"),
  plotWidth = 10,
  plotHeight = 10,
  rotateaxis = TRUE)

################################################################################
# SUPP  YAP LEVEL_01: UMAP CLUSTERID PHENOTYPE
################################################################################
umap_L1                    <- project.yap1.level01@embeddings$UMAP$df
umap_L1$ClustersID         <- project.yap1.level01$ClustersID
colnames(umap_L1)          <- c("dim1", "dim2","ClustersID")
cluster_colours     <- ggplotColours(length(unique(umap_L1$ClustersID)))

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersID)) +
  geom_point(size=2) + scale_color_manual(values=cluster_colours) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype")

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_ClusterID_phenotype.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersID)) +
  geom_point(size=2) + scale_color_manual(values=cluster_colours) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype") + NoLegend()

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_ClusterID_phenotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)

################################################################################
# SUPP  YAP LEVEL_01:  UMAP PHENOTYPE
################################################################################
umap_L1                    <- project.yap1.level01@embeddings$UMAP$df
umap_L1$ClustersPhenotype         <- project.yap1.level01$ClustersPhenotype
colnames(umap_L1)          <- c("dim1", "dim2","ClustersPhenotype")
phenotype_colours   <- ggplotColours(length(unique(umap_L1$ClustersPhenotype)))

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersPhenotype)) +
  geom_point(size=2) + scale_color_manual(values=phenotype_colours) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype")

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_phenotype.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersPhenotype)) +
  geom_point(size=2) + scale_color_manual(values=phenotype_colours) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype") + NoLegend()

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_phenotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)

################################################################################
# SUPP  YAP LEVEL_01:  UMAP GENOTYPE
################################################################################
umap_L1                   <- project.yap1.level01@embeddings$UMAP$df
umap_L1$Genotype         <- project.yap1.level01$Genotype
colnames(umap_L1)          <- c("dim1", "dim2","Genotype")
umap_L1$Genotype         <- factor(umap_L1$Genotype, levels=genotype_order)

#shuffle to avoid plotting everything on top of each other
set.seed(seed = 1997)
umap_L1 <- umap_L1[sample(x = 1:nrow(x = umap_L1)), ]


umap_genotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=Genotype)) +
  geom_point(size=2) + scale_color_manual(values=genotype_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Genotype")

umap_genotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_genotype.EPS', figure.dir.out)
ggsave(plot = umap_genotype_L1,
       filename =  umap_genotype_L1.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

umap_genotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=Genotype)) +
  geom_point(size=2) + scale_color_manual(values=genotype_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Genotype") + NoLegend()

umap_genotype_L1.export <- sprintf('%s/Yap1_Level_01_UMAP_genotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_genotype_L1,
       filename =  umap_genotype_L1.export,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)


################################################################################
# SUPP  YAP LEVEL_01:  BARPLOT GENOTYPE (PHENOTYPE)
################################################################################
# ---- Genotype || ClustersPhenotype  ----
meta_L1                    <- project.yap1.level01$ClustersPhenotype %>% as.data.frame()
meta_L1$Genotype         <- project.yap1.level01$Genotype
meta_L1$Genotype         <- factor(meta_L1$Genotype, levels=genotype_order)

colnames(meta_L1)          <- c("ClustersPhenotype", "Genotype")

summary_genotype <- meta_L1 %>%
  dplyr::group_by(., ClustersPhenotype) %>%
  dplyr::count(., Genotype)

barplot.genotype <- ggplot(summary_genotype, aes(fill=Genotype, y=n, x=ClustersPhenotype)) +
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values= genotype_cols) +
  ggtitle('Stacked barplot: Genotype') + theme_bw()

barplot_genotype_L1.export <- sprintf('%s/Yap1_Level_01_barplot_phenotype_genotype.EPS', figure.dir.out)
ggsave(filename = barplot_genotype_L1.export,
       plot = barplot.genotype,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

################################################################################
# YAP LEVEL_02:  UMAP PHENOTYPE
################################################################################
umap_L1                    <- project.yap1.level02@embeddings$UMAP$df
umap_L1$ClustersPhenotype         <- project.yap1.level02$Phenotype_Level_02
colnames(umap_L1)          <- c("dim1", "dim2","ClustersPhenotype")

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersPhenotype)) +
  geom_point(size=2) + scale_color_manual(values=ClustersPhenotype_ordered_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype")

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_phenotype.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersPhenotype)) +
  geom_point(size=2) + scale_color_manual(values=ClustersPhenotype_ordered_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype") + NoLegend()

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_phenotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)

################################################################################
# YAP LEVEL_02:  UMAP GENOTYPE
################################################################################
umap_L1                   <- project.yap1.level02@embeddings$UMAP$df
umap_L1$Genotype         <- project.yap1.level02$Genotype
colnames(umap_L1)          <- c("dim1", "dim2","Genotype")
umap_L1$Genotype         <- factor(umap_L1$Genotype, levels=genotype_order)

#shuffle to avoid plotting everything on top of each other
set.seed(seed = 1997)
umap_L1 <- umap_L1[sample(x = 1:nrow(x = umap_L1)), ]


umap_genotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=Genotype)) +
  geom_point(size=2) + scale_color_manual(values=genotype_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Genotype")

umap_genotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_genotype.EPS', figure.dir.out)
ggsave(plot = umap_genotype_L1,
       filename =  umap_genotype_L1.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)

umap_genotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=Genotype)) +
  geom_point(size=2) + scale_color_manual(values=genotype_cols) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Genotype") + NoLegend()

umap_genotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_genotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_genotype_L1,
       filename =  umap_genotype_L1.export,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)


################################################################################
# YAP LEVEL_02: UMAP MARKERS
################################################################################
plot_marker_genes_for_publication(
  ArchRProject=project.yap1.level02,
  MatrixToPlot="GeneScoreMatrix",
  geneList=marker_genes,
  pointSize=1,
  colours=c("#d9d9d9", "#7a0177"),
  legendPosition=c("top"),
  fileName="Yap1_Level_02_UMAP_features",
  filePath=figure.dir.out)

################################################################################
# YAP LEVEL_02: HIPPO ENRICHMENT
################################################################################
project.yap1.level02 <- ModuleScore(project = project.yap1.level02,
                       features = hippo_targets_nomeox,
                       name = 'Hippo',
                       seed = 1997,
                       reducedDims = 'InterativeLSI_Level_02' # to get imputed weights, if not available
)

p1 <- plotEmbedding(
  project.yap1.level02,
  colorBy = "cellColData",
  name = 'Hippo'
)

x <- ggplot_build(p1)$data[[1]]
enrichent.plot <- ggplot(x, aes(x=x, y=y)) + geom_point(aes(color=value), size=1) +
  theme_classic() + theme(legend.position='top') +
  scale_colour_gradientn(colours=c("#d9d9d9", "#7a0177"))

enrichment.export <- sprintf('%s/Yap1_Level_02_hippo_enrichment.EPS', figure.dir.out)
ggsave(plot = enrichent.plot,
       filename =  enrichment.export,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)


################################################################################
# SUPP  YAP LEVEL_02: DOTPLOT PHENOTYPE/CLUSTERID
################################################################################
plot <- make_dotplot(
  project=project.yap1.level02,
  groupBy=c("Phenotype_Level_02"),
  useGroups=c(unique(project.yap1.level02$Phenotype_Level_02) %>% sort()),
  groupName=c(unique(project.yap1.level02$Phenotype_Level_02) %>% sort()),
  markerGenes=dotplot_marker_genes_L2,
  exportFileName="Yap1_Level_02_Phenotype_markergenes",
  exportFilePath = figure.dir.out,
  accessibilityColour = c("#d9d9d9", "#7a0177"),
  plotWidth = 10,
  plotHeight = 10,
  rotateaxis = TRUE)

################################################################################
# SUPP  YAP LEVEL_02: DOTPLOT PHENOTYPE
################################################################################
plot <- make_dotplot(
  project=project.yap1.level02,
  groupBy=c("ClustersID_Level_02"),
  useGroups=ClustersID_ordered_Level_02_dotplot,
  groupName=ClustersID_ordered_Level_02_dotplot,
  markerGenes=dotplot_marker_genes_L2,
  exportFileName="Yap1_Level_02_ClusterID_markergenes",
  exportFilePath = figure.dir.out,
  accessibilityColour = c("#d9d9d9", "#7a0177"),
  plotWidth = 10,
  plotHeight = 9,
  rotateaxis = TRUE)

################################################################################
# SUPP  YAP LEVEL_02: UMAP CLUSTERID PHENOTYPE
################################################################################
umap_L1                    <- project.yap1.level02@embeddings$UMAP$df
umap_L1$ClustersID_Level_02         <- project.yap1.level02$ClustersID_Level_02
colnames(umap_L1)          <- c("dim1", "dim2","ClustersID_Level_02")
umap_L1$ClustersID_Level_02         <- factor(umap_L1$ClustersID_Level_02, levels=ClustersID_ordered_Level_02)

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersID_Level_02)) +
  geom_point(size=2) + scale_color_manual(values=ClustersID_ordered_cols_Level_02) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype")

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_ClusterID_phenotype.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=25,
       height=10,
       units="cm",
       scale=1)

umap_phenotype_L1 <- ggplot(umap_L1, aes(x=dim1, y=dim2, color=ClustersID_Level_02)) +
  geom_point(size=2) + scale_color_manual(values=ClustersID_ordered_cols_Level_02) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggtitle(label="Cell Phenotype") + NoLegend()

umap_phenotype_L1.export <- sprintf('%s/Yap1_Level_02_UMAP_ClusterID_phenotype_nolegend.EPS', figure.dir.out)
ggsave(plot = umap_phenotype_L1,
       filename =  umap_phenotype_L1.export,
       device="eps",
       width=25,
       height=15,
       units="cm",
       scale=1)


################################################################################
# SUPP  YAP LEVEL_02:  BARPLOT GENOTYPE (PHENOTYPE)
################################################################################
# ---- Genotype || ClustersID_Level_02  ----
meta_L1                    <- project.yap1.level02$ClustersID_Level_02 %>% as.data.frame()
meta_L1$Genotype         <- project.yap1.level02$Genotype
meta_L1$Genotype         <- factor(meta_L1$Genotype, levels=genotype_order)

colnames(meta_L1)          <- c("ClustersID_Level_02", "Genotype")

summary_genotype <- meta_L1 %>%
  dplyr::group_by(., ClustersID_Level_02) %>%
  dplyr::count(., Genotype)

barplot.genotype <- ggplot(summary_genotype, aes(fill=Genotype, y=n, x=ClustersID_Level_02)) +
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values= genotype_cols) +
  ggtitle('Stacked barplot: Genotype') + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

barplot_genotype_L1.export <- sprintf('%s/Yap1_Level_02_barplot_clusterid_phenotype_genotype.EPS', figure.dir.out)
ggsave(filename = barplot_genotype_L1.export,
       plot = barplot.genotype,
       device="eps",
       width=15,
       height=10,
       units="cm",
       scale=1)


################################################################################
# SUPP  YAP LEVEL_02: TRACKPLOT
################################################################################
# ---- variables ----
project <- project.yap1.level02
groupBy <- 'Phenotype_Genotype_Level_02'
useGroups <- unique(project$Phenotype_Genotype_Level_02) %>% sort()
name <- 'yap1_48hpf_Level02'

# ---- Make heatmap all meox1 peaks ----
heatmap_yap1_level02 <- make_meox_heatmaps(project,
                                           groupBy,
                                           useGroups,
                                           name)

# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 50000,
                                          down = 50000)


# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 9000)


# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 500)


################################################################################
# SUPP  YAP LEVEL_02: TRACKPLOT, VEC ONLY
################################################################################
# ---- variables ----
project <- project.yap1.level02
groupBy <- 'Phenotype_Genotype_Level_02'
useGroups <- c('VEC_Mutant', 'VEC_WT')
name <- 'yap1_48hpf_VECwtmt_Level02'

# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 50000,
                                          down = 50000)


# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 9000)


# ---- Make trackplot ----
plot_yap1_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 500)

################################################################################
# SUPP  YAP LEVEL_02: TRACKPLOT, VEC ONLY (marking all peaks)
################################################################################
project <- project.yap1.level02
groupBy <- 'Phenotype_Genotype_Level_02'
useGroups <- c('VEC_WT', 'VEC_Mutant')
name <- 'yap1_48hpf_VECwtmt_Level02'

p_meox1 <- plotBrowserTrack(ArchRProj = project,
                            geneSymbol = 'meox1',
                            groupBy = groupBy,
                            useGroups = useGroups,
                            upstream = 40000,
                            downstream = 9000)

plotName <- sprintf('%s_PlotTracksAllPeaks_down9000_up40000_meox1.pdf', name)
plotPDF(plotList = p_meox1,
        name = plotName,
        ArchRProj = project,
        addDOC = FALSE, width = 8, height = 5)



################################################################################
# SUPP  PROX1A LEVEL_02: TRACKPLOT, WT ONLY
################################################################################
# ---- variables ----
project <- project.4dpf.level02
groupBy <- 'Phenotype_Genotype'
useGroups <- unique(project$Phenotype_Genotype)[grepl(x = unique(project$Phenotype_Genotype), pattern = '_wt')] %>% sort()
name <- 'prox1a_4dpf_WT_Level02'
# ---- Make trackplot ----
plot_prox1a_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 50000,
                                          down = 50000)


# ---- Make trackplot ----
plot_prox1a_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 9000)


# ---- Make trackplot ----
plot_prox1a_level02 <- make_p3p4_trackplots(project,
                                          groupBy,
                                          useGroups,
                                          name,
                                          up = 40000,
                                          down = 500)

################################################################################
# SUPP  PROX1A LEVEL_02: TRACKPLOT, WT ONLY (all peaks)
################################################################################
project <- project.4dpf.level02
groupBy <- 'Phenotype_Genotype'
useGroups <- unique(project$Phenotype_Genotype)[grepl(x = unique(project$Phenotype_Genotype), pattern = '_wt')] %>% sort()
name <- 'prox1a_4dpf_WT_Level02'

p_meox1 <- plotBrowserTrack(ArchRProj = project,
                            geneSymbol = 'meox1',
                            groupBy = groupBy,
                            useGroups = useGroups,
                            upstream = 40000,
                            downstream = 9000)

plotName <- sprintf('%s_PlotTracksAllPeaks_down9000_up40000_meox1.pdf', name)
plotPDF(plotList = p_meox1,
        name = plotName,
        ArchRProj = project,
        addDOC = FALSE, width = 8, height = 5)


################################################################################
# SUPP  FIMO RESULTS HEATMAP TEAD MOTIFS
################################################################################
# ---- All meox1 peaks ----
# get count data for fimo results
fimo.count.results <- fimo.results %>%
  filter(., p.value < 0.001) %>%
  dplyr::group_by(., sequence_name) %>%
  dplyr::count(., motif_alt_id)

#need to add rows for a peak that had no hits
fimo.count.results.add.on <- data.frame(sequence_name = rep('12:27419091−27419591', 10),
                                        motif_alt_id = unique(fimo.count.results$motif_alt_id),
                                        n = rep(0, 10))

fimo.count.results.all <- rbind(fimo.count.results, fimo.count.results.add.on)
fimo.count.results.all$sequence_name <- factor(fimo.count.results.all$sequence_name, levels = order.meox1.peaks)

heatmap.delta <- ggplot(fimo.count.results.all, aes(x=sequence_name, y=motif_alt_id, fill=n)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Blues", direction=1, limits = c(1, 4)) +
  labs(title = "Heatmap, meox1 peaks", fill = '# TEAD motif match\npval<0.001') +
  xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)

heatmap.delta
FileName <- sprintf('%s/yap1_48hpf_Level02_meox1_peaks_TEADmotif_hits_heatmap.EPS', figure.dir.out)
ggsave(FileName,
       heatmap.delta,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)

# ---- selected peaks ----
fimo.count.results.selected <- fimo.count.results %>%
  filter(., sequence_name %in% order.meox1.peaks.short.tead)

fimo.count.results.selected$sequence_name <- factor(fimo.count.results.selected$sequence_name, levels = order.meox1.peaks.short.tead)

heatmap.delta <- ggplot(fimo.count.results.selected, aes(x=sequence_name, y=motif_alt_id, fill=n)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Blues", direction=1, limits = c(1, 4)) +
  labs(title = "Heatmap, meox1 peaks", fill = '# TEAD motif match\npval<0.001') +
  xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)

heatmap.delta
FileName <- sprintf('%s/yap1_48hpf_Level02_meox1_P3P4peaks_TEADmotif_hits_heatmap.EPS', figure.dir.out)
ggsave(FileName,
       heatmap.delta,
       device="eps",
       width=15,
       height=15,
       units="cm",
       scale=1)


################################################################################
# ADDITIONAL CODE AFTER REVISION
################################################################################

################################################################################
# UMAPS
################################################################################
# ---- Level 01: meox1 split genotype ----
gene <- 'meox1'
# WT
idxSample <- BiocGenerics::which(project.yap1.level01$Sample %in% "wt")
cellsSample <- project.yap1.level01$cellNames[idxSample]
project.yap1.level01_wt <-  project.yap1.level01[cellsSample,]
project.yap1.level01_wt <- addImputeWeights(project.yap1.level01_wt, reducedDims = 'IterativeLSI')

p1_wt <- plotEmbedding(
  ArchRProj = project.yap1.level01_wt,
  colorBy = "GeneScoreMatrix",
  name = gene,
  embedding = "UMAP",
  size = 1,
  imputeWeights = getImputeWeights(project.yap1.level01_wt))

x_wt <- ggplot_build(p1_wt)$data[[1]]
x_wt$Genotype = "WT"

# Mutant
idxSample <- BiocGenerics::which(project.yap1.level01$Sample %in% "mut")
cellsSample <- project.yap1.level01$cellNames[idxSample]
project.yap1.level01_mt <-  project.yap1.level01[cellsSample,]
project.yap1.level01_mt <- addImputeWeights(project.yap1.level01_mt, reducedDims = 'IterativeLSI')
p1_mt <- plotEmbedding(
  ArchRProj = project.yap1.level01_mt,
  colorBy = "GeneScoreMatrix",
  name = gene,
  embedding = "UMAP",
  size = 1,
  imputeWeights = getImputeWeights(project.yap1.level01_mt))

x_mt <- ggplot_build(p1_mt)$data[[1]]
x_mt$Genotype = "Mutant"

combined_df <- rbind(x_wt, x_mt)
gex.plot <- ggplot(combined_df, aes(x=x, y=y)) + geom_point(aes(color=value), size=1) +
  theme_classic() + theme(legend.position='top') +
  scale_colour_gradientn(colours=c("#d9d9d9", "#7a0177")) +
  facet_grid(cols = vars(Genotype)) + NoAxes()

# gex.plot

enrichment.export <- sprintf('%s_snATAC_yap1_L1_%s_split_genotype.pdf', figure.dir.out, gene)
ggsave(plot = gex.plot,
       filename =  enrichment.export,
       device="pdf",
       width=30,
       height=17,
       units="cm",
       scale=1)



# ---- Level 02: GEX ccn1 ----
cap <- plot_GEX_umaps(project = project.yap1.level02,
                      gene = "ccn1",
                      title_addition = "",
                      exportFilePath = figure.dir.out,
                      exportFileName = "snATAC_yap1_L2",
                      size_h = 6,
                      size_w = 10)


# ---- Level 02: GEX cdh6 ----
cap <- plot_GEX_umaps(project = project.yap1.level02,
                      gene = "cdh6",
                      title_addition = "",
                      exportFilePath = figure.dir.out,
                      exportFileName = "snATAC_yap1_L2",
                      size_h = 6,
                      size_w = 10)

# ---- Level 02: GEX ccn2a ----
cap <- plot_GEX_umaps(project = project.yap1.level02,
                      gene = "ccn2a",
                      title_addition = "",
                      exportFilePath = figure.dir.out,
                      exportFileName = "snATAC_yap1_L2",
                      size_h = 6,
                      size_w = 10)

# ---- Level 02: GEX meox1 ----
cap <- plot_GEX_umaps(project = project.yap1.level02,
                      gene = "meox1",
                      title_addition = "",
                      exportFilePath = figure.dir.out,
                      exportFileName = "snATAC_yap1_L2",
                      size_h = 6,
                      size_w = 10)

# ---- Level 02: cdh5 split genotype ----
gene <- 'cdh5'
# WT
idxSample <- BiocGenerics::which(project.yap1.level02$Sample %in% "wt")
cellsSample <- project.yap1.level02$cellNames[idxSample]
project.yap1.level02_wt <-  project.yap1.level02[cellsSample,]
project.yap1.level02_wt <- addImputeWeights(project.yap1.level02_wt, reducedDims = 'IterativeLSI_Level_02')

p1_wt <- plotEmbedding(
  ArchRProj = project.yap1.level02_wt,
  colorBy = "GeneScoreMatrix",
  name = gene,
  embedding = "UMAP",
  size = 1,
  imputeWeights = getImputeWeights(project.yap1.level02_wt))

x_wt <- ggplot_build(p1_wt)$data[[1]]
x_wt$Genotype = "WT"

# Mutant
idxSample <- BiocGenerics::which(project.yap1.level02$Sample %in% "mut")
cellsSample <- project.yap1.level02$cellNames[idxSample]
project.yap1.level02_mt <-  project.yap1.level02[cellsSample,]
project.yap1.level02_mt <- addImputeWeights(project.yap1.level02_mt, reducedDims = 'IterativeLSI_Level_02')
p1_mt <- plotEmbedding(
  ArchRProj = project.yap1.level02_mt,
  colorBy = "GeneScoreMatrix",
  name = gene,
  embedding = "UMAP",
  size = 1,
  imputeWeights = getImputeWeights(project.yap1.level02_mt))

x_mt <- ggplot_build(p1_mt)$data[[1]]
x_mt$Genotype = "Mutant"

combined_df <- rbind(x_wt, x_mt)
gex.plot <- ggplot(combined_df, aes(x=x, y=y)) + geom_point(aes(color=value), size=1) +
  theme_classic() + theme(legend.position='top') +
  scale_colour_gradientn(colours=c("#d9d9d9", "#7a0177")) +
  facet_grid(cols = vars(Genotype)) + NoAxes()

# gex.plot

enrichment.export <- sprintf('%ssnATAC_yap1_L2_%s_split_genotype.pdf', figure.dir.out, gene)
ggsave(plot = gex.plot,
       filename =  enrichment.export,
       device="pdf",
       width=30,
       height=12,
       units="cm",
       scale=1)



################################################################################
# DOTPLOTS
################################################################################
# ---- Level 01: meox1 phenotype_genotype ----
L1_dotplot <-         make_dotplot(project = project.yap1.level01,
                                   groupBy = 'Phenotype_Genotype',
                                   useGroups = project.yap1.level01$Phenotype_Genotype %>% unique() %>% sort(),
                                   groupName = project.yap1.level01$Phenotype_Genotype %>% unique() %>% sort(),
                                   markerGenes = 'meox1',
                                   accessibilityColour = c("#d9d9d9", "#7a0177"),
                                   rotateaxis = FALSE)


fileName <- sprintf('%s/DotPlot_snATAC_yap1_L1_meox1_Phenotype_Genotype.pdf', figure.dir.out)
ggsave(filename =fileName,
       plot =  L1_dotplot,
       device = 'pdf',
       width = 10,
       height = 10,
       units="cm",
       scale=1)

L1_dotplot <-         make_dotplot(project = project.yap1.level01,
                                   groupBy = 'Phenotype_Genotype',
                                   useGroups = project.yap1.level01$Phenotype_Genotype %>% unique() %>% sort(),
                                   groupName = project.yap1.level01$Phenotype_Genotype %>% unique() %>% sort(),
                                   markerGenes = 'meox1',
                                   accessibilityColour = c("#d9d9d9", "#7a0177"),
                                   rotateaxis = TRUE)


fileName <- sprintf('%s/DotPlot_snATAC_yap1_L1_meox1_Phenotype_Genotype_rotated.pdf', figure.dir.out)
ggsave(filename =fileName,
       plot =  L1_dotplot,
       device = 'pdf',
       width = 10,
       height = 10,
       units="cm",
       scale=1)

# ---- Level 01: meox1 phenotype, genotype ratio dotplot ----
df.plotting <- L1_dotplot$data %>% #inherited from above
  mutate(., genotype = ifelse(grepl('Mutant', Groups), 'mutant', 'wt')) %>%
  mutate(., simplename = stringi::stri_replace_all_regex(str = Groups, pattern = '_Mutant|_WT', replacement = ''))
rownames(df.plotting) <- NULL
df.plotting$simplename <- factor(df.plotting$simplename)
df.plotting$genotype <- factor(df.plotting$genotype, levels=c('wt', 'mutant'))
df.plotting.fc <- df.plotting %>%
  group_by(simplename) %>%
  summarise(., avg.exp.log2fc = log2(Score[2]/Score[1]), pct.exp.fc= Proportion[2]/Proportion[1]) #use score for unscaled

df.plotting.fc$simplename <- factor(df.plotting.fc$simplename, levels = Phenotype_ordered_Level_03)
df.plotting.fc$Delta <- as.factor('[WT]/[yap1-/-]')

dotplot.ratio		<- ggplot(df.plotting.fc) +
  geom_point(aes(x= simplename, y = Delta, color = avg.exp.log2fc, size = pct.exp.fc)) +
  scale_colour_gradient(low = '#332288',
                        high = '#CC6677',
                        na.value = "black",
                        guide = "colourbar",
                        aesthetics = "colour") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))  + xlab('') + ylab('') +
  ggtitle('Ratio of wt/yap1-/-\n% expressed and avg expression\nunscaled data')

fileName <- sprintf('%s/DotPlot_snATAC_yap1_L1_meox1_Phenotype_Genotype_meox1_ratio_reverse_colours.pdf', figure.dir.out)
ggsave(filename =fileName,
       plot =  dotplot.ratio,
       device = 'pdf',
       width = 10,
       height = 10,
       units="cm",
       scale=1)



################################################################################
# TRACKPLOTS
################################################################################
# ---- Level 02: prox1a ----
gene <- 'prox1a'
groupBy <- "Phenotype_Level_02"
p_prox1a <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                             geneSymbol = gene,
                             groupBy = groupBy,
                             useGroups = L2_Phenotype_ordered_Level_02)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
# plotName <- sprintf('%s/snATAC_yap1_L2_Trackplot_%s.pdf',figure.dir.out, gene)
plotPDF(plotList = p_prox1a,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)

# ---- Level 02: lyve1b ----
gene <- 'lyve1b'
groupBy <- "Phenotype_Level_02"
p_t <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                        geneSymbol = gene,
                        groupBy = groupBy,
                        useGroups = L2_Phenotype_ordered_Level_02)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
plotPDF(plotList = p_t,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)

# ---- Level 02: flt1 ----
gene <- 'flt1'
groupBy <- "Phenotype_Level_02"
p_t <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                        geneSymbol = gene,
                        groupBy = groupBy,
                        useGroups = L2_Phenotype_ordered_Level_02,
                        upstream = 80000,
                        downstream = 20000)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
plotPDF(plotList = p_t,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)

# ---- Level 02: hand2 ----
gene <- 'hand2'
groupBy <- "Phenotype_Level_02"
p_t <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                        geneSymbol = gene,
                        groupBy = groupBy,
                        useGroups = L2_Phenotype_ordered_Level_02)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
plotPDF(plotList = p_t,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)

# ---- Level 02: ccn2a ----
gene <- 'ccn2a'
groupBy <- "Phenotype_Level_02"
p_t <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                        geneSymbol = gene,
                        groupBy = groupBy,
                        useGroups = L2_Phenotype_ordered_Level_02)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
plotPDF(plotList = p_t,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)
# ---- Level 02: meox1 ----
gene <- 'meox1'
groupBy <- "Phenotype_Level_02"
p_t <- plotBrowserTrack(ArchRProj = project.yap1.level02,
                        geneSymbol = gene,
                        groupBy = groupBy,
                        useGroups = L2_Phenotype_ordered_Level_02)

plotName <- sprintf('snATAC_yap1_L2_Trackplot_%s_%s.pdf', gene, groupBy)
plotPDF(plotList = p_t,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 3)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out, overwrite = T)
# ---- EMBO paper Level 02: prox1a ----
groupBy <- 'Phenotype_Genotype'
gene <- 'prox1a'
useGroups <- unique(project.4dpf.level02$Phenotype_Genotype)[grepl(x = unique(project.4dpf.level02$Phenotype_Genotype), pattern = '_wt')] %>% sort()

p_prox1a <- plotBrowserTrack(ArchRProj = project.4dpf.level02,
                             geneSymbol = gene,
                             groupBy = groupBy,
                             useGroups = useGroups)

plotName <- sprintf('snATAC_EMBO_wtonly_L2_Trackplot_%s.pdf', gene)
# plotName <- sprintf('%s/snATAC_yap1_L2_Trackplot_%s.pdf',figure.dir.out, gene)
plotPDF(plotList = p_prox1a,
        name = plotName,
        ArchRProj = project.yap1.level02,
        addDOC = FALSE, width = 5, height = 5)

# grid::grid.newpage()
# grid::grid.draw(p_prox1a$prox1a)

file.copy(sprintf('%s/Plots/%s', getOutputDirectory(project.yap1.level02),plotName), figure.dir.out)


################################################################################
# SNAKEPLOT
################################################################################
#these plots were made here: code/snATACseq/analysis/Yap1_scatac_Level_02_DAP.R
path_allEC <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/240802_AllCells_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC_snakeplot.pdf"
file.copy(c(path_allEC), figure.dir.out)

