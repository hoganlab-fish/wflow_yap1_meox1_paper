#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited October 2024

#   Depends on R/4.2.0.Core + Renv 

################################################################################
# SET VARIABLES 
################################################################################
#Load variables
source('/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/Rev_figures_snATAC_yap1_variables.R')

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
#these plots were made here: /hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/yap1/Yap1_scatac_Level_02_DAP.R
#copying them across to revision folder
path_allEC <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/240802_AllCells_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC_snakeplot.pdf"
path_eachCT <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/240802_eachCT_sep_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC_snakeplot.pdf"
file.copy(c(path_allEC, path_eachCT), figure.dir.out)


