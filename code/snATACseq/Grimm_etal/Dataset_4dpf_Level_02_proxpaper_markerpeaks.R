#! /usr/bin/Rscript

#   First author : Michelle Meier


#   Depends on R/4.2.0.Core + Renv 

################################################################################
# SET VARIABLES 
################################################################################
#Load variables
source('/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/Dataset_4dpf_Level_02_proxpaper_markerpeaks_variables.R')

################################################################################
# FIND MARKER PEAKS 
################################################################################
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "Phenotype_Genotype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.1")
markerDF_list <- lapply(markerList, as.data.frame)

#combine into one big DF
complete_peaks_df <- purrr::map_df(markerDF_list, ~as.data.frame(.x), .id="cell_type") %>% 
  filter(., cell_type %in%  c('AEC_wt', 'VEC_wt', 'LEC_wt'))


# add peak info from peak sets
peak_set 					<-	data.frame(getPeakSet(ArchRProj=project.4dpf.level02)) %>%
      select(., -c("replicateScoreQuantile", "groupScoreQuantile","Reproducibility", "GroupReplicate", "N", "idx"))

complete_peaks_df <- full_join(complete_peaks_df, peak_set, by=c("seqnames", "start", "end"))

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_VEC_LEC_wt_fdr1_fc01.txt', save.dir)
write.table(x = complete_peaks_df, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)

#maybe just marker peaks isn't the ideal approach... We also have the mutant enriched cell type in there that we're not really interested in.

################################################################################
# MORE SPECIFIC DAP ANALYSIS
################################################################################
project.4dpf.level02$LECVEC_comb_pheno <- gsub(pattern = "LEC_wt|VEC_wt", replacement = "LEC_VEC_wt", project.4dpf.level02$Phenotype_Genotype)
# ---- lyve1b: VEC/LEC vs AEC ----
#add meta
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "LECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'LEC_VEC_wt',
  bgdGroups = "AEC_wt"
)

markerList_lyve <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_lyve <- lapply(markerList_lyve, as.data.frame)
complete_peaks_df_lyve <- markerDF_list_lyve[[1]] %>% 
  mutate(., comparison = "(LEC,VEC) vs AEC")


# ---- flt1: AEC vs VEC/LEC ----
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "LECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'AEC_wt',
  bgdGroups = "LEC_VEC_wt"
)

markerList_flt1 <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_flt1 <- lapply(markerList_flt1, as.data.frame)
complete_peaks_df_flt1 <- markerDF_list_flt1[[1]] %>% 
  mutate(., comparison = "AEC vs (LEC,VEC)")

# ---- prox1a+ccn2a: LEC vs VEC/AEC ----
project.4dpf.level02$AECVEC_comb_pheno <- gsub(pattern = "AEC_wt|VEC_wt", replacement = "AEC_VEC_wt", project.4dpf.level02$Phenotype_Genotype)
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "AECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'LEC_wt',
  bgdGroups = "AEC_VEC_wt"
)

markerList_prox <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_prox <- lapply(markerList_prox, as.data.frame)
complete_peaks_df_prox <- markerDF_list_prox[[1]] %>% 
  mutate(., comparison = "LEC vs (AEC,VEC)")


# ---- combine and export ----
total_df <- rbind(complete_peaks_df_lyve, complete_peaks_df_flt1, complete_peaks_df_prox)

# add peak info from peak sets
peak_set 					<-	data.frame(getPeakSet(ArchRProj=project.4dpf.level02)) %>%
      select(., -c("replicateScoreQuantile", "groupScoreQuantile","Reproducibility", "GroupReplicate", "N", "idx"))

total_df <- full_join(total_df, peak_set, by=c("seqnames", "start", "end"))

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_VEC_LEC_wt_fdr1_fc001_specific_comparisons.txt', save.dir)
write.table(x = total_df, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)

# ---- export for genes of interest only ----

total_df_sub <- total_df %>%
filter(., nearestGene %in% c("prox1a", "lyve1b", "ccn2a", "flt1"))

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_VEC_LEC_wt_fdr1_fc001_specific_comparisons_prox1a_lyve1b_ccn2a_flt1.txt', save.dir)
write.table(x = total_df_sub, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)

################################################################################
# MORE SPECIFIC DAP ANALYSIS, CONTINUED
################################################################################
# ---- ft1: AEC vs LEC ----
# for flt1, AEC vs LEC only
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "Phenotype_Genotype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'AEC_wt',
  bgdGroups = "LEC_wt"
)

markerList_flt1 <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_flt1 <- lapply(markerList_flt1, as.data.frame)
complete_peaks_df_flt1 <- markerDF_list_flt1[[1]] %>% 
  mutate(., comparison = "AEC_vs_LEC")

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_vs_LEC_wt_fdr1_fc001.txt', save.dir)
write.table(x = complete_peaks_df_flt1, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)



# ---- ft1: AEC vs VEC ----
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.4dpf.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "Phenotype_Genotype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'AEC_wt',
  bgdGroups = "VEC_wt"
)

markerList_flt1 <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
markerDF_list_flt1 <- lapply(markerList_flt1, as.data.frame)
complete_peaks_df_flt1 <- markerDF_list_flt1[[1]] %>% 
  mutate(., comparison = "AEC_vs_VEC")

complete_peaks_df_flt1 %>%
  filter(., seqnames == "chr24" & Log2FC > 0 & start > 21775000) %>% 
  View()


save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_vs_VEC_wt_fdr1_fc0.txt', save.dir)
write.table(x = complete_peaks_df_flt1, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)


################################################################################
# MAKE TRACKPLOTS FOR SELECTED GENES
################################################################################
selected_genes <- c('prox1a', 'lyve1b', 'hand2', 'flt1', 'ccn2a')

trackplots <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                            geneSymbol = selected_genes, 
                            groupBy = "Phenotype_Genotype", 
                            useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'))

plotName <- "MarkerPeaks_prox1a_lyve1b_hand2_flt1_ccn2a_Phenotype_Genotype.pdf"
plotPDF(plotList = trackplots, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- make extended trackplot for flt1 to capture whole gene body ----
flt1_long <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                               geneSymbol = 'flt1', 
                               groupBy = "Phenotype_Genotype", 
                               useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'), 
                              upstream = 80000,
                              downstream = 20000)
grid::grid.newpage()
grid::grid.draw(flt1_long$flt1)


plotName <- "MarkerPeaks_flt1_shifted_frame_Phenotype_Genotype.pdf"
plotPDF(plotList = flt1_long, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

#what does this look like for yap2 dataset?
flt_long_yapdataset <- plotBrowserTrack(ArchRProj = project.yap1.level02, 
                              geneSymbol = 'flt1', 
                              groupBy = "Phenotype_Level_02", 
                              upstream = 80000,
                              downstream = 20000)
grid::grid.newpage()
grid::grid.draw(flt_long_yapdataset$flt1)




################################################################################
# MAKE TRACKPLOTS FOR SELECTED GENES: HIGHLIGHT PEAKS
################################################################################
# get peak regions 
region <- getPeakSet(ArchRProj=project.4dpf.level02)

# ---- prox1a ----
region_prox1_1 <- region[region@seqnames == 'chr17' & region@ranges@start == 32867962]
#get P3 region
region_prox1_2 <- region[region@seqnames == 'chr17' & region@ranges@start == 32874572]
#combined regions 
region_filtered_prox1a <- c(region_prox1_1, region_prox1_2)

#make trackplot with only those peaks
prox1a_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                 geneSymbol = 'prox1a', 
                 groupBy = "Phenotype_Genotype", 
                 features = region_filtered_prox1a,
                 useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'))
#visualise 
grid::grid.newpage()
grid::grid.draw(prox1a_trackplot$prox1a)

#save
plotName <- "MarkerPeaks_prox1_Phenotype_Genotype_highlight_selected_peaks.pdf"
plotPDF(plotList = prox1a_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- ccn2a ----
region_ccn2a_1 <- region[region@seqnames == 'chr20' & region@ranges@start == 25340480]
#get P3 region
region_ccn2a_2 <- region[region@seqnames == 'chr20' & region@ranges@start == 25320605]
#combined regions 
region_filtered_ccn2a <- c(region_ccn2a_1, region_ccn2a_2)

#make trackplot with only those peaks
ccn2a_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                 geneSymbol = 'ccn2a', 
                 groupBy = "Phenotype_Genotype", 
                 features = region_filtered_ccn2a,
                 useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'))
#visualise 
#grid::grid.newpage()
#grid::grid.draw(ccn2a_trackplot$ccn2a)

#save
plotName <- "MarkerPeaks_ccn2a_Phenotype_Genotype_highlight_selected_peaks.pdf"
plotPDF(plotList = ccn2a_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- lyve1b ----
region_lyve1b_1 <- region[region@seqnames == 'chr18' & region@ranges@start == 16741305]
region_lyve1b_2 <- region[region@seqnames == 'chr18' & region@ranges@start == 16743452]
region_lyve1b_3 <- region[region@seqnames == 'chr18' & region@ranges@start == 16743985]
#combined regions 
region_filtered_lyve1b <- c(region_lyve1b_1, region_lyve1b_2,region_lyve1b_3)

#make trackplot with only those peaks
lyve1b_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                                     geneSymbol = 'lyve1b', 
                                     features = region_filtered_lyve1b,
                                     groupBy = "Phenotype_Genotype", 
                                     # ylim = c(0, 0.95),
                                     useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'))
#visualise 
grid::grid.newpage()
grid::grid.draw(lyve1b_trackplot$lyve1b)

#save
plotName <- "MarkerPeaks_lyve1b_Phenotype_Genotype_highlight_selected_peaks.pdf"
plotPDF(plotList = lyve1b_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- flt1: VECLEC combined vs AEC ----
region_flt1_1 <- region[region@seqnames == 'chr24' & region@ranges@start ==21775017 ]
region_flt1_2 <- region[region@seqnames == 'chr24' & region@ranges@start == 21785186]
region_flt1_3 <- region[region@seqnames == 'chr24' & region@ranges@start == 21800033]
#combined regions 
region_filtered_flt1 <- c(region_flt1_1, region_flt1_2,region_flt1_3)

#make trackplot with only those peaks
flt1_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                                     geneSymbol = 'flt1', 
                                     features = region_filtered_flt1,
                                     groupBy = "Phenotype_Genotype", 
                                     useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'), 
                                     upstream = 80000,
                                     downstream = 20000)
#visualise 
grid::grid.newpage()
grid::grid.draw(flt1_trackplot$flt1)

#save
plotName <- "MarkerPeaks_flt1_extendedrange_Phenotype_Genotype_highlight_selected_peaks_VECLEC_combined.pdf"
plotPDF(plotList = flt1_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- flt1: AEC marker ----
region_flt1_1 <- region[region@seqnames == 'chr24' & region@ranges@start ==21775017 ]
region_flt1_2 <- region[region@seqnames == 'chr24' & region@ranges@start == 21785186]
region_flt1_3 <- region[region@seqnames == 'chr24' & region@ranges@start == 21800033]
#combined regions 
region_filtered_flt1 <- c(region_flt1_1, region_flt1_2,region_flt1_3)

#make trackplot with only those peaks
flt1_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                                   geneSymbol = 'flt1', 
                                   features = region_filtered_flt1,
                                   groupBy = "Phenotype_Genotype", 
                                   useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'), 
                                   upstream = 80000,
                                   downstream = 20000)
#visualise 
grid::grid.newpage()
grid::grid.draw(flt1_trackplot$flt1)

#save
plotName <- "MarkerPeaks_flt1_extendedrange_Phenotype_Genotype_highlight_selected_peaks_VECLEC_combined.pdf"
plotPDF(plotList = flt1_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- flt1: AEC vs LEC ----
region_flt1_1 <- region[region@seqnames == 'chr24' & region@ranges@start ==21765485 ]
region_flt1_2 <- region[region@seqnames == 'chr24' & region@ranges@start == 21775017]
region_flt1_3 <- region[region@seqnames == 'chr24' & region@ranges@start == 21778283]
region_flt1_4 <- region[region@seqnames == 'chr24' & region@ranges@start == 21778811]
region_flt1_5 <- region[region@seqnames == 'chr24' & region@ranges@start == 21781432]
region_flt1_6 <- region[region@seqnames == 'chr24' & region@ranges@start == 21785186]
#combined regions 
region_filtered_flt1 <- c(region_flt1_1, region_flt1_2,region_flt1_3, region_flt1_4, region_flt1_5, region_flt1_6)

#make trackplot with only those peaks
flt1_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                                   geneSymbol = 'flt1', 
                                   # features = region_filtered_flt1,
                                   groupBy = "Phenotype_Genotype", 
                                   useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'), 
                                   upstream = 80000,
                                   downstream = 20000)
#visualise 
grid::grid.newpage()
grid::grid.draw(flt1_trackplot$flt1)

#save
plotName <- "MarkerPeaks_flt1_extendedrange_Phenotype_Genotype_highlight_selected_peaks_AEC_vs_LEC.pdf"
plotPDF(plotList = flt1_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)


# ---- flt1: AEC vs LEC ----
region_flt1_1 <- region[region@seqnames == 'chr24' & region@ranges@start ==21805248 ]
region_flt1_2 <- region[region@seqnames == 'chr24' & region@ranges@start == 21806771]
region_flt1_3 <- region[region@seqnames == 'chr24' & region@ranges@start == 21808118]
region_flt1_4 <- region[region@seqnames == 'chr24' & region@ranges@start == 21808627]
region_flt1_5 <- region[region@seqnames == 'chr24' & region@ranges@start == 21812254]

#combined regions 
region_filtered_flt1 <- c(region_flt1_1, region_flt1_2,region_flt1_3, region_flt1_4, region_flt1_5)

#make trackplot with only those peaks
flt1_trackplot <- plotBrowserTrack(ArchRProj = project.4dpf.level02, 
                                   geneSymbol = 'flt1', 
                                   features = region_filtered_flt1,
                                   groupBy = "Phenotype_Genotype", 
                                   useGroups = c("AEC_wt", 'LEC_wt', 'VEC_wt'), 
                                   upstream = 80000,
                                   downstream = 20000)
#visualise 
grid::grid.newpage()
grid::grid.draw(flt1_trackplot$flt1)

#save
plotName <- "MarkerPeaks_flt1_extendedrange_Phenotype_Genotype_highlight_selected_peaks_AEC_vs_VEC.pdf"
plotPDF(plotList = flt1_trackplot, 
        name = plotName, 
        ArchRProj = project.4dpf.level02, 
        addDOC = FALSE, width = 10, height = 5)



