#! /usr/bin/Rscript

#   First author : Michelle Meier


#   Depends on R/4.2.0.Core + Renv 

################################################################################
# SET VARIABLES 
################################################################################
#Load variables
source('/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/yap1/Yap1_scatac_Level_02_DAP_crosscheck_4dfp_peaks_variables.R')
region <- getPeakSet(ArchRProj=project.yap1.level02)
################################################################################
# SPECIFIC DAP ANALYSIS
################################################################################
#to cross check stat significance of selected genes
project.yap1.level02$LECVEC_comb_pheno <- gsub(pattern = "VEC|Specified_LEC", replacement = "LEC_VEC", project.yap1.level02$Phenotype_Level_02)

# ---- lyve1b: VEC/LEC vs AEC ----
#add meta
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.yap1.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "LECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'LEC_VEC',
  bgdGroups = "AEC"
)

markerList_lyve <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_lyve <- lapply(markerList_lyve, as.data.frame)
complete_peaks_df_lyve <- markerDF_list_lyve[[1]] %>% 
  mutate(., comparison = "(Specified_LEC,VEC) vs AEC")

#filter to facilitate search
complete_peaks_df_lyve %>% 
  filter(., seqnames == "chr18" & start > 16741200 & start < 16745000)


# ---- prox1a: LEC vs AEC/VEC ----
project.yap1.level02$AECVEC_comb_pheno <- gsub(pattern = "VEC|AEC", replacement = "AEC_VEC", project.yap1.level02$Phenotype_Level_02)

#add meta
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.yap1.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "AECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'Specified_LEC',
  bgdGroups = "AEC_VEC"
)

markerList_prox1 <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_prox1 <- lapply(markerList_prox1, as.data.frame)
complete_peaks_df_prox1 <- markerDF_list_prox1[[1]] %>% 
  mutate(., comparison = "Specified_LEC vs (AEC,VEC)")

#filter to facilitate search
complete_peaks_df_prox1 %>% 
  filter(., seqnames == "chr17" & start > 32867850 & start < 32874650)


# ---- flt1: AEC vs VEC/LEC ----
#add meta
markerPeaks <- getMarkerFeatures(
  ArchRProj = project.yap1.level02, 
  useMatrix = "PeakMatrix", 
  groupBy = "LECVEC_comb_pheno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = 'AEC',
  bgdGroups = "LEC_VEC"
)

markerList_flt1 <- getMarkers(markerPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.01")
markerDF_list_flt1 <- lapply(markerList_flt1, as.data.frame)
complete_peaks_df_flt1 <- markerDF_list_flt1[[1]] %>% 
  mutate(., comparison = "AEC vs (Specified_LEC,VEC)")

#filter to facilitate search
complete_peaks_df_flt1 %>% 
  filter(., seqnames == "chr24" & start > 21774090 & start < 21775200)
complete_peaks_df_flt1 %>% 
  filter(., seqnames == "chr24" & start > 21785050 & start < 21785250)
complete_peaks_df_flt1 %>% 
  filter(., seqnames == "chr24" & start > 21800000 & start < 21800150)

# ---- combine and export ----
total_df <- rbind(complete_peaks_df_lyve, complete_peaks_df_flt1, complete_peaks_df_prox1)

# add peak info from peak sets
peak_set 					<-	data.frame(getPeakSet(ArchRProj=project.yap1.level02)) %>%
      select(., -c("replicateScoreQuantile", "groupScoreQuantile","Reproducibility", "GroupReplicate", "N", "idx"))

total_df <- full_join(total_df, peak_set, by=c("seqnames", "start", "end"))

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_VEC_specifiedLEC_fdr1_fc001_specific_comparisons.txt', save.dir)
write.table(x = total_df, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)

# ---- export for genes of interest only ----

total_df_sub <- total_df %>%
filter(., nearestGene %in% c("prox1a", "lyve1b", "ccn2a", "flt1"))

save.dir.peaks <- sprintf('%s/MarkerPeaks_AEC_VEC_specifiedLEC_fdr1_fc001_specific_comparisons_prox1a_lyve1b_ccn2a_flt1.txt', save.dir)
write.table(x = total_df_sub, file = save.dir.peaks, quote = F, sep = '\t', col.names = T, row.names = T)



################################################################################
# MAKE TRACKPLOTS FOR SELECTED GENES: HIGHLIGHT PEAKS
################################################################################
# ---- lyve1b ----
#make trackplot and highlight peak to double check
region_lyve1b_1 <- region[region@seqnames == 'chr18' & region@ranges@start == 16744026]
region_lyve1b_2 <- region[region@seqnames == 'chr18' & region@ranges@start == 16741309]

region_filtered_lyve1b <- c(region_lyve1b_1, region_lyve1b_2)

#make trackplot with only those peaks
lyve1b_trackplot <- plotBrowserTrack(ArchRProj = project.yap1.level02, 
                                     geneSymbol = 'lyve1b', 
                                     features = region_filtered_lyve1b,
                                     groupBy = "Phenotype_Level_02", 
                                     useGroups = c('AEC', 'Specified_LEC','VEC'))
#visualise 
#grid::grid.newpage()
#grid::grid.draw(lyve1b_trackplot$lyve1b)

#save
plotName <- "MarkerPeaks_lyve1b_Phenotype_Level_02_highlight_selected_peaks.pdf"
plotPDF(plotList = lyve1b_trackplot, 
        name = plotName, 
        ArchRProj = project.yap1.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- prox1a ----
#make trackplot and highlight peak to double check
region_prox1_1 <- region[region@seqnames == 'chr17' & region@ranges@start == 32874574]
region_prox1_2 <- region[region@seqnames == 'chr17' & region@ranges@start == 32867978]

region_filtered_prox1 <- c(region_prox1_1, region_prox1_2)

#make trackplot with only those peaks
prox1a_trackplot <- plotBrowserTrack(ArchRProj = project.yap1.level02, 
                                     geneSymbol = 'prox1a', 
                                     features = region_filtered_prox1,
                                     groupBy = "Phenotype_Level_02", 
                                     useGroups = c('AEC', 'Specified_LEC','VEC'))
#visualise 
#grid::grid.newpage()
#grid::grid.draw(prox1a_trackplot$prox1a)

#save
plotName <- "MarkerPeaks_prox1a_Phenotype_Level_02_highlight_selected_peaks.pdf"
plotPDF(plotList = prox1a_trackplot, 
        name = plotName, 
        ArchRProj = project.yap1.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- ccn2a ----
region_ccn2a_1 <- region[region@seqnames == 'chr20' & region@ranges@start == 25340430]
region_ccn2a_2 <- region[region@seqnames == 'chr20' & region@ranges@start == 25320594]
#combined regions 
region_filtered_ccn2a <- c(region_ccn2a_1, region_ccn2a_2)

#make trackplot with only those peaks
ccn2a_trackplot <- plotBrowserTrack(ArchRProj = project.yap1.level02, 
                 geneSymbol = 'ccn2a', 
                 groupBy = "Phenotype_Level_02", 
                 features = region_filtered_ccn2a,
                 useGroups = c('AEC', 'Specified_LEC','VEC'))
#visualise 
#grid::grid.newpage()
#grid::grid.draw(ccn2a_trackplot$ccn2a)

#save
plotName <- "MarkerPeaks_ccn2a_Phenotype_Genotype_highlight_selected_peaks.pdf"
plotPDF(plotList = ccn2a_trackplot, 
        name = plotName, 
        ArchRProj = project.yap1.level02, 
        addDOC = FALSE, width = 10, height = 5)

# ---- flt1 ----
region_flt1_1 <- region[region@seqnames == 'chr24' & region@ranges@start ==21774996 ]
region_flt1_2 <- region[region@seqnames == 'chr24' & region@ranges@start == 21785155]

#combined regions 
region_filtered_flt1 <- c(region_flt1_1, region_flt1_2)

#make trackplot with only those peaks
flt1_trackplot <- plotBrowserTrack(ArchRProj = project.yap1.level02, 
                                   geneSymbol = 'flt1', 
                                   features = region_filtered_flt1,
                                   groupBy = "Phenotype_Level_02", 
                                   useGroups = c('AEC', 'Specified_LEC','VEC'), 
                                   upstream = 80000,
                                   downstream = 20000)
#visualise 
#grid::grid.newpage()
#grid::grid.draw(flt1_trackplot$flt1)

#save
plotName <- "MarkerPeaks_flt1_extendedrange_Phenotype_Level_02_highlight_selected_peaks.pdf"
plotPDF(plotList = flt1_trackplot, 
        name = plotName, 
        ArchRProj = project.yap1.level02, 
        addDOC = FALSE, width = 10, height = 5)
