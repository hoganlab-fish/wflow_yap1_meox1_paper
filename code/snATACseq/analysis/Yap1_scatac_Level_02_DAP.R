#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited August 2024

#   Depends on R/4.2.0.Core + Renv 

################################################################################
# SET VARIABLES 
################################################################################
#Load variables
source('./yap1/Yap1_scatac_Level_02_DAP_variables.R')

################################################################################
# RUN DAP ANALYSIS: EACH CELL TYPE
################################################################################
# ---- set general variables ----
groupBy					<-	"Phenotype_Genotype_Level_02"
FDRFilter 				<- 	FALSE
FDRThreshold 			<- 	1
RawPvalFilter			<- 	TRUE
RawPvalThreshold		<- 	1
Log2FCFilter			<- 	TRUE
Log2FCThreshold			<- 	0

# ---- set up lapply ----
list_all_results_each_CT <- lapply(project$Phenotype_Level_02 %>% unique(), function(x){
  group_1					<-sprintf("%s_Mutant", x)
  group_2_background		<-sprintf("%s_WT", x)
  dap_table 				<- 	dap_analysis(project = project, 
                                 useGroup = group_1, 
                                 bgdGroup = group_2_background, 
                                 groupBy = groupBy,
                                 exportFilePath = NULL,
                                 exportFileName = NULL,
                                 chrStatus = NULL,
                                 FDRFilter = FDRFilter,
                                 FDRThreshold = FDRThreshold,
                                 RawPvalFilter = RawPvalFilter,
                                 RawPvalThreshold = RawPvalThreshold,
                                 Log2FCFilter = Log2FCFilter,
                                 Log2FCThreshold = Log2FCThreshold)
  return(dap_table)
})
names(list_all_results_each_CT) <- project$Phenotype_Level_02 %>% unique()
################################################################################
# RUN DAP ANALYSIS: ALL CELLTYPES MUTANT VS WT
################################################################################
# ---- set general variables ----
group_1					<-"Mutant"
group_2_background		<-"WT"
groupBy					<-	"Genotype"
FDRFilter 				<- 	FALSE
FDRThreshold 			<- 	1
RawPvalFilter			<- 	TRUE
RawPvalThreshold		<- 	1
Log2FCFilter			<- 	TRUE
Log2FCThreshold			<- 	0

# ---- run DAP ----
dap_table 				<- 	dap_analysis(project = project, 
                                 useGroup = group_1, 
                                 bgdGroup = group_2_background, 
                                 groupBy = groupBy,
                                 exportFilePath = NULL,
                                 exportFileName = '240802_AllCells_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC.txt',
                                 chrStatus = NULL,
                                 FDRFilter = FDRFilter,
                                 FDRThreshold = FDRThreshold,
                                 RawPvalFilter = RawPvalFilter,
                                 RawPvalThreshold = RawPvalThreshold,
                                 Log2FCFilter = Log2FCFilter,
                                 Log2FCThreshold = Log2FCThreshold)


################################################################################
# SNAPKEPLOT FOR EACH ALL CELLS MUTANT VS WT 
################################################################################
snakeplot_all_ECs <- make_snakeplot(dap_table_genotype = dap_table, 
                                    CT = 'AllCells', 
                                    gene.list = hippo_targets,
                                    pval = 0.05)
snakeplot_all_ECs

################################################################################
# SNAPKEPLOTs FOR EACH CELLTYPE 
################################################################################
list_all_snakeplots <- lapply( 1:length(list_all_results_each_CT), function(x){
  snake <- make_snakeplot(dap_table_genotype = list_all_results_each_CT[[x]], 
                          CT =names(list_all_results_each_CT)[x], 
                          gene.list = hippo_targets,
                          pval = 0.05)
  return(snake)
})

all_snakes <- ggarrange(plotlist = list_all_snakeplots)
all_snakes

################################################################################
# SAVE PLOTS 
################################################################################
# ---- all endothelial cells ----
filename.snake <- sprintf('%s/240802_AllCells_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC_snakeplot.pdf', figure.dir.out)
ggsave(filename = filename.snake,
       plot = snakeplot_all_ECs,
       device = 'pdf',
       width = 8,
       height = 5)

# ---- all endothelial cells, each cell type separately ----
filename.snake <- sprintf('%s/240802_eachCT_sep_Mutant_vs_WT_all_DAPs_filteredBy_rawPval_Log2FC_snakeplot.pdf', figure.dir.out)
ggsave(filename = filename.snake,
       plot = all_snakes,
       device = 'pdf',
       width = 20,
       height = 5)
