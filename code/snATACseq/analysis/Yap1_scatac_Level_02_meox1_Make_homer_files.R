#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited November 2023

#   Depends on R/4.2.0.Core + Renv 

################################################################################
# SET VARIABLES 
################################################################################
#Load variables
source('/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/Yap1_scatac_Level_02_meox1_Make_homer_files_variables.R')

################################################################################
# MAKE HOMER FILES
################################################################################
# ---- Meox1 associated peaks ----
dap.table.all.filtered <- dap.table.all %>%
  filter(., nearestGene == 'meox1')

dap.table.all.filtered_peak_file  <-  data.frame(ID = rownames(dap.table.all.filtered),
                                                    chr = dap.table.all.filtered$Chr,
                                                    start = dap.table.all.filtered$Start,
                                                    end = dap.table.all.filtered$End,
                                                    strand = ".")
fileName <- sprintf('%s/%s_meox1_associated_peaks.txt', save.dir, project.name)
write.table(dap.table.all.filtered_peak_file,
            file = fileName,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

