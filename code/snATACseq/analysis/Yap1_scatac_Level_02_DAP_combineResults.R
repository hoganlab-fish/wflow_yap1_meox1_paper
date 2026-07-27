#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Depends on R/4.2.0.Core + Renv

################################################################################
# COMBINE DAP FOR RELEVANT CELL TYPES FOR EXTENDED DATA TABLES
################################################################################


# ---- set up: load libraries ----
library(here)
library(tidyverse)


# ---- set up: load data ----
path_all <- list.files(path = here("output/analysis/D10029_yap1_snATAC_dataset/DAP_analysis"), full.names = T)
list_all_res <- lapply(path_all, read_tsv)
names(list_all_res) <- gsub("240802_", "",gsub("_Mutant.*", "", basename(path_all)) )



# ---- count DAPs ----

lapply(1:length(list_all_res), function(INDEX){
  NAME =names(list_all_res)[INDEX]
  ENTRY = list_all_res[[INDEX]]
  ENTRY %>% filter(., RawPval<0.05) %>%
    count(Status_Group_01) %>%
    mutate(celltype=NAME)
}) %>% bind_rows()


# ---- combine + export for VEC, Specified LEC and AEC for extended data table ----
df_all_export <- list_all_res[c("AEC", "Specified_LEC", "VEC")] %>%
  bind_rows(., .id="celltype")
#actually this is probably too big to fit in one sheet. Add a tab for each






