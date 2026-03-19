#!/usr/bin/Rscript
# update of Meox1_Level_03_L3_Seurat_cluster_predicted_phenotype_Genotype_composition_for_prism.csv
library(dplyr)
library(qs2)

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_path <- "../Revision_analysis/queue_plots/Meox1_Level_03_L3_Seurat_cluster_predicted_phenotype_Genotype_composition_for_prism.csv"

data <- qs2::qs_read(infile_path)
table(data@meta.data[["L3_cluster_id"]])

for_prism <- data@meta.data %>%
    mutate(Genotype = recode(Genotype, 
        "wildtype" = "WT", 
        "meox1_mutant" = "Mutant"
    )) %>%
    group_by(Genotype, L3_cluster_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Genotype) %>%
    mutate(percent = count / sum(count) * 100) %>%
    ungroup()

write.csv(
    for_prism, outfile_path, row.names = FALSE
    )
