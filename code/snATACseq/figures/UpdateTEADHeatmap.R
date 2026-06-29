#! /usr/bin/Rscript

########## UPDATE FIMO RESULTS TEAD MOTIFS - ONLY SHOW RELEVANT PEAKS ##########

# LOAD LIBRARY ----
library(here)
library(tidyverse)

# LOAD DATA ----
path.fimo.results <- here("output/analysis/D10029_yap1_snATAC_dataset/MotifAnalysis/fimo.tsv")
fimo.results <- read_tsv(path.fimo.results)

# SET VARIABLES ----
order.meox1.peaks.short.tead <- c('12:27425914-27426414', #-36meox1
                                  '12:27429569-27430069') #-32meox1
save_dir <- here("output/figures")
name <- "yap1_ATAC_dataset"

# GET FIMO COUNTS ----
fimo.count.results <- fimo.results %>%
  filter(., `p-value` < 0.001) %>%
  dplyr::group_by(., sequence_name) %>%
  dplyr::count(., motif_alt_id)

fimo.count.results.selected <- fimo.count.results %>%
  filter(., sequence_name %in% order.meox1.peaks.short.tead)

fimo.count.results.selected$pretty_name <- ifelse(fimo.count.results.selected$sequence_name == "12:27425914-27426414", "-36meox1", "-32meox1")

# MAKE PLOT ----

heatmap.delta <- ggplot(fimo.count.results.selected, aes(x=motif_alt_id, y=pretty_name, fill=n)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Blues", direction=1, limits = c(1, 4)) +
  labs(title = "TEAD motifs", fill = '# TEAD motif match\npval<0.001') +
  xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

FileName <- sprintf('%s/%s_Level02_meox1_32_36peaks_TEADmotif_hits_heatmap.pdf', save_dir, name)
ggsave(FileName,
       heatmap.delta,
       device="pdf",
       width=15,
       height=8,
       units="cm",
       scale=1)


