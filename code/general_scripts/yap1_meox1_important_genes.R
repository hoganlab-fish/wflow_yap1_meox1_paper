#! /usr/bin/Rscript

################################################################################
#     COLLECTION OF GENE SETS AND LISTS USED IN THE MANUSCRIPTS                #
################################################################################

## ANNOTATION GENES ----
# using the same set of genes for annotation in both yap1 and meox1 dataset

level_01_annotation_genes <- c(
  "kdrl",
  "cdh6",
  "prox1a",
  "stab2",
  "lyve1b",
  "mrc1a",
  "flt1",
  "dll4",
  "esm1",
  "hlx1",
  "hand2",
  "fn1a",
  "gata2b",
  "runx1",
  "alas2",
  "blvrb",
  "spi1b",
  "mpx",
  "neurod4",
  "nova2",
  "cdh1",
  "epcam",
  "dlx5a",
  "dlx6a",
  "lum",
  "col1a1a",
  "cryba1b",
  "crybb1",
  "pmela",
  "dct",
  "ednrba",
  "pnp4a",
  "gch2",
  "aox5")


level_02_annotation_genes <- c('prox1a', #LEC
                               'cdh6', #VEC/LEC
                               'stab2', #LEC/VEC
                               'mrc1a', #LEC/VEC
                               'lyve1b', #VEC/VEC
                               'flt1',#AEC, mAEC
                               'dll4', #AEC, mAEC
                               'hlx1', #mAEC
                               'esm1', #mAEC
                               'pcna', #cycling
                               'mki67', #cycling
                               'hand2',#Endocardium
                               'fn1a' #Endocardium
)


# level_03_annotation_genes <- c(
#   "stab2", #cvpVEC
#   "stab1", #cvpVEC
#   "ifi30",#cvpVEC
#   "mafba", #cvpVEC
#   "cxcl12a", #msVEC
#   "esm1",  #prVEC ?
#   "hlx1", #prVEC ?
#   "mki67", #prVEC ?
#   "pcna", #prVEC ?
#   "il13ra1",
#   "il4r.1",
#   "nop53", #LEC
#   "ccn2a", #hmsVEC
#   "cdh6",#LEC +hmsVEC
#   "tbx1", #LEC
#   "prox1a" #LEC
#
#   )

level_03_annotation_genes <- c(
  "ifi30",#cvpVEC
  "mafba", #cvpVEC
  "cxcl12a", #msVEC
  "esm1",  #prVEC ?
  "hlx1", #prVEC ?
  "mki67", #prVEC ?
  "pcna", #prVEC ?
  "il13ra1",
  "il4r.1",
  "ccn2a", #hmsVEC
  "cdh6",#LEC +hmsVEC
  "nop53", #LEC
  "tbx1", #LEC
  "prox1a" #LEC

)

## GENE SETS FOR MODULE SCORES----
hippo_targets <- c(
  "yap1",
  "meox1",
  "wwtr1",
  "ccn2a",
  "ccn2b",
  "amotl2a",
  "amotl2b",
  "ccn1",
  "ccn1l2",
  "cavin1b",
  "cavin2a",
  "cavin2b",
  "bmp4")

hippo_targets_nomeox1 <- c(
  "yap1",
  "wwtr1",
  "ccn2a",
  "ccn2b",
  "amotl2a",
  "amotl2b",
  "ccn1",
  "ccn1l2",
  "cavin1b",
  "cavin2a",
  "cavin2b",
  "bmp4")

ap1_targets <- c(
  'fosb',
  'fosaa',
  'fosab',
  'fosl2',
  'fosl1a',
  'fosl1b',
  'jun',
  'jund',
  'junba',
  'junbb',
  'june')

## TRANSITION SCORE GENES (EMBO) ----

mapk_targets <- read.table(here("data/genesets/MAPK_genes.csv"), header = T)
wnt_targets <- read.table(here("data/genesets/Wnt_genes.csv"), header = T)

embo_paper_markers <- readRDS(here("data/genesets/EMBO_L2_VECLEC_markers_3_5dpf.RDS"))
lec_score <- embo_paper_markers$LEC %>%
  filter(., p_val_adj < 0.05) %>%
  top_n(., n = 25, wt = avg_log2FC) %>%
  rownames_to_column("Gene") %>% pull(Gene)

vec_score <- embo_paper_markers$VEC %>%
  filter(., p_val_adj < 0.05) %>%
  top_n(., n = 25, wt = avg_log2FC) %>%
  rownames_to_column("Gene") %>% pull(Gene)

## TRANSITION SCORE GENES (MEOX1, LEC) ----
# EMBO score is 3dpf-5dpf- we need something a bit younger, so try use the meox1 dataset instead
meox1_markers <- read_csv(here("output/analysis/D10051_meox1_dataset/Markers/Level_03_Markers_L3_celltype_fc0.00_minpct_0.01.csv"), show_col_types = F)

lec_dev_score <- meox1_markers %>%
  filter(., p_val_adj < 0.05 &  cluster == "LEC") %>%
  top_n(., n = 25, wt = avg_log2FC) %>%
  rownames_to_column("Gene") %>% pull(gene)




