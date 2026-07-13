#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Depends on R/4.0.2


################################################################################
# LOAD DATA
################################################################################
#Load variables
source('.Rev_EMBO_meox1_expression_variables.R')


################################################################################
# EMBO UMAP LEVEL 03
################################################################################
# ----- L3 phenotype ----
EMBO_Level_03_seuratObject_relevel <- EMBO_Level_03_seuratObject
EMBO_Level_03_seuratObject_relevel 		<- SetIdent(EMBO_Level_03_seuratObject_relevel, value="L2_Predicted_phenotype")
EMBO_Level_03_seuratObject_relevel$L2_Predicted_phenotype <- factor(EMBO_Level_03_seuratObject_relevel$L2_Predicted_phenotype, levels =c('LEC', 'LEC_prox1a_low','muLEC', 'VEC', 'VEC_preLEC') )
embo_colours <- c('#a1d99b',
                  '#74c476',
                  '#fe9929',
                  '#6baed6',
                  '#2171b5')
names(embo_colours) <- c('LEC', 'LEC_prox1a_low','muLEC', 'VEC', 'VEC_preLEC' )
umap_pheno <- plot_umap(seurat.object = EMBO_Level_03_seuratObject_relevel,
                        group = "L2_Predicted_phenotype",
                        save.dir = save.dir,
                        cols = embo_colours,
                        project.name = project.name.level.03_embo)


# ----- L3 phenotype ----
embo_stage_colours <- c('#e4c6d3','#be7794', '#a95175', '#823e59')
names(embo_colours) <- c( 'Stage_40hpf', 'Stage_3dpf', 'Stage_4dpf', 'Stage_5dpf')
umap_pheno <- plot_umap(seurat.object = EMBO_Level_03_seuratObject_relevel,
                        group = "Stage",
                        save.dir = save.dir,
                        cols = embo_stage_colours,
                        project.name = project.name.level.03_embo)


# ----- GEX meox1 one for each  stage ----
all_stages <- lapply(EMBO_Level_03_seuratObject_relevel$Stage %>% unique(), function(stage){
  EMBO_Level_03_seuratObject_relevel <- SetIdent(EMBO_Level_03_seuratObject_relevel, value = "Stage")
  tmp_seurat <- subset(EMBO_Level_03_seuratObject_relevel, idents = stage)
  plot_umap_feature(seurat.object = tmp_seurat,
                    feature = 'meox1',
                    save.dir = save.dir,
                    cols = expression_colours,
                    project.name = sprintf('%s_%s', project.name.level.03_embo, stage))
})
#
# umap <- FeaturePlot(EMBO_Level_03_seuratObject_relevel,
#                     features = 'meox1',
#                     cols = expression_colours, split.by = 'Stage')

# ----- GEX meox1 ----
EMBO_meox1_umap <- plot_umap_feature(seurat.object = EMBO_Level_03_seuratObject,
                  feature = 'meox1',
                  save.dir = save.dir,
                  cols = expression_colours,
                  project.name = project.name.level.03_embo)


