#!/usr/bin/Rscript
# regenerate fig 7d split by genotype EMBO LEC VEC score plots split by genotype
# ../../Saki_data/Revision_analysis/queue_plots/
#   queue_plots/embo_transition_scores_split_genotype_dotplot.pdf
#   queue_plots/embo_transition_scores_mut_LEC_strip.pdf
#   queue_plots/embo_transition_scores_mut_LEC.pdf
#   queue_plots/embo_transition_scores_mut_VEC_strip.pdf
#   queue_plots/embo_transition_scores_mut_VEC.pdf
#   queue_plots/embo_transition_scores_split_genotype_strip.pdf
#   queue_plots/embo_transition_scores_split_genotype.pdf
#   queue_plots/embo_transition_scores_wt_LEC_strip.pdf
#   queue_plots/embo_transition_scores_wt_LEC.pdf
#   queue_plots/embo_transition_scores_wt_VEC_strip.pdf
#   queue_plots/embo_transition_scores_wt_VEC.pdf
library(qs2)
library(Seurat)
library(SeuratObject)

infile_path <- "../paper/D10051_meox1_dataset_Level_03_annotated.qs2"
data <- qs_read(infile_path)

data$Genotype <- gsub("meox1_mutant", "meox1 mutant", data$Genotype)
data$Genotype <- gsub("wildtype", "WT", data$Genotype)
data$Genotype <- factor(
    data$Genotype, 
    levels = c("WT", "meox1 mutant")
)

midpoint <- mean(range(data@meta.data$EMBO_transient_score))

# combine genotype
data$celltype_genotype <- paste(data$L3_celltype, data$Genotype, sep = "___")

# calculate values only
base_plt <- DotPlot(data, features = "EMBO_transient_score", group.by = "celltype_genotype")
df <- base_plt$data

# reconstruct celltype genotype
df_split <- do.call(rbind, strsplit(as.character(df$id), "___"))
df$Celltype <- df_split[, 1]
df$Genotype <- df_split[, 2]

# retain order
cell_order <- c("LEC", "pre_muLEC", "hmVEC", "mVEC")
df$Celltype <- factor(df$Celltype, levels = rev(cell_order))
df$Genotype <- factor(df$Genotype, levels = c("WT", "meox1 mutant"))

# custom ggplot
plt <- ggplot(df, aes(x = Genotype, y = Celltype, size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(
        low = "#406880",
        high = "darkgreen",
        mid = "azure2",
        midpoint = midpoint
    ) +
    labs(
        x = "Genotype", 
        y = "Cell Type", 
        color = "VEC-to-LEC\nscore", 
        size = "% Cells\nExpressing"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line(color = "grey90")
    ) +
    theme_classic()

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype_dotplot.pdf",
    device = "pdf",
    height = 7,
    width = 5
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype_dotplot.png",
    device = "png",
    height = 7,
    width = 5
    )


wt <- data[,data@meta.data$Genotype == "WT"]
mut <- data[,data@meta.data$Genotype == "meox1 mutant"]

md <- data[[]]
coords <- Embeddings(data[["umap"]])
md <- cbind(md, coords)

midpoint <- mean(range(md$EMBO_transient_score))

# Order by absolute difference from the midpoint
md <- md[order(abs(md$EMBO_transient_score - midpoint)), ]

plt <- ggplot(md, aes(x = umap_1, y = umap_2, color = EMBO_transient_score)) +
    geom_point(size = 1) +
    scale_color_gradient2(low = "#406880", high = "darkgreen", mid = "azure2", midpoint = midpoint) +
    facet_wrap(~ Genotype) +
    theme_classic() +
    ggtitle("EMBO LEC and VEC transition scores split by genotype")

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype.pdf",
    device = "pdf",
    height = 5,
    width = 12
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype.png",
    device = "png",
    height = 5,
    width = 12
    )

plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype_strip.pdf",
    device = "pdf",
    height = 5,
    width = 10
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_split_genotype_strip.png",
    device = "png",
    height = 5,
    width = 10
    )

# split
scores <- c("EMBO_LEC_score_1", "EMBO_VEC_score_1")
lec_cols <- c("#d9d9d9", "darkgreen")
vec_cols <- c("#d9d9d9", "#406880")

# lec scores
## wt
plt <- FeaturePlot(wt, features = "EMBO_LEC_score_1", cols = lec_cols,
    pt.size = 1, min.cutoff = "q1", max.cutoff = "q99",
    order = TRUE) + ggtitle("wildtype EMBO LEC score")

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_LEC.pdf",
    device = "pdf",
    height = 7,
    width = 7
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_LEC.png",
    device = "png",
    height = 7,
    width = 7
    )

plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_LEC_strip.pdf",
    device = "pdf",
    height = 5,
    width = 5
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_LEC_strip.png",
    device = "png",
    height = 5,
    width = 5
    )

## mut
plt <- FeaturePlot(mut, features = "EMBO_LEC_score_1", cols = lec_cols,
    pt.size = 1, min.cutoff = "q1", max.cutoff = "q99",
    order = TRUE) + ggtitle("mutant EMBO LEC score")

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_LEC.pdf",
    device = "pdf",
    height = 7,
    width = 7
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_LEC.png",
    device = "png",
    height = 7,
    width = 7
    )

plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_LEC_strip.pdf",
    device = "pdf",
    height = 5,
    width = 5
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_LEC_strip.png",
    device = "png",
    height = 5,
    width = 5
    )

# vec scores
## wt
plt <- FeaturePlot(wt, features = "EMBO_VEC_score_1", cols = vec_cols,
    pt.size = 1, min.cutoff = "q1", max.cutoff = "q99",
    order = TRUE) + ggtitle("wildtype EMBO VEC score")

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_VEC.pdf",
    device = "pdf",
    height = 7,
    width = 7
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_VEC.png",
    device = "png",
    height = 7,
    width = 7
    )

plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_VEC_strip.pdf",
    device = "pdf",
    height = 5,
    width = 5
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_wt_VEC_strip.png",
    device = "png",
    height = 5,
    width = 5
    )    

## mut
plt <- FeaturePlot(mut, features = "EMBO_VEC_score_1", cols = vec_cols,
    pt.size = 1, min.cutoff = "q1", max.cutoff = "q99",
    order = TRUE) + ggtitle("mutant EMBO VEC score")

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_VEC.pdf",
    device = "pdf",
    height = 7,
    width = 7
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_VEC.png",
    device = "png",
    height = 7,
    width = 7
    )

plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())

ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_VEC_strip.pdf",
    device = "pdf",
    height = 5,
    width = 5
    )
ggsave(
    plot = plt,
    filename = "queue_plots/embo_scores/embo_transition_scores_mut_VEC_strip.png",
    device = "png",
    height = 5,
    width = 5
    )