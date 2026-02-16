#! /usr/bin/Rscript

################################################################################
#     COLLECTION OF PAPER SETTINGS                                             #
################################################################################

# COLOURS - GENERAL ----
genotype_colours <- c(
  "#dbe2c6", #wildtype
  "#657c95" #mutant
)

# names(genotype_colours) <- c("WT", "yap1_mutant") -> set name depending on dataset

expression_colours <- c("#d9d9d9", "#40004b")

cols_cellcycle <- c( '#ffff66', '#cc85ff','#a0d0e0')
names(cols_cellcycle) <- c("G1/G0", "S", "G2M")



# FUNCTIONS ----
make_dim_umap <- function(seurat,
                          group,
                          width,
                          height,
                          path,
                          name,
                          point_size = 0.5,
                          colours = NULL){
  plot <- DimPlot(seurat,
                  group.by = group,
                  shuffle = T,
                  cols = colours,
                  pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", group, ".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot + NoAxes() + NoLegend() + theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", group, ".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_dim_split_umap <- function(seurat,
                          group,
                          width,
                          height,
                          path,
                          name,
                          split,
                          point_size = 0.5,
                          colours = NULL){
  plot <- DimPlot(seurat,
                  group.by = group,
                  shuffle = T,
                  cols = colours,
                  split.by = split,
                  pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", group, "_split_", split ,".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot & NoAxes() & NoLegend() & theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", group, "_split_", split ,".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_feature_umap <- function(seurat,
                          feature,
                          width,
                          height,
                          path,
                          name,
                          point_size = 0.5,
                          colours = expression_colours){
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, ".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot + NoAxes() + NoLegend() + theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, ".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_feature_split_umap <- function(seurat,
                              feature,
                              width,
                              height,
                              path,
                              name,
                              split,
                              point_size = 0.5,
                              colours = expression_colours){
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      split.by = split,
                      pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, "_split_", split ,".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot & NoAxes() & NoLegend() & theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, "_split_", split ,".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_genescore_umap <- function(seurat,
                              feature,
                              width,
                              height,
                              path,
                              name,
                              point_size = 0.5,
                              colours = expression_colours){
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      pt.size = point_size,
                      min.cutoff = "q1",
                      max.cutoff =  "q99",
                      order = T)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, ".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot + NoAxes() + NoLegend() + theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, ".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_genescore_split_umap <- function(seurat,
                                feature,
                                width,
                                height,
                                path,
                                name,
                                point_size = 0.5,
                                split ,
                                colours = expression_colours){
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      split.by = split,
                      pt.size = point_size,
                      min.cutoff = "q1",
                      max.cutoff =  "q99",
                      order = T)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, "_split_", split ,".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot & NoAxes() & NoLegend() & theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature,"_split_", split ,".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_ratio_dotplots <- function(deg_input,
                                width,
                                height,
                                path,
                                name){
  dotplot <- deg_input %>%
    mutate(.,`percent ratio` = pct.2/pct.1,
           `log2FC` = -avg_log2FC) %>%
    ggplot(.) +
    geom_point(aes(x= group, y = Gene, colour = log2FC, size = `percent ratio`)) +
    scale_colour_gradient2(low ='#332288',
                         high ='#CC6677',
                         na.value = "black",
                         guide = "colourbar",
                         midpoint = 0,
                         aesthetics = "colour") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  + xlab('') + ylab('') +
    ggtitle('wt/mutant-/-\n% expressed and log2FC')

  filename <- sprintf('%s/%s_ratio_dotplot.pdf', save_dir, name)
  ggsave(filename,
         dotplot,
         device="pdf",
         width=width,
         height=height)

  return(dotplot)
}

