#!/usr/bin/Rscript
# create template axis
library(grid)

outfile_path <- "../Saki/Revision_analysis/queue_plots/mini_axis.pdf"

pdf(outfile_path, width=2, height=2)
grid.newpage()
grid.lines(c(0.1, 0.5), c(0.1, 0.1), 
           arrow=arrow(type="closed", length=unit(0.05, "inches"), angle=20), 
           gp=gpar(lwd=1, fill="black"))
grid.lines(c(0.1, 0.1), c(0.1, 0.5), 
           arrow=arrow(type="closed", length=unit(0.05, "inches"), angle=20), 
           gp=gpar(lwd=1, fill="black"))
grid.text("UMAP_1", x=0.3, y=0.05, gp=gpar(cex=0.8))
grid.text("UMAP_2", x=0.05, y=0.3, rot=90, gp=gpar(cex=0.8))
dev.off()