#!/usr/bin/Rscript
# slingshot meox1 vs control
# - versions tagged Strip intentionally have no legend/title/axes for use
# - versions tagged Label have legend/title/axes for reference
# - colour scheme, point size and order matches that of existing
# - note that umap coordinates used in all cases are the dataset (ie, umaps were not replotted per sample)
# D10051 - level 3 meox dataset
# S20200 - meox1 mutant
# S20201 - wildtype sample
# start_{celltype} - reflects the "start" point of the algorithm (which should reflect the most "primitive" celltype cluster). i plotted all just for reference (and partly as a sanity check).

# hmVEC
# /Revision_analysis/queue_plots/
    # meox1__FeaLabel_Pseudotime_Contours_start_hmVEC_lin1.D10051.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_hmVEC_lin1.S20200.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_hmVEC_lin1.S20201.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_hmVEC_lin1.D10051.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_hmVEC_lin1.S20200.pdf 
    # meox1__FeaStrip_Pseudotime_Contours_start_hmVEC_lin1.S20201.pdf

# mVEC
# /Revision_analysis/queue_plots/
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin1.D10051.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin1.S20200.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin1.S20201.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin1.D10051.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin1.S20200.pdf 
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin1.S20201.pdf
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin2.D10051.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin2.S20200.pdf 
    # meox1__FeaLabel_Pseudotime_Contours_start_mVEC_lin2.S20201.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin2.D10051.pdf
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin2.S20200.pdf 
    # meox1__FeaStrip_Pseudotime_Contours_start_mVEC_lin2.S20201.pdf    
    # meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_mVEC_lin1.D10051.pdf
    # meox1__BarStrip_Isolated_Zones_Stacked_Proportion_start_mVEC_lin1.D10051.pdf
    # meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_mVEC_lin2.D10051.pdf
    # meox1__BarStrip_Isolated_Zones_Stacked_Proportion_start_mVEC_lin2.D10051.pdf
    # meox1__BarLabel_Custom_mVEC_Lin1_Zones12_Stacked_Proportion.pdf
    # meox1__BarStrip_Custom_mVEC_Lin1_Zones12_Stacked_Proportion.pdf
    # meox1__BarLabel_Custom_mVEC_Lin1_Zones12_Dodged_Proportion.pdf
    # meox1__BarStrip_Custom_mVEC_Lin1_Zones12_Dodged_Proportion.pdf
    # meox1__BarLabel_Custom_mVEC_Lin2_Zones12_Stacked_Proportion.pdf
    # meox1__BarStrip_Custom_mVEC_Lin2_Zones12_Stacked_Proportion.pdf
    # meox1__BarLabel_Custom_mVEC_Lin2_Zones12_Dodged_Proportion.pdf
    # meox1__BarStrip_Custom_mVEC_Lin2_Zones12_Dodged_Proportion.pdf
    # meox1__DensLabel_Minima_start_mVEC_lin1.D10051.pdf 
    # meox1__DensLabel_Minima_start_mVEC_lin1.S20200.pdf 
    # meox1__DensLabel_Minima_start_mVEC_lin1.S20201.pdf 
    # meox1__DensLabel_Minima_start_mVEC_lin2.D10051.pdf 
    # meox1__DensLabel_Minima_start_mVEC_lin2.S20200.pdf 
    # meox1__DensLabel_Minima_start_mVEC_lin2.S20201.pdf

# /Revision_analysis/queue_plots/
    # meox1__BarStrip_Combined_Total_and_Cluster_Proportion.D10051.pdf 
    # meox1__BarStrip_Combined_Total_and_Cluster_Absolute.D10051.pdf 
    # meox1__BarLabel_Combined_Total_and_Cluster_Absolute.D10051.pdf 
    # meox1__BarLabel_Combined_Total_and_Cluster_Proportion.D10051.pdf
    # meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_hmVEC_lin1.D10051.pdf
    # meox1__BarStrip_Isolated_Zones_Stacked_Proportion_start_hmVEC_lin1.D10051.pdf

# /Revision_analysis/queue_plots/
    # meox1__BarLabel_Custom_hmVEC_Lin1_Zones123_Stacked_Proportion.pdf
    # meox1__BarStrip_Custom_hmVEC_Lin1_Zones123_Stacked_Proportion.pdf
    # meox1__BarLabel_Custom_hmVEC_Lin1_Zones123_Dodged_Proportion.pdf
    # meox1__BarStrip_Custom_hmVEC_Lin1_Zones123_Dodged_Proportion.pdf
    # meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_hmVEC_lin1.S20200.pdf 
    # meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_hmVEC_lin1.S20201.pdf

# /Revision_analysis/queue_plots/
    # meox1__DensLabel_Minima_start_hmVEC_lin1.D10051.pdf 
    # meox1__DensLabel_Minima_start_hmVEC_lin1.S20200.pdf 
    # meox1__DensLabel_Minima_start_hmVEC_lin1.S20201.pdf

# this is the set of plots
# /Revision_analysis/queue_plots/
    # meox1__DensLabel_WT_Anchored_Overlay_Direct_HeadToHead_start_hmVEC_lin1.pdf
    # meox1__BarLabel_WT_Anchored_Quantified_Proportions_start_hmVEC_lin1.pdf

# /Revision_analysis/queue_plots/
    # meox1__DensLabel_WT_Anchored_Overlay_Direct_HeadToHead_start_mVEC_lin1.pdf
    # meox1__BarLabel_WT_Anchored_Quantified_Proportions_start_mVEC_lin1.pdf
    # meox1__DensLabel_WT_Anchored_Overlay_Direct_HeadToHead_start_mVEC_lin2.pdf
    # meox1__BarLabel_WT_Anchored_Quantified_Proportions_start_mVEC_lin2.pdf

# /Revision_analysis/queue_plots/
#   meox1__DimLabel_WT_Anchored_Contours_Combined_start_hmVEC_lin1.pdf
#   meox1__DimLabel_WT_Anchored_Contours_Split_start_hmVEC_lin1.pdf

# /Revision_analysis/queue_plots/
#   meox1__DimLabel_WT_Anchored_Contours_Combined_start_mVEC_lin1.pdf
#   meox1__DimLabel_WT_Anchored_Contours_Split_start_mVEC_lin1.pdf
#   meox1__DimLabel_WT_Anchored_Contours_Combined_start_mVEC_lin2.pdf
#   meox1__DimLabel_WT_Anchored_Contours_Split_start_mVEC_lin2.pdf


# note that
# S20200 is mutant and
# S20201 is wildtype

libraries <- c(
    "data.table",
    "DelayedMatrixStats",
    "dplyr",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "qs2",
    "Seurat",
    "SingleCellExperiment",
    "slingshot",
    "tidyverse"
)

load_packages <- function(packages) {
    suppressPackageStartupMessages({
        cat(paste("Loading", packages, sep=" "), sep="\n")
        invisible(sapply(
            packages, library, character.only = TRUE, quietly = TRUE
            ))
    })
}
invisible(Sys.setenv(OPENBLAS_NUM_THREADS="1"))
load_packages(libraries)

# custom transition zone barplots
plot_custom_transition_zones <- function(
    master_df,         
    target_zones,      # e.g., c(1, 2, 3)
    custom_colors,     # named vector mapping colors to "Transition Zone X"
    title_prefix,      # e.g., "hmVEC Lin 1"
    file_prefix,       # e.g., "Custom_hmVEC_Lin1"
    outfile_dir,
    col_name           # automatically tracks D10051 vs S20200 vs S20201
) {
   
    if (nrow(master_df) == 0) { return() }
    zone_names <- paste0("Transition Zone ", target_zones)
    custom_df <- master_df %>% filter(Zone %in% zone_names)
    
    if (nrow(custom_df) == 0) { return() }

    cat(paste0("Generating custom zone plots for: ", title_prefix, " (Dataset: ", col_name, ")...\n"))
        
    # automatically use whatever samples are present in master_df for this run
    custom_df <- master_df %>%
        filter(Zone %in% zone_names) %>%
        group_by(Sample, Zone) %>% 
        tally(name = "Cell_Count") %>% 
        ungroup()
    
    if (nrow(custom_df) == 0) {
        cat("   -> No cells found in the targeted zones. Skipping plot.\n")
        return()
    }

    # recalculate percentages relative only to this specific isolated subset per sample
    custom_df <- custom_df %>%
        group_by(Sample) %>%
        mutate(Sub_Total = sum(Cell_Count), Proportion = (Cell_Count / Sub_Total) * 100) %>% 
        ungroup()
        
    custom_df$Zone <- factor(custom_df$Zone, levels = zone_names)
    
    custom_df <- custom_df %>%
        arrange(Sample, desc(Zone)) %>%
        group_by(Sample) %>% 
        mutate(label_y_prop = cumsum(Proportion) - (0.5 * Proportion)) %>% 
        ungroup()

    # stacked plot
    plt_stacked <- ggplot(custom_df, aes(x = Sample, y = Proportion, fill = Zone)) +
        geom_col(width = 0.5, color = "black", linewidth = 0.5) +
        geom_text_repel(aes(y = label_y_prop, label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), "")), 
                        nudge_x = 0.4, direction = "y", hjust = 0, segment.size = 0.3) +
        scale_fill_manual(values = custom_colors) +
        coord_cartesian(ylim = c(0, 100)) +
        labs(title = paste0(title_prefix, " - Zones ", paste(target_zones, collapse="-"), " %"), 
             x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
        theme_classic() + 
        theme(axis.text = element_text(color = "black", size = 11), 
              axis.title.y = element_text(face = "bold", size = 12), 
              plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

    ggsave(paste0(outfile_dir, "meox1__BarLabel_", file_prefix, "_Zones", paste(target_zones, collapse=""), "_Stacked_Proportion.", col_name, ".pdf"), plot = plt_stacked, height = 5, width = 6.5)
    ggsave(paste0(outfile_dir, "meox1__BarStrip_", file_prefix, "_Zones", paste(target_zones, collapse=""), "_Stacked_Proportion.", col_name, ".pdf"), plot = plt_stacked + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

    # dodged plot
    max_prop <- max(custom_df$Proportion, na.rm = TRUE)
    
    plt_dodged <- ggplot(custom_df, aes(x = Sample, y = Proportion, fill = Zone)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
        geom_text(aes(label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), ""), y = Proportion + (max_prop * 0.02)), 
                  position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
        scale_fill_manual(values = custom_colors) +
        coord_cartesian(ylim = c(0, min(100, max_prop * 1.15))) +
        labs(title = paste0(title_prefix, " - Zones ", paste(target_zones, collapse="-"), " Comparison"), 
             x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
        theme_classic() + 
        theme(axis.text = element_text(color = "black", size = 11), 
              axis.title.y = element_text(face = "bold", size = 12), 
              plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

    ggsave(paste0(outfile_dir, "meox1__BarLabel_", file_prefix, "_Zones", paste(target_zones, collapse=""), "_Dodged_Proportion.", col_name, ".pdf"), plot = plt_dodged, height = 5, width = 7)
    ggsave(paste0(outfile_dir, "meox1__BarStrip_", file_prefix, "_Zones", paste(target_zones, collapse=""), "_Dodged_Proportion.", col_name, ".pdf"), plot = plt_dodged + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)
}

# accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")

run_wt_anchored_slingshot <- function(data, col_map, col_order, outfile_dir) {
    sce <- as.SingleCellExperiment(data)
    celltype <- as.vector(unique(data@meta.data$L3_celltype))
    data$L3_celltype <- factor(data$L3_celltype, levels=col_order)
    col_name <- "WT_Anchored"

    for (start in celltype) {
        sce_run <- slingshot(
            sce, 
            clusterLabels = "L3_celltype", 
            reducedDim = "UMAP",           
            start.clus = start             
        )

        pseudotime_matrix <- slingPseudotime(sce_run)
        n_lineages <- ncol(pseudotime_matrix)

        for (lin in seq_len(n_lineages)) {
            
            meta <- paste0("Slingshot_Pseudotime_start_", start, "_lin", lin)
            pt_values <- pseudotime_matrix[, lin]
            
            cat("Processing WT-Anchored lineage:", meta, "\n")

            # 1. BUILD MASTER DATAFRAME FOR ALL CELLS
            master_zones_df <- data.frame(
                Barcode = colnames(data), 
                Pseudotime = pt_values,
                Cluster = data@meta.data$L3_celltype
            ) %>%
                mutate(Sample = case_when(
                    grepl("_1$", Barcode) ~ "meox1 -/- mutant", 
                    grepl("_2$", Barcode) ~ "WildType", 
                    TRUE ~ "Unknown"
                )) %>%
                filter(Sample != "Unknown")

            # Split profiles
            wt_pseudotime <- master_zones_df %>% filter(Sample == "WildType", !is.na(Pseudotime)) %>% pull(Pseudotime)
            mut_pseudotime <- master_zones_df %>% filter(Sample == "meox1 -/- mutant", !is.na(Pseudotime)) %>% pull(Pseudotime)
            full_pseudotime <- master_zones_df %>% filter(!is.na(Pseudotime)) %>% pull(Pseudotime)
            
            if(length(wt_pseudotime) < 10) {
                cat("   -> Insufficient WT cells for density calculations. Skipping.\n")
                next
            }

            # 2. CALCULATE ALL DENSITIES & EXTRACT VALLEYS
            dens_wt <- density(wt_pseudotime, na.rm = TRUE)
            valleys <- which(diff(sign(diff(dens_wt$y))) == 2) + 1
            valley_times <- dens_wt$x[valleys]
            
            dens_mut  <- density(mut_pseudotime, na.rm = TRUE)
            valleys_mut <- which(diff(sign(diff(dens_mut$y))) == 2) + 1
            valley_times_mut <- dens_mut$x[valleys_mut]

            dens_full <- density(full_pseudotime, na.rm = TRUE)

            accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")
            
            # Setup base dataframes for plotting
            df_wt   <- data.frame(X = dens_wt$x, Y = dens_wt$y)
            df_mut  <- data.frame(X = dens_mut$x, Y = dens_mut$y)
            df_full <- data.frame(X = dens_full$x, Y = dens_full$y)

            # --- PLOT 1: CLEAN WT ONLY (NO OVERLAYS) ---
# Setup base dataframes for plotting
            df_wt   <- data.frame(X = dens_wt$x, Y = dens_wt$y)
            df_mut  <- data.frame(X = dens_mut$x, Y = dens_mut$y)
            df_full <- data.frame(X = dens_full$x, Y = dens_full$y)

            # Create mapping dataframes for the zones to force ggplot to build a legend
            zones_wt_df <- data.frame(
                Zone = paste0("Transition Zone ", seq_along(valley_times)),
                xmin = valley_times - 0.5,
                xmax = valley_times + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
            zone_colors_wt <- setNames(accent_palette[((seq_along(valley_times) - 1) %% length(accent_palette)) + 1], zones_wt_df$Zone)

            # --- PLOT 1: CLEAN WT ONLY (NO OVERLAYS) ---
            p1 <- ggplot(df_wt, aes(x = X, y = Y)) +
                geom_rect(data = zones_wt_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15) +
                geom_line(color = "#377EB8", linewidth = 1) + 
                geom_vline(xintercept = valley_times, color = "red", linetype = "dashed") +
                scale_fill_manual(name = "Trajectory Regions", values = zone_colors_wt) +
                labs(title = paste0("WT-Only Baseline Profile\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Relative Probability Density") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Pure_WT_start_", start, "_lin", lin, ".pdf"), plot = p1, height = 5, width = 6.5)

            # --- PLOT 2: CLEAN MUTANT ONLY (NO OVERLAYS) ---
            p2 <- ggplot(df_mut, aes(x = X, y = Y)) +
                geom_rect(data = zones_wt_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15) +
                geom_line(color = "#E41A1C", linewidth = 1) + 
                geom_vline(xintercept = valley_times, color = "red", linetype = "dashed") +
                scale_fill_manual(name = "Trajectory Regions", values = zone_colors_wt) +
                labs(title = paste0("Mutant-Only Profile mapped to WT Roadmarks\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Relative Probability Density") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Pure_Mut_start_", start, "_lin", lin, ".pdf"), plot = p2, height = 5, width = 6.5)

            # --- PLOT 3: OVERLAY - COMBINED FULL DATASET VS WT BASELINE ---
            p3 <- ggplot(df_full, aes(x = X, y = Y)) +
                geom_rect(data = zones_wt_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15) +
                geom_line(aes(color = "Combined Dataset"), linewidth = 1) + 
                geom_line(data = df_wt, aes(x = X, y = Y, color = "WT Baseline"), linetype = "longdash", linewidth = 0.8) +
                geom_vline(xintercept = valley_times, color = "red") +
                scale_color_manual(name = "Subpopulation", values = c("Combined Dataset" = "black", "WT Baseline" = "#377EB8")) +
                scale_fill_manual(name = "Trajectory Regions", values = zone_colors_wt) +
                labs(title = paste0("Combined Pool vs WT Baseline Validation\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Overlay_Full_start_", start, "_lin", lin, ".pdf"), plot = p3, height = 5, width = 7.5)

            # --- PLOT 4: OVERLAY - DIRECT COMPARISON (FACETED SUBPLOTS) ---
            df_compare <- bind_rows(
                df_wt %>% mutate(Genotype = "WildType Baseline"),
                df_mut %>% mutate(Genotype = "meox1 -/- Mutant")
            )
            df_compare$Genotype <- factor(df_compare$Genotype, levels = c("WildType Baseline", "meox1 -/- Mutant"))

            p4 <- ggplot(df_compare, aes(x = X, y = Y, color = Genotype)) +
                geom_rect(data = zones_wt_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15, color = NA) +
                geom_line(linewidth = 1) + 
                geom_vline(xintercept = valley_times, color = "red", linetype = "solid") +
                scale_color_manual(name = "Genotype", values = c("WildType Baseline" = "#377EB8", "meox1 -/- Mutant" = "#E41A1C")) +
                scale_fill_manual(name = "Trajectory Regions", values = zone_colors_wt) +
                facet_wrap(~ Genotype, ncol = 1) + 
                labs(title = paste0("Faceted Head-to-Head Density Alignment\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density") + 
                theme_classic() + 
                theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                      strip.background = element_rect(fill = "grey90", color = "black"),
                      strip.text = element_text(face = "bold", size = 11))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Overlay_Direct_HeadToHead_start_", start, "_lin", lin, ".pdf"), plot = p4, height = 6, width = 7.5)

            # --- PLOT 5: FULL DENSITY OVERLAY WITH WT ZONES ONLY ---
            p5 <- ggplot(df_full, aes(x = X, y = Y)) +
                geom_rect(data = zones_wt_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15) +
                geom_line(linewidth = 1, color = "black") + 
                geom_vline(xintercept = valley_times, color = "#377EB8", linetype = "dashed", linewidth = 1) +
                scale_fill_manual(name = "Trajectory Regions", values = zone_colors_wt) +
                labs(title = paste0("Full Dataset Density with WT Zones\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Full_with_WT_Zones_start_", start, "_lin", lin, ".pdf"), plot = p5, height = 5, width = 6.5)

            # --- PLOT 6: FULL DENSITY OVERLAY WITH MUTANT ZONES ONLY ---
            zones_mut_df <- data.frame(
                Zone = paste0("Mutant Zone ", seq_along(valley_times_mut)),
                xmin = valley_times_mut - 0.5,
                xmax = valley_times_mut + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
            p6 <- ggplot(df_full, aes(x = X, y = Y)) +
                geom_rect(data = zones_mut_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), inherit.aes = FALSE, alpha = 0.15) +
                geom_line(linewidth = 1, color = "black") + 
                geom_vline(xintercept = valley_times_mut, color = "#E41A1C", linetype = "dashed", linewidth = 1) +
                scale_fill_manual(name = "Mutant Peak Regions", values = setNames(rep("#E41A1C", length(valley_times_mut)), zones_mut_df$Zone)) +
                labs(title = paste0("Full Dataset Density with Mutant Zones\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Full_with_Mutant_Zones_start_", start, "_lin", lin, ".pdf"), plot = p6, height = 5, width = 6.5)

            # --- PLOT 7: FULL DENSITY OVERLAY WITH BOTH ZONES ALIGNED ---
            wt_lines <- data.frame(val = valley_times, type = "WT Zone Boundary")
            mut_lines <- data.frame(val = valley_times_mut, type = "Mutant Zone Boundary")
            all_lines <- rbind(wt_lines, mut_lines)

            p7 <- ggplot(df_full, aes(x = X, y = Y)) +
                geom_line(linewidth = 1, color = "black") +
                geom_vline(data = all_lines, aes(xintercept = val, color = type), linetype = "dashed", linewidth = 1) +
                scale_color_manual(values = c("WT Zone Boundary" = "#377EB8", "Mutant Zone Boundary" = "#E41A1C")) +
                labs(title = paste0("Full Density with WT vs Mutant Zone Alignment\n(Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density", color = "Calculated Boundary") + 
                theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave(paste0(outfile_dir, "meox1__DensLabel_WT_Anchored_Full_with_Both_Zones_start_", start, "_lin", lin, ".pdf"), plot = p7, height = 5, width = 7.5)


            # 3. ASSIGN INTENSITY CLASSIFICATIONS BY ZONE BOUNDARIES
            master_zones_df <- master_zones_df %>% 
                mutate(Zone = case_when(
                    is.na(Pseudotime) ~ "Trajectory NA (Unmapped)", 
                    TRUE ~ "Core/Other Cells"
                ))
            
            master_zones_df$Sample <- factor(master_zones_df$Sample, levels = c("meox1 -/- mutant", "WildType"))
            
            # Using WT valley times to assign zones for the final barplots as this is the baseline
            for (i in seq_along(valley_times)) {
                v_time <- valley_times[i]
                is_in_zone <- !is.na(master_zones_df$Pseudotime) & (abs(master_zones_df$Pseudotime - v_time) < 0.5)
                master_zones_df$Zone[is_in_zone] <- paste0("Transition Zone ", i)
            }

            # ========================================================================
            # WT-ANCHORED UMAP CONTOURS (With Contour Density Legends)
            # ========================================================================
            cat("   -> Generating WT-Anchored UMAP Contours...\n")
            
            # 1. Prepare data
            umap_embed <- as.data.frame(Embeddings(data, reduction = "umap"))
            colnames(umap_embed)[1:2] <- c("umap_1", "umap_2")
            umap_embed$Barcode <- rownames(umap_embed)
            umap_df <- master_zones_df %>% left_join(umap_embed, by = "Barcode")
            max_pt <- max(umap_df$Pseudotime, na.rm = TRUE)            
            
            # Ensure Zone levels match the palette exactly
            umap_df$Zone <- factor(umap_df$Zone, levels = names(zone_colors_wt))
            
            # 2. Base UMAP Plot
            p_umap <- ggplot(umap_df %>% arrange(Pseudotime), aes(x = umap_1, y = umap_2)) +
                geom_point(aes(color = Pseudotime), size = 0.8, alpha = 0.6) +
                scale_color_viridis_c(option = "plasma", na.value = "grey90", limits = c(0, max_pt)) +
                theme_classic() +
                labs(title = paste0("WT-Anchored UMAP Contours\n(Start: ", start, " Lin: ", lin, ")"),
                     x = "UMAP 1", y = "UMAP 2") +
                theme(plot.title = element_text(hjust = 0.5, face = "bold"))

            # 3. Add Zones mapping fill to density level, and linetype to the Zone
            for (i in seq_along(valley_times)) {
                z_label <- paste0("Transition Zone ", i)
                zone_subdata <- umap_df %>% filter(Zone == z_label)
                
                if (nrow(zone_subdata) > 5) { 
                    p_umap <- p_umap + stat_density_2d(
                        data = zone_subdata, 
                        aes(x = umap_1, y = umap_2, fill = after_stat(level), linetype = Zone), 
                        geom = "polygon", alpha = 0.12, bins = 4, 
                        color = zone_colors_wt[z_label], linewidth = 0.8
                    )
                }
            }
            
            # 4. Integrate Legends & Scales matching your original formatting
            p_umap <- p_umap + 
                scale_fill_gradient(name = "Contour Density", low = "#ECEFF1", high = "#37474F") +
                scale_linetype_manual(name = "Trajectory Regions", values = setNames(rep("solid", length(valley_times)), names(zone_colors_wt)), drop = FALSE) +
                guides(
                    color = guide_colorbar(title = "Pseudotime", order = 1),
                    linetype = guide_legend(order = 2, override.aes = list(color = unname(zone_colors_wt), fill = NA, linewidth = 1.5)),
                    fill = guide_colorbar(order = 3)
                )

            # 5. Save Combined Plot
            ggsave(paste0(outfile_dir, "meox1__DimLabel_WT_Anchored_Contours_Combined_start_", start, "_lin", lin, ".pdf"), 
                   plot = p_umap, height = 7, width = 8.5)
            
            # 6. Save Split (Faceted) Plot
            p_umap_split <- p_umap + facet_wrap(~ Sample) +
                theme(strip.background = element_rect(fill = "grey90", color = "black"),
                      strip.text = element_text(face = "bold", size = 11))
            
            ggsave(paste0(outfile_dir, "meox1__DimLabel_WT_Anchored_Contours_Split_start_", start, "_lin", lin, ".pdf"), 
                   plot = p_umap_split, height = 6, width = 12.5)

            # 4. COMPUTE SUMMARIES AND PLOT UPDATED COUNTS + PERCENTAGE BARPLOTS
            zone_names <- c(paste0("Transition Zone ", seq_along(valley_times)), "Core/Other Cells", "Trajectory NA (Unmapped)")
            
            # Isolate target zones for the focused presentation
            target_zones_list <- paste0("Transition Zone ", seq_along(valley_times))
            
            barplot_data <- master_zones_df %>%
                filter(Zone %in% target_zones_list) %>%
                group_by(Sample, Zone) %>%
                summarise(Count = n(), .groups = 'drop') %>%
                group_by(Sample) %>%
                mutate(
                    TotalInSample = sum(Count),
                    Percentage = (Count / TotalInSample) * 100
                )

            if (nrow(barplot_data) > 0) {
                
                # Make dynamic palette for stacked bars
                zone_colors <- c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "#F57C00", "Transition Zone 3" = "#7B1FA2", "Transition Zone 4" = "#388E3C")
                
                bp <- ggplot(barplot_data, aes(x = Sample, y = Percentage, fill = Zone)) +
                    geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", size = 0.3) +
                    # Add string text mapping directly over stacks containing both % and count (n)
                    geom_text(
                        aes(label = paste0(sprintf("%.1f", Percentage), "%\n(n=", Count, ")")),
                        position = position_stack(vjust = 0.5), 
                        size = 3.2, 
                        fontface = "bold",
                        color = "white"
                    ) +
                    scale_fill_manual(values = zone_colors) +
                    labs(
                        title = paste0("WT-Anchored Zone Allocations\n(Start: ", start, " Lin: ", lin, ")"),
                        x = "Genotype Pool", 
                        y = "Proportional Composition (%)"
                    ) + 
                    theme_classic() +
                    theme(
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text = element_text(color = "black", size = 10)
                    )

                ggsave(paste0(outfile_dir, "meox1__BarLabel_WT_Anchored_Quantified_Proportions_start_", start, "_lin", lin, ".pdf"), plot = bp, height = 6, width = 5.5)
            }
            
            # --- HELPER FUNCTION EXECUTIONS ---
            has_wt  <- "WildType" %in% unique(master_zones_df$Sample)
            has_mut <- "meox1 -/- mutant" %in% unique(master_zones_df$Sample)
            
            # hmVEC Lin 1
            if (start == "hmVEC" && lin == 1) {
                if (has_wt && has_mut) {
                    plot_custom_transition_zones(
                        master_df      = master_zones_df,
                        target_zones   = c(1, 2, 3),
                        custom_colors  = c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "#F57C00", "Transition Zone 3" = "#7B1FA2"),
                        title_prefix   = "WT-Anchored hmVEC Lin 1 (Combined)",
                        file_prefix    = "Custom_WT_Anchored_hmVEC_Lin1",
                        col_name       = col_name,
                        outfile_dir    = outfile_dir
                    )
                }
            }

            # mVEC Lin 1
            if (start == "mVEC" && lin == 1) {
                if (has_wt && has_mut) {
                    plot_custom_transition_zones(
                        master_df      = master_zones_df,
                        target_zones   = c(1, 2),
                        custom_colors  = c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "#F57C00"),
                        title_prefix   = "WT-Anchored mVEC Lin 1 (Combined)",
                        file_prefix    = "Custom_WT_Anchored_mVEC_Lin1",
                        col_name       = col_name,
                        outfile_dir    = outfile_dir
                    )
                }
            }

            # mVEC Lin 2
            if (start == "mVEC" && lin == 2) {
                if (has_wt && has_mut) {
                    plot_custom_transition_zones(
                        master_df      = master_zones_df,
                        target_zones   = c(1, 2),
                        custom_colors  = c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "#F57C00"),
                        title_prefix   = "WT-Anchored mVEC Lin 2 (Combined)",
                        file_prefix    = "Custom_WT_Anchored_mVEC_Lin2",
                        col_name       = col_name,
                        outfile_dir    = outfile_dir
                    )
                }
            }
        }
    }
}

run_slingshot <- function(data, col_map, col_order, col_name) {
    sce <- as.SingleCellExperiment(data)
    celltype <- as.vector(unique(data@meta.data$L3_celltype))
    data$L3_celltype <- factor(data$L3_celltype, levels=col_order)
    
    # dataset-level composition barplots (unified scale)
    has_mutant <- any(grepl("_1$", colnames(data)))
    has_wt     <- any(grepl("_2$", colnames(data)))

    # if (has_mutant & has_wt) {
    if (nrow(data) > 0) {
        cat("Generating dataset-level composition barplots...\n")
        
        comp_df <- data.frame(Barcode = colnames(data), Cluster = data@meta.data$L3_celltype) %>%
            mutate(Sample = case_when(grepl("_1$", Barcode) ~ "meox1 -/- mutant", 
                                      grepl("_2$", Barcode) ~ "WildType", 
                                      TRUE ~ "Unknown")) %>%
            filter(Sample != "Unknown")

        comp_df$Sample <- factor(comp_df$Sample, levels = c("meox1 -/- mutant", "WildType"))
        sample_colors <- c("meox1 -/- mutant" = "#E41A1C", "WildType" = "#377EB8")

        # 1. calculate Total Composition
        tot_df <- comp_df %>% 
            group_by(Sample) %>% 
            tally(name = "Cell_Count") %>%
            mutate(Total = sum(Cell_Count), 
                   Proportion = (Cell_Count / Total) * 100, 
                   X_Category = "Total Dataset")

        # 2. calculate Cluster Composition
        clust_df <- comp_df %>% 
            group_by(Cluster, Sample) %>% 
            tally(name = "Cell_Count") %>%
            group_by(Cluster) %>% 
            mutate(Total = sum(Cell_Count), 
                   Proportion = (Cell_Count / Total) * 100,
                   X_Category = as.character(Cluster)) %>%
            ungroup() %>%
            select(-Cluster) 

        # 3. combine into a single dataframe
        combined_df <- bind_rows(tot_df, clust_df)
        
        # 4. factor X_Category so "Total Dataset" is first
        combined_df$X_Category <- factor(combined_df$X_Category, levels = c("Total Dataset", col_order))

        # 5. calculate manual Y positions for stacked labels
        combined_df <- combined_df %>%
            arrange(X_Category, desc(Sample)) %>%
            group_by(X_Category) %>%
            mutate(label_y = cumsum(Cell_Count) - (0.5 * Cell_Count))

        # single absolute
        p_combined_abs <- ggplot(combined_df, aes(x = X_Category, y = Cell_Count, fill = Sample)) +
            geom_col(width = 0.6, color = "black", linewidth = 0.5) +
            geom_text_repel(data = subset(combined_df, X_Category == "Total Dataset"),
                            aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), 
                            nudge_x = 0.60, direction = "y", hjust = 0, segment.size = 0.3) +
            geom_text_repel(data = subset(combined_df, X_Category != "Total Dataset"),
                            aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), 
                            nudge_x = 0.35, direction = "y", hjust = 0, segment.size = 0.3) +
            scale_fill_manual(values = sample_colors) +
            labs(title = "Cells per Cluster (Including Total)", x = "Cluster / Group", y = "Absolute Number of Cells") +
            theme_classic() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
                  axis.text.y = element_text(color = "black"), 
                  axis.title = element_text(face="bold"), 
                  plot.title = element_text(hjust=0.5, face="bold"))
        
        # single proportion
        p_combined_prop <- ggplot(combined_df, aes(x = X_Category, y = Proportion, fill = Sample)) +
            geom_col(width = 0.6, color = "black", linewidth = 0.5) +
            scale_fill_manual(values = sample_colors) +
            coord_cartesian(ylim = c(0, 100)) +
            labs(title = "Proportion per Cluster (Including Total)", x = "Cluster / Group", y = "% of Group") +
            theme_classic() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
                  axis.text.y = element_text(color = "black"), 
                  axis.title = element_text(face="bold"), 
                  plot.title = element_text(hjust=0.5, face="bold"))

        ggsave(paste0(outfile_dir, "meox1__BarLabel_Combined_Total_and_Cluster_Absolute.", col_name, ".pdf"), p_combined_abs, height=5, width=7)
        ggsave(paste0(outfile_dir, "meox1__BarLabel_Combined_Total_and_Cluster_Proportion.", col_name, ".pdf"), p_combined_prop, height=5, width=7)

        p_combined_abs_strip <- p_combined_abs + NoAxes() + NoLegend() + theme(plot.title=element_blank())
        ggsave(paste0(outfile_dir, "meox1__BarStrip_Combined_Total_and_Cluster_Absolute.", col_name, ".pdf"), p_combined_abs_strip, height=4, width=6)

        p_combined_prop_strip <- p_combined_prop + NoAxes() + NoLegend() + theme(plot.title=element_blank())
        ggsave(paste0(outfile_dir, "meox1__BarStrip_Combined_Total_and_Cluster_Proportion.", col_name, ".pdf"), p_combined_prop_strip, height=4, width=6)

        # per-sample unified composition plots (total + clusters)
        cat("Generating per-sample unified composition barplots...\n")
        
        for (current_sample in levels(comp_df$Sample)) {
            
            sample_df <- comp_df %>% filter(Sample == current_sample)
            
            if (nrow(sample_df) > 0) {
                
                # added group_by(sample) so the column isn't dropped by tally()
                tot_s <- sample_df %>% 
                    group_by(Sample) %>%
                    tally(name = "Cell_Count") %>% 
                    mutate(Total = Cell_Count, Proportion = 100, X_Category = paste0("Total (", current_sample, ")"))
                
                # added sample grouping, and added ungroup() before mutate so proportions calculate correctly
                clust_s <- sample_df %>% 
                    group_by(Cluster, Sample) %>% 
                    tally(name = "Cell_Count") %>% 
                    ungroup() %>%
                    mutate(Total = sum(Cell_Count), Proportion = (Cell_Count / Total) * 100, X_Category = as.character(Cluster)) %>% 
                    select(-Cluster)
                
                combined_s <- bind_rows(tot_s, clust_s)
                combined_s$X_Category <- factor(combined_s$X_Category, levels = c(paste0("Total (", current_sample, ")"), col_order))
                combined_s <- combined_s %>% mutate(label_y = Cell_Count * 0.5)
                safe_sample_name <- gsub("[ /]", "_", current_sample)
          
                p_samp_abs <- ggplot(combined_s, aes(x = X_Category, y = Cell_Count, fill = Sample)) +
                    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                    geom_text_repel(data = subset(combined_s, X_Category == paste0("Total (", current_sample, ")")), aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.35, direction = "y", hjust = 0, segment.size = 0.3) +
                    geom_text_repel(data = subset(combined_s, X_Category != paste0("Total (", current_sample, ")")), aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.35, direction = "y", hjust = 0, segment.size = 0.3) +
                    scale_fill_manual(values = sample_colors) +
                    labs(title = paste0(current_sample, " - Cells per Cluster (Including Total)"), x = "Cluster", y = "Absolute Number of Cells") +
                    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(face="bold"), plot.title = element_text(hjust=0.5, face="bold"), legend.position = "none") 
                
                p_samp_prop <- ggplot(combined_s, aes(x = X_Category, y = Proportion, fill = Sample)) +
                    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                    scale_fill_manual(values = sample_colors) +
                    coord_cartesian(ylim = c(0, 100)) +
                    labs(title = paste0(current_sample, " - Proportion per Cluster (Including Total)"), x = "Cluster", y = "% of Sample") +
                    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(face="bold"), plot.title = element_text(hjust=0.5, face="bold"), legend.position = "none")

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Sample_", safe_sample_name, "_Unified_Absolute.", col_name, ".pdf"), p_samp_abs, height=5, width=7)
                ggsave(paste0(outfile_dir, "meox1__BarLabel_Sample_", safe_sample_name, "_Unified_Proportion.", col_name, ".pdf"), p_samp_prop, height=5, width=7)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Sample_", safe_sample_name, "_Unified_Absolute.", col_name, ".pdf"), p_samp_abs + NoAxes() + NoLegend() + theme(plot.title=element_blank()), height=4, width=6)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Sample_", safe_sample_name, "_Unified_Proportion.", col_name, ".pdf"), p_samp_prop + NoAxes() + NoLegend() + theme(plot.title=element_blank()), height=4, width=6)
            }
        }
    }

    for (start in celltype) {
        sce <- slingshot(
            sce, 
            clusterLabels = "L3_celltype", 
            reducedDim = "UMAP",           
            start.clus = start             
        )

        pseudotime_matrix <- slingPseudotime(sce)
        n_lineages <- ncol(pseudotime_matrix)
        max_pseudotime_global <- max(pseudotime_matrix, na.rm = TRUE)

        curves <- slingCurves(sce, as.df = TRUE)
        curves <- curves %>%
            group_by(Lineage) %>%
            arrange(Order, .by_group = TRUE)

        weights_matrix <- slingCurveWeights(SlingshotDataSet(sce))

        for (lin in seq_len(n_lineages)) {
            
            meta <- paste0("Slingshot_Pseudotime_start_", start, "_lin", lin)
            data[[meta]] <- pseudotime_matrix[, lin]
            
            cat("Processing lineage:", meta, "\n")
            print(summary(data[[meta]]))

            lineage_curve <- curves %>% filter(Lineage == lin)
            dens <- density(data@meta.data[[meta]], na.rm = TRUE)

            valleys <- which(diff(sign(diff(dens$y))) == 2) + 1
            valley_times <- dens$x[valleys]

            meta_data <- data@meta.data
            cluster_transition <- paste0("Cluster_Boundaries_start_", start, "_lin", lin)
            meta_data[[cluster_transition]] <- "Core Cluster"

            for (v in valley_times) {
                meta_data[[cluster_transition]][
                    abs(meta_data[[meta]] - v) < 0.5
                ] <- "Transitional Boundary"
            }
            data[[cluster_transition]] <- meta_data[[cluster_transition]]

            umap_embed <- as.data.frame(Embeddings(data, reduction = "umap"))
            colnames(umap_embed)[1:2] <- c("umap_1", "umap_2")
            umap_embed$State <- meta_data[[cluster_transition]]

            transition_zones_df <- data.frame(
                umap_1     = umap_embed$umap_1,
                umap_2     = umap_embed$umap_2,
                Pseudotime = data@meta.data[[meta]]
            ) %>%
                mutate(Zone_ID = "Other")

            for (z_idx in seq_along(valley_times)) {
                v_time <- valley_times[z_idx]
                is_in_valley <- !is.na(transition_zones_df$Pseudotime) & (abs(transition_zones_df$Pseudotime - v_time) < 0.5)
                transition_zones_df$Zone_ID[is_in_valley] <- paste0("Transition Zone ", z_idx)
            }

            transition_zones_only <- transition_zones_df %>% filter(Zone_ID != "Other")

            accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")
            zone_names <- paste0("Transition Zone ", seq_along(valley_times))

            p_dens_df <- data.frame(
                pseudotime_axis = dens$x,
                density_axis    = dens$y
            )

            # minima profile 
            plt_minima <- ggplot(p_dens_df, aes(x = pseudotime_axis, y = density_axis))
            for (i in seq_along(valley_times)) {
                v <- valley_times[i]
                rect_color <- accent_palette[((i - 1) %% length(accent_palette)) + 1]
                plt_minima <- plt_minima + 
                    annotate("rect", xmin = v - 0.5, xmax = v + 0.5, 
                             ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_color)
            }
            plt_minima <- plt_minima +
                geom_line(linewidth = 0.8, color = "black") +
                geom_vline(xintercept = valley_times, color = "red", linewidth = 1, linetype = "solid") +
                labs(title = paste0("Local Minima Profile (Start: ", start, " Lin: ", lin, ")"), x = "Pseudotime", y = "Density") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                      axis.text = element_text(color = "black", size = 11),
                      axis.title = element_text(face = "bold", size = 12))

            ggsave(paste0(outfile_dir, "meox1__DensLabel_Minima_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_minima, height = 5, width = 6)
            ggsave(paste0(outfile_dir, "meox1__DensStrip_Minima_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_minima + NoAxes() + theme(plot.title = element_blank()), height = 4, width = 5)
            rm(p_dens_df, plt_minima)

            # dynamic lineage weights
            weights_col_name <- paste0("Lineage_Weights_start_", start, "_lin", lin)
            data[[weights_col_name]] <- weights_matrix[, lin]

            meta_data <- data@meta.data
            transition_col_name <- paste0("Transition_State_start_", start, "_lin", lin)

            meta_data <- meta_data %>%
                mutate(
                    "{transition_col_name}" := case_when(
                        is.na(.data[[meta]]) | .data[[weights_col_name]] == 0 ~ "Unrelated",
                        .data[[weights_col_name]] > 0.05 & .data[[weights_col_name]] < 0.95 ~ "Transitional",
                        .data[[weights_col_name]] >= 0.95 ~ "Mature/Core",
                        TRUE ~ "Unrelated"
                    )
                )
            data[[transition_col_name]] <- meta_data[[transition_col_name]]

            # dimplot celltype
            plt_dim <- DimPlot(data, group.by = "L3_celltype", cols=col_map)
            plt_dim_arr <- plt_dim + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            ggsave(paste0(outfile_dir, "meox1__DimLabel_Slingshot_Pseudotime_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr, height = 7, width = 7)
            ggsave(paste0(outfile_dir, "meox1__DimStrip_Slingshot_Pseudotime_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # dimplot celltype contours
            plt_dim_arr <- plt_dim + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            for (i in seq_along(valley_times)) {
                z_label <- paste0("Transition Zone ", i)
                zone_subdata <- transition_zones_only %>% filter(Zone_ID == z_label)
                plt_dim_arr <- plt_dim_arr + stat_density_2d(data = zone_subdata, aes(x = umap_1, y = umap_2, fill = after_stat(level), linetype = Zone_ID), geom = "polygon", alpha = 0.12, bins = 4, color = accent_palette[((i - 1) %% length(accent_palette)) + 1], linewidth = 1.0)
            }
            plt_dim_arr <- plt_dim_arr + scale_fill_gradient(name = "Contour Density", low = "#ECEFF1", high = "#37474F") + scale_linetype_manual(name = "Trajectory Regions", values = setNames(rep("solid", length(valley_times)), zone_names)) + guides(color = guide_legend(order = 1, override.aes = list(alpha = 1, fill = NA, size = 4, shape = 16, linetype = 0, stroke = 0)), fill = guide_colorbar(order = 3), linetype = guide_legend(order = 2, override.aes = list(color = accent_palette[1:length(valley_times)], linewidth = 1.5)))
            ggsave(paste0(outfile_dir, "meox1__DimLabel_Celltype_Contours_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr, height = 7, width = 8)
            ggsave(paste0(outfile_dir, "meox1__DimStrip_Celltype_Contours_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # dimplot cluster transition zones
            plt_dim_trans <- DimPlot(data, group.by = cluster_transition)
            plt_dim_arr <- plt_dim_trans + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            ggsave(paste0(outfile_dir, "meox1__DimLabel_Cluster_Transition_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr, height = 7, width = 7)
            ggsave(paste0(outfile_dir, "meox1__DimStrip_Cluster_Transition_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # dimplot transition state
            plt_dim_state <- DimPlot(data, group.by = transition_col_name, cols = c("Transitional" = "red", "Mature/Core" = "lightgrey", "Unrelated" = "whitesmoke"))
            plt_dim_arr <- plt_dim_state + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            ggsave(paste0(outfile_dir, "meox1__DimLabel_Transition_State_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr, height = 7, width = 7)
            ggsave(paste0(outfile_dir, "meox1__DimStrip_Transition_State_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dim_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # feaplot transition score weights
            plt_fea_weights <- FeaturePlot(data, features = weights_col_name) + scale_color_viridis_c(option = "plasma")
            plt_fea_arr <- plt_fea_weights + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            ggsave(paste0(outfile_dir, "meox1__FeaLabel_Transition_Scores_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_arr, height = 7, width = 7)
            ggsave(paste0(outfile_dir, "meox1__FeaStrip_Transition_Scores_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # feaplot pseudotime scores
            plt_fea_pseudo <- FeaturePlot(data, features = meta) + scale_color_viridis_c(option = "plasma", limits = c(0, max_pseudotime_global))
            plt_fea_arr <- plt_fea_pseudo + geom_path(data = lineage_curve, aes(x = umap_1, y = umap_2), size = 1.2, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"))
            ggsave(paste0(outfile_dir, "meox1__FeaLabel_Pseudotime_Scores_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_arr, height = 7, width = 7)
            ggsave(paste0(outfile_dir, "meox1__FeaStrip_Pseudotime_Scores_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # feaplot pseudotime contours
            plt_fea_contour <- plt_fea_pseudo
            for (i in seq_along(valley_times)) {
                z_label <- paste0("Transition Zone ", i)
                zone_subdata <- transition_zones_only %>% filter(Zone_ID == z_label)
                plt_fea_contour <- plt_fea_contour + stat_density_2d(data = zone_subdata, aes(x = umap_1, y = umap_2, fill = after_stat(level), linetype = Zone_ID), geom = "polygon", alpha = 0.12, bins = 4, color = accent_palette[((i - 1) %% length(accent_palette)) + 1], linewidth = 1.0)
            }
            plt_fea_contour_arr <- plt_fea_contour + scale_fill_gradient(name = "Contour Density", low = "#ECEFF1", high = "#37474F") + scale_linetype_manual(name = "Trajectory Regions", values = setNames(rep("solid", length(valley_times)), zone_names)) + guides(color = guide_colorbar(title = "Pseudotime", order = 1), fill = guide_colorbar(order = 3), linetype = guide_legend(order = 2, override.aes = list(color = accent_palette[1:length(valley_times)], linewidth = 1.5)))
            ggsave(paste0(outfile_dir, "meox1__FeaLabel_Pseudotime_Contours_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_contour_arr, height = 7, width = 8)
            ggsave(paste0(outfile_dir, "meox1__FeaStrip_Pseudotime_Contours_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_fea_contour_arr + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 5, width = 5)

            # barplots cell proportion (transition zones)
            # if (has_mutant & has_wt) {
            if (nrow(data) > 0) {
                
                master_zones_df <- data.frame(
                    Barcode = colnames(data), 
                    Pseudotime = data@meta.data[[meta]],
                    Cluster = data@meta.data$L3_celltype
                ) %>%
                    mutate(Sample = case_when(grepl("_1$", Barcode) ~ "meox1 -/- mutant", grepl("_2$", Barcode) ~ "WildType", TRUE ~ "Unknown"),
                           Zone = case_when(is.na(Pseudotime) ~ "Trajectory NA (Unmapped)", TRUE ~ "Core/Other Cells")) %>%
                    filter(Sample != "Unknown")

                master_zones_df$Sample <- factor(master_zones_df$Sample, levels = c("meox1 -/- mutant", "WildType"))
                sample_colors <- c("meox1 -/- mutant" = "#E41A1C", "WildType" = "#377EB8")

                for (i in seq_along(valley_times)) {
                    v_time <- valley_times[i]
                    is_in_zone <- !is.na(master_zones_df$Pseudotime) & (abs(master_zones_df$Pseudotime - v_time) < 0.5)
                    master_zones_df$Zone[is_in_zone] <- paste0("Transition Zone ", i)
                }

                sample_baselines <- master_zones_df %>% group_by(Sample) %>% tally(name = "Total_Sample_Cells")

                max_proportion_global <- 0
                for (i in seq_along(valley_times)) {
                    z_label <- paste0("Transition Zone ", i)
                    current_prop <- master_zones_df %>% filter(Zone == z_label) %>% group_by(Sample) %>% tally() %>% left_join(sample_baselines, by = "Sample") %>% mutate(Proportion = (n / Total_Sample_Cells) * 100) %>% pull(Proportion)
                    if (length(current_prop) > 0) { max_proportion_global <- max(max_proportion_global, max(current_prop)) }
                }
                max_proportion_global <- min(100, max_proportion_global * 1.05)

                for (i in seq_along(valley_times)) {
                    v_time <- valley_times[i]
                    z_label <- paste0("Transition Zone ", i)
                    
                    indiv_summary <- master_zones_df %>%
                        mutate(Indiv_Status = case_when(Zone == z_label ~ z_label, Zone == "Trajectory NA (Unmapped)" ~ "Trajectory NA (Unmapped)", TRUE ~ "Core/Other Cells")) %>%
                        group_by(Sample, Indiv_Status) %>% tally(name = "Cell_Count") %>% ungroup() %>%
                        left_join(sample_baselines, by = "Sample") %>% mutate(Proportion = (Cell_Count / Total_Sample_Cells) * 100)

                    isolated_zone <- indiv_summary %>% filter(Indiv_Status == z_label)

                    plt_isolated <- ggplot(isolated_zone, aes(x = Sample, y = Proportion, fill = Sample)) +
                        geom_col(width = 0.5, color = "black", linewidth = 0.5) +
                        scale_fill_manual(values = sample_colors) +
                        coord_cartesian(ylim = c(0, max_proportion_global)) + 
                        labs(title = paste0("Transition Zone ", i, " (Pseudotime ~ ", round(v_time, 2), ")"), x = NULL, y = "% of Sample Cells") +
                        theme_classic() + theme(legend.position = "none", axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Transition_Zone_", i, "_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_isolated, height = 5, width = 5)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Transition_Zone_", i, "_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_isolated + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                    indiv_summary$Indiv_Status <- factor(indiv_summary$Indiv_Status, levels = c("Trajectory NA (Unmapped)", "Core/Other Cells", z_label))
                    
                    indiv_summary <- indiv_summary %>% arrange(Sample, desc(Indiv_Status)) %>% group_by(Sample) %>% mutate(label_y = cumsum(Cell_Count) - 0.5 * Cell_Count)

                    plt_stacked_indiv <- ggplot(indiv_summary, aes(x = Sample, y = Cell_Count, fill = Indiv_Status)) +
                        geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                        geom_text_repel(aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.45, direction = "y", hjust = 0, segment.size = 0.3) +
                        scale_fill_manual(values = setNames(c("#D32F2F", "#CFD8DC", "#37474F"), c(z_label, "Core/Other Cells", "Trajectory NA (Unmapped)"))) +
                        labs(title = paste0("Dataset Alignment Fit: Transition Zone ", i), x = NULL, y = "Total Number of Cells", fill = "State Definition") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Transition_Zone_", i, "_StackedTotal_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_indiv, height = 5, width = 6)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Transition_Zone_", i, "_StackedTotal_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_indiv + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                    # cluster per zone
                    zone_df <- master_zones_df %>% filter(Zone == z_label)
                    
                    if (nrow(zone_df) > 0) {
                        tot_z <- zone_df %>% group_by(Sample) %>% tally(name = "Cell_Count") %>% mutate(Total = sum(Cell_Count), Proportion = (Cell_Count / Total) * 100, X_Category = paste0("Total (Zone ", i, ")"))
                        clust_z <- zone_df %>% group_by(Cluster, Sample) %>% tally(name = "Cell_Count") %>% group_by(Cluster) %>% mutate(Total = sum(Cell_Count), Proportion = (Cell_Count / Total) * 100, X_Category = as.character(Cluster)) %>% ungroup() %>% select(-Cluster)

                        combined_z <- bind_rows(tot_z, clust_z)
                        combined_z$X_Category <- factor(combined_z$X_Category, levels = c(paste0("Total (Zone ", i, ")"), col_order))

                        combined_z <- combined_z %>% arrange(X_Category, desc(Sample)) %>% group_by(X_Category) %>% mutate(label_y = cumsum(Cell_Count) - (0.5 * Cell_Count))

                        p_z_abs <- ggplot(combined_z, aes(x = X_Category, y = Cell_Count, fill = Sample)) +
                            geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                            geom_text_repel(data = subset(combined_z, X_Category == paste0("Total (Zone ", i, ")")), aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.35, direction = "y", hjust = 0, segment.size = 0.3) +
                            geom_text_repel(data = subset(combined_z, X_Category != paste0("Total (Zone ", i, ")")), aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.35, direction = "y", hjust = 0, segment.size = 0.3) +
                            scale_fill_manual(values = sample_colors) +
                            labs(title = paste0("Zone ", i, " Cells per Cluster (Including Total)"), x = "Cluster", y = "Absolute Number of Cells") +
                            theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(face="bold"), plot.title = element_text(hjust=0.5, face="bold"), legend.position = "none") 
                        
                        p_z_prop <- ggplot(combined_z, aes(x = X_Category, y = Proportion, fill = Sample)) +
                            geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                            scale_fill_manual(values = sample_colors) +
                            coord_cartesian(ylim = c(0, 100)) +
                            labs(title = paste0("Proportion per Cluster (Zone ", i, ")"), x = "Cluster", y = "% of Zone") +
                            theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(face="bold"), plot.title = element_text(hjust=0.5, face="bold"), legend.position = "none")

                        ggsave(paste0(outfile_dir, "meox1__BarLabel_Transition_Zone_", i, "_Unified_Absolute_start_", start, "_lin", lin, ".", col_name, ".pdf"), p_z_abs, height=5, width=7)
                        ggsave(paste0(outfile_dir, "meox1__BarLabel_Transition_Zone_", i, "_Unified_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), p_z_prop, height=5, width=7)
                        ggsave(paste0(outfile_dir, "meox1__BarStrip_Transition_Zone_", i, "_Unified_Absolute_start_", start, "_lin", lin, ".", col_name, ".pdf"), p_z_abs + NoAxes() + NoLegend() + theme(plot.title=element_blank()), height=4, width=6)
                        ggsave(paste0(outfile_dir, "meox1__BarStrip_Transition_Zone_", i, "_Unified_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), p_z_prop + NoAxes() + NoLegend() + theme(plot.title=element_blank()), height=4, width=6)
                    }
                }

                combined_summary <- master_zones_df %>% group_by(Sample, Zone) %>% tally(name = "Cell_Count") %>% ungroup()
                zone_levels <- c("Trajectory NA (Unmapped)", "Core/Other Cells", paste0("Transition Zone ", seq_along(valley_times)))
                combined_summary$Zone <- factor(combined_summary$Zone, levels = zone_levels)
                combined_summary <- combined_summary %>% arrange(Sample, desc(Zone)) %>% group_by(Sample) %>% mutate(label_y = cumsum(Cell_Count) - 0.5 * Cell_Count)

                zone_colors <- c("Trajectory NA (Unmapped)" = "#37474F", "Core/Other Cells" = "#CFD8DC")
                for (i in seq_along(valley_times)) { zone_colors[paste0("Transition Zone ", i)] <- accent_palette[((i - 1) %% length(accent_palette)) + 1] }

                plt_combined <- ggplot(combined_summary, aes(x = Sample, y = Cell_Count, fill = Zone)) +
                    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                    geom_text_repel(aes(y = label_y, label = ifelse(Cell_Count > 0, Cell_Count, "")), nudge_x = 0.45, direction = "y", hjust = 0, segment.size = 0.3) +
                    scale_fill_manual(values = zone_colors) +
                    labs(title = paste0("Transition Zones Combined (Start: ", start, " Lin: ", lin, ")"), x = NULL, y = "Total Number of Cells", fill = "Trajectory Regions") +
                    theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

                ggsave(paste0(outfile_dir, "meox1__BarLabel_All_Transition_Zones_Stacked_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_combined, height = 5, width = 6.5)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_All_Transition_Zones_Stacked_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_combined + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                # 1. isolate just the transition zones
                zones_only_summary <- master_zones_df %>% 
                    filter(grepl("Transition Zone", Zone)) %>% 
                    group_by(Sample, Zone) %>% 
                    tally(name = "Cell_Count") %>% 
                    ungroup()
                
                # calculate proportions based strictly on isolated zone totals per sample
                zones_only_summary <- zones_only_summary %>%
                    group_by(Sample) %>%
                    mutate(Zone_Total = sum(Cell_Count),
                           Proportion = (Cell_Count / Zone_Total) * 100) %>%
                    ungroup()
                
                # update factor levels to just the zones
                zone_only_levels <- paste0("Transition Zone ", seq_along(valley_times))
                zones_only_summary$Zone <- factor(zones_only_summary$Zone, levels = zone_only_levels)
                
                # calculate Y positions for stacked versions (both absolute and percentage)
                zones_only_summary <- zones_only_summary %>% 
                    arrange(Sample, desc(Zone)) %>% 
                    group_by(Sample) %>% 
                    mutate(label_y_abs = cumsum(Cell_Count) - (0.5 * Cell_Count),
                           label_y_prop = cumsum(Proportion) - (0.5 * Proportion)) %>%
                    ungroup()

                # plot a1 stacked zones only absolute
                plt_stacked_zones_only <- ggplot(zones_only_summary, aes(x = Sample, y = Cell_Count, fill = Zone)) +
                    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                    geom_text_repel(aes(y = label_y_abs, label = ifelse(Cell_Count > 0, Cell_Count, "")), 
                                    nudge_x = 0.45, direction = "y", hjust = 0, segment.size = 0.3) +
                    scale_fill_manual(values = zone_colors) +
                    labs(title = paste0("Transition Zones Only (Start: ", start, " Lin: ", lin, ")"), x = NULL, y = "Total Number of Cells", fill = "Trajectory Regions") +
                    theme_classic() + 
                    theme(axis.text = element_text(color = "black", size = 11), 
                          axis.title.y = element_text(face = "bold", size = 12), 
                          plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Isolated_Zones_Stacked_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_zones_only, height = 5, width = 6.5)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Isolated_Zones_Stacked_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_zones_only + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                # plot a2 stacked zones only percentage
                plt_stacked_zones_prop <- ggplot(zones_only_summary, aes(x = Sample, y = Proportion, fill = Zone)) +
                    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
                    geom_text_repel(aes(y = label_y_prop, label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), "")), 
                                    nudge_x = 0.45, direction = "y", hjust = 0, segment.size = 0.3) +
                    scale_fill_manual(values = zone_colors) +
                    coord_cartesian(ylim = c(0, 100)) +
                    labs(title = paste0("Transition Zones Composition % (Start: ", start, " Lin: ", lin, ")"), x = NULL, y = "% of Transition Zone Cells", fill = "Trajectory Regions") +
                    theme_classic() + 
                    theme(axis.text = element_text(color = "black", size = 11), 
                          axis.title.y = element_text(face = "bold", size = 12), 
                          plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Isolated_Zones_Stacked_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_zones_prop, height = 5, width = 6.5)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Isolated_Zones_Stacked_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_stacked_zones_prop + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                # plot b side by side zones only absolute
                max_count <- max(zones_only_summary$Cell_Count, na.rm = TRUE)
                
                plt_dodged_zones_only <- ggplot(zones_only_summary, aes(x = Sample, y = Cell_Count, fill = Zone)) +
                    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
                    geom_text(aes(label = ifelse(Cell_Count > 0, Cell_Count, ""), y = Cell_Count + (max_count * 0.02)), 
                              position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
                    scale_fill_manual(values = zone_colors) +
                    labs(title = paste0("Transition Zones by Sample (Start: ", start, " Lin: ", lin, ")"), x = NULL, y = "Total Number of Cells", fill = "Trajectory Regions") +
                    theme_classic() + 
                    theme(axis.text = element_text(color = "black", size = 11), 
                          axis.title.y = element_text(face = "bold", size = 12), 
                          plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Isolated_Zones_Dodged_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dodged_zones_only, height = 5, width = 7)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Isolated_Zones_Dodged_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dodged_zones_only + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)

                # plot b side by side zones only percentage
                max_prop <- max(zones_only_summary$Proportion, na.rm = TRUE)

                plt_dodged_zones_prop <- ggplot(zones_only_summary, aes(x = Sample, y = Proportion, fill = Zone)) +
                    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
                    geom_text(aes(label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), ""), y = Proportion + (max_prop * 0.02)), 
                              position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
                    scale_fill_manual(values = zone_colors) +
                    coord_cartesian(ylim = c(0, min(100, max_prop * 1.15))) + # dynamic ceiling with head room
                    labs(title = paste0("Transition Zones % by Sample (Start: ", start, " Lin: ", lin, ")"), x = NULL, y = "% of Transition Zone Cells", fill = "Trajectory Regions") +
                    theme_classic() + 
                    theme(axis.text = element_text(color = "black", size = 11), 
                          axis.title.y = element_text(face = "bold", size = 12), 
                          plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Isolated_Zones_Dodged_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dodged_zones_prop, height = 5, width = 7)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Isolated_Zones_Dodged_Proportion_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dodged_zones_prop + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)

                dataset_zones_only <- master_zones_df %>% 
                    filter(grepl("Transition Zone", Zone)) %>% 
                    group_by(Zone) %>% 
                    tally(name = "Cell_Count") %>% 
                    ungroup()
                
                dataset_zones_only$Zone <- factor(dataset_zones_only$Zone, levels = zone_only_levels)
                max_ds_count <- max(dataset_zones_only$Cell_Count, na.rm = TRUE)

                plt_dataset_zones <- ggplot(dataset_zones_only, aes(x = Zone, y = Cell_Count, fill = Zone)) +
                    geom_col(color = "black", linewidth = 0.5, width = 0.6) +
                    geom_text(aes(label = ifelse(Cell_Count > 0, Cell_Count, ""), y = Cell_Count + (max_ds_count * 0.02)), 
                              vjust = 0, size = 4) +
                    scale_fill_manual(values = zone_colors) +
                    labs(title = paste0("Dataset Cells per Transition Zone (", start, " Lin: ", lin, ")"), 
                         x = NULL, y = "Total Number of Cells") +
                    theme_classic() + 
                    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11),
                          axis.text.y = element_text(color = "black", size = 11), 
                          axis.title.y = element_text(face = "bold", size = 12), 
                          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                          legend.position = "none") # Legend omitted since x-axis explains it

                ggsave(paste0(outfile_dir, "meox1__BarLabel_Isolated_Zones_DatasetLevel_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dataset_zones, height = 5, width = 5)
                ggsave(paste0(outfile_dir, "meox1__BarStrip_Isolated_Zones_DatasetLevel_start_", start, "_lin", lin, ".", col_name, ".pdf"), plot = plt_dataset_zones + NoAxes() + theme(plot.title = element_blank()), height = 4, width = 4)

                # hmVEC lineage 1 only (zones 1, 2, 3)
                # 1. hmVEC lineage 1
                if (start == "hmVEC" && lin == 1) {
                    cat("Generating custom 3-zone percentage plots for hmVEC lineage 1...\n")

                    hmvec_3zones <- master_zones_df %>%
                        filter(Zone %in% c("Transition Zone 1", "Transition Zone 2", "Transition Zone 3")) %>%
                        group_by(Sample, Zone) %>% tally(name = "Cell_Count") %>% ungroup()

                    hmvec_3zones <- hmvec_3zones %>%
                        group_by(Sample) %>%
                        mutate(Sub_Total = sum(Cell_Count), Proportion = (Cell_Count / Sub_Total) * 100) %>% ungroup()

                    hmvec_3zones$Zone <- factor(hmvec_3zones$Zone, levels = c("Transition Zone 1", "Transition Zone 2", "Transition Zone 3"))

                    hmvec_3zones <- hmvec_3zones %>%
                        arrange(Sample, desc(Zone)) %>%
                        group_by(Sample) %>% mutate(label_y_prop = cumsum(Proportion) - (0.5 * Proportion)) %>% ungroup()

                    # override custom color mapping for this specific plot
                    hmvec_colors <- c(
                        "Transition Zone 1" = "#D32F2F", 
                        "Transition Zone 2" = "#F57C00", 
                        "Transition Zone 3" = "#7B1FA2"  
                    )
                    # accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")

                    plt_hmvec_stacked <- ggplot(hmvec_3zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(width = 0.5, color = "black", linewidth = 0.5) +
                        geom_text_repel(aes(y = label_y_prop, label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), "")), nudge_x = 0.4, direction = "y", hjust = 0, segment.size = 0.3) +
                        scale_fill_manual(values = hmvec_colors) +
                        coord_cartesian(ylim = c(0, 100)) +
                        labs(title = "hmVEC Lin 1 - Zones 1-3 Composition %", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_hmVEC_Lin1_Zones123_Stacked_Proportion.", col_name, ".pdf"), plot = plt_hmvec_stacked, height = 5, width = 6.5)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_hmVEC_Lin1_Zones123_Stacked_Proportion.", col_name, ".pdf"), plot = plt_hmvec_stacked + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                    max_hmvec_prop <- max(hmvec_3zones$Proportion, na.rm = TRUE)
                    plt_hmvec_dodged <- ggplot(hmvec_3zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
                        geom_text(aes(label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), ""), y = Proportion + (max_hmvec_prop * 0.02)), position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
                        scale_fill_manual(values = hmvec_colors) +
                        coord_cartesian(ylim = c(0, min(100, max_hmvec_prop * 1.15))) +
                        labs(title = "hmVEC Lin 1 - Zones 1-3 % Comparison", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_hmVEC_Lin1_Zones123_Dodged_Proportion.", col_name, ".pdf"), plot = plt_hmvec_dodged, height = 5, width = 7)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_hmVEC_Lin1_Zones123_Dodged_Proportion.", col_name, ".pdf"), plot = plt_hmvec_dodged + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)

                    rm(hmvec_3zones, plt_hmvec_stacked, plt_hmvec_dodged)
                }

                # mVEC lin1, reuse previous code
                if (start == "mVEC" && lin == 1) {
                    cat("Generating custom 2-zone percentage plots for mVEC lineage 1...\n")

                    mvec_2zones <- master_zones_df %>%
                        filter(Zone %in% c("Transition Zone 1", "Transition Zone 2")) %>%
                        group_by(Sample, Zone) %>% tally(name = "Cell_Count") %>% ungroup()

                    mvec_2zones <- mvec_2zones %>%
                        group_by(Sample) %>%
                        mutate(Sub_Total = sum(Cell_Count), Proportion = (Cell_Count / Sub_Total) * 100) %>% ungroup()

                    mvec_2zones$Zone <- factor(mvec_2zones$Zone, levels = c("Transition Zone 1", "Transition Zone 2"))

                    mvec_2zones <- mvec_2zones %>%
                        arrange(Sample, desc(Zone)) %>%
                        group_by(Sample) %>% mutate(label_y_prop = cumsum(Proportion) - (0.5 * Proportion)) %>% ungroup()

                    # override: custom color mapping for this specific plot
                    mvec_colors <- c(
                        "Transition Zone 1" = "#D32F2F", 
                        "Transition Zone 2" = "#F57C00" 
                    )
                    # accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")

                    plt_mvec_stacked <- ggplot(mvec_2zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(width = 0.5, color = "black", linewidth = 0.5) +
                        geom_text_repel(aes(y = label_y_prop, label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), "")), nudge_x = 0.4, direction = "y", hjust = 0, segment.size = 0.3) +
                        scale_fill_manual(values = mvec_colors) +
                        coord_cartesian(ylim = c(0, 100)) +
                        labs(title = "mVEC Lin 1 - Zones 1-2 Composition %", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_mVEC_Lin1_Zones12_Stacked_Proportion.", col_name, ".pdf"), plot = plt_mvec_stacked, height = 5, width = 6.5)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_mVEC_Lin1_Zones12_Stacked_Proportion.", col_name, ".pdf"), plot = plt_mvec_stacked + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                    max_mvec_prop <- max(mvec_2zones$Proportion, na.rm = TRUE)
                    plt_mvec_dodged <- ggplot(mvec_2zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
                        geom_text(aes(label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), ""), y = Proportion + (max_mvec_prop * 0.02)), position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
                        scale_fill_manual(values = mvec_colors) +
                        coord_cartesian(ylim = c(0, min(100, max_mvec_prop * 1.15))) +
                        labs(title = "mVEC Lin 1 - Zones 1-2 % Comparison", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_mVEC_Lin1_Zones12_Dodged_Proportion.", col_name, ".pdf"), plot = plt_mvec_dodged, height = 5, width = 7)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_mVEC_Lin1_Zones12_Dodged_Proportion.", col_name, ".pdf"), plot = plt_mvec_dodged + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)

                    rm(mvec_2zones, plt_mvec_stacked, plt_mvec_dodged)
                }                

                # mVEC lin1, reuse previous code
                if (start == "mVEC" && lin == 2) {
                    cat("Generating custom 2-zone percentage plots for mVEC lineage 2...\n")

                    mvec_2zones <- master_zones_df %>%
                        filter(Zone %in% c("Transition Zone 1", "Transition Zone 2")) %>%
                        group_by(Sample, Zone) %>% tally(name = "Cell_Count") %>% ungroup()

                    mvec_2zones <- mvec_2zones %>%
                        group_by(Sample) %>%
                        mutate(Sub_Total = sum(Cell_Count), Proportion = (Cell_Count / Sub_Total) * 100) %>% ungroup()

                    mvec_2zones$Zone <- factor(mvec_2zones$Zone, levels = c("Transition Zone 1", "Transition Zone 2"))

                    mvec_2zones <- mvec_2zones %>%
                        arrange(Sample, desc(Zone)) %>%
                        group_by(Sample) %>% mutate(label_y_prop = cumsum(Proportion) - (0.5 * Proportion)) %>% ungroup()

                    # override custom color mapping for this specific plot
                    mvec_colors <- c(
                        "Transition Zone 1" = "#D32F2F", 
                        "Transition Zone 2" = "#F57C00" 
                    )
                    # accent_palette <- c("#D32F2F", "#F57C00", "#7B1FA2", "#388E3C")

                    plt_mvec_stacked <- ggplot(mvec_2zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(width = 0.5, color = "black", linewidth = 0.5) +
                        geom_text_repel(aes(y = label_y_prop, label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), "")), nudge_x = 0.4, direction = "y", hjust = 0, segment.size = 0.3) +
                        scale_fill_manual(values = mvec_colors) +
                        coord_cartesian(ylim = c(0, 100)) +
                        labs(title = "mVEC Lin 2 - Zones 1-2 Composition %", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_mVEC_Lin2_Zones12_Stacked_Proportion.", col_name, ".pdf"), plot = plt_mvec_stacked, height = 5, width = 6.5)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_mVEC_Lin2_Zones12_Stacked_Proportion.", col_name, ".pdf"), plot = plt_mvec_stacked + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 3)

                    max_mvec_prop <- max(mvec_2zones$Proportion, na.rm = TRUE)
                    plt_mvec_dodged <- ggplot(mvec_2zones, aes(x = Sample, y = Proportion, fill = Zone)) +
                        geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
                        geom_text(aes(label = ifelse(Proportion > 0, paste0(round(Proportion, 1), "%"), ""), y = Proportion + (max_mvec_prop * 0.02)), position = position_dodge(width = 0.8), vjust = 0, size = 3.5) +
                        scale_fill_manual(values = mvec_colors) +
                        coord_cartesian(ylim = c(0, min(100, max_mvec_prop * 1.15))) +
                        labs(title = "mVEC Lin 2 - Zones 1-2 % Comparison", x = NULL, y = "% of Selected Zones", fill = "Trajectory Regions") +
                        theme_classic() + theme(axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

                    ggsave(paste0(outfile_dir, "meox1__BarLabel_Custom_mVEC_Lin2_Zones12_Dodged_Proportion.", col_name, ".pdf"), plot = plt_mvec_dodged, height = 5, width = 7)
                    ggsave(paste0(outfile_dir, "meox1__BarStrip_Custom_mVEC_Lin2_Zones12_Dodged_Proportion.", col_name, ".pdf"), plot = plt_mvec_dodged + NoAxes() + NoLegend() + theme(plot.title = element_blank()), height = 4, width = 4)

                    rm(mvec_2zones, plt_mvec_stacked, plt_mvec_dodged)
                }

                # master_zones_df$Sample holds exactly ONE level per run_slingshot() call:
                # "meox1 -/- mutant" when col_name == "S20200", or "WildType" when col_name == "S20201"
                has_wt  <- "WildType" %in% unique(master_zones_df$Sample)
                has_mut <- "meox1 -/- mutant" %in% unique(master_zones_df$Sample)

                if (start == "hmVEC" && lin == 1) {
                    if (has_wt) {
                        plot_custom_transition_zones(
                            master_df      = master_zones_df,
                            target_zones   = c(1, 2, 3),
                            custom_colors  = c(
                                "Transition Zone 1" = "#D32F2F",
                                "Transition Zone 2" = "#F57C00",
                                "Transition Zone 3" = "#7B1FA2"
                            ),
                            title_prefix   = "hmVEC Lin 1 (WT Only)",
                            file_prefix    = "Custom_hmVEC_Lin1_WT",
                            col_name       = col_name,
                            outfile_dir    = outfile_dir
                        )
                    }
                    if (has_mut) {
                        plot_custom_transition_zones(
                            master_df      = master_zones_df,
                            target_zones   = c(1, 2),
                            custom_colors  = c(
                                "Transition Zone 1" = "#D32F2F",
                                "Transition Zone 2" = "#F57C00"
                            ),
                            title_prefix   = "hmVEC Lin 1 (meox -/- Only)",
                            file_prefix    = "Custom_hmVEC_Lin1_Mut",
                            col_name       = col_name,
                            outfile_dir    = outfile_dir
                        )
                    }
                }

                # 1. hmVEC Lineage 1 (Zones 1, 2, 3)
                # if (start == "hmVEC" && lin == 1) {
                #     plot_custom_transition_zones(
                #         master_df      = master_zones_df,
                #         target_samples = c("meox1 -/- mutant", "WildType"), # Both samples
                #         target_zones   = c(1, 2, 3),                        # Specific zones
                #         custom_colors  = c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "lightgrey", "Transition Zone 3" = "lightgrey"),
                #         title_prefix   = "hmVEC Lin 1",
                #         file_prefix    = "Custom_hmVEC_Lin1",
                #         outfile_dir    = outfile_dir,
                #         col_name       = col_name
                #     )
                # }

                # 2. mVEC Lineage 2 (Zones 1 & 2)
                # if (start == "mVEC" && lin == 2) {
                #     plot_custom_transition_zones(
                #         master_df      = master_zones_df,
                #         target_samples = c("meox1 -/- mutant", "WildType"), 
                #         target_zones   = c(1, 2),                           
                #         custom_colors  = c("Transition Zone 1" = "#D32F2F", "Transition Zone 2" = "#F57C00"),
                #         title_prefix   = "mVEC Lin 2",
                #         file_prefix    = "Custom_mVEC_Lin2",
                #         outfile_dir    = outfile_dir,
                #         col_name       = col_name
                #     )
                # }

                # Clean up memory cleanly across all iterations
                suppressWarnings(rm(
                    master_zones_df, sample_baselines, indiv_summary, isolated_zone, plt_isolated, 
                    plt_stacked_indiv, zone_df, combined_z, p_z_abs, p_z_prop, combined_summary, 
                    plt_combined, zones_only_summary, plt_stacked_zones_only, plt_stacked_zones_prop, 
                    plt_dodged_zones_only, plt_dodged_zones_prop, dataset_zones_only, plt_dataset_zones, 
                    zone_colors
                ))

            }
            rm(transition_zones_only, transition_zones_df, zone_names)
        }
    }
}

run_featureplot <- function(data, features, outfile_dir) {
    for (i in features) {
        plt <- FeaturePlot(data, i, pt.size=1)
        plt_name <- paste0(outfile_dir, "/", i, ".pdf")
        ggsave(plt_name, plt, height=7, width=7)

        plt <- plt + NoAxes() + NoLegend() + theme(plot.title = element_blank())
        plt_name <- paste0(outfile_dir, "/", i, "_STRIP.pdf")
        ggsave(plt_name, plt, height=5, width=5)
    }
}

col_order <- c(
    "LEC", 
    "pre_muLEC", 
    "hmVEC",
    "mVEC"
)
# "pre_muLEC" = "#a1d99b",
col_map <- c(
    "LEC"       = "#41ae76",
    "pre_muLEC" = "#fec44f",
    "hmVEC"     = "#9ecae1",
    "mVEC"      = "#6baed6"
)

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/slingshot/"

data_D10051 <- qs_read(infile_path)
data_S20200 <- data_D10051[, grepl("_1$", colnames(data_D10051))]
data_S20201 <- data_D10051[, grepl("_2$", colnames(data_D10051))]

# run_slingshot(data_D10051, col_map, col_order, col_name="D10051")
# run_slingshot(data_S20200, col_map, col_order, col_name="S20200")
# run_slingshot(data_S20201, col_map, col_order, col_name="S20201")

features <- c("mki67", "pcna")
# run_featureplot(data_D10051, features, outfile_dir)

run_wt_anchored_slingshot(data_D10051, col_map, col_order, outfile_dir)