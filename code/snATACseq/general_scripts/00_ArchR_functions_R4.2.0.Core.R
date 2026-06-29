#! /usr/bin/Rscript

#	First Authors: Tyrone Chen, Oliver Yu, Lizzie Mason and Michelle Meier
#	Depends on ArchR R environment ....


#   Depends on R/4.2.0.Core

## ===========================================================================
##                             LOAD LIBRARIES
## ===========================================================================

require(ArchR)
#require(Gviz)

## ===========================================================================
##                             DEFINE FUNCTIONS
## ===========================================================================

#	Specifying the ggplot color scheme
ggplotColours <- function(n = 6, h = c(0, 360) + 15){

# Make the same spectrum of colours as GGPLOT

  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

## ===========================================================================
##                        Manipulating ArchRProject
## ===========================================================================

# Copy an existing ArchRProject to a new path
copy_ArchR_Project 				<- 	function(oldProjectPath,
							  				newProjectPath
							   				){

#		Arguments:
#		

#		Example:	project 	<-			copy_ArchR_Project(oldProjectPath = "/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_01",
#																newProjectPath = "/hogan_lab/Hogan_Lab_People/.....")
#							  									
#======================================================================================================

# load old ArchR Project
project_old				<-		loadArchRProject(path = oldProjectPath, force = FALSE, showLogo = FALSE)

# save project to new location
project_new 			<-	 	saveArchRProject(ArchRProj = project_old,
  												outputDirectory = newProjectPath,
  												overwrite = TRUE,
    											load = TRUE
  												)

return(project_new)
}


#======================================================================================================
#	Differential Accessibility Analysis
#======================================================================================================

#----- DAP analysis
dap_analysis			<- function(project,
		 								groupBy,
		 								useGroup,
										bgdGroup,
										exportFilePath = NULL,
										exportFileName = NULL,
										chrStatus = NULL,
										FDRFilter = TRUE,
										FDRThreshold = 0.05,
										RawPvalFilter = FALSE,
										RawPvalThreshold = 0.05,
										Log2FCFilter = FALSE,
										Log2FCThreshold = 1.25){
#		Arguments:
#		project = ArchrProject
#		groupBy = The column name in cellColData used for grouping cells together for marker feature identification (e.g. "ClustersPhenotype").
#		useGroup = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData (e.g. "LEC_wt").
#		bgdGroup = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData to be used for background calculations in marker feature identification.
#		exportFilePath = An alternative directory path to save DAP output to. Default is /DAP_analysis in outputDirectory of the ArchRProject. (DO NOT ADD "/" AT THE END)
#		exportFileName = A character vector that is used to label the output file (e.g. "Level_03_LEC_vs_VEC"). (no need for ".txt" or "open/closed")
#		chrStatus = either "open" or "closed", which is used to subset the peak set based on the Fold Change. Default will produce full table.
#		FDRFilter = A boolean indicating whether to use FDR to subset the peak set. Default is TRUE.
#		FDRThreshold = FDR threshold for subsetting the table. Default is 0.05.
#		RawPvalFilter = A boolean indicating whether to use raw Pval to subset the peak set. Default is FALSE
#		RawPvalThreshold = Raw P value threshold for subsetting the table. Default is 0.05.
#		Log2FCFilter = A boolean indicating whether to use raw Pval to subset the peak set. Default is FALSE
#		Log2FCThreshold = Log2FC threshold for subsetting the table. Default is 1.25.
#
#		Example: dap_LEC_vs_VEC_opn <- dap_analysis(project = ArchRProject, groupBy="Phenotype_Genoype", useGroup="LEC_wt", bgdGroup= "VEC_wt", exportFileName="Level_03_LECwt_vs_VECwt", chrStatus="open")

#======================================================================================================

# make output directory for DAP analysis
if (!dir.exists(file.path(getOutputDirectory(project),"DAP_analysis")))	{
	dir.create(file.path(getOutputDirectory(project),"DAP_analysis"))
}

# check chrStatus
if (!is.null(chrStatus)) {
	if (chrStatus %ni% c("open", "closed")) {
		chrStatus 			<-	NULL
		warning('chrStatus must be "open" / "closed" / NULL. Using default: All DAPs will be printed.',"\n")
	}
} 

# make raw dap table
dap_raw 					<- getMarkerFeatures(ArchRProj= project, 
                    					         useMatrix= "PeakMatrix", 
                    					         groupBy= groupBy,
                    					         testMethod = "wilcoxon",
                    					         bias = c("TSSEnrichment", "log10(nFrags)"),
                    					         useGroups = useGroup,   
                    					         bgdGroups = bgdGroup)

# extract peak info
cat("extracting peaks info ...\n")

dap_table					<- 	data.frame(assays(dap_raw)$Log2FC,
											assays(dap_raw)$FDR,
										  	assays(dap_raw)$Pval,
										  	assays(dap_raw)$AUC,
                    				      	assays(dap_raw)$Mean,
                  					      	assays(dap_raw)$MeanBGD,
                  					      	assays(dap_raw)$MeanDiff,
                  					      	rowData(dap_raw)$seqnames,
                  					     	rowData(dap_raw)$start,
                  					    	rowData(dap_raw)$end)

colnames(dap_table) 		<- 	c("Log2FC",
                    		    	"FDR",
                    		     	"RawPval",
                    		     	"AUC",
                    		     	"Group_01_mean",
                    		     	"Background_mean",
                    		     	"MeanDiff",
                    		     	"Chr",
                    		     	"Start",
                    		     	"End")

for (i in 1:9)	{
	dap_table$Chr[dap_table$Chr == paste0("chr",i)] <- paste0("chr",paste0("0",i))
}

dap_table					<-	dap_table[order(dap_table$Chr,dap_table$Start),]

# add peak info from peak sets
cat("Getting peak sets ...\n")
peak_set 					<-	data.frame(getPeakSet(ArchRProj=project))

peak_set					<-	peak_set[order(peak_set$seqnames,peak_set$start),]

dap_table 					<-	cbind(dap_table,
									 peak_set$score,
           					         peak_set$distToGeneStart,
           					         peak_set$nearestGene,
           					         peak_set$peakType,
           					         peak_set$distToTSS,
           					         peak_set$nearestTSS,
           					         peak_set$GC)

colnames(dap_table) 		<- c("Log2FC",
                    		     "FDR",
                    		     "RawPval",
                    		     "AUC",
                    		     "Group_01_mean",
                    		     "Background_mean",
                    		     "MeanDiff",
                    		     "Chr",
                    		     "Start",
                    		     "End",
                    		     "PeakScore",
                    		     "distToGeneStart",
                    		     "nearestGene",
                    		     "peakType",
                    		     "distToTSS",
           					     "nearestTSS",
           					     "GC")

dap_table$nearestGene 		<- 	gsub("NA_","",dap_table$nearestGene)

dap_table$Status_Group_01[dap_table$Log2FC > 0] <- "open"
dap_table$Status_Group_01[dap_table$Log2FC < 0] <- "closed"
dap_table$Status_Group_01[dap_table$Log2FC == 0] <- "NA"

dap_table$GeneDist 			<- 	paste(dap_table$nearestGene, dap_table$distToGeneStart, sep="_")

dap_table 					<- 	dap_table[!is.na(dap_table$FDR),]

dap_table$PeakID			<-	paste(dap_table$Chr,paste(dap_table$Start, dap_table$End, sep = "-"),sep = ":")

#======================================================================================================

# file name (filter)
if (FDRFilter == FALSE & RawPvalFilter == FALSE & Log2FCFilter == FALSE) {
	filename_filter			<- "raw"
}	else {
	filename_filter				<- "filteredBy"
}

# remove uncharacterised genes
#remove_genes 				<- dap_table$nearestGene[grep("^LOC", dap_table$nearestGene)]
#dap_table					<- dap_table[!dap_table$nearestGene %in% remove_genes,]

# subsetting with FDR
if (FDRFilter != FALSE) {
	cat("subsetting peaks with FDR < ", FDRThreshold,"\n")
	dap_table			<- subset(dap_table, FDR < FDRThreshold)
	filename_filter		<- paste0(filename_filter,"_FDR")
}

# subsetting with raw Pval
if (RawPvalFilter != FALSE) {
	cat("subsetting peaks with raw Pval < ", RawPvalThreshold,"\n")
	dap_table <- subset(dap_table, RawPval < RawPvalThreshold)
	filename_filter		<- paste0(filename_filter,"_rawPval")
}

# subsetting with Log2FC
if (Log2FCFilter != FALSE) {
	cat("subsetting peaks with abs(Log2FC) > ", Log2FCThreshold,"\n")
	dap_table <- subset(dap_table, abs(Log2FC) > Log2FCThreshold)
	filename_filter		<- paste0(filename_filter,"_Log2FC")
}

#======================================================================================================
# file name (chrStatus)
if (is.null(chrStatus)) {
	filename_chrStatus		<-	"all_DAPs"
}	else if (chrStatus == "closed") {
	filename_chrStatus		<-	"closed_DAPs"
}	else if	(chrStatus == "open") {
	filename_chrStatus		<-	"open_DAPs"
}

# default export file path
if (is.null(exportFilePath)) {
	exportFilePath <- file.path(getOutputDirectory(project),"DAP_analysis")
}	else if (!dir.exists(exportFilePath))	{
	warning("exportFilePath does not exist. Using default output path.", "\n")
	exportFilePath <- file.path(getOutputDirectory(project),"DAP_analysis")
}

exportFilePath		<-	paste0(exportFilePath,"/")

# default export file name
if (is.null(exportFileName)) {
	exportFileName	<-	paste(format(Sys.time(), "%y%m%d"),useGroup,"vs",bgdGroup,filename_chrStatus,filename_filter, sep = "_")
}

#======================================================================================================

# subsetting with chrStatus
if (is.null(chrStatus) == TRUE){
	
	dap_table <- dap_table[with(dap_table, order(Status_Group_01, abs(Log2FC), decreasing = TRUE)),]

	cat(nrow(dap_table[dap_table$Log2FC > 0,]),"peaks are more open in", useGroup, "compared to", bgdGroup,"\n")
	cat(nrow(dap_table[dap_table$Log2FC < 0,]),"peaks are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dap_table,
		file = paste(exportFilePath, exportFileName, ".txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, exportFileName, ".txt", sep=""), "\n")

} else if (chrStatus == "open") {
	dap_table <- subset(dap_table, Log2FC > 0)
	dap_table <- dap_table[order(-dap_table$Log2FC),]

	cat(nrow(dap_table),"peaks are more open in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dap_table,
		file = paste(exportFilePath, exportFileName, ".txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, exportFileName, ".txt", sep=""), "\n")

} else if (chrStatus == "closed") {
	dap_table <- subset(dap_table, Log2FC < 0)
	dap_table <- dap_table[order(dap_table$Log2FC),]
			
	cat(nrow(dap_table),"peaks are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dap_table,
		file = paste(exportFilePath, exportFileName, ".txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, exportFileName, ".txt", sep=""), "\n")

} else {
	dap_table <- dap_table[order(-dap_table$Log2FC),]

	cat(nrow(dap_table[dap_table$Log2FC > 0,]),"peaks are more open in", useGroup, "compared to", bgdGroup,"\n")
	cat(nrow(dap_table[dap_table$Log2FC < 0,]),"peaks are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dap_table,
		file = paste(exportFilePath, exportFileName, ".txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, exportFileName, ".txt", sep=""), "\n")
}

for (i in 1:9)	{
	dap_table$Chr[dap_table$Chr == paste0("chr",paste0("0",i))]	<-	paste0("chr",i)
}

rownames(dap_table) 		<-	paste(dap_table$Chr,paste(dap_table$Start, dap_table$End, sep = "-"),sep = ":")

return(dap_table)
}

#======================================================================================================
#======================================================================================================

#----- DAG analysis
dag_analysis 				<- function(project,
		 								groupBy,
		 								useGroup,
										bgdGroup,
										exportFilePath = NULL,
										exportFileName,
										chrStatus = NULL,
										FDRFilter = TRUE,
										FDRThreshold = 0.05,
										RawPvalFilter = FALSE,
										RawPvalThreshold = 0.05,
										Log2FCFilter = FALSE,
										Log2FCThreshold = 1.25){

#		Arguments:
#		project = ArchrProject
#		groupBy = The name of the column in cellColData to use for grouping cells together for marker feature identification (e.g. "ClustersPhenotype").
#		useGroup = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData (e.g. "LEC_wt").
#		bgdGroup = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData to be used for background calculations in marker feature identification.
#		exportFilePath = A directory path to save DAG output to. Default is outputDirectory of the ArchRProject.(DO NOT ADD "/" AT THE END)
#		exportFileName = A character vector that is used to label the output file (e.g. "Level_03_LEC_vs_VEC"). (no need for ".txt" or "open/closed")
#		chrStatus = either "open" or "closed", which is used to subset the table based on the Fold Change. Default will produce full table.
#		FDRFilter = A boolean indicating whether to use FDR to subset the gene table. Default is TRUE.
#		FDRThreshold = P value threshold for subsetting the table. Default is 0.05.
#		RawPvalFilter = A boolean indicating whether to use raw Pval to subset the gene table. Default is FALSE
#		RawPvalThreshold = Raw P value threshold for subsetting the table. Default is 0.05.
#		Log2FCFilter = A boolean indicating whether to use raw Pval to subset the gene table. Default is FALSE
#		Log2FCThreshold = Log2FC threshold for subsetting the table. Default is 1.25.
#
#
#		Example: dag_LEC_vs_VEC_opn <- dag_analysis(project = ArchRProject, groupBy="ClustersPhenotype", useGroup="LEC_wt", bgdGroup= "VEC_wt", exportFileName="Level_03_LECwt_vs_VECwt", chrStatus="open")

#======================================================================================================

# make output directory for DAP analysis
if (!dir.exists(file.path(getOutputDirectory(project),"DAG_analysis")))	{
	dir.create(file.path(getOutputDirectory(project),"DAG_analysis"))
}

# check chrStatus
if (!is.null(chrStatus)) {
	if (chrStatus %ni% c("open", "closed")) {
		chrStatus 			<-	NULL
		warning('chrStatus must be "open" / "closed" / NULL. Using default: All DAPs will be printed.',"\n")
	}
} 

# making raw dag table
dag_raw 					<- getMarkerFeatures(ArchRProj= project,
												useMatrix= "GeneScoreMatrix", 
												groupBy= groupBy,
												testMethod = "wilcoxon",
												bias = c("TSSEnrichment", "log10(nFrags)"),
												useGroups = useGroup,   
												bgdGroups = bgdGroup)

# extract info
dag_table 					<- data.frame(gsub("NA_","",make.names(rowData(dag_raw)$name, unique = TRUE)),
											assays(dag_raw)$Log2FC,
											assays(dag_raw)$FDR,
											assays(dag_raw)$Pval,
											assays(dag_raw)$AUC,
											assays(dag_raw)$Mean,
											assays(dag_raw)$MeanBGD,
											assays(dag_raw)$MeanDiff)

colnames(dag_table) 		<- c("Gene",
								 "Log2FC",
								 "FDR",
								 "RawPval",
								 "AUC",
								 "Group_01_mean",
								 "Background_mean",
								 "MeanDiff")

# remove unannotated genes
#remove_genes 				<- dag_table$Gene[grep("^LOC", dag_table$Gene)]
#dag_table					<- dag_table[!dag_table$Gene %in% remove_genes,]

# subsetting with raw Pval
if (RawPvalFilter != FALSE) {
cat("subsetting genes with raw Pval < ", RawPvalThreshold,"\n")
dag_table <- subset(dag_table, RawPval < RawPvalThreshold)
}

# subsetting with FDR
if (FDRFilter != FALSE) {
cat("subsetting genes with FDR < ", FDRThreshold,"\n")
dag_table <- subset(dag_table, FDR < FDRThreshold)
}

# subsetting with Log2FC
if (Log2FCFilter != FALSE) {
cat("subsetting genes with abs(Log2FC) > ", Log2FCThreshold,"\n")
dag_table <- subset(dag_table, abs(Log2FC) > Log2FCThreshold)
}

# default export file path
if (is.null(exportFilePath) == TRUE) {
	exportFilePath <- getOutputDirectory(project)
}

# filter and subset
if (is.null(chrStatus) == TRUE){

	dag_table$Status_Group_01[dag_table$Log2FC > 0] <- "open"
	dag_table$Status_Group_01[dag_table$Log2FC < 0] <- "closed"
	dag_table$Status_Group_01[dag_table$Log2FC == 0] <- "NA"
	dag_table <- dag_table[with(dag_table, order(Status_Group_01, abs(Log2FC), decreasing = TRUE)),]

	cat(nrow(dag_table[dag_table$Log2FC > 0,]),"genes are more open in", useGroup, "compared to", bgdGroup,"\n")
	cat(nrow(dag_table[dag_table$Log2FC < 0,]),"genes are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dag_table,
		file = paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_complete_table.txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_complete_table.txt", sep=""), "\n")

} else if (chrStatus == "open") {
	dag_table <- subset(dag_table, Log2FC > 0)
	dag_table <- dag_table[order(-dag_table$Log2FC),]

	cat(nrow(dag_table),"genes are more open in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dag_table,
		file = paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_open_table.txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_open_table.txt", sep=""), "\n")

} else if (chrStatus == "closed") {
	dag_table <- subset(dag_table, Log2FC < 0)
	dag_table <- dag_table[order(dag_table$Log2FC),]
			
	cat(nrow(dag_table),"genes are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dag_table,
		file = paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_closed_table.txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_closed_table.txt", sep=""), "\n")

} else {
	dag_table$Status_Group_01[dag_table$Log2FC > 0] <- "open"
	dag_table$Status_Group_01[dag_table$Log2FC < 0] <- "closed"
	dag_table <- dag_table[with(dag_table, order(Status_Group_01, abs(Log2FC), decreasing = TRUE)),]

	cat(nrow(dag_table[dag_table$Log2FC > 0,]),"genes are more open in", useGroup, "compared to", bgdGroup,"\n")
	cat(nrow(dag_table[dag_table$Log2FC < 0,]),"genes are more closed in", useGroup, "compared to", bgdGroup,"\n")

	write.table(dag_table,
		file = paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_complete_table.txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

	cat("Output File:", paste(exportFilePath, "/DAG_analysis/", exportFileName, "_DAG_complete_table.txt", sep=""), "\n")
}

return(dag_table)
}

#----- Get Marker Genes
makeMarkerGeneTable			<- function(project,
		 								groupBy,
										exportFilePath = NULL,
										exportFileName = NULL,
										FDRFilter = TRUE,
										FDRThreshold = 0.05,
										RawPvalFilter = FALSE,
										RawPvalThreshold = 0.05,
										Log2FCFilter = FALSE,
										Log2FCThreshold = 1.25){
#		Arguments:
#		project = ArchrProject
#		groupBy = The column name in cellColData used for grouping cells together for marker feature identification (e.g. "ClustersPhenotype").
#		exportFilePath = An alternative directory path to save DAP output to. Default is /DAP_analysis in outputDirectory of the ArchRProject. (DO NOT ADD "/" AT THE END)
#		exportFileName = An alternative file name that is used for the output file (e.g. "Level_03_LEC_vs_VEC"). Default is "DATE_MarkerGene_By_Phenotype_filteredBy_FDR"
#		FDRFilter = A boolean indicating whether to use FDR to subset the peak set. Default is TRUE.
#		FDRThreshold = FDR threshold for subsetting the table. Default is 0.05.
#		RawPvalFilter = A boolean indicating whether to use raw Pval to subset the peak set. Default is FALSE
#		RawPvalThreshold = Raw P value threshold for subsetting the table. Default is 0.05.
#		Log2FCFilter = A boolean indicating whether to use raw Pval to subset the peak set. Default is FALSE
#		Log2FCThreshold = Log2FC threshold for subsetting the table. Default is 1.25.
#
#		Example: markerGenes_ByPhenotype <- makeMarkerGeneTable(project = ArchRProject, groupBy="Phenotype_Genoype",FDRFilter = TRUE,FDRThreshold = 0.05)

#======================================================================================================

# make output directory for MarkerGenesTable
if (!dir.exists(file.path(getOutputDirectory(project),"MarkerGenes")))	{
	dir.create(file.path(getOutputDirectory(project),"MarkerGenes"))
}

# get gene score data
markersGS					<-	getMarkerFeatures(
    								ArchRProj = project, 
    								useMatrix = "GeneScoreMatrix", 
    								groupBy = groupBy,
    								bias = c("TSSEnrichment", "log10(nFrags)"),
    								testMethod = "wilcoxon"
									)

# extract gene score data
cat("extracting gene score data ...\n")

GS_table					<- 	data.frame(Group = character(),
											GeneName = character(),
											Chr = character(),
                  					     	Start = integer(),
                  					    	End = integer(),
											Log2FC = double(),
											FDR = double(),
										  	Pval = double(),
										  	AUC = double(),
                    				      	Mean = double(),
                  					      	MeanBGD = double(),
                  					      	MeanDiff = double()
                  					      	)

for (groups in colnames(markersGS)) {

	groupGS_table			<-	data.frame(rowData(markersGS)$name,
                  					      	rowData(markersGS)$seqnames,
                  					     	rowData(markersGS)$start,
                  					    	rowData(markersGS)$end,
                  					    	assays(markersGS)$Log2FC[[groups]],
											assays(markersGS)$FDR[[groups]],
										  	assays(markersGS)$Pval[[groups]],
										  	assays(markersGS)$AUC[[groups]],
                    				      	assays(markersGS)$Mean[[groups]],
                  					      	assays(markersGS)$MeanBGD[[groups]],
                  					      	assays(markersGS)$MeanDiff[[groups]]
                  					      	)

	groupGS_table			<-	cbind(Group = groups, groupGS_table)

	colnames(groupGS_table) 		<- c(
									"Group",
                    			    "GeneName",
                    			    "Chr",
    	                		    "Start",
        	            		    "End",
									"Log2FC",
    	                		    "FDR",
        	            		    "RawPval",
            	        		    "AUC",
                	    		    "Group_mean",
                    			    "Background_mean",
                    			    "MeanDiff"
            	        		    )

	groupGS_table			<-	groupGS_table[order(groupGS_table$FDR),]

	GS_table 				<-	rbind(GS_table,groupGS_table)
}

#======================================================================================================

# file name (filter)
if (FDRFilter == FALSE & RawPvalFilter == FALSE & Log2FCFilter == FALSE) {
	filename_filter			<- "raw"
}	else {
	filename_filter				<- "filteredBy"
}

# subsetting with FDR
if (FDRFilter != FALSE) {
	cat("filtering MarkerGenes with FDR < ", FDRThreshold,"\n")
	GS_table			<- subset(GS_table, FDR < FDRThreshold)
	filename_filter		<- paste0(filename_filter,"_FDR")
}

# subsetting with raw Pval
if (RawPvalFilter != FALSE) {
	cat("filtering MarkerGenes with raw Pval < ", RawPvalThreshold,"\n")
	GS_table 			<- subset(GS_table, RawPval < RawPvalThreshold)
	filename_filter		<- paste0(filename_filter,"_rawPval")
}

# subsetting with Log2FC
if (Log2FCFilter != FALSE) {
	cat("filtering MarkerGenes with abs(Log2FC) > ", Log2FCThreshold,"\n")
	GS_table 			<- subset(GS_table, abs(Log2FC) > Log2FCThreshold)
	filename_filter		<- paste0(filename_filter,"_Log2FC")
}

#======================================================================================================

# default export file path
if (is.null(exportFilePath)) {
	exportFilePath <- file.path(getOutputDirectory(project),"MarkerGenes")
}	else if (!dir.exists(exportFilePath))	{
	warning("exportFilePath does not exist. Using default output path.", "\n")
	exportFilePath <- file.path(getOutputDirectory(project),"MarkerGenes")
}

exportFilePath		<-	paste0(exportFilePath,"/")

# default export file name
if (is.null(exportFileName)) {
	exportFileName	<-	paste(format(Sys.time(), "%y%m%d"),"MarkerGenes_By",groupBy,filename_filter, sep = "_")
}

#======================================================================================================

# exporting marker genes table

MarkerGenes 		<- subset(GS_table, Log2FC > 0)

write.table(MarkerGenes,
		file = paste(exportFilePath, exportFileName, ".txt", sep=""),
		quote = FALSE,
		sep = "\t",
		row.names = FALSE)

cat("Output File:", paste(exportFilePath, exportFileName, ".txt", sep=""), "\n")

return(MarkerGenes)
}

#======================================================================================================
#	VIZUALISATION
#======================================================================================================

#----- Gene Score Dot Plots
make_dot_plot_ATAC 				<- 	function(project,
							  				groupBy,
							   				useGroups,
							   				groupName,
							   				markerGenes,
							   				exportFileName,
							   				exportFilePath = NULL,
							   				accessibilityColour = c("#d9d9d9","#d7301f"),
							   				plotWidth = 30,
							   				plotHeight = 12,
							   				rotateaxis = FALSE,
							  				scale.scores = TRUE
							   				){

#		Arguments:
#		project = ArchrProject
#		groupBy = The column name in cellColData used for grouping cells together for marker feature identification (e.g. "Level_02_Phenotype").
#		useGroups = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData (e.g. c("LEC","VEC","AEC")).
#		groupName = A character vector that contains axis labels in the dot plot (the order of the group name should be corresponding to the useGroups).
#		markerGenes = A character vector that contains a list of gene names.
#		exportFilePath = An alternative directory path to save the dot plot to. Default is ./Plots in outputDirectory of the ArchRProject.
#		exportFileName = A character vector that is used for the output file name (e.g. exportFileName = "Level_02_marker_genes", the output file name would be "DotPlot_Level_02_marker_genes.pdf") ).
#		accessibilityColour = A character vector that specifies the colour range for accessibility score.
#		plotWidth = A numeric vector that specifies the width of the dot plot. 
#		plotHeight = A numeric vector that specifies the height of the dot plot. 
#		rotateaxis = A boolean vector indicating whether to rotate the axis (default x axis is genes and y axis is groups).
#   scale = A boolean vector indicating whether to scale the gene prediction score across genes. Default is TRUE.

#		Example:	dotPlot 	<-			make_dot_plot_ATAC (project = project,
#							  									groupBy = "Level_02_Phenotype)",
#							   									useGroups = c("LEC","VEC","AEC"),
#							   									groupName = c("Lymphatics", "Veins", "Arterials"),
#							   									exportFileName = "Level_02_marker_genes"
#																)

#======================================================================================================

  # check if all the gene names are correct
  geneNamesError			<-	NULL
  
  for (genes in markerGenes) {
    if (genes %ni% project@geneAnnotation$genes$symbol) {
      geneNamesError		<-	paste(geneNamesError, genes, sep = ",")
    }
  }
  
  if (!is.null(geneNamesError)) {
    geneNamesError			<-	substring(geneNamesError, 2)
    stop(geneNamesError, " are not found in gene annotations.", "\n")
  }
  
  # check if all the useGroups are correct
  groupNamesError			<-	NULL
  
  for (groups in useGroups) {
    if (groups %ni% project@cellColData[[groupBy]]) {
      groupNamesError		<-	paste(groupNamesError, groups, sep = ",")
    }
  }
  
  if (!is.null(groupNamesError)) {
    groupNamesError			<-	substring(groupNamesError, 2)
    stop(groupNamesError, " are not found in ", groupBy, ".", "\n")
  }
  
  
  # === Start function ===
  # get matrix of accessibility score for each cell & gene
  gene_score_matrix 			<- 	getMatrixFromProject(ArchRProj = project, useMatrix = "GeneScoreMatrix")
  matrix <- assays(gene_score_matrix)$GeneScoreMatrix %>% as.matrix() # make sure it is a matrix
  rownames(matrix) <- rowData(gene_score_matrix)$name # use gene names as rownames
  
  # Get relevant metadata as a dataframe
  metadata_df <- colData(gene_score_matrix)[which(colData(gene_score_matrix)[[groupBy]] %in% useGroups), groupBy] %>%
    as.data.frame() %>% # we need ensure this is a dataframe because the line above creates a character vector
    dplyr::rename(., GroupBy = '.') %>% # rename the column to GroupBy
    mutate(., rownames = rownames(colData(gene_score_matrix))[which(colData(gene_score_matrix)[[groupBy]] %in% useGroups)]) #add in the barcodes as a column called rownames
  
  # # Get relevant gene expression
  # filtered_gex <- matrix[markerGenes,] %>% # only keep the genes we are interested in
  #   t() %>% as.data.frame() %>% #transpose matrix and then convert into dataframe because the functions bellow only work on dfs
  #   merge(., metadata_df, by.x = 0, by.y = 'rownames') %>% #merge metadata dataframe into this dataframe so we can access GroupBy information
  #   dplyr::select(., -`Row.names`) %>% # Remove "Row.names" column because we don't need it
  #   tidyr::gather(., key = 'Gene', 'GeneScore', -GroupBy) #convert dataframe from a wide to a long format

  if (length(markerGenes) == 1){ # for some reason, if only one gene the matrix is already transposed
    # Get relevant gene expression
    filtered_gex <- matrix[markerGenes,] %>% # only keep the genes we are interested in
      as.data.frame() %>% #transpose matrix and then convert into dataframe because the functions bellow only work on dfs
      merge(., metadata_df, by.x = 0, by.y = 'rownames') %>% #merge metadata dataframe into this dataframe so we can access GroupBy information
      dplyr::select(., -`Row.names`) %>% # Remove "Row.names" column because we don't need it
      tidyr::gather(., key = 'Gene', 'GeneScore', -GroupBy) %>% #convert dataframe from a wide to a long format
      dplyr::mutate(., Gene = markerGenes) # gene names are removed when there's only one column
    
  } else {
    # Get relevant gene expression
    filtered_gex <- matrix[markerGenes,] %>% # only keep the genes we are interested in
      t() %>% as.data.frame() %>% #transpose matrix and then convert into dataframe because the functions bellow only work on dfs
      merge(., metadata_df, by.x = 0, by.y = 'rownames') %>% #merge metadata dataframe into this dataframe so we can access GroupBy information
      dplyr::select(., -`Row.names`) %>% # Remove "Row.names" column because we don't need it
      tidyr::gather(., key = 'Gene', 'GeneScore', -GroupBy) #convert dataframe from a wide to a long format
  }
  
  # Get mean gene expression
  filtered_gex_means <- filtered_gex %>%
    dplyr::group_by(., GroupBy, Gene) %>% # group by GroupBy variable and Gene, allows us to then apply the code bellow on each grouping
    dplyr::summarise(., mean_genescore = mean(GeneScore)) #get mean gene score
  
  # Calculate proportions with > 0
  overall_n <- metadata_df %>% 
    dplyr::count(., GroupBy) %>% #count how many times each factor in GroupBy occurs
    dplyr::rename(., n_total = n)  #rename column to n_total
  
  percentages <- filtered_gex %>% 
    dplyr::group_by(., GroupBy, Gene) %>% # group by GroupBy variable and Gene, allows us to then apply the code bellow on each grouping
    dplyr::mutate(., gene_score_yes = ifelse(GeneScore == 0, 0, 1)) %>% #check if GeneScore is bigger than 0
    dplyr::filter(., gene_score_yes == 1) %>% # only keep the instances where GeneScore was > 0
    dplyr::count(., gene_score_yes) %>% # count how often GeneScore was > 0, will add in a column called n
    merge(., overall_n, by = 'GroupBy') %>% # merge with overall_n dataframe which gives us total number per group in GroupBy
    mutate(., percent = (n/n_total)*100) %>% # get percentages
    dplyr::select(., GroupBy, Gene, percent)
  
  # Make marker_accessibility_table from above 
  marker_accessibility_table <- filtered_gex_means %>%
    merge(., percentages, by = c('GroupBy', 'Gene')) %>%
    dplyr::rename(., Genes = Gene, Groups = GroupBy, Score = mean_genescore, Proportion=percent)
  
  if (isTRUE(scale.scores)){
    #scale expression score across genes
    marker_accessibility_table <- transform(marker_accessibility_table, norm = ave(Score, Genes, FUN = scale))

  }
  
  # change Group and Gene columns into factor
  if (isTRUE(rotateaxis)) {
    marker_accessibility_table$Groups 		<- 	factor(marker_accessibility_table$Groups, levels = groupName)
    marker_accessibility_table$Genes 		<- 	factor(marker_accessibility_table$Genes, levels = rev(markerGenes))
  }	else {
    marker_accessibility_table$Groups 		<- 	factor(marker_accessibility_table$Groups, levels = rev(groupName))
    marker_accessibility_table$Genes 		<- 	factor(marker_accessibility_table$Genes, levels = markerGenes)
  }
  
  # make dot plot
  p 		<- ggplot(marker_accessibility_table) +
    geom_point(aes(x= Genes, y = Groups, color = Score, size = Proportion)) + 
    scale_colour_gradient(low = accessibilityColour[1],
                          high = accessibilityColour[2],
                          na.value = "black",
                          guide = "colourbar",
                          aesthetics = "colour") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  # + theme(aspect.ratio=1/1)
  
  # rotate x and y axis if required
  if (isTRUE(rotateaxis)) {
    p <- p + coord_flip()
  }
  
  if (isTRUE(scale.scores)){
    #add it as a title
    p <- p + ggtitle('Values are scaled')
    
  }
  
  # export file path and names
  if (is.null(exportFilePath)) {
    exportFilePath 		<- 	paste(getOutputDirectory(project), "/Plots",sep = "")
  }	else if (!dir.exists(exportFilePath))	{
    warning("exportFilePath does not exist. Using default output path.", "\n")
    exportFilePath 		<- 	paste(getOutputDirectory(project), "/Plots",sep = "")
  }
  
  exportFileName			<-	paste0("DotPlot_",exportFileName,".pdf")
  
  # save dot plot into pdf
  ggsave(
    file.path(exportFilePath, exportFileName),
    plot = p,
    device = "pdf",
    scale = 1,
    width = plotWidth,
    height = plotHeight,
    units = "cm"
  )
  
  return(p)

}

#----- TRACK PLOTS
make_track_plot 				<- 	function(project,
																			groupBy,
																			useGroups,
																			geneName,
							   											upstream,
							   											downstream,
							   											peakTrack = TRUE,
							   											coAccessibilityTrack = TRUE,
							   											geneTrack = TRUE,
							   											trackSize,
							   											exportFileName,
							   											exportFilePath = NULL,
							   											pal = NULL,
							   											plotWidth = 30,
							   											plotHeight = 12
							   											){

#		Arguments:
#		project = ArchrProject
#		groupBy = The column name in cellColData used for grouping cells together for marker feature identification (e.g. "Level_02_Phenotype").
#		useGroups = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData (e.g. c("LEC","VEC","AEC")).
#		geneName = A character vector that contains a gene name.
#		upstream = An interger vector that indicates the distance upstream of the gene body.
#		downstream = An interger vector that indicates the distance downstream of the gene body.
#		trackSize = A numeric vector that determines the relative sizes of each track (make sure the number of trackSize matches the number of track type)
#		exportFilePath = An alternative directory path to save the dot plot to. Default is ./Plots in outputDirectory of the ArchRProject.
#		exportFileName = A character vector that is used for the output file name (e.g. exportFileName = "Level_02_prox1a", the output file name would be "TrackPlot_Level_02_marker_genes.pdf").
#		pal = A character vector that specifies the colour for each group. (Set to Null to use default pal)
#		plotWidth = A numeric vector that specifies the width of the dot plot. 
#		plotHeight = A numeric vector that specifies the height of the dot plot. 

#		Example:
#		track_plot					<-	make_track_plot	(project = project,
#																							groupBy = "Level_02_Phenotype",
#																							useGroups = c("LEC", "VEC", "AEC"),
#																							geneName = "prox1a",
#																							upstream = 10000,
#																							downstream = 30000,
#	   																					peakTrack = TRUE,
#	  																					coAccessibilityTrack = FALSE,
#	   																					geneTrack = TRUE,	  											
#																							trackSize = c(3,1,1),
#																							exportFileName = "Level_02_prox1a",
#																							pal = c("#502A83","#BD1622","#EF8F00"),
#																							plotWidth = 8,
#																							plotHeight = 4
#							   															)

#======================================================================================================

# check if the gene name is correct

if (geneName %ni% project@geneAnnotation$genes$symbol) {
	stop(geneName, " is not found in gene annotations.", "\n")
}

# check if all the useGroups are correct
groupNamesError				<-	NULL

for (groups in useGroups) {
	if (groups %ni% project@cellColData[[groupBy]]) {
		groupNamesError		<-		paste(groupNamesError, groups, sep = ",")
	}
}

if (!is.null(groupNamesError)) {
	groupNamesError			<-		substring(groupNamesError, 2)
	stop(groupNamesError, " are not found in ", groupBy, ".", "\n")
}

# add gene size to upstream or downstream

rowNumber							<-		which(project@geneAnnotation$genes$symbol %in% geneName)
geneSize							<-		abs(start(project@geneAnnotation$genes[rowNumber]) - end(project@geneAnnotation$genes[rowNumber]))

if (as.character(strand(project@geneAnnotation$genes[rowNumber])) == "+") {
	downstream					<-		downstream + geneSize
}	else {
	upstream						<-		upstream + geneSize
}

# select requested tracks
tracks 								<-		c("bulkTrack")

if (isTRUE(peakTrack)) {
	tracks 							<-		append(tracks, "featureTrack")
}

if (isTRUE(coAccessibilityTrack)) {
	tracks 							<-		append(tracks, "loopTrack")
}

if (isTRUE(geneTrack)) {
	tracks 							<-		append(tracks, "geneTrack")
}

# make track plot
p 										<-		plotBrowserTrack(ArchRProj = project,
  																						groupBy = groupBy,
  																						useGroups = useGroups,
  																						plotSummary = tracks,
  																						sizes = trackSize,
  																						geneSymbol = geneName,
  																						upstream = upstream,
  																						downstream = downstream,
  																						pal = pal
  																						)

# export file
exportFileName				<-		paste0("TrackPlot_",exportFileName,".pdf")

ArchR::plotPDF(
  p, name=exportFileName, ArchRProj=project, addDOC=FALSE, width=plotWidth, height=plotHeight
  )

# make output directory for DAP analysis
if (!dir.exists(file.path(getOutputDirectory(project),"Plots","TrackPlots")))	{
	dir.create(file.path(getOutputDirectory(project),"Plots","TrackPlots"))
}

# move output file to specified location
if (is.null(exportFilePath)) {
	file.rename(file.path(getOutputDirectory(project),"Plots",exportFileName),file.path(getOutputDirectory(project),"Plots","TrackPlots",exportFileName))
}	else if (!dir.exists(exportFilePath))	{
	warning("exportFilePath does not exist. Using default output path.", "\n")
	file.rename(file.path(getOutputDirectory(project),"Plots",exportFileName),file.path(getOutputDirectory(project),"Plots","TrackPlots",exportFileName))
} else	{
	file.rename(file.path(getOutputDirectory(project),"Plots",exportFileName),file.path(exportFilePath,exportFileName))
}
	
return(p)

}

#======================================================================================================
#	VISUALISATION OF DAP AND DAG RESULTS
#======================================================================================================

import_DAG_tables <- function(importFilePath, statusGroup_01){
	#	ATAC DAG lists
	DAGTable 						<- read.table(importFilePath, header=TRUE, sep=",")
	colnames(DAGTable)				<- c("Gene", "Log2FC", "Pval", "Group_01_mean", "Background_mean")
	keep_genes 						<- remove_unannotated_genes(DAGTable$Gene)
	DAGTable						<- DAGTable[DAGTable$Gene %in% keep_genes,]
	DAGTable$Status_Group_01 		<- c(rep(statusGroup_01, nrow(DAGTable)))

	return(DAGTable)

}


import_DAP_tables <- function(importFilePath){
	#	ATAC DAG lists
	DAPTable 						<- read.table(importFilePath, header=TRUE, sep="\t")

	colNames <- c("Log2FC",
	"FDR",
	"RawPval",
	"AUC",
	"Group_01_mean",
	"Background_mean",
	"MeanDiff",
	"Chr",
	"Start",
	"End",
	"PeakScore",
	"distToGeneStart",
	"nearestGene",
	"peakType",
	"distTSS",
	"nearestTSS",
	"GC",
	"Status_Group_01")

	colnames(DAPTable)				<- colNames
	DAPTable 						<- DAPTable[!is.na(DAPTable$FDR),]

	DAPTable$PeakID 				<- paste(DAPTable$nearestGene, DAPTable$distTSS, DAPTable$peakType, sep="_")

	return(DAPTable)

}

make_pieChart_ATAC_peakType <- function(DAPTable,exportFilePath, exportFileName){

	peakType_df 			<- as.data.frame(table(DAPTable$peakType))
	colnames(peakType_df) 	<- c("Class", "Number")
	peakType_df$Class 		<- factor(peakType_df$Class, levels=c("Distal", "Intronic", "Promoter", "Exonic"))

	bp 						<- ggplot(peakType_df, aes(x="", y=Number, fill=Class)) + geom_bar(width = 1, stat = "identity")
	pie 					<- bp + coord_polar("y", start=0)
	pie 					<- pie + theme_minimal()

	peakType_pie_Filename 	<- paste(exportFilePath, exportFileName, "_peakType_pieChart.EPS", sep="")
	ggsave(peakType_pie_Filename, pie, device="eps", width=35, height=16, units = "cm")
}

make_pieChart_ATAC_DAG 	<- function(DAGTable, geneListColors, exportFilePath, exportFileName){

gene_df 				<- as.data.frame(table(DAGTable$barplotColors))
colnames(gene_df) 		<- c("Class", "Number")

bp 						<- ggplot(gene_df, aes(x="", y=Number, fill=Class)) + geom_bar(width = 1, stat = "identity")
pie 					<- bp + coord_polar("y", start=0)
pie 					<- pie + scale_fill_manual(values = c("lightgray", geneListColors))
pie 					<- pie + theme_minimal()

DAG_pie_Filename 		<- paste(exportFilePath, exportFileName, "_DAG_pieChart.EPS", sep="")

ggsave(DAG_pie_Filename, pie, device="eps", width=35, height=16, units = "cm")

}


make_pieChart_ATAC_DAP 	<- function(DAGTable, geneListColors, exportFilePath, exportFileName){

gene_df 				<- as.data.frame(table(DAGTable$barplotColors))
colnames(gene_df) 		<- c("Class", "Number")

bp 						<- ggplot(gene_df, aes(x="", y=Number, fill=Class)) + geom_bar(width = 1, stat = "identity")
pie 					<- bp + coord_polar("y", start=0)
pie 					<- pie + scale_fill_manual(values = c("lightgray", geneListColors))
pie 					<- pie + theme_minimal()

DAG_pie_Filename 		<- paste(exportFilePath, exportFileName, "_DAG_pieChart.EPS", sep="")

ggsave(DAG_pie_Filename, pie, device="eps", width=35, height=16, units = "cm")

}


make_FoldChange_BarPlot_ATAC_peaks 	<- function(DAPTable, geneLists, geneListColors, exportFilePath, exportFileName){

#	Merge the DEG information and plotting metadata and order the data frame
DAPTable 						<- DAPTable[!duplicated(DAPTable$PeakID),]
DAPTable 						<- DAPTable[with(DAPTable, order(Status_Group_01, -Log2FC)),]

barplotColors 					<- vector()

for ( gene in as.character(DAPTable$nearestGene) ){

	if ( gene %in% geneLists[[1]] == TRUE ){

		barplotColors <- c(barplotColors, names(geneLists[1]))

	} else if ( gene %in% geneLists[[2]] == TRUE ) {

		barplotColors <- c(barplotColors, names(geneLists[2]))

	} else {

		barplotColors <- c(barplotColors, "background")
	}
}

DAPTable$PeakID 				<- factor(DAPTable$PeakID, levels=as.character(DAPTable$PeakID))
DAPTable$barplotColors 			<- factor(barplotColors, levels=c("background", names(geneLists)))

DAP_barplotFilename 			<- paste(exportFilePath, exportFileName, "_DAP_Log2FC_barplot.EPS", sep="")
DAP_barplot_table_filename		<- paste(exportFilePath, exportFileName, "_DAP_Log2FC_barplot_table.txt", sep="")

p 								<- ggplot(data=DAPTable, aes(x=PeakID, y=Log2FC)) +
										geom_col(aes(fill=barplotColors)) +
										scale_fill_manual(values = c("lightgray", geneListColors)) +
										theme_classic() +
										theme(axis.title.x=element_blank(),
        									 axis.text.x=element_blank(),
        									 axis.ticks.x=element_blank())
ggsave(DAP_barplotFilename, p, device="eps", width=35, height=15, units = "cm")

write.table(DAPTable,
	file=DAP_barplot_table_filename,
	quote=FALSE,
	sep="\t",
	row.names=FALSE)

return(DAPTable)

}

make_FoldChange_BarPlot_ATAC_peaks_N_3 	<- function(DAPTable, geneLists, geneListColors, exportFilePath, exportFileName){

#	Merge the DEG information and plotting metadata and order the data frame
DAPTable 						<- DAPTable[with(DAPTable, order(Status_Group_01, -Log2FC)),]

barplotColors 					<- vector()

for ( gene in as.character(DAPTable$nearestGene) ){

	if ( gene %in% geneLists[[1]] == TRUE ){

		barplotColors <- c(barplotColors, names(geneLists[1]))

	} else if ( gene %in% geneLists[[2]] == TRUE ) {

		barplotColors <- c(barplotColors, names(geneLists[2]))

	} else if ( gene %in% geneLists[[3]] == TRUE ) {

			barplotColors <- c(barplotColors, names(geneLists[3]))

	} else {

		barplotColors <- c(barplotColors, "background")
	}
}

DAPTable$PeakID 				<- factor(DAPTable$PeakID, levels=as.character(DAPTable$PeakID))
DAPTable$barplotColors 			<- factor(barplotColors, levels=c("background", names(geneLists)))

DAP_barplotFilename 			<- paste(exportFilePath, exportFileName, "_DAP_Log2FC_barplot.EPS", sep="")
DAP_barplot_table_filename		<- paste(exportFilePath, exportFileName, "_DAP_Log2FC_barplot_table.txt", sep="")

p 								<- ggplot(data=DAPTable, aes(x=PeakID, y=Log2FC)) +
										geom_col(aes(fill=barplotColors)) +
										scale_fill_manual(values = c("lightgray", geneListColors)) +
										theme_classic() +
										theme(axis.title.x=element_blank(),
        									 axis.text.x=element_blank(),
        									 axis.ticks.x=element_blank())
ggsave(DAP_barplotFilename, p, device="eps", width=35, height=15, units = "cm")

write.table(DAPTable,
	file=DAP_barplot_table_filename,
	quote=FALSE,
	sep="\t",
	row.names=FALSE)

return(DAPTable)

}


make_FoldChange_BarPlot_ATAC 	<- function(DAGTable, geneLists, geneListColors, exportFilePath, exportFileName){

#	Merge the DEG information and plotting metadata and order the data frame
DAGTable 						<- DAGTable[!duplicated(DAGTable$Gene),]
DAGTable 						<- DAGTable[with(DAGTable, order(Status_Group_01, -Log2FC)),]

barplotColors 					<- vector()

for ( gene in as.character(DAGTable$Gene) ){

	if ( gene %in% geneLists[[1]] == TRUE ){

		barplotColors <- c(barplotColors, names(geneLists[1]))

	} else if ( gene %in% geneLists[[2]] == TRUE ) {

		barplotColors <- c(barplotColors, names(geneLists[2]))

	} else {

		barplotColors <- c(barplotColors, "background")
	}
}

DAGTable$Gene 					<- factor(DAGTable$Gene, levels=as.character(DAGTable$Gene))
DAGTable$barplotColors 			<- factor(barplotColors, levels=c("background", names(geneLists)))

DAG_barplotFilename 			<- paste(exportFilePath, exportFileName, "_DAG_Log2FC_barplot.EPS", sep="")
DAG_barplot_table_filename		<- paste(exportFilePath, exportFileName, "_DAG_Log2FC_barplot_table.txt", sep="")

p 								<- ggplot(data=DAGTable, aes(x=Gene, y=Log2FC)) +
										geom_col(aes(fill=barplotColors)) +
										scale_fill_manual(values = c("lightgray", geneListColors)) +
										theme_classic() +
										theme(axis.title.x=element_blank(),
        									 axis.text.x=element_blank(),
        									 axis.ticks.x=element_blank())
ggsave(DAG_barplotFilename, p, device="eps", width=35, height=15, units = "cm")

write.table(DAGTable,
	file=DAG_barplot_table_filename,
	quote=FALSE,
	sep="\t",
	row.names=FALSE)

return(DAGTable)

}


#======================================================================================================
#	PLOT MARKER GENES
#======================================================================================================

plot_marker_genes     <- function(
              ArchRProject,
              MatrixToPlot,
              geneList,
              pointSize,
              colourScheme,
              filePath,
              fileName){

#   This function will remove the legend and use inbuilt ArchR colour palettes

for ( gene in geneList ) {

file_save <- paste(filePath,"/",fileName,"_gene_accessibility_",gene,".EPS",sep="")

  p1 <- plotEmbedding(
        ArchRProj = ArchRProject, 
        colorBy = "GeneScoreMatrix", 
        name = gene, 
        continuousSet = "horizonExtra",
        embedding = "UMAP",
        size = 1,
        imputeWeights = getImputeWeights(ArchRProject))

    x <- ggplot_build(p1)$data[[1]]
    ggplot(x, aes(x=x, y=y)) + geom_point(aes(color=value), size=pointSize) + 
      theme_classic() + theme(legend.position="none") + 
      scale_colour_gradientn(colours=ArchRPalettes$horizonExtra)
    
    ggsave(file_save,
     device="eps",
     width=25,
     height=15,
     units="cm",
     scale=1)

  }
}


plot_marker_genes_for_publication     <- function(
              ArchRProject,
              MatrixToPlot,
              geneList,
              pointSize,
              colours,
              legendPosition,
              filePath,
              fileName,
              width.plot = 25,
              height.plot = 15){

#   This function will allow you to keep/remove the legend and select your own colours

for ( gene in geneList ) {

file_save <- paste(filePath,"/",fileName,"_gene_accessibility_",gene,".EPS",sep="")

  p1 <- plotEmbedding(
        ArchRProj = ArchRProject, 
        colorBy = "GeneScoreMatrix", 
        name = gene,
        embedding = "UMAP",
        size = 1,
        imputeWeights = getImputeWeights(ArchRProject))

    x <- ggplot_build(p1)$data[[1]]
    ggplot(x, aes(x=x, y=y)) + geom_point(aes(color=value), size=pointSize) + 
      theme_classic() + theme(legend.position=legendPosition) + 
      scale_colour_gradientn(colours=c("#d9d9d9", "#7a0177"))
    
    ggsave(file_save,
     device="eps",
     width=width.plot,
     height=height.plot,
     units="cm",
     scale=1)

  }
}

#======================================================================================================
#	RNA-SEQ INTEGRATION
#======================================================================================================

define_scRNA_genotype_from_barcodes <- function(genotypeDataTable, RNAseq_sampleLabel){

# Arguments:
#   genotypeDataTable  = 
#   RNAseq_sampleLabel = character vector corresponding to sample genotype

  barcodes      <- strsplit(genotypeDataTable$RNA_barcode, "-")
  
  RNA_genotype    <- vector()
  
  for ( cell in 1:length(barcodes) ){
    temp      <- barcodes[[cell]][2]
    RNA_genotype  <- c(RNA_genotype, temp)
  }

  genotypeDataTable           <- cbind(genotypeDataTable, RNA_genotype)
  genotypeDataTable$RNA_genotype    <- factor(genotypeDataTable$RNA_genotype, labels=RNAseq_sampleLabel)

  genotypeDataTable$ATAC_RNA_genotypes <- paste(genotypeDataTable$ATAC_genotype, genotypeDataTable$RNA_genotype, sep="_")
  genotypeDataTable$ATAC_RNA_genotypes <- factor(genotypeDataTable$ATAC_RNA_genotypes)

  return(genotypeDataTable)
}


make_pieChart_Genotype  <- function(
              genotypeDataTable,
              exportFilePath,
              exportFileName){

  match_df      <- as.data.frame(table(genotypeDataTable$ATAC_RNA_genotypes))
  colnames(match_df)  <- c("Class", "Number")

  bp          <- ggplot(match_df, aes(x="", y=Number, fill=Class)) + geom_bar(width = 1, stat = "identity")
  pie         <- bp + coord_polar("y", start=0)
  pie         <- pie + scale_fill_manual(values = c("#31a354", "#de2d26", "#fee0d2", "#e5f5e0"))
  pie         <- pie + theme_minimal()

  pie_Filename    <- paste(exportFilePath, "/Plots/ATAC_RNA_Integration_Genotype_match_pieChart_",exportFileName,".EPS", sep="")

  ggsave(pie_Filename, pie, device="eps", width=35, height=16, units = "cm")

}



#======================================================================================================
#	PREPARE FOR CELLXGENE
#======================================================================================================

prep_for_cellxgene <- function(
	
  project_path, use_matrix, maximum_score=0, outfile_dir="./"
){
  
#	Author: Tyrone Chen

  #	Arguments:
  #		project_path = path to ArchR project
  #   use_matrix = ArchR matrix: must be in ArchR::getAvailableMatrices(project)
  #   maximum_score = set maximum score as filter
  #		outfile_dir = character vector containing the dataset directory
  #     "GeneIntegrationMatrix" "GeneScoreMatrix"       "MotifMatrix"
  #     "PeakMatrix"            "TileMatrix"
  project <- loadArchRProject(path=project_path, showLogo=FALSE)

  print("Available matrices:")
  print(ArchR::getAvailableMatrices(project))

  if (!use_matrix %in% ArchR::getAvailableMatrices(project)) {
    stop("Matrix is not present, check ArchR::getAvailableMatrices(project)!")
  }

	# get matrix
	project_matrix <- getMatrixFromProject(
    ArchRProj = project, useMatrix = use_matrix
  )

	#	Cell metadata
	cell_meta_data <- project_matrix@colData
	cell_meta_data <- cell_meta_data[order(row.names(cell_meta_data)),]

	#	UMAP coordinates
	umap_coordinates <- project@embeddings$UMAP$df
	umap_coordinates <- umap_coordinates[order(row.names(umap_coordinates)),]
	colnames(umap_coordinates) <- c("UMAP_1","UMAP_2")
	combined_data <- cbind(umap_coordinates, cell_meta_data)

	meta_out <- paste(outfile_dir,"/CellXGene_MetaData_",use_matrix,".txt",sep="")
	write.table(combined_data, file=meta_out, sep="\t", quote=FALSE)

	#	Expression data
	expression_data <- as.data.frame(
    as.matrix(assay(project_matrix)), row.names = rowData(project_matrix)$name
  )
	expression_data <- replace(
    expression_data, expression_data > maximum_score, maximum_score
  )
	expression_data <- expression_data[,order(colnames(expression_data))]

	expn_out <- paste(outfile_dir,"/CellXGene_ExpnData_",use_matrix,".txt",sep="")
	write.table(expression_data, file=expn_out, sep="\t", quote=FALSE)
}

parse_argv <- function() {
  p <- argparser::arg_parser(
    "Take ArchR project, matrix id, generate matrices of data and metadata."
  )
  # Add command line arguments
  p <- argparser::add_argument(
    p, "--project_path", type="character", nargs=1, default=NA,
    help="path to project"
  )
  p <- argparser::add_argument(
    p, "--use_matrix", type="character", nargs=1, default=NA,
    help="matrix identity, must be in ArchR::getAvailableMatrices(project)!"
  )
  p <- argparser::add_argument(
    p, "--maximum_score", type="numeric", nargs=1, default=Inf,
    help="score threshold"
  )
  p <- argparser::add_argument(
    p, "--outfile_dir", type="character", nargs=1, default=NA,
    help="path to output directory"
  )
  argv <- argparser::parse_args(p)
  return(argv)
}

main <- function() {
  argv <- parse_argv()
  project_path <- argv$project_path
  use_matrix <- argv$use_matrix
  maximum_score <- argv$maximum_score
  outfile_dir <- argv$outfile_dir
  prep_for_cellxgene(project_path, use_matrix, maximum_score, outfile_dir)
}