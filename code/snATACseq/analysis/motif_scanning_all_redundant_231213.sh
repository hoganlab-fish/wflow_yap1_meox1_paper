#!/bin/bash
#SBATCH --job-name="meox_enhancer_all_motif_redundant"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.err"
#SBATCH --partition=prod_med
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Oliver.Yu@petermac.org

#   Author: Oliver Yu
#   Last updated: Dec 2023

#======================================================================================

outputDir="/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq"

#======================================================================================

module load bedtools/2.27.1

reference_genome="/team_folders/hogan_lab/genomes/DR_GRCz11.97_scATAC/fasta/genome.fa" 
bedFile="/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/input/Motif_analysis_files_202311/Yap1_atac_Level_02_meox1_associated_peaks.bed"
fastaOutputFile="/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/input/Motif_analysis_files_202311/Yap1_atac_Level_02_meox1_associated_peaks.fasta"

bedtools getfasta \
  -fi ${reference_genome} \
  -bed ${bedFile} \
  -fo ${fastaOutputFile}

#======================================================================================

module load meme/5.4.1

fastaFile=${fastaOutputFile}
motifOutputDir="${outputDir}/Motif_analysis/MEME/231213_meox1_enhancers_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}