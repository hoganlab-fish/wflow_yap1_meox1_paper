#!/bin/bash
#SBATCH --job-name="meox_enhancer_all_motif_redundant_-32meox1_vertebrates"
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

fastaFileDir="/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/input/Motif_analysis_-32meox1_vertebrates_240119"

#======================================================================================

module load meme/5.4.1

fastaFile="${fastaFileDir}/-32meox1_barbel.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_barbel_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_carp.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_carp_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_chicken.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_chicken_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_goldfish.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_goldfish_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_human.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_human_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_mouse.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_mouse_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}

fastaFile="${fastaFileDir}/-32meox1_zebrafish.fasta"
motifOutputDir="${outputDir}/Motif_analysis/MEME/240119_-32meox1/240119_-32meox1_zebrafish_all_motifs_redundant_1e2"
motifFile="/team_folders/hogan_lab/Hogan_Lab_People/Oliver_Yu/motif_files/JASPAR_2024/JASPAR2024_CORE_redundant_pfms_meme.txt"
PvalThresh=1e-2

fimo --o ${motifOutputDir} --thresh ${PvalThresh} ${motifFile} ${fastaFile}