#! /bin/sh
#
#SBATCH --job-name="Make_CellXGene"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=1-00
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.err"
#SBATCH --partition=prod
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Elizabeth.Mason@petermac.org

#	Make CellXGene file for ATAC data
#	Test Michelle's script

cd /team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset/CellXGene


tmux new -s ATACproject2cellxgene

/team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/ATAC_pipeline/ATAC_project2cellxgene.sh \
/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset/ \
GeneScoreMatrix \
999999999 \
/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset/CellXGene &> ATAC_project2cellxgene_out.txt