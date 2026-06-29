#!/bin/sh
# convert archr project to cellxgene file

matrix="GeneScoreMatrix"

Rscript /hogan_lab/Hogan_Lab_People/Lizzie_Mason/scripts/scripts_hogan_lab/scripts_projects/R_snATAC_seq_functions/archr_to_matrix.r \
  --project_path "/team_folders/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/Level_02_Endothelial_cells/" \
  --use_matrix ${matrix} \
  --maximum_score 999999999 \
  --outfile_dir "/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/CellXGene"

python /hogan_lab/Hogan_Lab_People/Lizzie_Mason/scripts/scripts_hogan_lab/scripts_projects/R_snATAC_seq_functions/matrix_to_cellxgene.py \
  "CellXGene_ExpnData_${matrix}.txt" \
  "CellXGene_MetaData_${matrix}.txt" \
  "CellXGene_${matrix}.h5ad"
