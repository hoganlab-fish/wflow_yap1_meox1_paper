#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=14:00:00
#SBATCH --mail-user=michelle.meier@petermac.org
#SBATCH --mail-type=END
#SBATCH --partition=prod_med
#SBATCH --output="logs/ArchR_makearrow_%j.out"
#SBATCH --mem-per-cpu=4062

# -----------------------
module load R/4.2.0.Core

# ==============================================================================
# GET VARIABLES FROM INPUT ARGUMENTS
# ==============================================================================
#set variables, one trailing argument
csv_file=$1 

#get sample paths and names
infile_paths=$(tail -n +2 $csv_file | cut -d ',' -f2 | tr '\n' ',') 
sample_names=$(tail -n +2 $csv_file | cut -d ',' -f1 | tr '\n' ',') 
arrow_path=$(sed -n '2p' $csv_file | cut -d ',' -f3)

mkdir -p ${arrow_path}
cd ${arrow_path}

#make a log folder to capture log outputs
mkdir -p logs

# ==============================================================================
# ENSURING RENV IS INSTALLED
# ==============================================================================
# check if renv/library exists --> this would indicate that the environment was built previously
renv_dir="./renv/library"
if [ -d "$renv_dir" ]; then
echo "Using found r environment ${renv_dir}"

else
echo "Could not find a r environment in current directory"
echo "Start building ATAC renv"
#make execution folder by copying renv and move there
cp -r /team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_env/4.2.0.Core_atac .
#annoying clean up (can we avoid this?)
mv 4.2.0.Core_atac/renv .
mv 4.2.0.Core_atac/renv.lock .
mv 4.2.0.Core_atac/.Renviron .
mv 4.2.0.Core_atac/.Rprofile . #moving without wildcard to ensure full transparency of the files we need
rmdir 4.2.0.Core_atac

# intitate R to bootstrap and install renv
R -e 'renv::restore()'
fi

# ==============================================================================
# RUNNING ARROW FILE SCRIPT
# ==============================================================================

Rscript /team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/ATAC_pipeline/src/ArchR_makearrow.R \
--arrow_path=$arrow_path \
--infile_paths=$infile_paths \
--ncpus=16 \
--sample_names=$sample_names