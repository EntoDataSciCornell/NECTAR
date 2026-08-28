#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=2-10:00:00
#SBATCH --partition=regular
#SBATCH --chdir=/local/storage/NECTAR/Code/12_Analyses/12_04_SDM_Evaluation
#SBATCH --job-name=SDM_CV
#SBATCH --mail-user=ctc225@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-200%18
#SBATCH --output=logs/SDM_CV_%A_%a.out
#SBATCH --error=logs/SDM_CV_%A_%a.err

module purge
unset R_HOME
module load R/4.6.1

export R_LIBS_USER="/local/storage/NECTAR/Code/12_Analyses/12_04_SDM_Evaluation/R_libs/sdm_cv_R_4.6.1"
export SDM_ROOT="/local/storage/NECTAR"

Rscript 01_sdm_spatial_cv.R $SLURM_ARRAY_TASK_ID
