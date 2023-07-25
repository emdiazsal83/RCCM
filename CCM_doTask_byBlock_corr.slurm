#!/bin/bash
#SBATCH --array=1%16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --exclude=nodo16,nodo17,nodo18,nodo19,nodo20
#SBATCH --job-name=R_CCM_job
#SBATCH --output=log/R_CCM_job_%A_%a.out

# Para usar el R del sistema
# module load R

# Run command
Rscript CCM_getCorr_byBlock.R
