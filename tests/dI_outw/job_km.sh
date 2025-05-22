#!/bin/sh 
#
#SBATCH --job-name=dI_KMCaLs_outw
#SBATCH --output=data/output/dI_KMCaLs_outw-%j.txt
#SBATCH --error=data/output/dI_KMCaLs_outw--%j.err
#SBATCH --partition=batch
#
#SBATCH --time=2:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_outw_km.jl $i_batch_KM