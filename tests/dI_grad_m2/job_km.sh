#!/bin/sh 
#
#SBATCH --job-name=dI_KMCaLs_m2
#SBATCH --output=data/output/dI_KMCaLs_m2-%j.txt
#SBATCH --error=data/output/dI_KMCaLs_m2--%j.err
#SBATCH --partition=batch
#
#SBATCH --time=2:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_grad_f_repl_by_s_km.jl $i_pCaLs $i_batch_KM