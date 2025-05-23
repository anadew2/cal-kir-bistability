#!/bin/sh 
#
#SBATCH --job-name=dI_mc
#SBATCH --output=data/output/dI_mc-%j.txt
#SBATCH --error=data/output/dI_mc-%j.err
#SBATCH --partition=batch
#
#SBATCH --time=3:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_var_mc_kir.jl $i_batch_mc