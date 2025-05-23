#!/bin/sh 
#
#SBATCH --job-name=dI_KirCaLs_inw
#SBATCH --output=data/output/dI_KirCaLs_inw-%j.txt
#SBATCH --error=data/output/dI_KirCaLs_inw--%j.err
#SBATCH --partition=batch
#
#SBATCH --time=2:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_inw_kir.jl $i_batch_Kir