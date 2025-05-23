#!/bin/sh 
#
#SBATCH --job-name=dI_KirCaLs
#SBATCH --output=data/output/dI_KirCaLs-%j.txt
#SBATCH --error=data/output/dI_KirCaLs--%j.err
#SBATCH --partition=batch
#
#SBATCH --time=2:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_grad.jl $i_pCaLs $i_batch_Kir