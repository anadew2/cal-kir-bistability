#!/bin/sh 
#
#SBATCH --error=data/output/kir_fI-%j.err
#SBATCH --partition=batch
#
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_fI_curve.jl $i_pCaLs $i_Kir $i_batch