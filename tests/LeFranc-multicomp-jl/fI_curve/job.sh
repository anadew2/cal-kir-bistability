#!/bin/sh 
#
#SBATCH --partition=batch
#
#SBATCH --time=7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_fI_curve.jl $i_pCaLs_d $i_Kir_s