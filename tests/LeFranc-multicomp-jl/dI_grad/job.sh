#!/bin/sh 
#
#SBATCH --partition=batch
#
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_grad_Kir_s.jl $i_pCaLs_d $i_Kir_s