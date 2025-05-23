#!/bin/sh 
#
#SBATCH --partition=batch
#
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

julia script_dI_noise_cyc_kir.jl $i_std $i_batch