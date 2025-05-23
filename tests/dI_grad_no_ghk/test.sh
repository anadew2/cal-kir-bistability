#!/bin/sh 

set -o allexport

sta_i_pCaLs=1
end_i_pCaLs=16

sta_i_batch_Kir=1
end_i_batch_Kir=6
n_batch_Kir=6

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do 
        echo " $i_pCaLs"
        for i_batch_Kir in $(eval echo "{$sta_i_batch_Kir..$end_i_batch_Kir..$n_batch_Kir}") 
            do 
                echo " $i_batch_Kir"
                sbatch --partition=batch --job-name dI_KirCaLs_no_ghk --export=ALL job.sh
            done
    done