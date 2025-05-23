#!/bin/sh 

set -o allexport

sta_i_pCaLs=1
end_i_pCaLs=21

sta_i_batch_KM=1
end_i_batch_KM=16
n_batch_KM=6

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do 
        echo " $i_pCaLs"
        for i_batch_KM in $(eval echo "{$sta_i_batch_KM..$end_i_batch_KM..$n_batch_KM}") 
            do 
                echo " $i_batch_KM"
                sbatch --partition=batch --job-name dI_KMCaLs_no_ghk --export=ALL job_km.sh
            done
    done