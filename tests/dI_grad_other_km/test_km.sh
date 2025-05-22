#!/bin/sh 

set -o allexport

sta_i_pCaLs=24
end_i_pCaLs=24 

sta_i_batch_KM=1
end_i_batch_KM=17 
n_batch_KM=1

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do 
        echo " $i_pCaLs"
        for i_batch_KM in $(eval echo "{$sta_i_batch_KM..$end_i_batch_KM..$n_batch_KM}") 
            do 
                echo " $i_batch_KM"
                sbatch --partition=batch --job-name dI_KMCaLs --export=ALL job_km.sh
            done
    done