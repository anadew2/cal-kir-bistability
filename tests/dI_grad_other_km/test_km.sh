#!/bin/sh 

set -o allexport

sta_i_pCaLs=1
end_i_pCaLs=24 

sta_i_batch_KM=1
end_i_batch_KM=16 
n_batch_KM=1

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do 
        echo " $i_pCaLs"
        for i_batch_KM in $(eval echo "{$sta_i_batch_KM..$end_i_batch_KM..$n_batch_KM}") 
            do 
                echo " $i_batch_KM"
                sbatch --partition=batch --job-name dI_oth_KM --output data/output/dI_oth_KM-cals-${i_pCaLs}-km-"$i_batch_KM"-%j.txt --error data/output/dI_oth_KM-cals-${i_pCaLs}-km-"$i_batch_KM"-%j.err --export=ALL job_km.sh
            done
    done