#!/bin/sh 

set -o allexport

sta_N_i_pCaLs=1
end_N_i_pCaLs=3

sta_N_i_KM=1
end_N_i_KM=3

sta_N_i_batch=1
end_N_i_batch=1

for i_pCaLs in $(eval echo "{$sta_N_i_pCaLs..$end_N_i_pCaLs}") 
    do
        echo "$i_pCaLs"
        for i_KM in $(eval echo "{$sta_N_i_KM..$end_N_i_KM}") 
            do 
                echo " $i_KM"
                for i_batch in $(eval echo "{$sta_N_i_batch..$end_N_i_batch}") 
                    do 
                        echo " $i_batch"               
                        sbatch --partition=batch --job-name km_fI --output data/output/km_fI-cals-${i_pCaLs}-km-"$i_KM"-batch-"$i_batch"-%j.txt --error data/output/km_fI-cals-${i_pCaLs}-km-"$i_KM"-batch-"$i_batch"-%j.err --export=ALL job_km.sh
                    done
            done
    done