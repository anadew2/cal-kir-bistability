#!/bin/sh 

set -o allexport

sta_i_batch_KM=1
end_i_batch_KM=51
n_batch_KM=5

for i_batch_KM in $(eval echo "{$sta_i_batch_KM..$end_i_batch_KM..$n_batch_KM}") 
    do 
        echo " $i_batch_KM"
        sbatch --partition=batch --job-name dI_KMCaLs_outw --export=ALL job_km.sh
    done