#!/bin/sh 

set -o allexport

sta_i_batch_mc=1
end_i_batch_mc=300
n_batch_mc=5

for i_batch_mc in $(eval echo "{$sta_i_batch_mc..$end_i_batch_mc..$n_batch_mc}") 
    do 
        echo " $i_batch_mc"
        sbatch --partition=batch --job-name dI_mc_km --export=ALL job_km.sh
    done
