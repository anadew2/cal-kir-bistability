#!/bin/sh 

set -o allexport

sta_N_std=1
end_N_std=24 #24 #26 for I0~I1

sta_N_batch=1
end_N_batch=20 #20

for i_std in $(eval echo "{$sta_N_std..$end_N_std}") 
    do
        echo "$i_std"
        for i_batch in $(eval echo "{$sta_N_batch..$end_N_batch}") 
            do 
                echo " $i_batch"
                sbatch --partition=batch --job-name kir_noise --output data/output/kir_noise-sigma-${i_std}-batch-"$i_batch"-%j.txt --error data/output/kir_noise-sigma-${i_std}-batch-"$i_batch"-%j.err --export=ALL job_st.sh
            done
    done