#!/bin/sh 

set -o allexport

sta_i_batch_Kir=1
end_i_batch_Kir=51
n_batch_Kir=5

for i_batch_Kir in $(eval echo "{$sta_i_batch_Kir..$end_i_batch_Kir..$n_batch_Kir}") 
    do 
        echo " $i_batch_Kir"
        sbatch --partition=batch --job-name dI_KirCaLs_inw --export=ALL job.sh
    done