#!/bin/sh 

set -o allexport

sta_N_i_pCaLs=1
end_N_i_pCaLs=3

sta_N_i_Kir=1
end_N_i_Kir=3

sta_N_i_batch=1
end_N_i_batch=1

for i_pCaLs in $(eval echo "{$sta_N_i_pCaLs..$end_N_i_pCaLs}") 
    do
        echo "$i_pCaLs"
        for i_Kir in $(eval echo "{$sta_N_i_Kir..$end_N_i_Kir}") 
            do 
                echo " $i_Kir"
                for i_batch in $(eval echo "{$sta_N_i_batch..$end_N_i_batch}") 
                    do 
                        echo " $i_batch"
                        sbatch --partition=batch --job-name kir_fI --output data/output/kir_fI-cals-${i_pCaLs}-kir-"$i_Kir"-batch-"$i_batch"-%j.txt --error data/output/kir_fI-cals-${i_pCaLs}-kir-"$i_Kir"-batch-"$i_batch"-%j.err --export=ALL job_kir.sh
                    done
            done
    done