#!/bin/sh 

set -o allexport

sta_i_pCaLs=1
end_i_pCaLs=5

sta_i_Kir=1
end_i_Kir=31

sta_i_batch_KM=1
end_i_batch_KM=31
n_batch_KM=1

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do 
        echo " $i_pCaLs"
        for i_Kir in $(eval echo "{$sta_i_Kir..$end_i_Kir}") 
            do 
                echo " $i_Kir"
                for i_batch_KM in $(eval echo "{$sta_i_batch_KM..$end_i_batch_KM..$n_batch_KM}") 
                    do 
                        echo " $i_batch_KM"
                        sbatch --partition=batch --job-name dI_altgt --output data/output/dI_altgt-ipCaLs-${i_pCaLs}-iKir-${i_Kir}-iKM-${i_batch_KM}-batch-"$i_batch"-%j.txt --error data/output/dI_altgt-ipCaLs-${i_pCaLs}-iKir-${i_Kir}-iKM-${i_batch_KM}-batch-"$i_batch"-%j.err --export=ALL job.sh
                    done
            done
    done