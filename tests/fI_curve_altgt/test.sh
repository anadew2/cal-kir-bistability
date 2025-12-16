#!/bin/sh 

set -o allexport

sta_i_pCaLs=1
end_i_pCaLs=2

sta_i_Kir=1
end_i_Kir=2

sta_i_KM=1
end_i_KM=2

sta_N_i_batch=1
end_N_i_batch=1

for i_pCaLs in $(eval echo "{$sta_i_pCaLs..$end_i_pCaLs}") 
    do
        echo "$i_pCaLs"
        for i_Kir in $(eval echo "{$sta_i_Kir..$end_i_Kir}") 
            do 
                echo " $i_Kir"
                for i_KM in $(eval echo "{$sta_i_KM..$end_i_KM}") 
                    do 
                    echo " $i_KM"
                        for i_batch in $(eval echo "{$sta_N_i_batch..$end_N_i_batch}") 
                            do 
                                echo " $i_batch"
                                sbatch --partition=batch --job-name fI_altgt --output data/output/dI_altgt-ipCaLs-${i_pCaLs}-iKir-${i_Kir}-iKM-${i_KM}-batch-"$i_batch"-%j.txt --error data/output/dI_altgt-ipCaLs-${i_pCaLs}-iKir-${i_Kir}-iKM-${i_KM}-batch-"$i_batch"-%j.err --export=ALL job.sh
                            done
                    done
            done
    done