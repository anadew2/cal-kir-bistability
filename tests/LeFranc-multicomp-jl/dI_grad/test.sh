#!/bin/sh 

set -o allexport

sta_i_pCaLs_d=1
end_i_pCaLs_d=31

sta_i_Kir_s=1
end_i_Kir_s=41
n_batch_Kir_s=1


for i_pCaLs_d in $(eval echo "{$sta_i_pCaLs_d..$end_i_pCaLs_d}") 
    do
        echo "$i_pCaLs_d"
        for i_Kir_s in $(eval echo "{$sta_i_Kir_s..$end_i_Kir_s..$n_batch_Kir_s}") 
            do 
                echo " $i_Kir_s"
                sbatch --partition=batch --job-name dI_comp --output data/output_J_I1down/dI_comp-ipCaLs_d-${i_pCaLs_d}-iKir_s-${i_Kir_s}-%j.txt --error data/output_J_I1down/dI_comp-ipCaLs_d-${i_pCaLs_d}-iKir_s-${i_Kir_s}-%j.err --export=ALL job.sh
            done
    done