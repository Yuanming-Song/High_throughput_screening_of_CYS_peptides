#!/bin/bash

# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")
amino_acids=("G")
amino_acids1=("V")
# Loop through each amino acid for amino1
for amino1 in "${amino_acids[@]}"; do
    # Loop through each amino acid for amino2
    for amino2 in "${amino_acids1[@]}"; do
        # Check if neither amino1 nor amino2 is H; if so, skip to the next iteration
        #if [[ "$amino1" != "H" && "$amino2" != "H" ]]; then
        #    continue  # Skip this loop iteration
        #fi
        
        # Submit the job using sbatch
        sbatch --job-name="${amino1}_C_${amino2}_mon" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
               --account=dtobias_lab \
               --partition=free-gpu \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=8 \
               --time=6:00:00 \
               --out=out/${amino1}_C_${amino2}.out \
               --gres=gpu:A30:1 \
               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Cys_unstapled_tripeptide_SIRAH_run_GPU.sh ${amino1} C ${amino2} "
    done
done

#sbatch --job-name="SASA_calc" \
#               --mail-user=yuanmis1@uci.edu --mail-type=ALL \
#               --account=dtobias_lab \
#               --partition=standard \
#               --nodes=1 \
#               --ntasks-per-node=1 \
#               --cpus-per-task=10 \
#               --time=24:00:00 \
#               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SASAcalc.sh"
