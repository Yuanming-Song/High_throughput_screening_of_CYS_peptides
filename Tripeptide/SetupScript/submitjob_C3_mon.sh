#!/bin/bash

# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")

# Loop through each amino acid for amino1
for amino1 in "${amino_acids[@]}"; do
    # Loop through each amino acid for amino2
    for amino2 in "${amino_acids[@]}"; do
        # Check if neither amino1 nor amino2 is H; if so, skip to the next iteration
#        if [[ "$amino1" != "H" && "$amino2" != "H" ]]; then
#            continue  # Skip this loop iteration
#        fi
        
        # Submit the job using sbatch
        sbatch --job-name="${amino1}_${amino2}_C" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
               --account=dtobias_lab \
               --partition=standard \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=40 \
               --time=3:00:00 \
               --out=out/${amino1}_${amino2}_C.out \
               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SetupScript/Cys_unstapled_tripeptide_MARTINI_run.sh ${amino1} ${amino2} C"
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
