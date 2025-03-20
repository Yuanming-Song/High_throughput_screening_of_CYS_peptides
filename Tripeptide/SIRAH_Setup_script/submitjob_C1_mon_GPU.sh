#!/bin/bash
# Function to check if an amino acid is in a given array
is_in_array() {
    local element="$1"
    shift
    for item; do
        if [[ "$item" == "$element" ]]; then
            return 0
        fi
    done
    return 1
}
# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")
# Define arrays for positively and negatively charged residues
positive_charged=("K" "R")
negative_charged=("D" "E")
amino_acids1=("T")
# Loop through each amino acid for amino1
for amino1 in "${amino_acids[@]}"; do
    # Loop through each amino acid for amino2
    for amino2 in "${amino_acids1[@]}"; do
        # Check if neither amino1 nor amino2 is H; if so, skip to the next iteration
        #if [[ "$amino1" != "H" && "$amino2" != "H" ]]; then
        #    continue  # Skip this loop iteration
        #fi
        # Determine the GPU type based on the charge of amino1 and amino2
        
        if (is_in_array "$amino1" "${positive_charged[@]}" || is_in_array "$amino2" "${positive_charged[@]}") ||
            (is_in_array "$amino1" "${negative_charged[@]}" || is_in_array "$amino2" "${negative_charged[@]}"); then
            gpu_type="A100:1"
        else
            gpu_type="A30:1"
        fi

        # Submit the job using sbatch
        sbatch --job-name="C_${amino1}_${amino2}" \
            --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
            --account=dtobias_lab \
            --partition=free-gpu \
            --nodes=1 \
            --ntasks-per-node=1 \
            --cpus-per-task=8 \
            --time=5:00:00 \
            --gres=gpu:$gpu_type \
            --out=out/C_${amino1}_${amino2}.out \
            --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide/SIRAH/Cys_unstapled_tripeptide_SIRAH_run_GPU.sh C ${amino1} ${amino2} "
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
