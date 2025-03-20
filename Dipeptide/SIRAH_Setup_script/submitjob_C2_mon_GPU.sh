#!/bin/bash

# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")
# Define arrays for positively and negatively charged residues
positive_charged=("K" "R")
negative_charged=("D" "E")
# File containing job IDs
job_id_file="/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Tripeptide_mon_C3_redo/submitted_job_ids.txt"

# Read the last line of the file for dependency job IDs
dependency_jobs=$(tail -n 1 "$job_id_file")

# Ensure the job IDs are available
if [[ -z "$dependency_jobs" ]]; then
    echo "No job IDs found in the last line of $job_id_file."
    #exit 1
fi

# Convert the dependency jobs to an array
read -ra job_ids <<<"$dependency_jobs"

# Initialize a counter for job dependencies
counter=0

# Define the file to store job IDs
job_id_file="submitted_job_ids.txt"

# Clear the file if it exists
current_run_jobs=""

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

# Loop through each amino acid for amino2
for amino2 in "${amino_acids[@]}"; do
    # Determine the GPU type based on the charge of amino1 and amino2
    if (is_in_array "$amino2" "${positive_charged[@]}" || is_in_array "$amino2" "${negative_charged[@]}"); then
        gpu_type="A100:1"
    else
        gpu_type="A30:1"
    fi
    # Use the counter to pick a dependency job
    dependency_job="${job_ids[$counter]}"
    if [ -d "${amino2}_C" ]; then
        continue
    fi
    #            --dependency=afterany:$dependency_job \

    # Submit the job using sbatch
    new_job_id=$(
        sbatch --parsable \
         --job-name="${amino2}_C_mon" \
            --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
            --account=dtobias_lab_GPU \
            --partition=free-gpu \
            --nodes=1 \
            --ntasks-per-node=1 \
            --cpus-per-task=8 \
            --time=5:00:00 \
            --out=out/${amino2}_C.out \
            --gres=gpu:${gpu_type} \
            --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Cys_unstapled_dipeptide_MARTINI_run.sh ${amino2} C "
    )
    # Check if the job submission was successful
    if [ $? -eq 0 ]; then
        # Append the job ID to the file
        current_run_jobs+="$new_job_id "
    else
        echo "Failed to submit job for ${amino2}_C_mon"
    fi

    # Increment the counter
    counter=$((counter + 1))

    # Reset the counter if it exceeds the number of dependency jobs
    if [[ $counter -ge ${#job_ids[@]} ]]; then
        counter=0
    fi
done
echo "${current_run_jobs% }" >>"$job_id_file"
#sbatch --job-name="SASA_calc" \
#               --mail-user=yuanmis1@uci.edu --mail-type=ALL \
#               --account=dtobias_lab \
#               --partition=standard \
#               --nodes=1 \
#               --ntasks-per-node=1 \
#               --cpus-per-task=10 \
#               --time=24:00:00 \
#               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SASAcalc.sh"
