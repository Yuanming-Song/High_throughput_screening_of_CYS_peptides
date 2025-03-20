#!/bin/bash


amino_acids1=("N")
# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y")

# Define arrays for positively and negatively charged residues
positive_charged=("K" "R")
negative_charged=("D" "E")

# Define the file to store job IDs
job_id_file="submitted_job_ids.txt"

# Clear the file if it exists
current_run_jobs=""

# File containing job IDs
job_id_file="/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Tripeptide_dis_C3/submitted_job_ids.txt"

# Read the last line of the file for dependency job IDs
dependency_jobs=$(tail -n 1 "$job_id_file")

# Ensure the job IDs are available
if [[ -z "$dependency_jobs" ]]; then
    echo "No job IDs found in the last line of $job_id_file."
    exit 1
fi

# Convert the dependency jobs to an array
read -ra job_ids <<<"$dependency_jobs"

# Initialize a counter for job dependencies
counter=0
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

# Loop through each amino acid for amino1
for amino1 in "${amino_acids1[@]}"; do
    # Loop through each amino acid for amino2
    for amino2 in "${amino_acids[@]}"; do
        # Determine the GPU type based on the charge of amino1 and amino2
        if (is_in_array "$amino1" "${positive_charged[@]}" || is_in_array "$amino2" "${positive_charged[@]}") ||
            (is_in_array "$amino1" "${negative_charged[@]}" || is_in_array "$amino2" "${negative_charged[@]}"); then
            gpu_type="A100:1"
        else
            gpu_type="A30:1"
        fi
    # Use the counter to pick a dependency job
        dependency_job="${job_ids[$counter]}"
		    #--dependency=afterany:$dependency_job \
                #--deadline=6:00 \

        # Submit the job using sbatch
        job_id=$(
            sbatch --parsable \
            --job-name="${amino1}_${amino2}_C_mon" \
                --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
                --account=dtobias_lab \
                --partition=free_gpu \
                --nodes=1 \
                --ntasks-per-node=1 \
                --cpus-per-task=8 \
                --time=5:00:00 \
                --out=out/${amino1}_${amino2}_C.out \
                --gres=gpu:$gpu_type \
                --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide/SIRAH/Cys_unstapled_tripeptide_SIRAH_run_GPU.sh ${amino1} ${amino2} C "
        )
        # Check if the job submission was successful
        if [ $? -eq 0 ]; then
            # Append the job ID to the file
            current_run_jobs+="$job_id "
        else
            echo "Failed to submit job for ${amino1}_${amino2}_C_mon"
        fi
        # Increment the counter
        counter=$((counter + 1))

        # Reset the counter if it exceeds the number of dependency jobs
        if [[ $counter -ge ${#job_ids[@]} ]]; then
            counter=0
        fi
    done
done
# Remove the trailing space and append the collected job IDs as a new line in the file
echo "${current_run_jobs% }" >>"$job_id_file"
#--begin=23:00 --deadline=5:00 --dependency=after:34323389 34323497 34323496 34323495 34323494 34323493 34323492 34323491 34323490 34323489 34323488 34323487 34323486 34323485 34323484 34323483 34323482 34323481 34323480 34323479 \
