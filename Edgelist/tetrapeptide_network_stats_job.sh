#!/bin/bash
# Define arrays for positions and first non-C residues (inpres)
positions=(2 3 4)
inpres=("A" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "E" "D")

# Set single_node parameter (TRUE by default)
single_node="TRUE"

# Set paths
script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Rscript/network_stats_tetrapeptide.R"
log_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/log"
err_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/err"

# Create output directories if they don't exist
mkdir -p "${log_dir}" "${err_dir}"

# Create array to store job IDs
declare -a job_ids

# Loop over positions and inpres combinations
for pos in "${positions[@]}"; do
    for inpres_val in "${inpres[@]}"; do
        jobname="C${pos}_${inpres_val}_network_stats"
        out_log="${log_dir}/${jobname}.out"
        err_log="${err_dir}/${jobname}.err"
        
        # Submit job and capture job ID
        job_id=$(sbatch --parsable \
                       --job-name="${jobname}" \
                       --mail-user=yuanmis1@uci.edu \
                       --mail-type=FAIL \
                       --account=dtobias_lab \
                       --partition=standard \
                       --nodes=1 \
                       --ntasks-per-node=1 \
                       --cpus-per-task=40 \
                       --time=24:00:00 \
                       --output="${out_log}" \
                       --error="${err_log}" \
                       --wrap="module load R; R -f ${script} --args ${pos} ${inpres_val} ${single_node}")
        
        # Store job ID
        job_ids+=($job_id)
    done
done

# Create a comma-separated list of job IDs
dependency_list=$(IFS=,; echo "${job_ids[*]}")

# Export the dependency list to a file for the combine script to use
echo "${dependency_list}" > "${log_dir}/network_stats_job_ids.txt"

echo "All jobs submitted. Job IDs saved to ${log_dir}/network_stats_job_ids.txt" 
