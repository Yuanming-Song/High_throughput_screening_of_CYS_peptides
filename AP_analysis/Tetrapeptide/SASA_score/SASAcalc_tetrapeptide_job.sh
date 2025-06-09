#!/bin/bash
# Define arrays for positions and first non-C residues
positions=(3)
inpres=("A" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "E" "D")

# Set paths
script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Script/SASAcalc_tetrapeptide.sh"
script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Script/SASAcalc_tetrapeptide_mon.sh"
log_dir="log/"
err_dir="err/"

# Create log and error directories if they don't exist
mkdir -p "$log_dir" "$err_dir"

# Loop over positions and inpres combinations
for pos in "${positions[@]}"; do
    for inpres_val in "${inpres[@]}"; do
        jobname="SASA_C${pos}_${inpres_val}_tetra"
        out_log="${log_dir}/${jobname}_mon.out"
        err_log="${err_dir}/${jobname}_mon.err"
        
        # Submit the job
        sbatch --job-name="${jobname}" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
               --account=dtobias_lab \
               --partition=standard \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=1 \
               --time=24:00:00 \
               --output="${out_log}" \
               --error="${err_log}" \
               --wrap="module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1; ${script} ${pos} ${inpres_val}"
    done
done 
