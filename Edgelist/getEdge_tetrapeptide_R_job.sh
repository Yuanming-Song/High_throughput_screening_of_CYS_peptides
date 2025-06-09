#!/bin/bash
# Define arrays for positions and first non-C residues (inpres)
positions=(2 3 4)
#2 3 4)
inpres=("A" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "E" "D" "C")
# Set paths
script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Rscript/getEdge_tetrapeptide_martini.R"
log_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/log"
err_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/err"

# Loop over positions and inpres combinations
for pos in "${positions[@]}"; do
    for inpres_val in "${inpres[@]}"; do
        jobname="${pos}_${inpres_val}_edge_tetra"
        out_log="${log_dir}/${jobname}.out"
        err_log="${err_dir}/${jobname}.err"
        sbatch --job-name="${jobname}" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
               --account=dtobias_lab \
               --partition=standard \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=40 \
               --time=12:00:00 \
               --output="${out_log}" \
               --error="${err_log}" \
               --wrap="module load R; R -f ${script} --args ${pos} ${inpres_val}"
    done
done
