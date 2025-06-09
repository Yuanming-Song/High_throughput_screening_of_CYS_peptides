#!/bin/bash
positions=(1)
inpres=("A" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "E" "D")

script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/gyrT/Rscript/getGyrT_tetrapeptide.R"
log_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/gyrT/log"
err_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/gyrT/err"

mkdir -p "$log_dir" "$err_dir"

for pos in "${positions[@]}"; do
    for inpres_val in "${inpres[@]}"; do
        jobname="C${pos}_${inpres_val}_gyrT_tetra"
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