#!/bin/bash

# Set paths
script="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Rscript/combine_network_stats_tetrapeptide.R"
log_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/log"
err_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/err"

# Read dependency list from file
dependency_list=$(cat "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/log/network_stats_job_ids.txt")

# Submit combine job for each position (1-4)
for pos in {2..4}
do
    jobname="combine_network_stats_C${pos}"
    out_log="${log_dir}/${jobname}.out"
    err_log="${err_dir}/${jobname}.err"

    sbatch --job-name="${jobname}" \
           --mail-user=yuanmis1@uci.edu \
           --mail-type=END,FAIL \
           --account=dtobias_lab \
           --partition=standard \
           --nodes=1 \
           --ntasks-per-node=1 \
           --cpus-per-task=20 \
           --time=2:00:00 \
           --output="${out_log}" \
           --error="${err_log}" \
           --wrap="module load R; R -f ${script} --args ${pos} 1 TRUE ${pos}"
           #--dependency=afterok:${dependency_list} 

done

echo "Combine jobs submitted for positions 1-4 with dependency on jobs: ${dependency_list}" 
