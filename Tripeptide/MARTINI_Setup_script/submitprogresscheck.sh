sbatch --job-name="CG Sim Check" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
               --account=dtobias_lab \
               --partition=standard \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=1 \
               --time=24:00:00 \
               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/checkprogress.sh"