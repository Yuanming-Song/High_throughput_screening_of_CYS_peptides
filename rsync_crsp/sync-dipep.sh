source="/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/MARTINI21"
target="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/."
jobname="sync_dipep"
sbatch --job-name="$jobname" \
       --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
       --account=dtobias_lab \
       --partition=free \
       --nodes=1 \
       --ntasks-per-node=1 \
       --cpus-per-task=40 \
       --time=24:00:00 \
       --out=${jobname}.log \
       --wrap="rsync -avzh --progress --stats ${source} ${target}"
