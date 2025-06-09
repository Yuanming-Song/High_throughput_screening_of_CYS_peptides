sbatch --job-name="sync_csh-aa_remove" \
       --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
       --account=dtobias_lab \
       --partition=free \
       --nodes=1 \
       --ntasks-per-node=1 \
       --cpus-per-task=20 \
       --time=24:00:00 \
       --out=sync_csh-aa_remove.log \
       --wrap="rsync -avzh --progress --stats --remove-source-files /dfs9/tw/yuanmis1/mrsec/CSH-CSSC/aa yuanmis1@gplogin3.ps.uci.edu:/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/. && find /dfs9/tw/yuanmis1/mrsec/CSH-CSSC/aa -type d -empty -delete"
