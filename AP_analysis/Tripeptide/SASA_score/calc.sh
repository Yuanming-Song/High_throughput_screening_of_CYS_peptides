sbatch --job-name="SASA_calc" \
               --mail-user=yuanmis1@uci.edu --mail-type=FAIL\
               --account=dtobias_lab \
               --partition=free \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=5 \
               --time=24:00:00 \
               --out=SASA_cal_dim.out \
               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Script/SASAcalc.sh"
