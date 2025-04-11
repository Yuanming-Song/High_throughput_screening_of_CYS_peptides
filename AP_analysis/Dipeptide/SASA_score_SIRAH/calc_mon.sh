sbatch --job-name="mon_di_SASA_calc_SIRAH_mon" \
		       --mail-user=yuanmis1@uci.edu --mail-type=FAIL\
               --account=dtobias_lab \
               --partition=standard \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=5 \
               --time=24:00:00 \
               --out=SASA_calc_SIRAH_mon.log \
               --wrap="bash /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Script/SASAcalc_SIRAH_mon_dipeptide.sh "
