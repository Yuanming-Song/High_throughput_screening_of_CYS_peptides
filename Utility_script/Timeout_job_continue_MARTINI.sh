#!/usr/bin/env bash
set -euo pipefail

# Toggle debug mode: 
#  - when DEBUG=true, only print the first sbatch command
#  - when DEBUG=false, actually submit all jobs
DEBUG=false

workdir="$PWD"
out_dir="$workdir/out"
avoid_file="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/nodes_to_avoid_gromacs"

# build comma‑separated exclude list
if [[ -r "$avoid_file" ]]; then
  exclude_nodes=$(paste -sd, "$avoid_file")
  exclude_arg=(--exclude="$exclude_nodes")
else
  exclude_arg=()
fi

base_dir=$(basename "$workdir")          # e.g. Tetrapeptide_mon_C3
sys="${base_dir#Tetrapeptide_}"          # e.g. mon_C3

first_printed=false

for logfile in "$out_dir"/*.out; do
  # only consider time‑limit cancellations
  if tail -n1 "$logfile" | grep -q "CANCELLED.*DUE TO TIME LIMIT"; then
    seqname=$(basename "$logfile" .out)
    subdir="$workdir/$seqname"

# determine whether to resume or start fresh by picking the latest .cpt if any exist
latest_cpt=$(ls -1t "${subdir}/${seqname}_CG_md"*cpt 2>/dev/null | head -n1 || true)
if [[ -f "$latest_cpt" ]]; then
  wrap_cmd="module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1; \
cd \"$subdir\"; \
cp \"${latest_cpt}\" \"${seqname}_CG_md.cpt\"; \
gmx_mpi mdrun -deffnm \"${seqname}_CG_md\" -cpi \"${seqname}_CG_md.cpt\" -append > MD_conti.log 2>&1"
else
      wrap_cmd="module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1; \
cd \"$subdir\"; \
gmx_mpi mdrun -deffnm \"${seqname}_CG_md\" > MD.log 2>&1"
    fi

    sbatch_cmd=( sbatch
      --job-name="${seqname}_${sys}_conti"
      --mail-user=yuanmis1@uci.edu
      --mail-type=FAIL
      --account=dtobias_lab
      --partition=standard
      --nodes=1
      --ntasks-per-node=1
      --cpus-per-task=40
      --time=24:00:00
      --output="out/${seqname}.out"
      "${exclude_arg[@]}"
      --wrap="$wrap_cmd"
    )

    if $DEBUG; then
      if ! $first_printed; then
        printf "Would run:\n  %s\n" "${sbatch_cmd[*]}"
        first_printed=true
        exit 0
      fi
      # do not actually sbatch
    else
      "${sbatch_cmd[@]}"
    fi
  fi
done