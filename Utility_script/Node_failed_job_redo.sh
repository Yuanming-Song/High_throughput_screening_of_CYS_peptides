#!/usr/bin/env bash
set -euo pipefail

DEBUG=false
workdir="$PWD"
out_dir="$workdir/out"
avoid_file="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/nodes_to_avoid_gromacs"

# build exclude list
if [[ -r "$avoid_file" ]]; then
  exclude_nodes=$(paste -sd, "$avoid_file")
  exclude_arg=( --exclude="$exclude_nodes" )
else
  exclude_arg=()
fi

base=$(basename "$workdir")
if [[ "$base" == *mon* ]]; then
  wrapper="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/Cys_unstapled_tetrapeptide_MARTINI_run.sh"
elif [[ "$base" == *dis* ]]; then
  wrapper="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/Cys_stapled_tetrapeptide_MARTINI_run.sh"
else
  echo "ERROR: Cannot determine state from '$base'" >&2
  exit 1
fi

first_printed=false

for logfile in "$out_dir"/*.out; do
  if tail -n2 "$logfile" | grep -q "ERROR: MD simulation failed"; then
    seqname=$(basename "$logfile" .out)
    IFS='_' read -r aa1 aa2 aa3 aa4 <<< "$seqname"

    # build wrap string: remove old dir, then run wrapper
    wrap_cmd="rm -rf \"$workdir/$seqname\"; bash \"$wrapper\" $aa1 $aa2 $aa3 $aa4"

    sbatch_cmd=(
      sbatch
      --job-name="${seqname}_retry"
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
    else
      "${sbatch_cmd[@]}"
    fi
  fi
done