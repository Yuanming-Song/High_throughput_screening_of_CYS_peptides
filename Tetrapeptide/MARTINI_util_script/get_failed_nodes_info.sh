#!/usr/bin/env bash
set -euo pipefail

out_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI22/Tetrapeptide_mon_C3/out"
output_file="./nodes_to_avoid_gromacs"

declare -A bad_nodes

# Loop over all .out files
for logfile in "$out_dir"/*.out; do
  # Check if the last 20 lines contain “error” (case‐insensitive)
  if tail -n 20 "$logfile" | grep -qi "error"; then
    # From those lines, extract the first occurrence of a “node” entry
    # e.g. a line like “Node: nid1234” or “Running on node nid1234”
    node_line=$(tail -n 20 "$logfile" | grep -i "node" | head -n1)
    # Pull out the node name (assumes it’s the word after “node” or “Node:”)
    node=$(echo "$node_line" | awk '{
      for(i=1;i<NF;i++){
        if(tolower($i)=="node" || tolower($i)=="node:"){
          print $(i+1)
          exit
        }
      }
    }')
    # If we found something, record it
    if [[ -n "$node" ]]; then
      bad_nodes["$node"]=1
    fi
  fi
done

# Write unique node names to file
touch "$output_file"
for n in "${!bad_nodes[@]}"; do
  echo "$n"
done >> "$output_file"

# Remove redundancy: sort and unique
sort -u "$output_file" -o "$output_file"