#!/bin/bash

# Remote and local paths
REMOTE_USER="yuanmis1"
REMOTE_HOST="gplogin3.ps.uci.edu"
REMOTE_DIR="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/analysis/cssc_50mM_cg"
LOCAL_DIR="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/CG"

# Create local directory if it doesn't exist
mkdir -p "$LOCAL_DIR"

# List of files to copy
FILES=(
  "getEdgesSmallsirah.tcl"
  "data/CSSC_50mM_cg_sol_ion.psf.dat.nodes"
  "data/CSSC_50mM_cg_sol_ion.psf.nodes"
  "../analysis/getDistAveAllv3.tcl"
  "../analysis/myMeasureContacts.tcl"
)

for FILE in "${FILES[@]}"; do
  BASENAME=$(basename "$FILE")
  if [ -f "$LOCAL_DIR/$BASENAME" ]; then
    echo "$BASENAME already exists in CG, skipping copy."
  else
    # Determine the correct remote directory
    if [[ "$FILE" == ../analysis/* ]]; then
      REMOTE_FILE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/analysis/${FILE##../analysis/}"
    else
      REMOTE_FILE="$REMOTE_DIR/$FILE"
    fi
    scp "$REMOTE_USER@$REMOTE_HOST:$REMOTE_FILE" "$LOCAL_DIR/"
  fi

done 