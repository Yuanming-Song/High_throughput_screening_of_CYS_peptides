#!/bin/bash

# This script copies required .rda files from the remote server to the local AA and CG directories.
# Directories: /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/AA and CG
# Update the README after running this script.

REMOTE_USER="yuanmis1"
REMOTE_HOST="gplogin3.ps.uci.edu"
REMOTE_BASE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC"

# Local directories
LOCAL_BASE="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC"
AA_DIR="$LOCAL_BASE/AA"
CG_DIR="$LOCAL_BASE/CG"

# Create directories if they do not exist
mkdir -p "$AA_DIR"
mkdir -p "$CG_DIR"

# AA files
scp "$REMOTE_USER@$REMOTE_HOST:$REMOTE_BASE/aa/analysis/csh32_50mM/data/csh32_50mM_every10.edgel.stack.rda" "$AA_DIR/"
scp "$REMOTE_USER@$REMOTE_HOST:$REMOTE_BASE/aa/analysis/cssc16_50mM/data/cssc16_50mM_every10.edgel.stack.rda" "$AA_DIR/"

# CG files
scp "$REMOTE_USER@$REMOTE_HOST:$REMOTE_BASE/cg/analysis/csh_50mM_cg/data/csh_50mM_cg_big_every100_100to407.edgel.stack.rda" "$CG_DIR/"
scp "$REMOTE_USER@$REMOTE_HOST:$REMOTE_BASE/cg/analysis/cssc_50mM_cg/data/cssc_50mM_cg_big_every100_1to359.edgel.stack.rda" "$CG_DIR/"

# Reminder
echo "All files copied. Please update the README to document this data transfer and script usage." 