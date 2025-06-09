#!/bin/bash

# Set variables for remote and local paths
REMOTE_USER="yuanmis1"
REMOTE_HOST="gplogin3.ps.uci.edu"

# CSH paths
CSH_INIT_REMOTE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/csh_50mM_cg/sirah_setup/CSH_50mM_cg_sol_ion.gro"
CSH_FINAL_REMOTE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/csh_50mM_cg/sirah_run/CSH_50mM_cg_md.part0048.gro"
# CSSC paths
CSSC_INIT_REMOTE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/cssc_50mM_cg/sirah_setup/CSSC_50mM_cg_sol_ion.gro"
CSSC_FINAL_REMOTE="/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/cssc_50mM_cg/sirah_run/CSSC_50mM_cg_md.part0058.gro"

# Local destination
LOCAL_DIR="$(dirname "$0")"

# Copy .gro files if not already present
if [ ! -f "$LOCAL_DIR/CSH_50mM_cg_sol_ion.gro" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:$CSH_INIT_REMOTE" "$LOCAL_DIR/CSH_50mM_cg_sol_ion.gro"
else
  echo "CSH_50mM_cg_sol_ion.gro already exists, skipping copy."
fi
if [ ! -f "$LOCAL_DIR/CSH_50mM_cg_md.part0048.gro" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:$CSH_FINAL_REMOTE" "$LOCAL_DIR/CSH_50mM_cg_md.part0048.gro"
else
  echo "CSH_50mM_cg_md.part0048.gro already exists, skipping copy."
fi
if [ ! -f "$LOCAL_DIR/CSSC_50mM_cg_sol_ion.gro" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:$CSSC_INIT_REMOTE" "$LOCAL_DIR/CSSC_50mM_cg_sol_ion.gro"
else
  echo "CSSC_50mM_cg_sol_ion.gro already exists, skipping copy."
fi
if [ ! -f "$LOCAL_DIR/CSSC_50mM_cg_md.part0058.gro" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:$CSSC_FINAL_REMOTE" "$LOCAL_DIR/CSSC_50mM_cg_md.part0058.gro"
else
  echo "CSSC_50mM_cg_md.part0058.gro already exists, skipping copy."
fi

# Copy vdwradii.dat if not already present
if [ ! -f "$LOCAL_DIR/vdwradii.dat" ]; then
  cp /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_SIRAH/vdwradii.dat "$LOCAL_DIR/vdwradii.dat"
else
  echo "vdwradii.dat already exists, skipping copy."
fi

# Copy CSSC md.xtc file if not already present
if [ ! -f "$LOCAL_DIR/CSSC_50mM_cg_md.xtc" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/cssc_50mM_cg/sirah_run/CSSC_50mM_cg_md.xtc" "$LOCAL_DIR/CSSC_50mM_cg_md.xtc"
else
  echo "CSSC_50mM_cg_md.xtc already exists, skipping copy."
fi

# Load GROMACS module (if on cluster)
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1

# Check if CSSC_50mM_cg_md.tpr exists locally, if not, copy it from the remote server
if [ ! -f "$LOCAL_DIR/CSSC_50mM_cg_md.tpr" ]; then
  scp "$REMOTE_USER@$REMOTE_HOST:/DFS-B/DATA/tobias/yuanmis1/mrsec/CSH-CSSC/cg/cssc_50mM_cg/sirah_run/CSSC_50mM_cg_md.tpr" "$LOCAL_DIR/CSSC_50mM_cg_md.tpr"
fi

# If .tpr is still not present, use the initial .gro as a fallback for trjconv
if [ -f "$LOCAL_DIR/CSSC_50mM_cg_md.tpr" ]; then
  gmx trjconv -f CSSC_50mM_cg_md.xtc -s CSSC_50mM_cg_md.tpr -o CSSC_50mM_cg_md_500ns.gro -dump 500000 <<EOF
0
EOF
else
  gmx trjconv -f CSSC_50mM_cg_md.xtc -s CSSC_50mM_cg_sol_ion.gro -o CSSC_50mM_cg_md_500ns.gro -dump 500000 <<EOF
0
EOF
fi

# Run gmx sasa for each .gro file (using itself as both -f and -s)
echo 2 | gmx sasa -f CSH_50mM_cg_sol_ion.gro -s CSH_50mM_cg_sol_ion.gro -o CSH_ini.xvg > CSH_ini_sasa.log 2>&1
echo 2 | gmx sasa -f CSH_50mM_cg_md.part0048.gro -s CSH_50mM_cg_md.part0048.gro -o CSH_final.xvg > CSH_final_sasa.log 2>&1
echo 2 | gmx sasa -f CSSC_50mM_cg_sol_ion.gro -s CSSC_50mM_cg_sol_ion.gro -o CSSC_ini.xvg > CSSC_ini_sasa.log 2>&1
echo 2 | gmx sasa -f CSSC_50mM_cg_md_500ns.gro -s CSSC_50mM_cg_md_500ns.gro -o CSSC_final.xvg > CSSC_final_sasa.log 2>&1
#echo 2 | gmx sasa -f CSSC_50mM_cg_md.part0058.gro -s CSSC_50mM_cg_md.part0058.gro -o CSSC_final.xvg > CSSC_final_sasa.log 2>&1
# Extract SASA values (last row, second column)
CSH_INI=$(awk '/^[^@#]/{val=$2} END{print val}' CSH_ini.xvg)
CSH_FINAL=$(awk '/^[^@#]/{val=$2} END{print val}' CSH_final.xvg)
CSSC_INI=$(awk '/^[^@#]/{val=$2} END{print val}' CSSC_ini.xvg)
CSSC_FINAL=$(awk '/^[^@#]/{val=$2} END{print val}' CSSC_final.xvg)

# Calculate and print ratios
CSH_RATIO=$(echo "$CSH_INI / $CSH_FINAL " | bc -l)
CSSC_RATIO=$(echo "$CSH_INI * 320 / 321 / $CSSC_FINAL " | bc -l) 

echo "CSH final/ini SASA ratio: $CSH_RATIO"
echo "CSSC final/ini SASA ratio: $CSSC_RATIO"

echo "CSH $CSH_RATIO $CSSC_RATIO" > CSH_AP.txt 