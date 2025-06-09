# CSH-CSSC SASA Analysis

## Script: copy_and_sasa.sh

This script automates the process of copying the required .gro files and `vdwradii.dat`, running SASA analysis for CSH and CSSC, and calculating the final/initial SASA ratios for each system.

### What it does:
- Copies initial and final .gro files for both CSH and CSSC from the remote server.
- Copies `vdwradii.dat` from the Dipeptide SASA_score_SIRAH directory.
- Runs SASA analysis for both initial and final structures using `SASAcalc_SIRAH_mon_dipeptide.sh`.
- Calculates and prints the ratio of final to initial SASA for both CSH and CSSC.

### Usage

```bash
bash copy_and_sasa.sh
```

The script will print the SASA ratios to the terminal.

---

**Note:**
- Ensure you have SSH access to the remote server and the required permissions to copy files.
- The SASA analysis script (`SASAcalc_SIRAH_mon_dipeptide.sh`) must be present in `../Script/` relative to this directory.
- The script assumes `vdwradii.dat` is available at `/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_SIRAH/vdwradii.dat`. 