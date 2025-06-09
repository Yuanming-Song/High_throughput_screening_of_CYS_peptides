# ML-MD-Peptide

This repository contains scripts and analysis tools for coarse-grained molecular dynamics (CGMD) simulations of peptide systems.

Each subdirectory contains its own README with detailed information about the scripts, file formats, and analysis procedures specific to that directory.

## Directory Structure

### Simulation Setup Scripts
- `Dipeptide/`, `Tripeptide/`, `Tetrapeptide/`: Setup and slurm job submission scripts for CGMD simulations
  - Covers all single cysteine-containing sequences
  - Both monomer and disulfide-linked dimer simulations
  - Note: Only setup and analysis scripts are tracked in git, not simulation data
  - See README in each directory for specific setup procedures

### Analysis Tools
- `AP_analysis/`: Aggregation propensity score calculation
  - SASA calculation scripts for different peptide lengths
  - Analysis tools for different force fields (MARTINI 2.2, 2.1, 3, SIRAH)
  - Detailed documentation in directory's README

- `Csizedist/`: Cluster size distribution analysis
  - RDA files format: first column is cluster size, other columns are density
  - Column names indicate peptide sequence
  - Analysis scripts for different force fields
  - See directory's README for detailed data format specifications

- `csize_SIRAH/`: Cluster size distribution analysis specific to SIRAH force field
  - Analysis scripts customized for SIRAH force field simulations
  - Similar format to Csizedist but with SIRAH-specific adaptations
  - See directory's README for SIRAH-specific analysis details

- `Edgelist/`: Network analysis tools
  - Scripts for generating and analyzing peptide interaction networks
  - Edge list generation for each peptide system
  - Refer to directory's README for network analysis procedures

- `gyrT/`: Gyration tensor analysis
  - Tools for analyzing aggregate shapes
  - Gyration tensor calculation scripts
  - Check directory's README for calculation methods

Note: While simulation directories contain extensive data files, this repository only tracks setup and analysis scripts. Large data files, simulation results, and directories like DL_for_Peptide, Rebinned_Data, rsync_crsp, Judred, and ML_MD_project_base_R are not included in version control.

## Dipeptide/SIRAH/Box_size/submitjob_boxsize.sh

This script automates the setup and submission of dipeptide simulations for different box sizes. For each box size (currently placeholders: 10, 12, 14), it creates a directory (e.g., box10), then creates 'mon' (monomer) and 'dis' (dimer) subdirectories. In each subdirectory, it runs both FC and CF sequences using the updated SIRAH setup scripts, passing the box size as an argument.

**Usage:**
```bash
bash submitjob_boxsize.sh
```

This will:
- Create directories for each box size
- Run both monomer and dimer simulations for FC and CF in each box size
- Use the correct SIRAH setup scripts for each case

Edit the box_sizes array in the script to change the box sizes as needed.
