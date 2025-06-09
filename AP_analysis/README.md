# AP Analysis Directory

This directory contains scripts and data for analyzing Aggregation Propensity (AP) scores of peptides with different force fields and configurations. The AP score is calculated as the ratio of Solvent Accessible Surface Area (SASA) at the beginning and end of each simulation, providing a measure of how much the peptide's surface exposure changes during aggregation.

## Directory Structure

```
AP_analysis/
├── Script/                 # Analysis scripts
├── Dipeptide/             # Data for dipeptide analysis
│   ├── SASA_score_M21/    # MARTINI 2.1 force field results
│   ├── SASA_score_SIRAH/  # SIRAH force field results
├── Tripeptide/            # Data for tripeptide analysis
│   ├── SASA_score/        # Default MARTINI 2.2 results
│   ├── SASA_score_M3/     # MARTINI 3 results
│   ├── SASA_score_SIRAH/  # SIRAH force field results
└── Tetrapeptide/          # Data for tetrapeptide analysis
    ├── SASA_score/        # Default MARTINI 2.2 results
      └── SubData/          # Parallel processing results by C position
```

## Usage Notes

1. Before running any MARTINI SASA calculations:
   - Run `vdwradiigenerator.sh` first
   - This generates required van der Waals radii parameters

2. For SASA calculations:
   - Use appropriate script based on peptide length (di/tri/tetra)
   - Use `_mon` scripts for monomers
   - Use non-`_mon` scripts for disulfide-linked dimers

3. For common ini SASA AP score calculations:
   - Use appropriate `SASA_cal_mon_ini_*.sh` script based on peptide length
   - These scripts use the initial monomer SASA as reference
   - They calculate both monomer AP (initial/final monomer SASA) and dimer AP (initial monomer/final dimer SASA) 

## Script Directory

The Script directory contains various SASA calculation scripts for different peptide lengths and force fields:

### SASA Calculation Scripts

#### MARTINI 2.2 Force Field
- `SASAcalc.sh`: Tripeptide disulfide linked dimer SASA calculation script
- `SASAcalc_mon.sh`: tripeptide monomer SASA calculation script
- `SASAcalc_dipeptide.sh`: For disulfide-linked dipeptide dimers
- `SASAcalc_dipeptide_mon.sh`: For dipeptide monomers
- `SASAcalc_tripeptide.sh`: For disulfide-linked tripeptide dimers
- `SASAcalc_tripeptide_mon.sh`: For tripeptide monomers
- `SASAcalc_tetrapeptide.sh`: For disulfide-linked tetrapeptide dimers
- `SASAcalc_tetrapeptide_mon.sh`: For tetrapeptide monomers

#### SIRAH Force Field
- `SASAcalc_SIRAH.sh`: Generic SIRAH force field SASA calculation
- `SASAcalc_SIRAH_mon.sh`: Generic SIRAH force field monomer SASA calculation
- `SASAcalc_SIRAH_dipeptide.sh`: For SIRAH force field dipeptide dimers
- `SASAcalc_SIRAH_mon_dipeptide.sh`: For SIRAH force field dipeptide monomers

### AP Score Calculation Scripts

These scripts calculate the Aggregation Propensity score using only the initial SASA of the monomer as reference:
- `SASA_cal_mon_ini_dipeptide.sh`: AP calculation for dipeptides
- `SASA_cal_mon_ini_tripeptide.sh`: AP calculation for tripeptides
- `SASA_cal_mon_ini_tetrapeptide.sh`: AP calculation for tetrapeptides

The AP score is calculated as:
- For monomers: AP = initial SASA / final SASA
- For dimers with monomer reference: AP = initial monomer SASA / final dimer SASA

### Utility Scripts

- `vdwradiigenerator.sh`: Generates van der Waals radii file for MARTINI force field. **Essential** to run this before any MARTINI SASA calculations.

## Data Directories

### Force Field Naming Convention

- Default (no suffix): MARTINI 2.2 force field
- `M21`: MARTINI 2.1 force field
- `SIRAH`: SIRAH force field

### Directory Contents

Each peptide directory (Di/Tri/Tetrapeptide) contains:
- Analysis results organized by force field:
  - `SASA_score/`: Results using MARTINI 2.2 (default)
  - `SASA_score_M21/`: Results using MARTINI 2.1
  - `SASA_score_M3/`: Results using MARTINI 3
  - `SASA_score_SIRAH/`: Results using SIRAH force field
- Force field specific subdirectories for simulation data
- Tetrapeptide only:
`SubData/`: Contains intermediate SASA calculation results (.dat files)
  - Files are organized by cysteine (C) position and the next non-cysteine residue
  - Format: `SASA_result_Tetrapeptide_C{position}_{next_residue}.dat`
  - Enables parallel processing of different sequence combinations
  - Each .dat file contains SASA ratios for sequences with:
    - C at the specified position
    - The specified residue immediately after C
    - All possible combinations of other positions

### File Naming Convention

- Files without `_mon` suffix: For disulfide-linked dimers
- Files with `_mon` suffix: For monomers
- Files with `_ini` suffix: using common initial SASA (monomer simulation setup)

