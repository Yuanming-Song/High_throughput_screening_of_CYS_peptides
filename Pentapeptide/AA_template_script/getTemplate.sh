#!/bin/bash
# Set minimization input file location (modify as needed)
mininp="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/AA_template_script/min.namd.tcl"
tcldir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/AA_template_script"

# Load NAMD module
module load namd/2.14b2/gcc.8.4.0-cuda.10.1.243


# Loop through indices 1 to 4 for Tetrapeptide
for ind in {1..4}; do
    echo "Processing Tetrapeptide with index: $ind"
    # Create directory for this simulation
    mkdir -p AA_template_dimer_C${ind}
    cd AA_template_dimer_C${ind}
    mkdir -p minimization/

    cd minimization/
    
    # Run the initial psfgen TCL script with the index as argument
    /dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none -e ${tcldir}/psfgen_Dimer_Tetrapeptide_template.tcl -args $ind
    
    # Copy the minimization input file into the current directory
    cp $mininp .
    
    # Run minimization using NAMD (using the basename of the minimization file)
    namd2 min.namd.tcl
    
    # Load minimized configuration and write new pdb using VMD TCL script
    /dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none *psf min.coor -e ${tcldir}/getPDBtemplate.tcl
    
    # Return to the Tetrapeptide directory
    cd ../../
done
