#!/bin/bash

# Define the full path to the ITP file
ITP_FILE="/dfs9/tw/yuanmis1/mrsec/FFssFF/CG_single/martini/martini_v2.1.itp"
OUTPUT_FILE="vdwradii.dat"

# VDW radii based on bead type
R_RADIUS=0.264
S_RADIUS=0.230
T_RADIUS=0.191

# Initialize the output file
echo -e ";resname\tatom_name\tvdw_radius" >$OUTPUT_FILE

# Flag to indicate if in [ atomtypes ] section
in_atomtypes_section=false

# Read the ITP file line by line
while IFS= read -r line; do
    # Trim leading and trailing whitespace
    line=$(echo "$line" | xargs)

    # Check for the start of the [ atomtypes ] section
    if [[ $line =~ ^\[.*atomtypes.*\] ]]; then
        in_atomtypes_section=true
        continue
    fi

    # If another section starts, exit the atomtypes section
    if [[ $line =~ ^\[.*\] ]]; then
        in_atomtypes_section=false
    fi

    # If in the atomtypes section, process lines
    if $in_atomtypes_section; then
        # Ignore comments and empty lines
        if [[ -z $line || $line =~ ^\; ]]; then
            continue
        fi

        # Extract the atom name (first column)
        atom_name=$(echo $line | awk '{print $1}')

        # Determine the VDW radius based on atom name
        if [[ $atom_name =~ ^S ]]; then
            vdw_radius=$S_RADIUS
        elif [[ $atom_name =~ ^T ]]; then
            vdw_radius=$T_RADIUS
        else
            vdw_radius=$R_RADIUS
        fi

        # Write to the output file
        echo -e "???\t$atom_name\t$vdw_radius" >>$OUTPUT_FILE
    fi
done <"$ITP_FILE"

echo -e "???\tBB\t$R_RADIUS" >>$OUTPUT_FILE

echo "VDW data written to $OUTPUT_FILE"
