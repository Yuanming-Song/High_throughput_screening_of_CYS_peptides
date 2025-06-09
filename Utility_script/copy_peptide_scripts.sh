#!/bin/bash

# Define peptide lengths and their names
declare -A peptide_names=(
    [5]="Penta"
   # [6]="Hexa"
   # [7]="Hepta"
   # [8]="Octa"
   # [9]="Nona"
   # [10]="Deca"
)

# Base directory
BASE_DIR="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide"
TETRA_DIR="${BASE_DIR}/Tetrapeptide"

for length in {5}; do
    name="${peptide_names[$length]}peptide"
    target_dir="${BASE_DIR}/${name}"
    
    echo "Processing ${name}..."
    
    # Copy AA_template_script directory
    cp -r "${TETRA_DIR}/AA_template_script" "${target_dir}/"
    
    # Copy MARTINI_Setup_script directory
    cp -r "${TETRA_DIR}/MARTINI_Setup_script" "${target_dir}/"
    
    # Copy MARTINI_util_script directory
    cp -r "${TETRA_DIR}/MARTINI_util_script" "${target_dir}/"
    
    # Modify psfgen template script
    sed -i.bak "s/Tetrapeptide/${name}/g" "${target_dir}/AA_template_script/psfgen_Dimer_Tetrapeptide_template.tcl"
    
    # Update the secondary structure string length in MARTINI scripts
    sec_struct=$(printf 'E%.0s' $(seq 1 $length))
    for script in "${target_dir}/MARTINI_Setup_script/"*peptide_MARTINI_run.sh; do
        if [ -f "$script" ]; then
            sed -i.bak "s/secstructure=\"EEEE\"/secstructure=\"${sec_struct}\"/g" "$script"
            # Rename the script to match the new peptide length
            new_script=$(echo "$script" | sed "s/tetrapeptide/${name,,}/g")
            mv "$script" "$new_script"
        fi
    done
    
    # Clean up backup files
    find "${target_dir}" -name "*.bak" -delete
    
    echo "Completed ${name}"
done

echo "All peptide script directories have been created and modified." 