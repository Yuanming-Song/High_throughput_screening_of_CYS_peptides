#!/bin/bash
#########################################################################
# Script for setting up and running MARTINI coarse-grained MD simulations
# of stapled peptides.
#
# Usage examples:
#   For tripeptide with C at position 1:
#     ./script.sh C A D
#     -> Creates directory: C_A_D
#     -> Uses template: .../Tripeptide/AA_template_dimer_C1/
#
#   For pentapeptide with C at position 3:
#     ./script.sh A D C E F
#     -> Creates directory: A_D_C_E_F
#     -> Uses template: .../Pentapeptide/AA_template_dimer_C3/
#
# Input format:
#   - Takes 2-10 single-letter amino acid codes as arguments
#   - Must include exactly one Cysteine (C) for stapling
#   - Directory name will be created by joining all residues with underscore
#   - Example: A B C D -> directory name: A_B_C_D
#
# Secondary structure:
#   - Automatically generates secondary structure string
#   - Length matches peptide length
#   - Example: for 4 residues -> "EEEE"
#########################################################################

# Base directory for all tools and resources
BaseDir="/dfs9/tw/yuanmis1"
# box size and simulation parameters
MartiniDir="${BaseDir}/mrsec/FFssFF/CG_single/ubiquitin/martini/"
MartiniDir_mdp_file="${BaseDir}/mrsec/ML-MD-Peptide/Tripeptide/MARTINI_Setup_script/dynamic.mdp"
Martini2psf="${BaseDir}/mrsec/ML-MD-Peptide/Tripeptide/MARTINI_Setup_script/getMARTINIpsf.sh"
makeAAtemplate="${BaseDir}/mrsec/ML-MD-Peptide/Tetrapeptide/AA_template_script/mutate_template.tcl"
vmd_exc="${BaseDir}/vmd/bin/vmd"
martinize_exc="${BaseDir}/tools/martinize.py/martinize.py"
box_size_nm=14.3
ncopy=150

# Print node information
echo "Running on node: $(hostname)"

# Check argument count (2-10 residues allowed)
if [ "$#" -lt 2 ] || [ "$#" -gt 10 ]; then
    echo "Usage: $0 <residue1> [residue2 ... residue10]"
    echo "Examples:"
    echo "  $0 C A D        -> Creates C_A_D (tripeptide)"
    echo "  $0 A D C E F    -> Creates A_D_C_E_F (pentapeptide)"
    echo "Note: Must include exactly one Cysteine (C) for stapling"
    exit 1
fi

# Count Cysteines in input
cys_count=0
cys_pos=-1
pos=1
for arg in "$@"; do
    if [ "$arg" == "C" ]; then
        ((cys_count++))
        cys_pos=$pos
    fi
    ((pos++))
done

# Verify exactly one Cysteine
if [ "$cys_count" -ne 1 ]; then
    echo "Error: Must have exactly one Cysteine (C) for stapling"
    echo "Found $cys_count Cysteine(s) in input sequence"
    exit 1
fi

# Create directory name by joining residues with underscore
# Example: A B C D -> A_B_C_D
dir_name=""
sequence=""
for arg in "$@"; do
    if [ -z "$dir_name" ]; then
        dir_name="$arg"
        sequence="$arg"
    else
        dir_name="${dir_name}_${arg}"
        sequence="${sequence}${arg}"
    fi
done

# Generate secondary structure string (all E's, matching peptide length)
# Example: 4 residues -> "EEEE"
secstructure=""
for ((i=1; i<=$#; i++)); do
    secstructure="${secstructure}E"
done

echo "Processing sequence: $sequence"
echo "Secondary structure: $secstructure"
echo "Cysteine position: $cys_pos"
echo "Directory name: $dir_name"

# Check if the directory exists
if [ -d "$dir_name" ]; then
    echo "Directory '$dir_name' already exists, submit manually."
    exit 0
fi

# Create the directory and enter it
mkdir "$dir_name"
cd "$dir_name"

# Use template to create initial structure
# Example: for "A B C D" with C at pos 3, calls:
# vmd -dispdev none -e script.tcl -args A B C D
${vmd_exc} -dispdev none -e ${makeAAtemplate} -args "$@" >mutate_${dir_name}.log 2>&1

# Check if any of the arguments is H (Histidine)
if [[ "$*" == *"H"* ]]; then
    pdb_file="${dir_name}.pdb"
    if [[ -f "$pdb_file" ]]; then
        sed -i 's/HSD/HIS/g' "$pdb_file"
        echo "Replaced HSD with HIS in $pdb_file"
    else
        echo "File $pdb_file does not exist."
    fi
fi

# Martinize the structure
module load python/2.7.17
python2.7 ${martinize_exc} -f ${dir_name}.pdb -o ${dir_name}_CG.top -x ${dir_name}_CG.pdb -ff martini22 -v -cys auto -name ${dir_name} -ss ${secstructure} >${dir_name}_CG_martinize.log 2>&1

module load gromacs/2022.1/gcc.8.4.0-cuda.11.7.1

# Insert multiple copies of the molecule
gmx insert-molecules -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -nmol ${ncopy} -ci ${dir_name}_CG.pdb -radius 0.3 -o ${dir_name}_CG_Lattice.gro >${dir_name}_CG_multi_copy_setup.log 2>&1

echo "Converting pdb to gro file is done, check editconf.log file"

# Solvate the system
gmx solvate -cp ${dir_name}_CG_Lattice.gro -cs ${MartiniDir}water.gro -radius 0.21 -o ${dir_name}_CG_solvated.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} >solvate.log 2>&1
echo "Solvating solute in water is done, check solvate.log file"

# Process topology file and add water
if [ -f "solvate.log" ]; then
    # Get number of water molecules
    watnum=$(grep 'Number of solvent molecules:' solvate.log | awk '{print $5}')

    # Update topology with number of copies and water
    sed -i "/${dir_name}_A+${dir_name}_B/s/1/${ncopy}/" ${dir_name}_CG.top
    echo -e "\nW ${watnum}" >>${dir_name}_CG.top

    # Add required includes
    sed -i "1i $lines_to_insert" ${dir_name}_CG.top
    sed -i '/#include "martini.itp"/d' ${dir_name}_CG.top
    echo "Number of water molecules added: $watnum"
else
    echo "Error: solvate.log not found or does not contain expected output."
    exit 1
fi

# Run minimization
gmx grompp -p ${dir_name}_CG.top -c ${dir_name}_CG_solvated.gro -f ${MartiniDir}minimization.mdp -o minimization.tpr >gromppMin.log 2>&1
echo "grompp minimization is done, check gromppMin.log file"

# Add ions for neutralization
echo 13 | gmx genion -s minimization.tpr -o ${dir_name}_CG_solvated_ions.gro -neutral -rmin 0.5 -p ${dir_name}_CG.top -nname CL- -pname NA+ >genion.log 2>&1
echo "Adding ion is done, check genion.log file"

# Run minimization after adding ions
gmx grompp -p ${dir_name}_CG.top -c ${dir_name}_CG_solvated_ions.gro -f ${MartiniDir}minimization.mdp -o minimization.tpr -maxwarn 1 >gromppMinIon.log 2>&1
echo "grompp minization of neutralized is done, check gromppMinIon.log file"

gmx mdrun -deffnm minimization -v >minimization.log 2>&1
echo "Minimization is done, check minimization.log file"

# Run equilibration
gmx grompp -p ${dir_name}_CG.top -c minimization.gro -f ${MartiniDir}equilibration.mdp -o equilibration.tpr -maxwarn 2 >gromppEqu.log 2>&1
echo "grompp equalibration is done, check gromppEqu.log file"

gmx mdrun -deffnm equilibration -v >equilibration.log 2>&1
echo "Equalibration is done, check equilibration.log file"

# Setup and run production MD
gmx grompp -p ${dir_name}_CG.top -c equilibration.gro -f ${MartiniDir_mdp_file} -o ${dir_name}_CG_md.tpr -maxwarn 1 >gromppMD.log 2>&1
echo "grompp MD production is done, check gromppMD.log file"

# Generate visualization files
${Martini2psf} -b ${dir_name}_CG
echo "MD production run is going, check ${dir_name}_CG_md.log"

# Run production MD
gmx mdrun -nt $SLURM_CPUS_PER_TASK -deffnm ${dir_name}_CG_md >${dir_name}_CG_md.log 2>&1

# Verify simulation completed
if ! ls *md*.gro >/dev/null 2>&1; then
    echo "ERROR: MD simulation failed - no *md*.gro file found"
    echo "Node: $(hostname)"
    exit 1
fi
