#!/bin/bash
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
# Lines to insert at the beginning of the .top file, for itp files 
lines_to_insert='#include "'${BaseDir}'/mrsec/FFssFF/CG_single/martini/martini_v2.2.itp"\n#include "'${BaseDir}'/mrsec/ff/martini/martini_v2.0_ions.itp"'

# Check if exactly three arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <arg1> <arg2> <arg3> <arg4>"
    exit 1
fi

# Assign arguments to variables
arg1=$1
arg2=$2
arg3=$3
arg4=$4
echo "${arg1}_${arg2}_${arg3}_${arg4}"
# Define the directory name
dir_name="${arg1}_${arg2}_${arg3}_${arg4}"

# Check if the directory exists
if [ -d "$dir_name" ]; then
    echo "Directory '$dir_name' already exists, submit manually."
    #exit 1  # Stop the script with a non-zero exit status
else
    # Create the directory if it does not exist
    mkdir "$dir_name"
    echo "Directory '$dir_name' created successfully."
fi
cd ${arg1}_${arg2}_${arg3}_${arg4}

# Check each argument and update IndexC with the 1-based index
if [ "$arg1" == "C" ]; then
    IndexC=1  
elif [ "$arg2" == "C" ]; then
    IndexC=2
elif [ "$arg3" == "C" ]; then
    IndexC=3
elif [ "$arg4" == "C" ]; then
    IndexC=4
fi

# Use a linear peptide template, mutate residue based on input
${vmd_exc} -dispdev none -e ${makeAAtemplate} -args ${arg1} ${arg2} ${arg3} ${arg4} AA_template_dimer_C$IndexC >mutate_${arg1}_${arg2}_${arg3}_${arg4}.log 2>&1
# Check if any of the arguments is H
if [[ "$arg1" == "H" || "$arg2" == "H" || "$arg3" == "H" || "$arg4" == "H" ]]; then
    # Define the pdb file name based on the arguments
    pdb_file="${arg1}_${arg2}_${arg3}_${arg4}.pdb"

    # Check if the file exists
    if [[ -f "$pdb_file" ]]; then
        # Replace HSD with HIS in the file
        sed -i 's/HSD/HIS/g' "$pdb_file"
        echo "Replaced HSD with HIS in $pdb_file"
    else
        echo "File $pdb_file does not exist."
    fi
fi
# Martinize
module load python/2.7.17
python2.7 ${martinize_exc} -f ${arg1}_${arg2}_${arg3}_${arg4}.pdb -o ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -x ${arg1}_${arg2}_${arg3}_${arg4}_CG.pdb -ff martini22 -v -cys auto -name ${arg1}_${arg2}_${arg3}_${arg4} -ss EEEEEEEE >${arg1}_${arg2}_${arg3}_${arg4}_CG_martinize.log 2>&1
module load gromacs/2022.1/gcc.8.4.0-cuda.11.7.1
#Load gromacs
#Add something to topology file?
# Modify setup script
gmx insert-molecules -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -nmol 150 -ci ${arg1}_${arg2}_${arg3}_${arg4}_CG.pdb -radius 0.3 -o ${arg1}_${arg2}_${arg3}_${arg4}_CG_Lattice.gro >${arg1}_${arg2}_${arg3}_${arg4}_CG_multi_copy_setup.log 2>&1

#gmx editconf -f ${arg1}_${arg2}_${arg3}_${arg4}_CG_Lattice.pdb -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -bt cubic -o ${arg1}_${arg2}_${arg3}_${arg4}_CG.gro >editconf.log 2>&1
echo "Converting pdb to gro file is done, check editconf.log file"
#solvate molecule
gmx solvate -cp ${arg1}_${arg2}_${arg3}_${arg4}_CG_Lattice.gro -cs ${MartiniDir}water.gro -radius 0.21 -o ${arg1}_${arg2}_${arg3}_${arg4}_CG_solvated.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} >solvate.log 2>&1
echo "Solvating solute in water is done, check solvate.log file"
#in solvate.log, find the number comes after the Number of solvent molecules:
# Check if solvate.log exists and contains the required information
if [ -f "solvate.log" ]; then
    # Extract number of solvent molecules from solvate.log
    watnum=$(grep 'Number of solvent molecules:' solvate.log | awk '{print $5}')

    # Replace placeholder in ${arg1}_${arg2}_${arg3}_${arg4}_CG.top with watnum

    # Use sed to find the line, modify the number, and add a new line after it
    sed -i "/${arg1}_${arg2}_${arg3}_${arg4}_A+${arg1}_${arg2}_${arg3}_${arg4}_B/s/1/150/" ${arg1}_${arg2}_${arg3}_${arg4}_CG.top
    echo -e "\nW ${watnum}" >>${arg1}_${arg2}_${arg3}_${arg4}_CG.top

    # Use sed to insert lines at the beginning of the file
    sed -i "1i $lines_to_insert" ${arg1}_${arg2}_${arg3}_${arg4}_CG.top
    sed -i '/#include "martini.itp"/d' ${arg1}_${arg2}_${arg3}_${arg4}_CG.top
    echo "Number of water molecules added: $watnum"
else
    echo "Error: solvate.log not found or does not contain expected output."
    exit 1
fi

#Generate minimization file
gmx grompp -p ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -c ${arg1}_${arg2}_${arg3}_${arg4}_CG_solvated.gro -f ${MartiniDir}minimization.mdp -o minimization.tpr >gromppMin.log 2>&1
echo "grompp minimization is done, check gromppMin.log file"
#Add ion to neutralize system
echo 13 | gmx genion -s minimization.tpr -o ${arg1}_${arg2}_${arg3}_${arg4}_CG_solvated_ions.gro -neutral -rmin 0.5 -p ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -nname CL- -pname NA+ >genion.log 2>&1
echo "Adding ion is done, check genion.log file"
#Generate minimization file again for neutralized system
gmx grompp -p ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -c ${arg1}_${arg2}_${arg3}_${arg4}_CG_solvated_ions.gro -f ${MartiniDir}minimization.mdp -o minimization.tpr -maxwarn 1 >gromppMinIon.log 2>&1
echo "grompp minization of neutralized is done, check gromppMinIon.log file"
#Run minimization
gmx mdrun -deffnm minimization -v >minimization.log 2>&1
echo "Minimization is done, check minimization.log file"
#Generate Equilibration tpr
gmx grompp -p ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -c minimization.gro -f ${MartiniDir}equilibration.mdp -o equilibration.tpr -maxwarn 2 >gromppEqu.log 2>&1
echo "grompp equalibration is done, check gromppEqu.log file"
#Run equilibration
gmx mdrun -deffnm equilibration -v >equilibration.log 2>&1
echo "Equalibration is done, check equilibration.log file"
#Generate first dynamics tpr file
gmx grompp -p ${arg1}_${arg2}_${arg3}_${arg4}_CG.top -c equilibration.gro -f ${MartiniDir_mdp_file} -o ${arg1}_${arg2}_${arg3}_${arg4}_CG_md.tpr -maxwarn 1 >gromppMD.log 2>&1
echo "grompp MD production is done, check gromppMD.log file"
#generate psf file for easy viewing in VMD
${Martini2psf} -b ${arg1}_${arg2}_${arg3}_${arg4}_CG
echo "MD production run is going, check ${arg1}_${arg2}_${arg3}_${arg4}_CG_md.log"
#Run the dynamic
gmx mdrun -nt $SLURM_CPUS_PER_TASK -deffnm ${arg1}_${arg2}_${arg3}_${arg4}_CG_md >${arg1}_${arg2}_${arg3}_${arg4}_CG_md.log 2>&1
