#!/bin/bash
# Define base directory variable
BASEDIR="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide"

#location for mdp files 
mdpdir="$BASEDIR/Tripeptide/SIRAH/Setup_test/mdp_files_longer/"
mdpdir="$BASEDIR/Tripeptide/SIRAH_Setup_script/mdp_files/"
#radius for insert molecule
radius=0.5
# box size
inibox_size=11.5
box_size_nm=11
# how many copies of peptide to put in
ncopy=300
#SIRAH CG script
SIRAH="$BASEDIR/../ff/sirah.ff/tools/CGCONV/cgconv.pl"
#ff file location
SIRAHff="$BASEDIR/../ff/sirah.ff/"
#Load gromacs
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1

# Check if exactly two or three arguments are provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <arg1> <arg2> [box_size_nm]"
    exit 1
fi

# Assign arguments to variables
arg1=$1
arg2=$2
if [ "$#" -eq 3 ]; then
    box_size_nm=$3
    inibox_size=$3
else
    box_size_nm=11
    inibox_size=11.5
fi
echo "${arg1}_${arg2}"
# Define the directory name
dir_name="${arg1}_${arg2}"

# Check if the directory exists
if [ -d "$dir_name" ]; then
    echo "Directory '$dir_name' already exists, submit manually."
    rm $dir_name/*
    #mv $dir_name Mon_wt_CYX/
    mkdir "$dir_name"
    #exit 1  # Stop the script with a non-zero exit status
else
    # Create the directory if it does not exist
    mkdir "$dir_name"
    echo "Directory '$dir_name' created successfully."
fi
cd ${arg1}_${arg2}

# Check each argument and update IndexC with the 1-based index
if [ "$arg1" == "C" ]; then
    IndexC=1
elif [ "$arg2" == "C" ]; then
    IndexC=2
fi

# Use a linear peptide template, mutate residue based on input
/dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none -e $BASEDIR/Dipeptide/AA_template_script/mutate_template_mon_dipeptide.tcl -args ${arg1} ${arg2}  AA_template_dimer_C$IndexC > mutate_${arg1}_${arg2}.log 2>&1
# Check if any of the arguments is H
if [[ "$arg1" == "H" || "$arg2" == "H"  ]]; then
    # Define the pdb file name based on the arguments
    pdb_file="${arg1}_${arg2}.pdb"
    
    # Check if the file exists
    if [[ -f "$pdb_file" ]]; then
        # Replace HSD with HIS in the file
        sed -i 's/HSD/HIS/g' "$pdb_file"
        echo "Replaced HSD with HIS in $pdb_file"
    else
        echo "File $pdb_file does not exist."
    fi
fi
ln -s $SIRAHff sirah.ff
ln -s sirah.ff/residuetypes.dat
ln -s sirah.ff/specbond.dat
ln -s sirah.ff/vdwradii.dat 
#replace CYS with CYX
#sed -i  "s/CYS/CYX/g" "${arg1}_${arg2}.pdb"

$SIRAH -i ${arg1}_${arg2}.pdb -o ${arg1}_${arg2}_CG.pdb
echo 1 | gmx_mpi pdb2gmx -f ${arg1}_${arg2}_CG.pdb -o ${arg1}_${arg2}_CG.gro -ff 'sirah' -p ${arg1}_${arg2}_SIRAH_single.top -merge all >${arg1}_${arg2}_pdb2gmx.log 2>&1 

##select backbone restraint
echo -e "a GN GO\n\nq" | gmx_mpi make_ndx -f ${arg1}_${arg2}_CG.gro -o ${arg1}_${arg2}_CG.ndx >make_ndx_single.log 2>&1
#make restrain file 
echo -e "GN_GO" | gmx_mpi genrestr -f ${arg1}_${arg2}_CG.gro -n ${arg1}_${arg2}_CG.ndx -o bkbres.itp >genrestr.log 2>&1
echo -e "GN_GO" | gmx_mpi genrestr -f ${arg1}_${arg2}_CG.gro -n ${arg1}_${arg2}_CG.ndx -o bkbres_soft.itp -fc 100 100 100 >genrestr_soft.log 2>&1



#Find the last ATOM line in ${arg1}_${arg2}.pdb and get the peplength
peplength=$(awk '/^ATOM/ {if ($2 > max) max = $2} END {print max}' "${arg1}_${arg2}_CG.pdb")
# Process ${arg1}_${arg2}_CG.pdb, inserting TER lines as specified
output_file=${arg1}_${arg2}_CG_redo.pdb
> "$output_file"  # Clear or create the output file

# Loop through each line of ${arg1}_${arg2}_CG.pdb
awk -v peplength="$peplength" -v output_file="$output_file" '
  BEGIN { rescounter = 1 }
  /^ATOM/ {
    print $0 >> output_file  # Print the ATOM line to the output file
    target_index = rescounter * peplength 
    if ($2 == target_index) {
      print "TER" >> output_file  # Insert TER only after the matching ATOM line
      rescounter++  # Increment rescounter only when a match is found
    }
    next
  }
  { print $0 >> output_file }  # Print any non-ATOM lines to the output file
' ${arg1}_${arg2}_CG.pdb


# generate new top with disulfide bond
echo 1 | gmx_mpi pdb2gmx -f ${arg1}_${arg2}_CG_redo.pdb -o ${arg1}_${arg2}_CG_redo.gro -ff 'sirah' -p ${arg1}_${arg2}_SIRAH.top -merge all >${arg1}_${arg2}_Lattice_pdb2gmx.log 2>&1 

#include these itp files
sed -i '/; Include water topology/i\
; Backbone restraints\n#ifdef GN_GO\n#include "bkbres.itp"\n#endif\n; Backbone soft restrains\n#ifdef GN_GO_SOFT\n#include "bkbres_soft.itp"\n#endif
' "${arg1}_${arg2}_SIRAH.top" 
#insert molecules
gmx_mpi insert-molecules -box ${inibox_size} ${inibox_size} ${inibox_size} -nmol ${ncopy} -ci ${arg1}_${arg2}_CG_redo.gro -radius ${radius} -o ${arg1}_${arg2}_CG_Lattice.gro -rot xyz >${arg1}_${arg2}_CG_multi_copy_setup.log 2>&1

#change the number of peptide in the top file manually
sed -i -e "s/\(Protein_chain_A\s\+\)1/\1${ncopy}/" "${arg1}_${arg2}_SIRAH.top"

#since box size for insert peptides and solvation are different
center_value=$(echo "${box_size_nm}/2" | bc -l)
center_value=$(echo "($box_size_nm - $inibox_size) / 2" | bc -l)

#shift the center to the center of solvation
gmx_mpi editconf -f ${arg1}_${arg2}_CG_Lattice.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -bt cubic -o ${arg1}_${arg2}_CG_Lattice.gro -d 1 -translate $center_value $center_value $center_value > editconf.log 2>&1

#solvate molecule
gmx_mpi solvate -cp ${arg1}_${arg2}_CG_Lattice.gro -cs ${SIRAHff}wt416.gro -o ${arg1}_${arg2}_CG_solvated.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -radius 0.1 -p ${arg1}_${arg2}_SIRAH.top > solvate.log 2>&1
echo "Solvating solute in water is done, check solvate.log file"

#make ndx file
echo q | gmx_mpi make_ndx -f ${arg1}_${arg2}_CG_solvated.gro -o ${arg1}_${arg2}_CG_solvated.ndx > make_ndx.log 2>&1

#generate the tpr file to remove water
gmx_mpi grompp -f ${mdpdir}em1_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po delete1.mdp -c ${arg1}_${arg2}_CG_solvated.gro -r ${arg1}_${arg2}_CG_solvated.gro -o ${arg1}_${arg2}_CG_solvated_em1.tpr -maxwarn 2 > grompp_em1.log 2>&1

#make the selection
gmx_mpi select -f ${arg1}_${arg2}_CG_solvated.gro -s ${arg1}_${arg2}_CG_solvated_em1.tpr -n ${arg1}_${arg2}_CG_solvated.ndx -on rm_close_wt4.ndx -select 'not (same residue as (resname WT4 and within 0.3 of group Protein))'  > select_wt.log 2>&1

#keep only the selection
gmx_mpi editconf -f ${arg1}_${arg2}_CG_solvated.gro -o ${arg1}_${arg2}_CG_solvated2.gro -n rm_close_wt4.ndx  > rm_wat.log 2>&1 

# Count WP1 in both files
count1=$(grep -c WP1 "${arg1}_${arg2}_CG_solvated.gro")
count2=$(grep -c WP1 "${arg1}_${arg2}_CG_solvated2.gro")

# Check if WT4 line has only one occurrence and replace the count
if grep -q "^WT4" "${arg1}_${arg2}_SIRAH.top"; then
    sed -i "s/^\(WT4\s\+\)$count1/\1$count2/" "${arg1}_${arg2}_SIRAH.top"
else
    echo "Error: WT4 line not found or multiple occurrences."
    #exit 1
fi

# Generate the tpr file
gmx_mpi grompp -f ${mdpdir}em1_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po delete2.mdp -c ${arg1}_${arg2}_CG_solvated2.gro -r ${arg1}_${arg2}_CG_solvated2.gro -o em1.tpr -maxwarn 2 > grompp_em1_2.log 2>&1

#Add ion to neutralize system
echo -e "resname WT4\n" | gmx_mpi genion -s em1.tpr -o ${arg1}_${arg2}_CG_solvated_ions.gro -rmin 0.5 -p ${arg1}_${arg2}_SIRAH.top -nname ClW -pname NaW -neutral -np 0 > genion.log 2>&1
echo "Adding ion is done, check genion.log file"

##make backbone restraint ndx file
echo -e "a GN GO\n\nq" | gmx_mpi make_ndx -f ${arg1}_${arg2}_CG_solvated_ions.gro -o ${arg1}_${arg2}_CG_solvated2.ndx > make_ndx2.log 2>&1

count1=$(grep -c WP1 "${arg1}_${arg2}_CG_solvated2.gro")
count2=$(grep -c WP1 "${arg1}_${arg2}_CG_solvated_ions.gro")

# Check if WT4 line has only one occurrence and replace the count
if grep -q "^WT4" "${arg1}_${arg2}_SIRAH.top"; then
    sed -i "s/^\(WT4\s\+\)$count1/\1$count2/" "${arg1}_${arg2}_SIRAH.top"
else
    echo "Error: WT4 line not found or multiple occurrences."
    #exit 1
fi

##genrate psf
./sirah.ff/tools/g_top2psf.pl -i ${arg1}_${arg2}_SIRAH.top -o ${arg1}_${arg2}_SIRAH.psf

#Generate minimization file again for neutralized system
gmx_mpi grompp -p ${arg1}_${arg2}_SIRAH.top -c ${arg1}_${arg2}_CG_solvated_ions.gro -f ${mdpdir}em1_CGPROT.mdp -o ${arg1}_${arg2}_SIRAH_em1.tpr -n ${arg1}_${arg2}_CG_solvated2.ndx -r ${arg1}_${arg2}_CG_solvated_ions.gro -po em1.mdp -maxwarn 1 > grompp_em1_CGPROT.log 2>&1

#Run minimization
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_em1 -v > em1_CGPROT.log 2>&1
echo "Side chain minimization is done, check em1_CGPROT.log file"

#grompp tpr file for em2
gmx_mpi grompp -f ${mdpdir}em2_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po em2.mdp -n ${arg1}_${arg2}_CG_solvated2.ndx -c ${arg1}_${arg2}_SIRAH_em1.gro -o ${arg1}_${arg2}_SIRAH_em2.tpr > grompp_em2_CGPROT.log 2>&1

#Run minimization
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_em2 -v >  EM2.log 2>&1
echo "All atom minimization done, check em2_CGPROT.log file"

#grompp tpr file for solvent condensation, otherwise the system will collapse due to lack of WT4
gmx_mpi grompp -f ${mdpdir}NPT_solvent_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po NPT_solvent.mdp -n ${arg1}_${arg2}_CG_solvated2.ndx -c ${arg1}_${arg2}_SIRAH_em2.gro -o ${arg1}_${arg2}_SIRAH_NPT_solvent.tpr -r ${arg1}_${arg2}_CG_solvated_ions.gro > grompp_NPT_solvent.log 2>&1
#Run npt equalib
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_NPT_solvent -v >  NPT_solvent.log 2>&1
echo "NPT equal done, check NPT_solvent.log file"

#grompp tpr file
gmx_mpi grompp -f ${mdpdir}eq1_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po eq1.mdp -n ${arg1}_${arg2}_CG_solvated2.ndx  -o ${arg1}_${arg2}_SIRAH_eq1.tpr -c ${arg1}_${arg2}_SIRAH_NPT_solvent.gro -r ${arg1}_${arg2}_SIRAH_NPT_solvent.gro > grompp_eq1_CGPROT.log 2>&1
#Run NVT equal
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_eq1  >  EQ1.log 2>&1
echo "eq1 done"

#grompp tpr file
gmx_mpi grompp -f ${mdpdir}eq2_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po eq2.mdp -n ${arg1}_${arg2}_CG_solvated2.ndx -c ${arg1}_${arg2}_SIRAH_eq1.gro -r ${arg1}_${arg2}_SIRAH_eq1.gro -o ${arg1}_${arg2}_SIRAH_eq2.tpr > grompp_eq2_CGPROT.log 2>&1

#more NVT equal
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_eq2  >  EQ2.log 2>&1
echo "eq2 done"

#grompp tpr file
gmx_mpi grompp -f ${mdpdir}md_CGPROT.mdp -p ${arg1}_${arg2}_SIRAH.top -po md.mdp -n ${arg1}_${arg2}_CG_solvated2.ndx -c ${arg1}_${arg2}_SIRAH_eq2.gro -o ${arg1}_${arg2}_SIRAH_md.tpr > grompp_md_CGPROT.log 2>&1

#production run
gmx_mpi mdrun -deffnm ${arg1}_${arg2}_SIRAH_md  >  MD.log 2>&1