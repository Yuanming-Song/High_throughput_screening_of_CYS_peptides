#!/bin/bash

#location for mdp files 
mdpdir=/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/mdp_files_longer/

mdpdir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/mdp_files/"
#radius for insert molecule
radius=0.5
# box size
inibox_size=12.5
box_size_nm=13
# how many copies of peptide to put in
ncopy=300
#SIRAH CG script
SIRAH="/dfs9/tw/yuanmis1/mrsec/ff/sirah.ff/tools/CGCONV/cgconv.pl"
#ff file location
SIRAHff="/dfs9/tw/yuanmis1/mrsec/ff/sirah.ff/"
#Load gromacs
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1

# Check if exactly three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <arg1> <arg2> <arg3> for tripeptide"
    exit 1
fi

# Assign arguments to variables
arg1=$1
arg2=$2
arg3=$3
echo "peptide ${arg1}_${arg2}_${arg3}"
# Define the directory name
dir_name="${arg1}_${arg2}_${arg3}"

# Check if the directory exists
if [ -d "$dir_name" ]; then
    echo "Directory '$dir_name' already exists, submit manually."
    #rm ${arg1}_${arg2}_${arg3}/* 
    #exit 1  # Stop the script with a non-zero exit status
else
    # Create the directory if it does not exist
    mkdir "$dir_name"
    echo "Directory '$dir_name' created successfully."
fi
cd ${arg1}_${arg2}_${arg3}

# Check each argument and update IndexC with the 1-based index
if [ "$arg1" == "C" ]; then
    IndexC=1
elif [ "$arg2" == "C" ]; then
    IndexC=2
elif [ "$arg3" == "C" ]; then
    IndexC=3
fi
# Location for exisiting peptide pdb
MartiniDir="/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide_mon_C${IndexC}/${arg1}_${arg2}_${arg3}/"
#remove files just in case
rm *trr* *edr* *ndx* *log* *gro* *tpr* *xtc* *itp* *top* *cpt*
#building short cuts
#cp ${MartiniDir}/${arg1}_${arg2}_${arg3}.pdb .
/dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none -e /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SetupScript/mutate_template_mon.tcl -args ${arg1} ${arg2} ${arg3} AA_template_dimer_C$IndexC > mutate_${arg1}_${arg2}_${arg3}.log 2>&1
# Check if any of the arguments is H
if [[ "$arg1" == "H" || "$arg2" == "H" || "$arg3" == "H" ]]; then
    # Define the pdb file name based on the arguments
    pdb_file="${arg1}_${arg2}_${arg3}.pdb"
    
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

# Input and output file names
input_file="${arg1}_${arg2}_${arg3}.pdb"
output_file="${arg1}_${arg2}_${arg3}_filtered.pdb"

# Filter lines and save to the output file
awk '$5 != "B"' "$input_file" > "$output_file"

#get CG rep
$SIRAH -i $output_file -o ${arg1}_${arg2}_${arg3}_CG.pdb

#get top and gro files
echo 1 | gmx pdb2gmx -f ${arg1}_${arg2}_${arg3}_CG.pdb -o ${arg1}_${arg2}_${arg3}_CG.gro -ff 'sirah' -p ${arg1}_${arg2}_${arg3}_SIRAH.top -merge all > ${arg1}_${arg2}_${arg3}_pdb2gmx.log 2>&1 

##select backbone restraint
echo -e "a GN GO\n\nq" | gmx make_ndx -f ${arg1}_${arg2}_${arg3}_CG.gro -o ${arg1}_${arg2}_${arg3}_CG.ndx > make_ndx_single.log 2>&1

#make restrain itp file 
echo -e "GN_GO" | gmx genrestr -f ${arg1}_${arg2}_${arg3}_CG.gro -n ${arg1}_${arg2}_${arg3}_CG.ndx -o bkbres.itp > genrestr.log 2>&1
echo -e "GN_GO" | gmx genrestr -f ${arg1}_${arg2}_${arg3}_CG.gro -n ${arg1}_${arg2}_${arg3}_CG.ndx -o bkbres_soft.itp -fc 100 100 100 > genrestr_soft.log 2>&1

#include these itp files
sed -i '/; Include water topology/i\
; Backbone restraints\n#ifdef GN_GO\n#include "bkbres.itp"\n#endif\n; Backbone soft restrains\n#ifdef GN_GO_SOFT\n#include "bkbres_soft.itp"\n#endif
' "${arg1}_${arg2}_${arg3}_SIRAH.top" 

#insert molecules
gmx insert-molecules -box ${inibox_size} ${inibox_size} ${inibox_size} -nmol ${ncopy} -ci ${arg1}_${arg2}_${arg3}_CG.gro -radius ${radius} -o ${arg1}_${arg2}_${arg3}_CG_Lattice.gro -rot xyz > ${arg1}_${arg2}_${arg3}_CG_multi_copy_setup.log 2>&1

#change the number of peptide in the top file manually
sed -i -e "s/\(Protein_chain_A\s\+\)1/\1${ncopy}/" "${arg1}_${arg2}_${arg3}_SIRAH.top"

#since box size for insert peptides and solvation are different
center_value=$(echo "${box_size_nm}/2" | bc -l)
center_value=$(echo "($box_size_nm - $inibox_size) / 2" | bc -l)

#shift the center to the center of solvation
gmx editconf -f ${arg1}_${arg2}_${arg3}_CG_Lattice.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -bt cubic -o ${arg1}_${arg2}_${arg3}_CG_Lattice.gro -d 1 -translate $center_value $center_value $center_value > editconf.log 2>&1

#solvate molecule
gmx solvate -cp ${arg1}_${arg2}_${arg3}_CG_Lattice.gro -cs ${SIRAHff}wt416.gro -o ${arg1}_${arg2}_${arg3}_CG_solvated.gro -box ${box_size_nm} ${box_size_nm} ${box_size_nm} -radius 0.1 -p ${arg1}_${arg2}_${arg3}_SIRAH.top > solvate.log 2>&1
echo "Solvating solute in water is done, check solvate.log file"

#make ndx file
echo q | gmx make_ndx -f ${arg1}_${arg2}_${arg3}_CG_solvated.gro -o ${arg1}_${arg2}_${arg3}_CG_solvated.ndx > make_ndx.log 2>&1

#generate the tpr file to remove water
gmx grompp -f ${mdpdir}em1_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po delete1.mdp -c ${arg1}_${arg2}_${arg3}_CG_solvated.gro -r ${arg1}_${arg2}_${arg3}_CG_solvated.gro -o ${arg1}_${arg2}_${arg3}_CG_solvated_em1.tpr -maxwarn 2 > grompp_em1.log 2>&1

#make the selection
gmx select -f ${arg1}_${arg2}_${arg3}_CG_solvated.gro -s ${arg1}_${arg2}_${arg3}_CG_solvated_em1.tpr -n ${arg1}_${arg2}_${arg3}_CG_solvated.ndx -on rm_close_wt4.ndx -select 'not (same residue as (resname WT4 and within 0.3 of group Protein))'  > select_wt.log 2>&1

#keep only the selection
gmx editconf -f ${arg1}_${arg2}_${arg3}_CG_solvated.gro -o ${arg1}_${arg2}_${arg3}_CG_solvated2.gro -n rm_close_wt4.ndx  > rm_wat.log 2>&1 

# Count WP1 in both files
count1=$(grep -c WP1 "${arg1}_${arg2}_${arg3}_CG_solvated.gro")
count2=$(grep -c WP1 "${arg1}_${arg2}_${arg3}_CG_solvated2.gro")

# Check if WT4 line has only one occurrence and replace the count
if grep -q "^WT4" "${arg1}_${arg2}_${arg3}_SIRAH.top"; then
    sed -i "s/^\(WT4\s\+\)$count1/\1$count2/" "${arg1}_${arg2}_${arg3}_SIRAH.top"
else
    echo "Error: WT4 line not found or multiple occurrences."
    #exit 1
fi

# Generate the tpr file
gmx grompp -f ${mdpdir}em1_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po delete2.mdp -c ${arg1}_${arg2}_${arg3}_CG_solvated2.gro -r ${arg1}_${arg2}_${arg3}_CG_solvated2.gro -o em1.tpr -maxwarn 2 > grompp_em1_2.log 2>&1

#Add ion to neutralize system
echo -e "resname WT4\n" | gmx genion -s em1.tpr -o ${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro -neutral -rmin 0.5 -p ${arg1}_${arg2}_${arg3}_SIRAH.top -nname ClW -pname NaW > genion.log 2>&1
echo "Adding ion is done, check genion.log file"

##make backbone restraint ndx file
echo -e "a GN GO\n\nq" | gmx make_ndx -f ${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro -o ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx > make_ndx2.log 2>&1

count1=$(grep -c WP1 "${arg1}_${arg2}_${arg3}_CG_solvated2.gro")
count2=$(grep -c WP1 "${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro")

# Check if WT4 line has only one occurrence and replace the count
if grep -q "^WT4" "${arg1}_${arg2}_${arg3}_SIRAH.top"; then
    sed -i "s/^\(WT4\s\+\)$count1/\1$count2/" "${arg1}_${arg2}_${arg3}_SIRAH.top"
else
    echo "Error: WT4 line not found or multiple occurrences."
    #exit 1
fi

##genrate psf
./sirah.ff/tools/g_top2psf.pl -i ${arg1}_${arg2}_${arg3}_SIRAH.top -o ${arg1}_${arg2}_${arg3}_SIRAH.psf

#Generate minimization file again for neutralized system
gmx grompp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -c ${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro -f ${mdpdir}em1_CGPROT.mdp -o ${arg1}_${arg2}_${arg3}_SIRAH_em1.tpr -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx -r ${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro -po em1.mdp -maxwarn 1 > grompp_em1_CGPROT.log 2>&1

#Run minimization
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_em1 -v > em1_CGPROT.log 2>&1
echo "Side chain minimization is done, check em1_CGPROT.log file"

#grompp tpr file for em2
gmx grompp -f ${mdpdir}em2_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po em2.mdp -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx -c ${arg1}_${arg2}_${arg3}_SIRAH_em1.gro -o ${arg1}_${arg2}_${arg3}_SIRAH_em2.tpr > grompp_em2_CGPROT.log 2>&1

#Run minimization
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_em2 -v >  EM2.log 2>&1
echo "All atom minimization done, check em2_CGPROT.log file"

#grompp tpr file for solvent condensation, otherwise the system will collapse due to lack of WT4
gmx grompp -f ${mdpdir}NPT_solvent_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po NPT_solvent.mdp -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx -c ${arg1}_${arg2}_${arg3}_SIRAH_em2.gro -o ${arg1}_${arg2}_${arg3}_SIRAH_NPT_solvent.tpr -r ${arg1}_${arg2}_${arg3}_CG_solvated_ions.gro > grompp_NPT_solvent.log 2>&1
#Run npt equalib
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_NPT_solvent -v >  NPT_solvent.log 2>&1
echo "NPT equal done, check NPT_solvent.log file"

#grompp tpr file
gmx grompp -f ${mdpdir}eq1_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po eq1.mdp -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx  -o ${arg1}_${arg2}_${arg3}_SIRAH_eq1.tpr -c ${arg1}_${arg2}_${arg3}_SIRAH_NPT_solvent.gro -r ${arg1}_${arg2}_${arg3}_SIRAH_NPT_solvent.gro > grompp_eq1_CGPROT.log 2>&1
#Run NVT equal
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_eq1 -v >  EQ1.log 2>&1
echo "eq1 done"

#grompp tpr file
gmx grompp -f ${mdpdir}eq2_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po eq2.mdp -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx -c ${arg1}_${arg2}_${arg3}_SIRAH_eq1.gro -r ${arg1}_${arg2}_${arg3}_SIRAH_eq1.gro -o ${arg1}_${arg2}_${arg3}_SIRAH_eq2.tpr > grompp_eq2_CGPROT.log 2>&1

#more NVT equal
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_eq2 -v >  EQ2.log 2>&1
echo "eq2 done"

#grompp tpr file
gmx grompp -f ${mdpdir}md_CGPROT.mdp -p ${arg1}_${arg2}_${arg3}_SIRAH.top -po md.mdp -n ${arg1}_${arg2}_${arg3}_CG_solvated2.ndx -c ${arg1}_${arg2}_${arg3}_SIRAH_eq2.gro -o ${arg1}_${arg2}_${arg3}_SIRAH_md.tpr > grompp_md_CGPROT.log 2>&1

#production run
gmx mdrun -deffnm ${arg1}_${arg2}_${arg3}_SIRAH_md -v >  MD.log 2>&1 