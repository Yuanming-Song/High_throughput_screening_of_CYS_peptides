#!/bin/bash
Cindex=2
#load file
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1
# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "C")
#in case of debugging
#amino_acids=("A" "C")
rm */*mon*
for Cindex in 1 2; do
    dt=1000
    outfile="SASA_result_SIRAH_C${Cindex}_mon.dat"
    rm $outfile
    # Set main directory for all files
    #maindir="/share/crsp/lab//dtobias/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide_dis_C${Cindex}/" # Replace with your actual main directory path
    maindir="/share/crsp/lab//dtobias/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Dipeptide_mon_C${Cindex}/" # Replace with your actual main directory path
    # Loop through each amino acid for amino1
    for amino1 in "${amino_acids[@]}"; do
        # Loop through each amino acid for amino2
        for amino2 in "${amino_acids[@]}"; do
            # Set the working directory
            workdir="${maindir}/${amino1}_${amino2}"

            if [[ "$Cindex" -eq 1 && "$amino1" != "C" ]]; then
                continue
            elif [[ "$Cindex" -eq 2 && "$amino2" != "C" ]]; then
                continue
            fi
            if [[ "$Cindex" -eq 1 && ("$amino2" == "C") ]]; then
                continue
            elif [[ "$Cindex" -eq 2 && ("$amino1" == "C") ]]; then
                continue
            elif [[ "$Cindex" -eq 3 && ("$amino1" == "C" || "$amino2" == "C") ]]; then
                continue
            elif [[ "$amino1" == "C" && "$amino2" == "C" ]]; then
                continue
            fi
            # Check if the directory exists
            if [ ! -d "$workdir" ]; then
                echo "nodir $workdir"
                continue
            fi

            # Check if the .xtc file exists
            xtc_file="${workdir}/${amino1}_${amino2}_SIRAH_md.xtc"
            if [ ! -f "$xtc_file" ]; then
                echo "job not started ${amino1}_${amino2}"
                continue
            fi

            # Check if the .log file exists
            log_file="${workdir}/${amino1}_${amino2}_SIRAH_md.log"
            if [ ! -f "$log_file" ]; then
                echo "no log file $log_file"
                continue
            fi

            # Check if the log file contains "Finished mdrun" in the last 10 lines
            if ! tail -n 10 "$log_file" | grep -q "^Finished mdrun"; then
                echo "job not finished ${amino1}_${amino2}"
                continue
            fi

            # Run gmx sasa
            gro_file="${workdir}/${amino1}_${amino2}_CG_solvated_ions.gro"
            output_xvg="out/${amino1}_${amino2}_mon.xvg"
            echo 1 | gmx sasa -f "$xtc_file" -s "$gro_file" -o "$output_xvg" -dt $dt >log/SASA_${amino1}_${amino2}_mon.log 2>&1
            output_xvgi="out/${amino1}_${amino2}_mon_ini.xvg"
            echo 1 | gmx sasa -f "$gro_file" -s "$gro_file" -o "$output_xvgi" >log/SASA_${amino1}_${amino2}_mon_ini.log 2>&1

            # Extract the first row, second column, and the last row, second column, and calculate the ratio
            first_value=$(awk '/^[^@#]/ {print $2; exit}' "$output_xvgi")
            last_value=$(awk 'END {print $2}' "$output_xvg")

            if [ -n "$first_value" ] && [ -n "$last_value" ]; then
                ratio=$(echo "scale=2; $first_value / $last_value" | bc)
                echo "${amino1} ${amino2} $ratio" >>$outfile
            fi
        done
    done
done
