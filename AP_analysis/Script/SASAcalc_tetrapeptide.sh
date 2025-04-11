#!/bin/bash
# Set main directory for all files

#load file
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <C_position> <first_non_C_residue>"
    echo "Example: $0 1 A"
    exit 1
fi

# Get command line arguments
Cindex=$1
first_res=$2

# Define the list of natural amino acids excluding C
amino_acids=("A" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "C")

# Clean up previous files
#rm -f */*dim*.??? log/*?_?.???

# Set output files
outfile="SubData/SASA_result_Tetrapeptide_C${Cindex}_${first_res}.dat"
missing_vdw_log="missing_vdw_dim_C${Cindex}_${first_res}.log"
rm -f $outfile

maindir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI22/Tetrapeptide_dis_C${Cindex}/"

# Loop through each amino acid for amino1
for amino1 in "${amino_acids[@]}"; do
    # Skip if C position is 1 and amino1 is not C
    if [[ "$Cindex" -eq 1 && "$amino1" != "C" ]]; then
        continue
    fi
    # Skip if C position is 4 and amino1 is not first_res
    if [[ "$Cindex" -eq 4 && "$amino1" != "$first_res" ]]; then
        continue
    fi
    # Skip if C position is not 1 but amino1 is C
    if [[ "$Cindex" -ne 1 && "$amino1" == "C" ]]; then
        continue
    fi

    # Loop through each amino acid for amino2
    for amino2 in "${amino_acids[@]}"; do
        # Skip if C position is 2 and amino2 is not C
        if [[ "$Cindex" -eq 2 && "$amino2" != "C" ]]; then
            continue
        fi
        # Skip if C position is 1 and amino2 is not first_res
        if [[ "$Cindex" -eq 1 && "$amino2" != "$first_res" ]]; then
            continue
        fi
        # Skip if C position is not 2 but amino2 is C
        if [[ "$Cindex" -ne 2 && "$amino2" == "C" ]]; then
            continue
        fi

        # Loop through each amino acid for amino3
        for amino3 in "${amino_acids[@]}"; do
            # Skip if C position is 3 and amino3 is not C
            if [[ "$Cindex" -eq 3 && "$amino3" != "C" ]]; then
                continue
            fi
            # Skip if C position is 2 and amino3 is not first_res
            if [[ "$Cindex" -eq 2 && "$amino3" != "$first_res" ]]; then
                continue
            fi
            # Skip if C position is not 3 but amino3 is C
            if [[ "$Cindex" -ne 3 && "$amino3" == "C" ]]; then
                continue
            fi

            # Loop through each amino acid for amino4
            for amino4 in "${amino_acids[@]}"; do
                # Skip if C position is 4 and amino4 is not C
                if [[ "$Cindex" -eq 4 && "$amino4" != "C" ]]; then
                    continue
                fi
                # Skip if C position is 3 and amino4 is not first_res
                if [[ "$Cindex" -eq 3 && "$amino4" != "$first_res" ]]; then
                    continue
                fi
                # Skip if C position is not 4 but amino4 is C
                if [[ "$Cindex" -ne 4 && "$amino4" == "C" ]]; then
                    continue
                fi

                # Set the working directory
                workdir="${maindir}/${amino1}_${amino2}_${amino3}_${amino4}"

                # Check if the directory exists
                if [ ! -d "$workdir" ]; then
                    echo "no dir ${amino1}_${amino2}_${amino3}_${amino4}" 
                    continue
                fi 

                # Check if the .xtc file exists
                xtc_file="${workdir}/${amino1}_${amino2}_${amino3}_${amino4}_CG_md.xtc"
                if [ ! -f "$xtc_file" ]; then
                    echo "job not started ${amino1}_${amino2}_${amino3}_${amino4}"
                    continue
                fi

                # Check if the .log file exists
                log_file="${workdir}/${amino1}_${amino2}_${amino3}_${amino4}_CG_md.log"
                if [ ! -f "$log_file" ]; then
                    echo "no log file for ${amino1}_${amino2}_${amino3}_${amino4}"
                    continue
                fi

                # Check if the log file contains "Finished mdrun" in the last 10 lines
                if ! tail -n 10 "$log_file" | grep -q "^Finished mdrun"; then
                    echo "job not finished ${amino1}_${amino2}_${amino3}_${amino4}"
                    continue
                fi

                # Run gmx sasa
                gro_file="${workdir}/${amino1}_${amino2}_${amino3}_${amino4}_CG_solvated_ions.gro"
                output_xvg="out/${amino1}_${amino2}_${amino3}_${amino4}.xvg"
                echo 1 | gmx sasa -f "$xtc_file" -s "$gro_file" -o "$output_xvg" >log/SASA_${amino1}_${amino2}_${amino3}_${amino4}.log 2>&1
                output_xvgi="out/${amino1}_${amino2}_${amino3}_${amino4}_dim_ini.xvg"
                echo 1 | gmx sasa -f "$gro_file" -s "$gro_file" -o "$output_xvgi" >log/SASA_${amino1}_${amino2}_${amino3}_${amino4}_dim_ini.log 2>&1
                if grep -q "WARNING: could not" log/SASA_${amino1}_${amino2}_${amino3}_${amino4}_dim_ini.log; then
                    echo "Warning in ${amino1}_${amino2}_${amino3}_${amino4}: Skipping this loop." >>$missing_vdw_log
                    continue
                fi
                # Extract the first row, second column, and the last row, second column, and calculate the ratio
                first_value=$(awk '/^[^@#]/ {print $2; exit}' "$output_xvgi")
                last_value=$(awk 'END {print $2}' "$output_xvg")

                if [ -n "$first_value" ] && [ -n "$last_value" ]; then
                    ratio=$(echo "scale=2; $first_value / $last_value" | bc)
                    echo "${amino1} ${amino2} ${amino3} ${amino4} $ratio" >>$outfile
                fi
            done
        done
    done
done
