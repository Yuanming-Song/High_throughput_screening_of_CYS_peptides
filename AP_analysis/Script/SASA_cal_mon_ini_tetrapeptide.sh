#!/bin/bash

# Directory containing the .xvg files
dir="out/"
for Cindex in  3; do
    # Output file name based on Cindex
    output_file="SASA_result_with_common_mon_ini_tetra_C${Cindex}.txt"
    >"$output_file" # Clear the output file

    # Search for files matching the pattern #_#_#_#.xvg
    for file in "$dir"*_*_*_*.xvg; do
        if [[ $file =~ ([A-Z])_([A-Z])_([A-Z])_([A-Z]).xvg ]]; then
            # Extract letters from the filename
            letter1=${BASH_REMATCH[1]}
            letter2=${BASH_REMATCH[2]}
            letter3=${BASH_REMATCH[3]}
            letter4=${BASH_REMATCH[4]}

            # Check Cindex condition
            if [[ $Cindex -eq 1 && $letter1 != "C" ]]; then
                continue # Skip files where letter1 is not "C"
            elif [[ $Cindex -eq 2 && $letter2 != "C" ]]; then
                continue # Skip files where letter2 is not "C"
            elif [[ $Cindex -eq 3 && $letter3 != "C" ]]; then
                continue # Skip files where letter3 is not "C"
            elif [[ $Cindex -eq 4 && $letter4 != "C" ]]; then
                continue # Skip files where letter4 is not "C"
            fi

            # Construct corresponding _mon.xvg filename
            mon_file="${dir}${letter1}_${letter2}_${letter3}_${letter4}_mon_ini.xvg"
            # Check if the corresponding _mon.xvg file exists
            if [[ -f "$mon_file" ]]; then
                # Read first non-comment line from _mon.xvg
                mon_in=$(awk '!/^#|@/ {print $2; exit}' "$mon_file")

                mon_file="${dir}${letter1}_${letter2}_${letter3}_${letter4}_mon.xvg"
                if [[ -f "$mon_file" ]]; then
                    # Get the second column from the last line of _mon.xvg
                    mon_fin=$(awk '!/^#|@/ {last=$2} END {print last}' "$mon_file")
                fi
                # Get the second column from the last line of the original .xvg file
                dim_fin=$(awk '!/^#|@/ {last=$2} END {print last}' "$file")

                # Calculate the required ratios
                ratio_mon_fin=$(echo "scale=4; $mon_in/$mon_fin" | bc)
                ratio_dim_fin=$(echo "scale=4; $mon_in/$dim_fin" | bc)

                # Write to the output file
                echo "$letter1 $letter2 $letter3 $letter4 $ratio_mon_fin $ratio_dim_fin" >>"$output_file"
            fi
        fi
    done
done
echo "Results written to $output_file" 
