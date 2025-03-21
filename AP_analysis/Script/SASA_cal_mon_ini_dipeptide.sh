#!/bin/bash

# Cindex variable (set to 1, 2, or 3 as needed)
#Cindex=3

# Directory containing the .xvg files
dir="out/"

# Output file name based on Cindex
output_file="SASA_result_with_common_mon_ini.txt"
>"$output_file" # Clear the output file

# Search for files matching the pattern #_#_#.xvg
for file in "$dir"*_*.xvg; do
    if [[ $file =~ ([A-Z])_([A-Z]).xvg ]]; then
        # Extract letters from the filename
        letter1=${BASH_REMATCH[1]}
        letter2=${BASH_REMATCH[2]}

        # Check Cindex condition
        # if [[ $Cindex -eq 1 && $letter1 != "C" ]]; then
        #     continue  # Skip files where letter1 is not "C"
        # elif [[ $Cindex -eq 2 && $letter2 != "C" ]]; then
        #     continue  # Skip files where letter2 is not "C"
        # elif [[ $Cindex -eq 3 && $letter3 != "C" ]]; then
        #     continue  # Skip files where letter3 is not "C"
        # fi

        # Construct corresponding _mon.xvg filename
        mon_file="${dir}${letter1}_${letter2}_mon.xvg"

        mon_file_ini="${dir}${letter1}_${letter2}_mon_ini.xvg"

        # Check if the corresponding _mon.xvg file exists
        if [[ -f "$mon_file_ini" ]]; then
            # Read first non-comment line from _mon.xvg
            mon_in=$(awk '!/^#|@/ {print $2; exit}' "$mon_file_ini")

            # Get the second column from the last line of _mon.xvg
            mon_fin=$(awk '!/^#|@/ {last=$2} END {print last}' "$mon_file")

            # Get the second column from the last line of the original .xvg file
            dim_fin=$(awk '!/^#|@/ {last=$2} END {print last}' "$file")
            if [[ -f "$file" ]]; then
                # Calculate the required ratios
                ratio_mon_fin=$(echo "scale=4; $mon_in/$mon_fin" | bc)
                ratio_dim_fin=$(echo "scale=4; $mon_in/$dim_fin" | bc)

                # Write to the output file
                echo "$letter1 $letter2  $ratio_mon_fin $ratio_dim_fin" >>"$output_file"
            fi
        fi
    fi
done

#echo "Results written to $output_file"
