#!/bin/bash

# Define the main directory
maindir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI22"

# Define the list of amino acids (excluding Cys for non-C positions)
amino_acids=("A" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "E" "D" "C")

# Generate all possible tetrapeptide sequences
for aa1 in "${amino_acids[@]}"; do
    for aa2 in "${amino_acids[@]}"; do
        for aa3 in "${amino_acids[@]}"; do
            for aa4 in "${amino_acids[@]}"; do
                sequence="${aa1}${aa2}${aa3}${aa4}"
                jobname="${aa1}_${aa2}_${aa3}_${aa4}"
                # Count occurrences of cysteine
                c_count=$(grep -o "C" <<<"$sequence" | wc -l)

                # Skip if there is not exactly **one** cysteine
                if [[ "$c_count" -ne 1 ]]; then
                    continue
                fi

                # Determine the cysteine index (1-based)
                cindex=0
                for ((i = 0; i < ${#sequence}; i++)); do
                    if [[ "${sequence:$i:1}" == "C" ]]; then
                        cindex=$((i + 1))
                        break
                    fi
                done

                # Define directory path
                cindex_dir="$maindir/Tetrapeptide_dis_C${cindex}"
                peptide_dir="$cindex_dir/${jobname}"

                # Check if peptide folder exists
                if [[ ! -d "$peptide_dir" ]]; then
                    echo "No ${jobname}"
                else
                    shopt -s nullglob # Ensure globbing doesn't return literal patterns if no match
                    files=("$peptide_dir"/*md*gro)

                    if [ ${#files[@]} -gt 0 ]; then
                        continue
                    else
                        echo "No *md*gro file found for ${jobname}"
                        # delete the dir
                        rm -rf "$peptide_dir"
                    fi
                fi
                # Go to the directory associated with cindex
                cd "$cindex_dir" || continue
                #echo "now submitting job"
                #continue
                # Resubmit job
                sbatch --job-name="${jobname}" \
                    --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
                    --account=dtobias_lab \
                    --partition=standard \
                    --nodes=1 \
                    --ntasks-per-node=1 \
                    --cpus-per-task=40 \
                    --time=24:00:00 \
                    --out=out/${jobname}.out \
                    --wrap="bash //dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/Cys_stapled_tetrapeptide_MARTINI_run.sh ${aa1} ${aa2} ${aa3} ${aa4}"
                # Return to the original directory
                cd "$maindir" || exit

            done
        done
    done
done
