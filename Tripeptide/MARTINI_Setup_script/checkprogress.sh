#there are tripeptide and dipeptide simulation
#each simulation has a dimer and monomer version
#tripeptide has 3 sub dir corresponding to C being at different position
#each of these sub dirs has a pair of dimer and monomer subdir
#dipeptide is similar but only has 2 sub dir 
# Here are the names of these subdir 
#Dipeptide_dis_C1              Tripeptide_dis_C3 Dipeptide_dis_C2   Tripeptide_mon_C1 Dipeptide_mon_C1  Tripeptide_dis_C1   Tripeptide_mon_C2 Dipeptide_mon_C2  Tripeptide_dis_C2   Tripeptide_mon_C3

#dis means disulfide stapled dimer, mon means monomer
#under each subdir, there are peptide sequence, they always contain one C, like C_A_A under Tripeptide_mon_C1 and Tripeptide_dis_C1
#in each dir, there is a *md*log file, if it doesn't exit, if it doesn't contain "Writing checkpoint, step 2000000 " , then the job is not complete or failed, if it has error in it, then it failed otherwise incomplete.

#i want to loop through all the dipeptide and tripeptide that contains one and only one C, and check whether which sequence failed or incomplete

#!/bin/bash

# Define log file
log_file="simulation_status.log"
echo "Simulation Status Log" > "$log_file"
maindir="/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/"
# Define directories
dirs=("Dipeptide_dis_C1" "Dipeptide_dis_C2" "Dipeptide_mon_C1" "Dipeptide_mon_C2" \
      "Tripeptide_dis_C1" "Tripeptide_dis_C2" "Tripeptide_dis_C3" \
      "Tripeptide_mon_C1" "Tripeptide_mon_C2" "Tripeptide_mon_C3")

# Define possible amino acids excluding C
amino_acids=("A" "G" "L" "F" "M" "N" "P" "Q" "R" "S" "T" "V" "Y" "W" "H" "I" "K")

# Loop through each directory
for dir in "${dirs[@]}"; do
  dir=${maindir}/$dir
  # Determine sequence type based on directory name
  if [[ "$dir" == *Tripeptide_dis_C1* || "$dir" == *Tripeptide_mon_C1* ]]; then
    for a1 in "${amino_acids[@]}"; do
      for a2 in "${amino_acids[@]}"; do
        sequence="C_${a1}_${a2}"
        seq_dir="$dir/$sequence"

        if [[ -d "$seq_dir" ]]; then
          md_log_file=$(find "$seq_dir" -name "${sequence}*md*log" -print -quit)

          if [[ -z "$md_log_file" ]]; then
            echo "$sequence in $seq_dir- Incomplete (No md log file found)" >> "$log_file"
          elif grep -q "Error" "$md_log_file"; then
            echo "$sequence in $seq_dir- Failed (Error found)" >> "$log_file"
          elif ! grep -q "Writing checkpoint, step 2000000" "$md_log_file"; then
            echo "$sequence in $seq_dir- Incomplete (Did not reach final checkpoint)" >> "$log_file"
          
          fi
        else
          echo "$sequence in $seq_dir- Directory does not exist" >> "$log_file"
        fi
      done
    done

  elif [[ "$dir" == *Tripeptide_dis_C2* || "$dir" == *Tripeptide_mon_C2* ]]; then
    for a1 in "${amino_acids[@]}"; do
      for a2 in "${amino_acids[@]}"; do
        sequence="${a1}_C_${a2}"
        seq_dir="$dir/$sequence"

        if [[ -d "$seq_dir" ]]; then
          md_log_file=$(find "$seq_dir" -name "${sequence}*md*log" -print -quit)

          if [[ -z "$md_log_file" ]]; then
            echo "$sequence in $seq_dir- Incomplete (No md log file found)" >> "$log_file"
          elif grep -q "Error" "$md_log_file"; then
            echo "$sequence in $seq_dir- Failed (Error found)" >> "$log_file"
          elif ! grep -q "Writing checkpoint, step 2000000" "$md_log_file"; then
            echo "$sequence in $seq_dir- Incomplete (Did not reach final checkpoint)" >> "$log_file"
          fi
        else
          echo "$sequence in $seq_dir- Directory does not exist" >> "$log_file"
        fi
      done
    done

  elif [[ "$dir" == *Tripeptide_dis_C3* || "$dir" == *Tripeptide_mon_C3* ]]; then
    for a1 in "${amino_acids[@]}"; do
      for a2 in "${amino_acids[@]}"; do
        sequence="${a1}_${a2}_C"
        seq_dir="$dir/$sequence"

        if [[ -d "$seq_dir" ]]; then
          md_log_file=$(find "$seq_dir" -name "${sequence}*md*log" -print -quit)

          if [[ -z "$md_log_file" ]]; then
            echo "$sequence in $seq_dir- Incomplete (No md log file found)" >> "$log_file"
          elif grep -q "Error" "$md_log_file"; then
            echo "$sequence in $seq_dir- Failed (Error found)" >> "$log_file"
          elif ! grep -q "Writing checkpoint, step 2000000" "$md_log_file"; then
            echo "$sequence in $seq_dir- Incomplete (Did not reach final checkpoint)" >> "$log_file"
          fi
        else
          echo "$sequence in $seq_dir- Directory does not exist" >> "$log_file"
        fi
      done
    done

  elif [[ "$dir" == *Dipeptide_dis_C1* || "$dir" == *Dipeptide_mon_C1* ]]; then
    for a1 in "${amino_acids[@]}"; do
      sequence="C_${a1}"
      seq_dir="$dir/$sequence"

      if [[ -d "$seq_dir" ]]; then
        md_log_file=$(find "$seq_dir" -name "${sequence}*md*log" -print -quit)

        if [[ -z "$md_log_file" ]]; then
          echo "$sequence in $seq_dir- Incomplete (No md log file found)" >> "$log_file"
        elif grep -q "Error" "$md_log_file"; then
          echo "$sequence in $seq_dir- Failed (Error found)" >> "$log_file"
        elif ! grep -q "Writing checkpoint, step 2000000" "$md_log_file"; then
          echo "$sequence in $seq_dir- Incomplete (Did not reach final checkpoint)" >> "$log_file"
        fi
      else
        echo "$sequence in $seq_dir- Directory does not exist" >> "$log_file"
      fi
    done

  elif [[ "$dir" == *Dipeptide_dis_C2* || "$dir" == *Dipeptide_mon_C2* ]]; then
    for a1 in "${amino_acids[@]}"; do
      sequence="${a1}_C"
      seq_dir="$dir/$sequence"

      if [[ -d "$seq_dir" ]]; then
        md_log_file=$(find "$seq_dir" -name "${sequence}*md*log" -print -quit)

        if [[ -z "$md_log_file" ]]; then
          echo "$sequence in $seq_dir- Incomplete (No md log file found)" >> "$log_file"
        elif grep -q "Error" "$md_log_file"; then
          echo "$sequence in $seq_dir- Failed (Error found)" >> "$log_file"
        elif ! grep -q "Writing checkpoint, step 2000000" "$md_log_file"; then
          echo "$sequence in $seq_dir- Incomplete (Did not reach final checkpoint)" >> "$log_file"
        fi
      else
        echo "$sequence in $seq_dir- Directory does not exist" >> "$log_file"
      fi
    done
  fi
done