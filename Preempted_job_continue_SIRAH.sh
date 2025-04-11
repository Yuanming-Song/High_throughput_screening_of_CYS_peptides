#!/bin/bash

# Directories to check
dirs=$(ls -d */ | sed 's|/$||')

# Files without logs or without "Finished mdrun" status
no_log_file=()
no_finished_mdrun=()

for dir in $dirs; do
  log_file=$(ls "$dir"/${dir}_SIRAH_md*log 2>/dev/null | head -n 1)

  if [[ -z $log_file ]]; then
    no_log_file+=("$dir")
  else
    if ! tail -n 10 "$log_file" | grep -q "Finished mdrun"; then
      no_finished_mdrun+=("$dir")
    fi
  fi
done

echo "Directories without log files:"
printf "%s\n" "${no_log_file[@]}"

echo "Directories with log files missing 'Finished mdrun':"
printf "%s\n" "${no_finished_mdrun[@]}"

# Directories with missing "Finished mdrun" status

# Base directory for output logs
mkdir -p out
sys="mon"
for subdir in "${no_finished_mdrun[@]}"; do
  #if [[ "$subdir" == "N_W_C" ]]; then
  #  continue
  #fi

  # Check if the checkpoint file exists
  if [[ -f "${subdir}/${subdir}_SIRAH_md.cpt" ]]; then
    # Submit job to continue from the last checkpoint
    sbatch --job-name="${subdir}_${sys}_conti" \
      --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
      --account=dtobias_lab_gpu \
      --partition=free-gpu \
      --nodes=1 \
      --ntasks-per-node=1 \
      --cpus-per-task=8 \
      --time=6:00:00 \
      --output="out/${subdir}.out" \
      --gres=gpu:A30:1 \
      --wrap="module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1; cd ${subdir}; gmx_mpi mdrun -deffnm ${subdir}_SIRAH_md -cpi ${subdir}_SIRAH_md.cpt -append > MD_conti.log 2>&1"
  else
    # Submit job to start a new run
    sbatch --job-name="${subdir}_${sys}_conti" \
      --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
      --account=dtobias_lab_gpu \
      --partition=free-gpu \
      --nodes=1 \
      --ntasks-per-node=1 \
      --cpus-per-task=8 \
      --time=6:00:00 \
      --output="out/${subdir}.out" \
      --gres=gpu:A30:1 \
      --wrap="module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1; cd ${subdir}; gmx_mpi mdrun -deffnm ${subdir}_SIRAH_md > MD.log 2>&1"
  fi
done

for subdir in "${no_log_file[@]}"; do
  if [[ "$subdir" == "out" ]]; then
    continue
  fi
  read -p "Are you sure you want to delete the directory '$subdir'? (y/n): " confirm
  if [[ "$confirm" == "y" || "$confirm" == "Y" ]]; then
    rm -r "$subdir" &
    echo "Directory '$subdir' has been deleted."
  else
    echo "Deletion of '$subdir' canceled."
  fi
done
