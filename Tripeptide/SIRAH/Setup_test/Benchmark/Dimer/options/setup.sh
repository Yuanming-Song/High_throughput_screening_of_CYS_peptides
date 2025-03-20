#!/bin/bash

# Define GROMACS versions and corresponding module names
declare -A gmx_versions=(
    ["gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"]="gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"
)

# Define core count for benchmarking
core_num=8

# Base directory for benchmarks
base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/options/"

# List of options to benchmark
options=("-nb gpu" "-pme gpu" "-bonded gpu" "-update gpu")
option_names=("nb" "pme" "bonded" "update")

# Loop through each GROMACS version
for version_key in "${!gmx_versions[@]}"; do
    gmx_module="${gmx_versions[$version_key]}"
    
    # Create a directory for this GROMACS version
    version_dir="${base_dir}/${version_key}"
    mkdir -p "$version_dir"

    # Determine the appropriate GROMACS command
            if [[ "$gmx_module" == *"openmpi"* ]]; then
                gmx_cmd="gmx_mpi"
            else
                gmx_cmd="gmx"
            fi

    # Benchmark for each option by removing it
    for i in "${!options[@]}"; do
        current_option="${options[$i]}"
        option_name="${option_names[$i]}"
        
        # Create a subdirectory for this option
        option_dir="${version_dir}/${option_name}"
        mkdir -p "$option_dir"
        
        for replica in {1..3}; do
            # Define job name and output file
            job_name="${version_key}_${option_name}_replica${replica}"
            output_file="${option_dir}/${job_name}.out"
            
            # Build the options list without the current option
            remaining_options=$(printf "%s " "${options[@]}" | sed "s|${current_option}||")
            
            # Create the SLURM job script
            job_script="${option_dir}/run_benchmark_replica${replica}.sh"
            cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --job-name="$job_name"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$core_num
#SBATCH --time=0:10:00
#SBATCH --output=$output_file
#SBATCH --gres=gpu:A30:1

module load $gmx_module
cd $option_dir
$gmx_cmd mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica${replica} $remaining_options
EOL
            
            # Make the job script executable
            chmod +x "$job_script"
            
            # Submit the job
            sbatch "$job_script"
        done
    done
    
    # Run a blank option benchmark (no GPU options)
    blank_dir="${version_dir}/blank"
    mkdir -p "$blank_dir"
    
    for replica in {1..3}; do
        # Define job name and output file
        job_name="${version_key}_blank_replica${replica}"
        output_file="${blank_dir}/${job_name}.out"
        
        # Create the SLURM job script
        job_script="${blank_dir}/run_benchmark_replica${replica}.sh"
        cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --job-name="$job_name"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$core_num
#SBATCH --time=0:10:00
#SBATCH --output=$output_file
#SBATCH --gres=gpu:A30:1

module load $gmx_module
cd $blank_dir
$gmx_cmd mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica${replica}
EOL
        
        # Make the job script executable
        chmod +x "$job_script"
        
        # Submit the job
        sbatch "$job_script"
    done
done