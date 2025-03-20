#!/bin/bash

# Define GROMACS versions and corresponding module names
declare -A gmx_versions=(
    ["gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"]="gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"
)

# Define GPU types and quantities
declare -A gpu_configs=(
    ["A100"]="1 2"
    ["A30"]="1 2 3 4"
)

# Base directory for benchmarks
base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/gpu/"

# Loop through each GROMACS version
for version_key in "${!gmx_versions[@]}"; do
    gmx_module="${gmx_versions[$version_key]}"
    
    # Create a directory for this GROMACS version
    version_dir="${base_dir}/${version_key}"
    mkdir -p "$version_dir"
    
    # Loop through each GPU configuration
    for gpu_type in "${!gpu_configs[@]}"; do
        # Create a directory for each GPU type under the version directory
        gpu_type_dir="${version_dir}/${gpu_type}"
        mkdir -p "$gpu_type_dir"
        
        for gpu_qty in ${gpu_configs[$gpu_type]}; do
            # Create a subdirectory for each GPU quantity under the GPU type directory
            gpu_qty_dir="${gpu_type_dir}/${gpu_qty}gpus"
            mkdir -p "$gpu_qty_dir"
            
            for replica in {1..3}; do
                # Define job name and output file
                job_name="${version_key}_${gpu_type}_${gpu_qty}gpus_replica${replica}"
                output_file="${gpu_qty_dir}/${job_name}.out"
                
                # Determine the appropriate GROMACS command
                if [[ "$gmx_module" == *"openmpi"* ]]; then
                    gmx_cmd="gmx_mpi"
                else
                    gmx_cmd="gmx"
                fi
                
                # Create the SLURM job script
                job_script="${gpu_qty_dir}/run_benchmark_${gpu_type}_${gpu_qty}gpus_replica${replica}.sh"
                cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --job-name="$job_name"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40  # Adjust if needed
#SBATCH --time=0:30:00
#SBATCH --output=$output_file
#SBATCH --gres=gpu:${gpu_type}:${gpu_qty} 

module load $gmx_module
cd $gpu_qty_dir
$gmx_cmd mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_${gpu_type}_${gpu_qty}gpus_replica${replica} -nb gpu -pme gpu -bonded gpu -update gpu
EOL
                
                # Make the job script executable
                chmod +x "$job_script"
                
                # Submit the job
                sbatch "$job_script"
            done
        done
    done
done