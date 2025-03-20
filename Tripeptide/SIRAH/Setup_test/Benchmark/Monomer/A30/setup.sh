#!/bin/bash

# Define GROMACS versions and corresponding module names
declare -A gmx_versions=(
    #["gcc.11.2.0"]="gromacs/2024.2/gcc.11.2.0"
    #["gcc.11.2.0-cuda.11.7.1"]="gromacs/2024.2/gcc.11.2.0-cuda.11.7.1"
    ["gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"]="gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1"
    #["gcc.11.2.0-openmpi.5.0.1"]="gromacs/2024.2/gcc.11.2.0-openmpi.5.0.1"
)

# Define core counts for benchmarking
core_counts=(1 4 8 16 32 40 64)

# Base directory for benchmarks
base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/A30/"

# Loop through each GROMACS version
for version_key in "${!gmx_versions[@]}"; do
    gmx_module="${gmx_versions[$version_key]}"
    
    # Create a directory for this GROMACS version
    version_dir="${base_dir}/${version_key}"
    mkdir -p "$version_dir"
    
    # Loop through each core count
    for core_num in "${core_counts[@]}"; do
        # Create a subdirectory for this core count
        bench_dir="${version_dir}/${core_num}cores"
        mkdir -p "$bench_dir"
        
                for replica in {1..3}; do
            # Define job name and output file
            job_name="${version_key}_${core_num}cores_replica${replica}"
            output_file="${bench_dir}/${job_name}.out"
            
            # Determine the appropriate GROMACS command
            if [[ "$gmx_module" == *"openmpi"* ]]; then
                gmx_cmd="gmx_mpi"
            else
                gmx_cmd="gmx"
            fi
            
            # Create the SLURM job script
            job_script="${bench_dir}/run_benchmark_replica${replica}.sh"
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
cd $bench_dir
$gmx_cmd mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica${replica} -nb gpu -pme gpu -bonded gpu -update gpu
EOL
            
            # Make the job script executable
            chmod +x "$job_script"
            
            # Submit the job
            sbatch "$job_script"
        done
    done
done