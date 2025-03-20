#!/bin/bash
#SBATCH --job-name="gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_12cores_replica2"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0:30:00
#SBATCH --output=/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/A100//gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1/12cores/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_12cores_replica2.out
#SBATCH --gres=gpu:A100:1 

module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1
cd /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/A100//gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1/12cores
gmx_mpi mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Dimer/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica2 -nb gpu -pme gpu -bonded gpu -update gpu
