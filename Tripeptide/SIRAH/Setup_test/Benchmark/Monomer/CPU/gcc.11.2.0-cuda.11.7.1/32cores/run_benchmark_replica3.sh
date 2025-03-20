#!/bin/bash
#SBATCH --job-name="gcc.11.2.0-cuda.11.7.1_32cores_replica3"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --output=/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/gcc.11.2.0-cuda.11.7.1/32cores/gcc.11.2.0-cuda.11.7.1_32cores_replica3.out

module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1
cd /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/gcc.11.2.0-cuda.11.7.1/32cores
gmx mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica3
