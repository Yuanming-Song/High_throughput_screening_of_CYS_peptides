#!/bin/bash
#SBATCH --job-name="gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_16cores_replica1"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1/16cores/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_16cores_replica1.out

module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1
cd /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1/16cores
gmx_mpi mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/C_A_A_SIRAH_md.tpr -deffnm C_A_A_SIRAH_md_replica1
