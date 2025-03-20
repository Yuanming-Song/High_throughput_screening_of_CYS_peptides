#!/bin/bash
#SBATCH --job-name="gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_100nlist_replica2"
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=dtobias_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:10:00
#SBATCH --gres=gpu:A30:1
#SBATCH --output=/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Monomer/nstlist/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1/100/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1_100nlist_replica2.out

module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1
cd /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Monomer/nstlist/gcc.11.2.0-cuda.11.7.1.openmpi.5.0.1//100
gmx_mpi mdrun -s /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/SIRAH/Setup_test/Benchmark/Monomer/nstlist/C_A_A_SIRAH_md_100.tpr -deffnm C_A_A_SIRAH_md_replica2
