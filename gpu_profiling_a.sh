#!/bin/bash

#SBATCH --job-name=profiling_a
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --gpus=a100-40:1
#SBATCH --constraint=xgph
#SBATCH --time=01:00:00
#SBATCH --output=output_%N_%j_profiling_h.slurmlog
#SBATCH --error=error_%N_%j_profiling_h.slurmlog

echo "Job is running on $(hostname), started at $(date)"

# Get some output about GPU status
nvidia-smi 

# Compile
make clean
make

# Make directory for profiling results
mkdir -p nsys
mkdir -p ncu

# Run the implementation
echo -e "\n====> Running...\n"
for i in 10 20 25 30
do
    # nsys profile -o nsys/kernel_${i} ./matcher tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.out
    # nsys profile -o nsys/benchmark_${i} ./bench-h100 tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.ans
    ncu --set full --clock-control none -o ncu/kernel_${i} ./matcher tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.out
    ncu --set full --clock-control none -o ncu/benchmark_${i} ./bench-a100 tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.ans
done

echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"