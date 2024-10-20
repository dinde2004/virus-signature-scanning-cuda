#!/bin/bash

#SBATCH --job-name=my_job_h
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --gpus=h100-96:1
#SBATCH --constraint=xgpi
#SBATCH --time=01:00:00
#SBATCH --output=output_%N_%j_job_h.slurmlog
#SBATCH --error=error_%N_%j_job_h.slurmlog

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
    nsys profile -t cuda -o nsys/job_${i}.qdrep ./matcher tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.out
    nsys profile -t cuda -o nsys/benchmark_${i}.qdrep ./bench-h100 tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.ans
    # ncu --set full --clock-control none -o ncu/job_5_${i}.ncu-report ./matcher tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.out
    # ncu --set full --clock-control none -o ncu/benchmark_${i}.ncu-report ./bench-h100 tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.ans
done

echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"