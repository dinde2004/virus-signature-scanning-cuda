#!/bin/bash

#SBATCH --job-name=benchmark_job_a
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --gpus=a100-40:1
#SBATCH --constraint=xgph
#SBATCH --time=01:00:00
#SBATCH --output=output_%N_%j_benchmark_a.slurmlog
#SBATCH --error=error_%N_%j.slurmlog

echo "Job is running on $(hostname), started at $(date)"

# Get some output about GPU status
nvidia-smi 

# Run the benchmark
echo -e "\n====> Running...\n"
for i in {0..9}
do
    ./bench-a100 tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.ans
done

echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"
