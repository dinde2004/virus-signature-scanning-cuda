#!/bin/bash

#SBATCH --job-name=gen_job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --gpus=h100-96:1
#SBATCH --constraint=xgpi
#SBATCH --time=01:00:00
#SBATCH --output=output_%N_%j_gen_job.slurmlog
#SBATCH --error=error_%N_%j_gen_job.slurmlog

echo "Job is running on $(hostname), started at $(date)"

# Get some output about GPU status
nvidia-smi 

# Compile
make clean
make

# Run the test generator
echo -e "\n====> Running...\n"
# No wildcards tests
for i in {0..4}
do
    ./gen_sig 1000 3000 10000 0.0 > tests/sig_${i}.fasta
    ./gen_sample tests/sig_${i}.fasta 2000 20 1 2 100000 200000 10 30 0.0 > tests/samp_${i}.fastq
done
# With wildcards tests
for i in {5..9}
do
    ./gen_sig 1000 3000 10000 0.1 > tests/sig_${i}.fasta
    ./gen_sample tests/sig_${i}.fasta 2000 20 1 2 100000 200000 10 30 0.1 > tests/samp_${i}.fastq
done

echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"