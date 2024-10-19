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

# Run the implementation
echo -e "\n====> Running...\n"
for i in {0..9}
do
    ./matcher tests/samp_${i}.fastq tests/sig_${i}.fasta > tests/out_${i}.out
done

echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"