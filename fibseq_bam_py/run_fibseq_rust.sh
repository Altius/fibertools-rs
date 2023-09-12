#!/bin/bash -fx
#SBATCH --partition=pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --get-user-env

  echo slurm node: $SLURMD_NODENAME , jobid: $SLURM_JOB_ID
  module load fiberseq-rs
  python fibseq_rust.py $1 $2
