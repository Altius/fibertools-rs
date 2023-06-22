#!/bin/bash -fx
#SBATCH --partition=pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --get-user-env

  echo slurm node: $SLURMD_NODENAME , jobid: $SLURM_JOB_ID
  module add rust
  python ../target/debug/ft extract $1 -v -m 0 -a ./test_all.tsv -s
