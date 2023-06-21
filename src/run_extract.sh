#!/bin/bash -fx
#SBATCH --partition=pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --get-user-env

  echo slurm node: $SLURMD_NODENAME , jobid: $SLURM_JOB_ID
  module add rust
  python ft extract $1 -v -m 1 --m6a ./test_m6a.tsv --cpg ./test_cpg.tsv -a ./test_all.tsv -f -q -s
