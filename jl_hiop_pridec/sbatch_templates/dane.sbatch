#!/bin/bash
#SBATCH --job-name=$CASE
#SBATCH --time=$TIME:00:00
#SBATCH --ntasks=$NTASKS
#SBATCH -A hiop
#SBATCH -p pbatch

$HIOP_SH
