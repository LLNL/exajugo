#!/bin/bash
#SBATCH --job-name=$CASE
#SBATCH --time=$TIME:00:00
#SBATCH --ntasks=$NTASKS

export CONTINGENCY_FILE=$CONT_FILE
export MAX_ITER=$MAX_ITER
export JULIA_SRC_FILE=./julia_src/hiop.jl

mkdir -p output
export OUTPUT_DIR=output/$CASE_$(date +%Y%m%d_%H%M%S)

# Run the executable using srun

srun --output=$OUTPUT_DIR/rank_%t/log_$CASE_%j_%t.out \
     --error=$OUTPUT_DIR/rank_%t/log_$CASE_%j_%t.err \
     bash -c 'export OUTPUT_DIR_RANK="'$OUTPUT_DIR'/rank_${SLURM_PROCID}"/; ./build/jl_NlpPriDec.exe $CASE'