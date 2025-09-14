#!/bin/bash

export CONTINGENCY_FILE=$CONT_FILE
export MAX_ITER=$MAX_ITER
export JULIA_SRC_FILE=./julia_src/hiop.jl

mkdir -p output
export OUTPUT_DIR=output/$CASE_$(date +%Y%m%d_%H%M%S)

mkdir -p "${OUTPUT_DIR}"

export OUTPUT_LOGS=$OUTPUT_DIR/rank_logs
mkdir -p "${OUTPUT_LOGS}"

export OUTPUT_ITER=$OUTPUT_DIR/iterations
mkdir -p "${OUTPUT_ITER}"


{
  echo "CONTINGENCY_FILE=${CONTINGENCY_FILE}"
  echo "MAX_ITER=${MAX_ITER}"
  echo "JULIA_SRC_FILE=${JULIA_SRC_FILE}"
  echo "PATH_TO_EXAJUGO=${PATH_TO_EXAJUGO}"
  echo "PATH_TO_INSTANCES=${PATH_TO_INSTANCES}"
  echo "PATH_TO_HSLLIB=${PATH_TO_HSLLIB}"
  echo "HIOP_INSTALL_DIR=${HIOP_INSTALL_DIR}"
  echo "PATH_TO_TSSLOPE=${PATH_TO_TSSLOPE}"
  echo "NORMALIZE_X=${NORMALIZE_X}"
  echo "GRAD_MULTIPLIER=${GRAD_MULTIPLIER}"
  echo "BASE_CASE_OPTIONS=${BASE_CASE_OPTIONS}"
  echo "CONTINGENCY_CASE_OPTIONS=${CONTINGENCY_CASE_OPTIONS}"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
} > "${OUTPUT_DIR}/environment_vars.txt"

# Copy submitted batch file to output directory

cp "${BATCH_FILE}" "${OUTPUT_DIR}"
options_file=hiop_pridec.options

if [ -f "${options_file}" ]; then
    cp "${options_file}" "${OUTPUT_DIR}"
fi

if [ -n "${BASE_CASE_OPTIONS:-}" ]; then
    options_file="$BASE_CASE_OPTIONS"
    if [ -f "$options_file" ]; then
        cp "$options_file" "$OUTPUT_DIR"
    fi
fi

if [ -n "${CONTINGENCY_CASE_OPTIONS:-}" ]; then
    options_file="$CONTINGENCY_CASE_OPTIONS"
    if [ -f "$options_file" ]; then
        cp "$options_file" "$OUTPUT_DIR"
    fi
fi

# Run the executable using srun

srun --output=$OUTPUT_LOGS/rank_%t/log_$CASE_%j_%t.out \
     --error=$OUTPUT_LOGS/rank_%t/log_$CASE_%j_%t.err \
     bash -c 'export OUTPUT_DIR_RANK="'$OUTPUT_LOGS'/rank_${SLURM_PROCID}"/; ./build/jl_NlpPriDec.exe $CASE'
