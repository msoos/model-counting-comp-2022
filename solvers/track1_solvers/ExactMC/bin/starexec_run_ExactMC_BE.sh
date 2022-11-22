#!/bin/bash

echo "c o This script is for Track 1 in model counting competition 2022"

file=$1
clean_file=$(mktemp XXXXXX.cnf)
preprocessed_file=$(mktemp XXXXXX.cnf)

mc=$(grep "^c t " $file)
echo "c o found header: $mc"
echo "c o file: ${file} clean_file: ${clean_file} preprocessed_file: ${preprocessed_file}"

timeout_total=${STAREXEC_WALLCLOCK_LIMIT}

time_begin=$(date +%s)
timeout_be=$((timeout_total/12))
grep -v "^c" ${file} > ${clean_file}
echo "c o Running B+E with timeout: ${timeout_be} seconds"
./doalarm ${timeout_be} ./B+E_linux -cpu-lim=${timeout_be} ${clean_file} > ${preprocessed_file}
preprocessing_status=$(grep "^p cnf" ${preprocessed_file})
if [[ ${preprocessing_status} == *"p cnf"* ]]; then
   echo "c o B+E succeeded"
   grep -v "^c" ${preprocessed_file} > ${clean_file}
else
   echo "c o WARNING! B+E did NOT succeed"
fi
time_end=$(date +%s)

timeout_mc=$((timeout_total + time_begin - time_end))
total_mem_gb=$(echo "scale=3; ${STAREXEC_MAX_MEM}/1024*0.75" | bc)
echo "c o Running ExactMC, timeleft: ${timeout_mc} seconds, memo: ${total_mem_gb} GB"
./KCBox ExactMC --competition --memo ${total_mem_gb} --quiet ${clean_file}

