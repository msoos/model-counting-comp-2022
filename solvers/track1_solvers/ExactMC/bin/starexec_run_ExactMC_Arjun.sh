#!/bin/bash

echo "c o This script is for Track 1 in model counting competition 2022"

file=$1
clean_file=$(mktemp XXXXXX.cnf)
preprocessed_file=$(mktemp XXXXXX.cnf)
solving_file=$(mktemp XXXXXX.out)

mc=$(grep "^c t " $file)
echo "c o found header: $mc"

echo "c o file: ${file} clean_file: ${clean_file} preprocessed_file: ${preprocessed_file} solving_file: ${solving_file}"

timeout_total=${STAREXEC_WALLCLOCK_LIMIT}

time_begin=$(date +%s)
timeout_be=$((timeout_total/12))
grep -v "^c" ${file} > ${clean_file}
echo "c o Running Arjun with timeout: ${timeout_be}"
timeout ${timeout_be} ./arjun --recomp 1 --backbone 1 ${file} --elimtofile ${preprocessed_file} >/dev/null 2>&1
preprocessing_status=$(grep "^p cnf" ${preprocessed_file})
if [[ ${preprocessing_status} == *"p cnf"* ]]; then
   echo "c o OK, Arjun succeeded"
   grep -v "^c" ${preprocessed_file} > ${clean_file}
   multi=`grep "^c MUST MUTIPLY BY" ${preprocessed_file} | sed "s/2\*\*//" | awk '{print $5}'`
else
   echo "c o WARNING Arjun did NOT succeed"
   multi=0
fi
echo "c o MULTI will be 2**${multi}"
time_end=$(date +%s)

timeout_mc=$((timeout_total + time_begin - time_end))
total_mem_gb=$(echo "scale=3; ${STAREXEC_MAX_MEM}/1024*0.75" | bc)
echo "c o Running ExactMC, timeleft: ${timeout_mc} seconds, memo: ${total_mem_gb} GB"
./KCBox ExactMC --memo ${total_mem_gb} --quiet ${clean_file} > ${solving_file}
count=$(grep "c o Number of models:" ${solving_file} | awk '{print $6}')

export BC_LINE_LENGTH=1000000000
tuned_count=$(echo "${count}*(2^${multi})" | bc -l)
if [[ ${tuned_count} == "0" ]]; then log_10_count="-inf"
else log_10_count=$(echo "scale=15; l($tuned_count)/l(10)" | bc -l)
fi
cat ${solving_file}
echo "c s type mc"
echo "c s log10-estimate ${log_10_count}"
echo "c s exact arb int ${tuned_count}"


