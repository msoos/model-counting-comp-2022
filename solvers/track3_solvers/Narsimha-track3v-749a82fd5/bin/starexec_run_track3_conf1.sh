#!/bin/bash

file=$1
echo "c o This script is for projected  model counting"

solfile=$(mktemp)
indfile=$(mktemp)
cleanfile2=$(mktemp)
cleanfile=$(mktemp)
preprocessed_cnf_file=$(mktemp)
cache_size=8000
tout_ganak=1200
tout_be=210

grep "c p show" $file | sed -E "s/c p show (.*)/c ind \1 0/" > $indfile
grep -v "^c" $file > $cleanfile2
cat $cleanfile2 $indfile > $cleanfile

./doalarm ${tout_be} ./arjun --backbone 1 $cleanfile --elimtofile $preprocessed_cnf_file | sed "s/^/c o /"
found=`grep "^p cnf" $preprocessed_cnf_file`
if [[ $found == *"p cnf"* ]]; then
   echo "c o OK, Arjun succeeded"
   cp $preprocessed_cnf_file $cleanfile
   multi=`grep "^c MUST MUTIPLY BY" $preprocessed_cnf_file| sed "s/2\*\*//" | awk '{print $5}'`
else
   echo "c o WARNING Arjun did NOT succeed"
   multi=0
fi

./doalarm ${tout_ganak} ./ganak -cs ${cache_size} -t ${tout_ganak} $cleanfile > $solfile
solved_by_ganak=`grep "^s .*SATISFIABLE" $solfile`
if [[ $solved_by_ganak == *"SATISFIABLE"* ]]; then
    sed -E "s/^(.)/c o \1/" $solfile
    sat=`grep "^s .*SATISFIABLE" $solfile`
    count=`grep "^s .*mc" $solfile | awk '{print $3}'`

    export BC_LINE_LENGTH=99999000000
    if [[ "$count" == "0" ]]; then
        log_10_count="-inf"
    else
        count=`echo "$count*(2^$multi)" | bc -l`
        log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `
    fi

    echo $sat
    echo "c s type pmc"
    echo "c s log10-estimate $log_10_count"
    echo "c s exact arb int $count"
    exit 0
else
    ./approxmc --epsilon 0.9 $cleanfile | tee $solfile | sed "s/^/c o /"
    solved_by_approxmc=`grep "^s .*SATISFIABLE" $solfile`
    if [[ $solved_by_approxmc == *"SATISFIABLE"* ]]; then
        sat=`grep "^s .*SATISFIABLE" $solfile`
        count=`grep "^s .*mc" $solfile | awk '{print $3}'`
        export BC_LINE_LENGTH=99999000000
        if [[ "$count" == "0" ]]; then
            log_10_count="-inf"
        else
            count=`echo "$count*(2^$multi)" | bc -l`
            log_10_count=`echo "scale=15; l($count)/l(10)" | bc -l `
        fi
    
        echo $sat
        echo "c s type pmc"
        echo "c s log10-estimate $log_10_count"
        echo "c s approx arb int $count"
        exit 0
    else
        echo "c o ApproxMC did NOT work"
        exit -1
    fi
fi
