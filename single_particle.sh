#!/bin/bash
#full trj data file
TRJ=$1
#max number of particles to parse
NMAX=$2
for (( i=1; i <= $NMAX; i++ ))
do
    head -1 ${TRJ} > ${TRJ}_${i}
    awk -v pid=${i} '{if($1 == pid) print $0}' ${TRJ} >> ${TRJ}_${i}
done