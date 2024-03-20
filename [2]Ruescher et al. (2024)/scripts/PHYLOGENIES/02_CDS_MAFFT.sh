#!/bin/bash

wd=$(pwd)
dataPath=${wd}/../data/PHYLOGENIES

threads=7
# prepare list of sample and extension

input=${dataPath}/01_CDS/01_SEQUENCES

cd ${input}

if [ -f "samples.txt" ]; then
    rm samples.txt
    ls *.fa >> samples.txt
else
    ls *.fa >> samples.txt
fi

for i in $(<samples.txt);do

echo ${i}

input=${dataPath}/01_CDS/01_SEQUENCES
output=${dataPath}/02_MAFFT

mkdir -p ${output}

/usr/local/bin/mafft --thread ${threads} --auto ${input}/${i} > ${output}/${i}

input=${output}
output=${dataPath}/03_TRIMAL/MAFFT

mkdir -p ${output}

/usr/local/bin/trimal\
 -in ${input}/${i}\
 -out ${output}/${i}\
 -automated1

done

exit
