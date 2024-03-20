#!/bin/bash


wd=$(pwd)
dataPath=${wd}/../data/PHYLOGENIES

threads=7

# prepare list of samples

input=${dataPath}/03_TRIMAL/MAFFT

cd ${input}

if [ -f "samples.txt" ]; then
    rm samples.txt
    ls *.fa >> samples.txt
else
    ls *.fa >> samples.txt
fi

for i in $(<samples.txt);do

  echo ${i}

  input=${dataPath}/03_TRIMAL/MAFFT
  output=${dataPath}/02_MAFFT/FILTERED

  mkdir -p ${output}

  /usr/local/bin/mafft --thread ${threads} --localpair ${input}/${i} > ${output}/${i}

  input=${output}
  output=${dataPath}/03_TRIMAL/MAFFT/FILTERED

  mkdir -p ${output}

  /usr/local/bin/trimal\
   -in ${input}/${i}\
   -out ${output}/${i}\
   -automated1

  input=${output}
  output=${dataPath}/04_IQTREE/MAFFT

  mkdir -p ${output}

  /usr/local/bin/iqtree -s ${input}/${i} -alrt 1000 -bb 1000 -nt AUTO -pre ${output}/${i}

done

exit
