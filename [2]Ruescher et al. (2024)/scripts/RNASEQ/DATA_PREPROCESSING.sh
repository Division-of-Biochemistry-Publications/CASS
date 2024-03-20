#!/bin/bash

WD=$(pwd) #working directory
DD=${WD}/01_RAW_DATA # raw data folder
MD=${WD}/00_META_DATA # meta data folder
GD=${MD}/Mesculenta # genome folder
threads=18 # change according to your work station

echo start preprocessing

# prepare list of sample and extension

cd ${DD}

rm samples.txt

extension=$(ls *_1*| head -1)
extension="${extension#*.}"

for i in *_1*;do
echo "${i%_1*}" >> samples.txt
done

cd ${WD}

# QC using FASTQC

echo Raw data FASTQC

input=${DD}
output=./02_FASTQC

mkdir -p ${output}/MULTIQC

fastqc --outdir ${output} -t ${threads} ${input}/*.${extension}
multiqc --outdir ${output}/MULTIQC ${output}

# trimming using BBDUK

input=${DD}
output=./03_TRIMMED_DATA

mkdir -p ${output}

for i in $(<${DD}/samples.txt);do

echo $i adapter and quality trimming

echo trimming ${i}

bbduk\
 k=23\
 mink=11\
 ktrim=r\
 hdist=1\
 ref=${MD}/adapters.fa\
 qtrim=rl\
 trimq=30\
 minlen=35\
 in1=${input}/${i}_1.${extension}\
 in2=${input}/${i}_2.${extension}\
 out1=${output}/${i}_1.${extension}\
 out2=${output}/${i}_2.${extension}\
 t=${threads}

done

# QC using FASTQC

echo Trimmed data FASTQC

input=${output}
output=./04_TRIMMED_DATA_FASTQC

mkdir -p ${output}/MULTIQC

fastqc --outdir ${output} -t ${threads} ${input}/*.${extension}
multiqc --outdir ${output}/MULTIQC ${output}

# Mapping

input=./03_TRIMMED_DATA
output=./05_MAPPED_DATA


echo generate STAR genome

STAR\
 --runMode genomeGenerate\
 --genomeSAindexNbases=13\
 --runThreadN=${threads}\
 --genomeDir=${GD}/STAR_GENOME\
 --genomeFastaFiles=${GD}/assembly/Mesculenta_671_v8.0_MT_PT.fa\
 --sjdbGTFfile=${GD}/annotation/Mesculenta_671_v8.1.gene_exons_MT_PT.gtf
 --sjdbOverhang=149

mkdir -p ${output}

echo star
for i in $(<${DD}/samples.txt);do

STAR\
 --runThreadN=${threads}\
 --genomeDir=${GD}/STAR_GENOME\
 --readFilesIn=${input}/${i}_1.${extension} ${input}/${i}_2.${extension}\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --alignIntronMax=30000\
 --alignMatesGapMax=30000\
 --outReadsUnmapped=Fastx\
 --quantMode GeneCounts\
 --outFilterMultimapNmax=1000\
 --outFileNamePrefix=${output}/${i}_

done

echo featureCounts

input=./05_MAPPED_DATA
output=./07_FEATURECOUNTS

mkdir -p ${output}

featureCounts\
 -t gene\
 -g ID\
 -T ${threads}\
 -p --countReadPairs\
 -a ${GD}/annotation/Mesculenta_671_v8.1.gene_exons_MT_PT.gff3\
 -o ${output}/FEATURECOUNTS_gene_UNIQUE_${i}.txt\
 ${input}/*Aligned.sortedByCoord.out.bam
