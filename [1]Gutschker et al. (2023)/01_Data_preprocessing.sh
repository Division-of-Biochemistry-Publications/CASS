#!/bin/bash

WD=./
DD=${WD}/01_RAW_DATA
mkdir ${WD}/02_FASTQC_OUTPUT
mkdir ${WD}/02_FASTQC_OUTPUT/MULTIQC

##############################################################################
QC
##############################################################################
fastqc --outdir ${WD}/02_FASTQC_OUTPUT -t 11 ${DD}/*.fq.gz
multiqc --outdir ${WD}/02_FASTQC_OUTPUT/MULTIQC ${WD}/02_FASTQC_OUTPUT

mkdir ${WD}/03_TRIMMED_DATA
mkdir ${WD}/03_TRIMMED_DATA/tmp
mkdir ${WD}/04_TRIMMED_DATA_FASTQC
mkdir ${WD}/04_TRIMMED_DATA_FASTQC/MULTIQC

##############################################################################
TRIMMING
##############################################################################
cd ${DD}
for i in *_1.fq.gz;do
  SAMPLE=$(echo ${i} | sed "s/_1\.fq\.gz//")
bbduk\
 k=23\
 mink=11\
 ktrim=r\
 hdist=1\
 ref=${WD}/00_META_DATA/adapters.fa\
 in1=${SAMPLE}_1.fq.gz\
 in2=${SAMPLE}_2.fq.gz\
 out1=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_1_ktrimr.fq.gz\
 out2=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_2_ktrimr.fq.gz\
 tpe=t\
 tbo=t\
 t=11
bbduk\
 k=23\
 mink=11\
 ktrim=l\
 hdist=1\
 ref=${WD}/00_META_DATA/adapters.fa\
 in1=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_1_ktrimr.fq.gz\
 in2=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_2_ktrimr.fq.gz\
 out1=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_1_ktriml.fq.gz\
 out2=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_2_ktriml.fq.gz\
 tpe=t\
 tbo=t\
 t=11
bbduk\
 qtrim=rl\
 trimq=30\
 maq=30\
 minlen=35\
 in1=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_1_ktriml.fq.gz\
 in2=${WD}/03_TRIMMED_DATA/tmp/${SAMPLE}_2_ktriml.fq.gz\
 out1=${WD}/03_TRIMMED_DATA/${SAMPLE}_1.fq.gz\
 out2=${WD}/03_TRIMMED_DATA/${SAMPLE}_2.fq.gz\
 t=11
done

rm -r ${WD}/03_TRIMMED_DATA/tmp

##############################################################################
QC
##############################################################################
fastq --outdir=${WD}/04_TRIMMED_DATA_FASTQC -t 11 ${WD}/03_TRIMMED_DATA/*.fq.gz;
multiqc --outdir ${WD}/04_TRIMMED_DATA_FASTQC/MULTIQC ${WD}/04_TRIMMED_DATA_FASTQC


##############################################################################
MAPPING
##############################################################################
DD=${WD}/03_TRIMMED_DATA
mkdir ${WD}/05_MAPPED_DATA
mkdir ${WD}/05_MAPPED_DATA/tmp
mkdir ${WD}/05_MAPPED_DATA
cd ${DD}

STAR\
 --runMode genomeGenerate\
 --genomeSAindexNbases=13\
 --runThreadN=11\
 --genomeDir=${WD}/00_META_DATA/Mesculenta/v8.1/STAR_GENOME\
 --genomeFastaFiles=${WD}/00_META_DATA/Mesculenta/v8.1/assembly/Mesculenta_671_v8.0.fa\
 --sjdbGTFfile=${WD}/00_META_DATA/Mesculenta/v8.1/annotation/Mesculenta_671_v8.1.gene_exons.gff3\
 --sjdbOverhang=149\
 --sjdbGTFtagExonParentTranscript=ID\
 --sjdbGTFfeatureExon=exon\
 --sjdbGTFtagExonParentGene=Parent

for i in *_1.fq.gz;do
  SAMPLE=$(echo ${i} | sed "s/_1\.fq\.gz//")
STAR\
 --runThreadN=11\
 --genomeDir=${WD}/00_META_DATA/Mesculenta/v8.1/STAR_GENOME\
 --readFilesIn=${WD}/03_TRIMMED_DATA/${SAMPLE}_1.fq.gz ${WD}/03_TRIMMED_DATA/${SAMPLE}_2.fq.gz\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --alignIntronMax=30000\
 --alignMatesGapMax=30000\
 --outReadsUnmapped=Fastx\
 --outFilterMultimapNmax=1000\
 --outFileNamePrefix=${WD}/05_MAPPED_DATA/tmp/${SAMPLE}_
done

for i in *_1.fq.gz;do
  SAMPLE=$(echo ${i} | sed "s/_1\.fq\.gz//")
STAR\
 --runThreadN=11\
 --genomeDir=${WD}/00_META_DATA/Mesculenta/v8.1/STAR_GENOME\
 --readFilesIn=${WD}/03_TRIMMED_DATA/${SAMPLE}_1.fq.gz ${WD}/03_TRIMMED_DATA/${SAMPLE}_2.fq.gz\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --alignIntronMax=30000\
 --alignMatesGapMax=30000\
 --outReadsUnmapped=Fastx\
 --outFilterMultimapNmax=1000\
 --outMultimapperOrder=Random\
 --runRNGseed=1234\
 --outFileNamePrefix=${WD}/05_MAPPED_DATA/${SAMPLE}_\
 --sjdbFileChrStartEnd ${WD}/05_MAPPED_DATA/tmp/*_SJ.out.tab
done

rm -r ${WD}/05_MAPPED_DATA/tmp

##############################################################################
FEATURECOUNTS
##############################################################################
mkdir ${WD}/07_FEATURECOUNTS

for feature in gene CDS exon;do
featureCounts\
 -t ${feature}\
 -g ID\
 -T 11\
 -M\
 --fraction\
 -p\
 -a ${WD}/00_META_DATA/Mesculenta/v8.1/annotation/Mesculenta_671_v8.1.gene_exons.gff3\
 -o ${WD}/07_FEATURECOUNTS/FEATURECOUNTS_${feature}_FRACTION.txt\
 ${WD}/05_MAPPED_DATA/*_Aligned.sortedByCoord.out.bam
 
featureCounts\
 -t ${feature}\
 -g ID\
 -T 11\
 -M\
 --primary\
 -p\
 -a ${WD}/00_META_DATA/Mesculenta/v8.1/annotation/Mesculenta_671_v8.1.gene_exons.gff3\
 -o ${WD}/07_FEATURECOUNTS/FEATURECOUNTS_${feature}_PRIMARY.txt\
 ${WD}/05_MAPPED_DATA/*_Aligned.sortedByCoord.out.bam

featureCounts\
 -t ${feature}\
 -g ID\
 -T 11\
 -p --countReadPairs\
 -a ${WD}/00_META_DATA/Mesculenta/v8.1/annotation/Mesculenta_671_v8.1.gene_exons.gff3\
 -o ${WD}/07_FEATURECOUNTS/FEATURECOUNTS_${feature}_UNIQUE.txt\
 ${WD}/05_MAPPED_DATA/*_Aligned.sortedByCoord.out.bam
done

##############################################################################
VARIANT CALLING
##############################################################################
mkdir ${WD}/15_VARIANT_CALLING

bcftools mpileup \
--fasta-ref ${WD}/00_META_DATA/Mesculenta/v8.1/assembly/Mesculenta_671_v8.0.fa \
--threads 15 \
--skip-indels \
--output ${WD}/15_VARIANT_CALLING/PYT52.bcf \
--output-type b \
${WD}/05_MAPPED_DATA/*_Aligned.sortedByCoord.out.bam \

bcftools call \
-m \
--output ${WD}/15_VARIANT_CALLING/PYT52.vcf \
--output-type v \
--threads 15 \
--variants-only \
${WD}/15_VARIANT_CALLING/PYT52.bcf


 
