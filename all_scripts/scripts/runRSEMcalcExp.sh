#!/bin/bash

#need to make sure this is the correct path where star software is installed
star_path='/usr/local/bin/'

#default locations for reference files - need to set these in case you use default:
gtffile='/singlecellcenter/etai/ReferenceData/Gencode/gencode.v25.primary_assembly.ERCC92.annotation.gtf'
rsemrefpath='/singlecellcenter/etai/ReferenceData/RSEM/human/human'
refgenome='/singlecellcenter/etai/ReferenceData/Gencode/GRCh38.primary_assembly.genome.ERCC92.fa'
    
args=( "$@" )

#input params:
#1. Paired end (2) or not (1)
#2. fastq1
#3. fastq2
#4. outprefix
#5. [gtffile]
#6. [rsemrefpath]
#7. [refgenome]

if [ "$#" -ge 2 ]
then
  echo "Start execution.."
else
  echo "Illegal number of parameters"
  echo "2/1 fastq1 fastq2/ outprefix [gtffile rsemrefpath refgenome]"
  exit 0
fi

echo "Having: --------- $#"
if [ $((${args[0]} + 0)) -eq 2 ]
then
  rf1=${args[1]}
  rf2=${args[2]}
  out=${args[3]}
  pe='--paired-end'
  if  [ "$#" -eq 7 ]
  then
    echo "Input includes genome reference"
    gtffile=${args[4]}
    rsemrefpath=${args[5]}
    refgenome=${args[6]}
  else
    echo "Using default reference - human v25 - make sure that this is what you wish to use"
  fi
  
  #outdir=${args[3]}
else
  rf1=${args[1]}
  rf2=''
  out=${args[2]}
  pe=''
  if  [ "$#" -eq 5 ]
  then
    echo "Input includes genome reference"
    gtffile=${args[3]}
    rsemrefpath=${args[4]}
    refgenome=${args[5]}
  else
  #need to make sure this is the default path to your reference files
    echo "Using default reference - human v25 - make sure that this is what you wish to use"
  fi
  #outdir=${args[2]}
fi


cpus=4

rsem-calculate-expression -p ${cpus} ${pe} \
--star --star-path ${star_path} \
--star-gzipped-read-file \
--estimate-rspd \
--append-names \
${rf1} ${rf2} \
${rsemrefpath} ${out}.rsem

#--output-genome-bam \
#bash ~/WORK/RNA/Workflows/runRSEMcalcExp.sh 2 170119_M01209_0014_000000000-AWW4H_Tri12_S2_L001_R1_001.fastq.gz 170119_M01209_0014_000000000-AWW4H_Tri12_S2_L001_R2_001.fastq.gz