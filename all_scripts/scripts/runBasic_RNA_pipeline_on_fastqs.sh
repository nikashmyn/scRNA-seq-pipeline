#!/bin/bash

#Written by Etai Jacob, PhD

#runBasic_scRNA_MiSeq_pipeline.sh
#bash runBasic_RNA_pipeline_on_fastqs.sh /czlab/etai/ExperimentsData/RPE_bulk/Selwin_Trisomies/ 3 *.1.fastq.gz >& SelwinTris.basicRun.stdouterr &


dstDir=$1 
cpus=$2 
#pattern=$3 #fastq file names input pattern e.g. *.R1.fastq.gz
scriptsdir=$3 #path to the folder which includes all auxiliary scripts
fastqctool="fastqc" #need to set this up for your fastqc tool path


echo ${dstDir}
echo ${cpus}
#echo ${pattern}
echo ${scriptsdir}

#in this version the genome ref and anno dir is defaulted to gencode/hg38
#TODO: add parameter in this script to set also 

echo "Running alignments and basic analysis on the following files:"
echo $PWD
cd ${dstDir}
#ls ${dstDir}/fastqs/${pattern}

#echo "Running alignments and QC.."
#  mkdir -p ${dstDir}/star
#STAR
#echo "Running STAR..."
#  Rscript ${scriptsdir}/generateRun_STARcmds.R ${dstDir}/fastqs ${dstDir}/star ${scriptsdir} > ${dstDir}/star.cmds
#  echo  ${dstDir}/star.cmds
#  parallel --jobs ${cpus} < star.cmds &> star.cmds.stdouterr
#  mkdir -p ${dstDir}/fastqc
  # this next command sues fastqctool and there is a display error with this part. Doesnt affect run. 
#  find ${dstDir}/star -name '*.Aligned.sortedByCoord.out.bam' | xargs ${fastqctool} -o ${dstDir}/fastqc -f bam &>  ${dstDir}/fastqc.stdouterr &

#RSEM
#echo "Runing RSEM..."
#  mkdir -p ${dstDir}/rsem
#  Rscript ${scriptsdir}/generateRSEMcalcExpCmds.R ${dstDir}/fastqs ${dstDir}/rsem ${scriptsdir} > ${dstDir}/rsem.cmds
#  parallel --jobs ${cpus} < rsem.cmds &> rsem.cmds.stdouterr

#GATK
echo "Prepare files for GATK runs..."
  mkdir -p ${dstDir}/gatk
  Rscript ${scriptsdir}/generatePrepareSTARBamOutputForGATK_Genotyping_cmds.R ${dstDir}/star ${dstDir}/gatk *.Aligned.sortedByCoord.out.bam$ ${scriptsdir} > ${dstDir}/prepareGATK.cmds
  parallel --jobs ${cpus} < prepareGATK.cmds &> prepareGATK.cmds.stdouterr
  
#echo "Executing Picard QC.."
#echo ${scriptsdir}
#parallel --jobs ${cpus} bash ${scriptsdir}/CollectMultipleMetricsForMultipleBams.sh {}  ::: ${dstDir}/star/*.Aligned.sortedByCoord.out.bam &


