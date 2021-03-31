#!/bin/bash

#/czlab/etai/ExperimentsData/Stamatis/run.v2/star/
inputfile=$1 #/czlab/etai/ExperimentsData/Stamatis/run.v2/star/170216_A4.R1.fastq.gz.Aligned.sortedByCoord.out.bam
outputdir=$2 #/czlab/etai/ExperimentsData/Stamatis/run.v2/gatk
SM=`basename ${inputfile} | awk -F'.' '{print $1}'`;

#make sure to set these file names:
genome_ref="/singlecellcenter/etai/ReferenceData/Gencode/GRCh38.primary_assembly.genome.ERCC92.fa"
GATK="/czlab/Software/GATK/current/GenomeAnalysisTK.jar"
picard="/czlab/Software/picard/current/bin/picard.jar"

java -Xmx4g -jar ${picard} AddOrReplaceReadGroups INPUT=${inputfile} \
RGID=${SM} RGSM=${SM} RGLB=${SM} RGPL=illumina RGPU=MiOrHiSeq OUTPUT=${inputfile}.withRG.bam \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
1> ${outputdir}/${SM}.AddOrReplaceReadGroups.stdouterr 2&>1


java -jar ${picard} BuildBamIndex I=${inputfile}.withRG.bam

java -Xmx16g -jar ${GATK} -T SplitNCigarReads -R ${genome_ref} -I ${inputfile}.withRG.bam \
-o ${outputdir}/${SM}.gatk.bam \
#-rf ReassignOneMappingQuality \
#-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS \
-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM \
1> ${outputdir}/${SM}.SplitNCigarReads.stdouterr 2&>1

# java -Xmx4g -jar ${picard} AddOrReplaceReadGroups INPUT=${outputdir}/${SM}.splitN.bam \
# RGID=${SM} RGSM=${SM} RGLB=${SM} RGPL=illumina RGPU=MiOrHiSeq OUTPUT=${outputdir}/${SM}.gatk.bam \
# CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
# 1> ${outputdir}/${SM}.AddOrReplaceReadGroups.stdouterr 2>&1

#rm ${outputdir}/${SM}.splitN.ba*

echo completed ${inputfile}
