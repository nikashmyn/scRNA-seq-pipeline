#!/bin/bash

args=$@

#mind that you might be required to change the following fasta reference file:
refFasta="/singlecellcenter/etai/ReferenceData/Gencode/GRCh38.primary_assembly.genome.ERCC92.fa"
#you will need also to change this path to the location of your picard tools:
picardpath="/czlab/Software/picard/current/bin/picard.jar"

for bamfile in ${args[@]}
do
  echo "Executing Picard QC pipline on ${bamfile}..."
  java -jar ${picardpath} CollectMultipleMetrics I=${bamfile} R=${refFasta} O=${bamfile}.MultipleMetrics.txt
  
  
done

# ls *.fastq.gz.Aligned.sortedByCoord.out.bam | xargs bash ~/WORK/SEQ2MAT/CollectMultipleMetricsForMultipleBams.sh
