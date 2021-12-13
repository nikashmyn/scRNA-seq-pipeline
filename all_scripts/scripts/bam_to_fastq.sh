inputdir=$1 #input bam or fastq files are located e.g. "/papathanasiou/SIS1025d/"
outputdir=$2 #folder (folder should exist) where project data and all results should be stored (path should exist). e.g. "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/"
numofcpus=$3 #number of cpus to run the workflow

#trasforming bam files to fastqs:
#--------------------------------
echo "Transforming bams to fastqs..."
mkdir -p ${outputdir}
cd ${outputdir}
parallel -j ${numofcpus} samtools sort -n {} -o ${PWD}/{/}.qsort ::: ${inputdir}/*.bam # Maybe add back & at end of this line if it is broken
parallel -j ${numofcpus} bedtools bamtofastq -i {} -fq ${PWD}/{/}.R1.fastq -fq2 ${PWD}/{/}.R2.fastq  ::: ls ${outputdir}/*.bam.qsort
parallel -j ${numofcpus} gzip {}  ::: ls ${PWD}/*.fastq

