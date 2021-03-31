

inputdir=$1 #input bam or fastq files are located e.g. "/papathanasiou/SIS1025d/"
outdir=$2 #folder (folder should exist) where project data and all results should be stored (path should exist). e.g. "/singlecellcenter/etai/ExperimentsData/Stamatis/SIS1025d/"
bamorfastq=$3 #fastq.gz - for fastq.gz input: R1 and R2 .fastq.gz. bam - for bam input
numofcpus=$4 #number of cpus to run the workflow
scriptsdir=$5 #path to the folder which includes all auxiliary scripts
runid=$6 #a name/id to be associated with this execution

##################################################################
#Bash basic genotype pipeline cmds (from fastq to TPM and SNP RDs):
##################################################################

#trasforming bam files to fastqs:
#--------------------------------
if [ $bamorfastq == "bam" ]
then
  echo "Transforming bams to fastqs..."
  mkdir ${outputdir}/bams
  cd ${outputdir}/bams
  parallel -j ${numofcpus} samtools sort -n {} ${PWD}/{/}.qsort ::: ${inputdir}/*.bam &
  mkdir ${outputdir}/fastqs
  cd ${outputdir}/fastqs
  parallel -j ${numofcpus} bedtools bamtofastq -i {} -fq ${PWD}/{/}.R1.fastq -fq2 ${PWD}/{/}.R2.fastq  ::: ls ${outputdir}/bams/*.bam
  parallel -j ${numofcpus} gzip {}  ::: ls ${PWD}/*.fastq
else
  mkdir ${outputdir}/fastqs
  cd ${outputdir}/fastqs
  cp -s ${inputdir}/*.fastq.gz ./ #linux command only - link fastq files to run directory
fi

#Running Pipeline:
#-----------------
echo "Running alignment pipeline..."
cd ${outputdir}/
bash ${scriptsdir}/runBasic_RNA_pipeline_on_fastqs.sh ${outputdir}/ ${numofcpus} *bam.R1.fastq.gz > ${runid}.runBasic_RNA_pipeline_on_fastqs.stdouterr 2>&1

#Running genotyping:
echo "Running genotyping..."
cd ${outputdir}
mkdir ${outputdir}/vcf
cd ${outputdir}/vcf
ls ${outputdir}/gatk/*.bam > bamfiles.list
bash ${scriptsdir}/RPE-1_GRCh38_Genotype_etai.sh ${PWD}/bamfiles.list ${outputdir}/vcf/${runid}.gatkBamFiles 1
parallel --jobs ${numofcpus} < ${runid}.gatkBamFiles_RPE_hets_GT.UG_jobs.list > ${runid}.gatkBamFiles_RPE_hets_GT.UG_jobs.stdouterr 2>&1
Rscript ${scriptsdir}/prepareGenotypeDataFromRNAvcfs.R 14 ${outputdir}/vcf ${outputdir}/vcf/${runid}.gatkBamFiles_RPE_hets.GT_UG.%s.vcf

