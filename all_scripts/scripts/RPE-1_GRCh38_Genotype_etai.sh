#!/bin/bash

GATK='/czlab/Software/GATK/current/GenomeAnalysisTK.jar';
picard='/czlab/Software/picard/current/bin/picard.jar';

# GRCh38 reference
reference='/czlab/References/GRCh38/GRCh38.primary_assembly.genome.fa'; 

#CNV12.1 - not in use: (a file copy also exists on google drive)
#genotypes_dir='/singlecellcenter/RPE-1/Genotyping/v.2'
#genotype_file=${genotypes_dir}/RPE-1.hets.SNP.vcf.gz #RPE-1.hets.chr1-X.BA.SNP.sites.vcf.gz;

#CNV12.0: (a file copy also exists on google drive - you will need to specifiy a location for this file)
genotype_file=/czlab/Data/RPE-1_Genotyping_GRCh38/v.1/RPE-1.hets.vcf.gz #/czlab/Data/RPE-1_Genotyping_GRCh38/v.0/RPE-1.hets.vcf #/czlab/Data/RPE-1_Genotyping_GRCh38/RPE-1.hets.vcf
output_type=RPE_hets;

# Input bams list
while read line  
do
#echo ${line}
input_file_list="${input_file_list} -I ${line}"	
done < $1 
echo ${input_file_list}
NAME=$2
NTHREADS=$3

tmp=".tmp";
run_GATK="java -Xmx32000m -jar ${GATK} --reference_sequence ${reference} --num_cpu_threads_per_data_thread ${NTHREADS} --output_mode EMIT_ALL_SITES --standard_min_confidence_threshold_for_calling 0"; # -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar";

Genotype_UG="${run_GATK} --genotyping_mode GENOTYPE_GIVEN_ALLELES --analysis_type UnifiedGenotyper --num_threads ${NTHREADS}"
Genotype_HC="${run_GATK} --genotyping_mode GENOTYPE_GIVEN_ALLELES --analysis_type HaplotypeCaller --min_mapping_quality_score 30"

HC_cmd="${Genotype_HC} ${input_file_list} --alleles ${genotype_file}" 
UG_cmd="${Genotype_UG} ${input_file_list} --alleles ${genotype_file}"

Discovery_cmd="${run_GATK} ${input_file_list} --genotyping_mode DISCOVERY --analysis_type HaplotypeCaller --min_mapping_quality_score 30"

HC_output=${NAME}_${output_type}_GT.HC_jobs.list
UG_output=${NAME}_${output_type}_GT.UG_jobs.list
Discovery_output=${NAME}_HC_Discovery_jobs.list

output_type=${output_type}.GT
output_prefix="$NAME";


for i in {1..22}
do
	echo "$HC_cmd -L chr${i} --out ${output_prefix}_${output_type}_HC.${i}.vcf" >> ${HC_output}
	echo "$UG_cmd -L chr${i} --out ${output_prefix}_${output_type}_UG.${i}.vcf">> ${UG_output}
	echo "$Discovery_cmd -L chr${i} --out ${output_prefix}_HC_Discovery.${i}.vcf" >> ${Discovery_output}
done

echo "$HC_cmd -L chrX --out ${output_prefix}_${output_type}_HC.X.vcf" >> ${HC_output}
echo "$UG_cmd -L chrX --out ${output_prefix}_${output_type}_UG.X.vcf">> ${UG_output}
echo "$Discovery_cmd -L chrX --out ${output_prefix}_HC_Discovery.X.vcf" >> ${Discovery_output}
echo "$HC_cmd -L chrY --out ${output_prefix}_${output_type}_HC.Y.vcf" >> ${HC_output}
echo "$UG_cmd -L chrY --out ${output_prefix}_${output_type}_UG.Y.vcf">> ${UG_output}
echo "$Discovery_cmd -L chrY --out ${output_prefix}_HC_Discovery.Y.vcf" >> ${Discovery_output}
echo "$HC_cmd -L chrM --out ${output_prefix}_${output_type}_HC.MT.vcf" >> ${HC_output}
echo "$UG_cmd -L chrM --out ${output_prefix}_${output_type}_UG.MT.vcf" >> ${UG_output}
echo "$Discovery_cmd -L chrM --out ${output_prefix}_HC_Discovery.MT.vcf" >> ${Discovery_output}
