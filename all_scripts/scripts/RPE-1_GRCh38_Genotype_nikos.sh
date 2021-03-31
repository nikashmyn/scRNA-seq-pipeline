#!/bin/bash

# Input bams list
while read line  
do
#echo ${line}
input_file_list="${input_file_list} -I ${line}"	
done < $1 
echo ${input_file_list}
NAME=$2
NTHREADS=$3
reference=$4
genotype_file=$5
GATK=$6
picard=$7

output_type=RPE_hets;

tmp=".tmp";
run_GATK_UG="java -Xmx32000m -jar ${GATK} --reference_sequence ${reference} --num_cpu_threads_per_data_thread ${NTHREADS} --output_mode EMIT_ALL_SITES --standard_min_confidence_threshold_for_calling 0"; # -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar";
run_GATK_HC="java -Xmx32000m -jar ${GATK} --reference_sequence ${reference} --num_cpu_threads_per_data_thread ${NTHREADS} --output_mode EMIT_ALL_SITES --standard_min_confidence_threshold_for_calling 0"; # -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar";

Genotype_UG="${run_GATK_UG} --genotyping_mode GENOTYPE_GIVEN_ALLELES --analysis_type UnifiedGenotyper --num_threads ${NTHREADS}"
Genotype_HC="${run_GATK_HC} --genotyping_mode GENOTYPE_GIVEN_ALLELES --analysis_type HaplotypeCaller --min_mapping_quality_score 30"

HC_cmd="${Genotype_HC} ${input_file_list} --alleles ${genotype_file}" 
UG_cmd="${Genotype_UG} ${input_file_list} --alleles ${genotype_file}"

Discovery_cmd="${run_GATK_HC} ${input_file_list} --genotyping_mode DISCOVERY --min_mapping_quality_score 30"

HC_output=${NAME}_${output_type}_GT.HC_jobs
UG_output=${NAME}_${output_type}_GT.UG_jobs
Discovery_output=${NAME}_HC_Discovery_jobs

output_type=${output_type}.GT
output_prefix="$NAME";


for i in {1..22}
do
	echo "$HC_cmd -L chr${i} --out ${output_prefix}_${output_type}_HC.${i}.vcf" > ${HC_output}.chr${i}.sh
	echo "$UG_cmd -L chr${i} --out ${output_prefix}_${output_type}_UG.${i}.vcf" > ${UG_output}.chr${i}.sh
	echo "$Discovery_cmd -L chr${i} --out ${output_prefix}_HC_Discovery.${i}.vcf" > ${Discovery_output}.chr${i}.sh
done

echo "$HC_cmd -L chrX --out ${output_prefix}_${output_type}_HC.X.vcf" > ${HC_output}.chrX.sh
echo "$UG_cmd -L chrX --out ${output_prefix}_${output_type}_UG.X.vcf"> ${UG_output}.chrX.sh
echo "$Discovery_cmd -L chrX --out ${output_prefix}_HC_Discovery.X.vcf" > ${Discovery_output}.chrX.sh
echo "$HC_cmd -L chrY --out ${output_prefix}_${output_type}_HC.Y.vcf" > ${HC_output}.chrY.sh
echo "$UG_cmd -L chrY --out ${output_prefix}_${output_type}_UG.Y.vcf"> ${UG_output}.chrY.sh
echo "$Discovery_cmd -L chrY --out ${output_prefix}_HC_Discovery.Y.vcf" > ${Discovery_output}.chrY.sh
echo "$HC_cmd -L chrM --out ${output_prefix}_${output_type}_HC.MT.vcf" > ${HC_output}.chrM.sh
echo "$UG_cmd -L chrM --out ${output_prefix}_${output_type}_UG.MT.vcf" > ${UG_output}.chrM.sh
echo "$Discovery_cmd -L chrM --out ${output_prefix}_HC_Discovery.MT.vcf" > ${Discovery_output}.chrM.sh
