#!/usr/bin/env bash

# Path to hg38.fa is needed

# Default values
sample_name="sample"
sample_sex="M"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -name)
            sample_name="$2"
            shift 2
            ;;
        -sex)
            sample_sex="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            echo "Usage: $0 -name <sample_name> -sex <M/F>"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$sample_name" ]; then
    echo "Error: Must specify sample name using -name parameter"
    echo "Usage: $0 -name <sample_name> -sex <M/F>"
    exit 1
fi

if [ -z "$sample_sex" ] || { [ "$sample_sex" != "M" ] && [ "$sample_sex" != "F" ]; }; then
    echo "Error: Must specify sample sex using -sex parameter (M or F)"
    echo "Usage: $0 -name <sample_name> -sex <M/F>"
    exit 1
fi

echo "Processing sample: $sample_name, sex: $sample_sex"

# Mapping to reference genome using bwa
# Parameters can be configured to the appropriate values:
# -t 16
# -R "@RG\tID:"${sample_name}"\tSM:"${sample_name}"\tPL:DIPSEQ\tLB:"${sample_name}"\tPU:"${sample_name}
# -@16

bwa mem -M -t 16 -R "@RG\tID:"${sample_name}"\tSM:"${sample_name}"\tPL:DIPSEQ\tLB:"${sample_name}"\tPU:"${sample_name} /path/to/hg38.fa ${sample_name}_1_clean.fq.gz ${sample_name}_2_clean.fq.gz | samtools sort -@16 -o ${sample_name}.bam -


# Tandem repeat genotyping using GangSTR

GangSTR --bam-samps ${sample_name} --samp-sex ${sample_sex} --bam ${sample_name}.bam --ref /path/to/hg38.fa --regions hg38_TR_Atlas.bed --out ${sample_name}_GangSTR


# GangSTR genotyping results filtering
# Parameters can be configured to the appropriate values:
# --gangstr-min-call-DP 10 
# --gangstr-min-call-Q 0.9 
# --gangstr-max-call-DP 1000 
# --gangstr-filter-spanbound-only 
# --gangstr-filter-badCI 
# --filter-regions hg38_segdup.sorted.bed.gz 
# --filter-regions-names SEGDUP

dumpSTR --vcf ${sample_name}_GangSTR.vcf --out ${sample_name}_GangSTR_dumpSTR --gangstr-min-call-DP 10 --gangstr-min-call-Q 0.9 --gangstr-max-call-DP 1000 --gangstr-filter-spanbound-only --gangstr-filter-badCI --filter-regions hg38_segdup.sorted.bed.gz --filter-regions-names SEGDUP
bcftools view -f 'PASS' ${sample_name}_GangSTR_dumpSTR.vcf -o ${sample_name}_GangSTR_dumpSTR_PASS.vcf


# Tandem repeat genotyping using ExpansionHunter

ExpansionHunter --reads ${sample_name}.bam --reference /path/to/hg38.fa --variant-catalog variant_catalog.json --output-prefix ${sample_name}_EH -m streaming -n 16


# ExpansionHunter genotyping results filtering
# Parameters can be configured to the appropriate values:
# --eh-min-call-LC 10 
# --filter-regions hg38_segdup.sorted.bed.gz 
# --filter-regions-names SEGDUP

dumpSTR --vcf ${sample_name}_EH.vcf --out ${sample_name}_EH_dumpSTR --eh-min-call-LC 10 --filter-regions hg38_segdup.sorted.bed.gz --filter-regions-names SEGDUP
bcftools view -f 'PASS' ${sample_name}_EH_dumpSTR.vcf -o ${sample_name}_EH_dumpSTR_PASS.vcf


# Genotyping results of GangSTR and ExpansionHunter were integrated using ensembleTR.
EnsembleTR --out ${sample_name}_ensembletr.vcf --ref /path/to/hg38.fa --vcfs ${sample_name}_GangSTR_dumpSTR_PASS.vcf,${sample_name}_EH_dumpSTR_PASS.vcf