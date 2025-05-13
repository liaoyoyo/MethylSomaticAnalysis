#!/bin/bash

thread_numbers=64

toolsPath="/big8_disk/liaoyoyo2001/HCC1395_methylation_analysis/tools"
ref_fasta="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
somatic_tagging_longphase_path="/big8_disk/liaoyoyo2001/HCC1395_methylation_analysis/longphase-somatic-v1.73"

###-------------------HCC1395 methylation-------------------
hcc1395_tumor_bam="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
hcc1395_normal_bam="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
hcc1395_clairs_v040_pileup_vcf="/big8_disk/liaoyoyo2001/data/vcf/ClairS_ss/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS.vcf"
hcc1395_clairs_v040_normal_phased_vcf="/big8_disk/liaoyoyo2001/data/vcf/HCC1395BL_methyl_phase.vcf.gz"
somatic_out_path="/big8_disk/liaoyoyo2001/longphase_somatic_output_ss"
somatic_snp_benchmark_result_path="/big8_disk/liaoyoyo2001/data/longphase_somatic_result_ss/"
consider_bed="true"
hcc1395_seqc_high_conf_vcf="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf"
#hcc1395_seqc_super_set_vcf="/big8_disk/zhenyu112/HCC1395/SEQC2/sSNV.MSDUKT.superSet.v1.2.vcf"
hcc1395_seqc_high_conf_bed="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"

bash ${toolsPath}/somatic_benchmark.sh \
    --toolsPath ${toolsPath} \
    --tumor_bam ${hcc1395_tumor_bam} \
    --tumor_vcf ${hcc1395_clairs_v040_pileup_vcf} \
    --normal_bam ${hcc1395_normal_bam} \
    --normal_vcf ${hcc1395_clairs_v040_normal_phased_vcf} \
    --reference_fasta ${ref_fasta} \
    --thread_numbers ${thread_numbers} \
    --consider_bed ${consider_bed} \
    --high_conf_vcf ${hcc1395_seqc_high_conf_vcf} \
    --high_conf_bed ${hcc1395_seqc_high_conf_bed} \
    --somatic_tag_out_dir ${somatic_out_path} \
    --read_log_path ${somatic_out_path} \
    --somatic_snp_benchmark_result_path ${somatic_snp_benchmark_result_path} \
    --data_set_name "HCC1395_methylation"
