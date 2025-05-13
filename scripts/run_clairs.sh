# Parameters


#load data path
# TUMOR_BAM_DIR="/big8_disk/mingen112/test_data/HCC1395/ONT/orig_bam/tumor"
# NORMAL_BAM_DIR="/big8_disk/mingen112/test_data/HCC1395/ONT/orig_bam/normal"
TUMOR_BAM_DIR="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG"
NORMAL_BAM_DIR="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG"
REF_DIR="/big8_disk/ref"
HIGH_CON_DIR="/big8_disk/data/HCC1395/SEQC2"


REF="GRCh38_no_alt_analysis_set.fasta"
# NORMAL_BAM="alignment-sort-hcc1395bl.bam"
# TUMOR_BAM="alignment-sort-hcc1395.bam"
NORMAL_BAM="HCC1395BL.bam"
# TUMOR_BAM=("merged_t10_n40.bam" "merged_t20_n30.bam" "merged_t30_n20.bam" "merged_t40_n10.bam")
TUMOR_BAM=("HCC1395.bam")


#benchmark
BASELINE_VCF_FILE_PATH="high-confidence_sSNV_in_HC_regions_v1.2.1.vcf"
BASELINE_BED_FILE_PATH="High-Confidence_Regions_v1.2.bed"

OUTPUT_VCF_FILE_PATH="output.vcf.gz"



paths=(
    "${TUMOR_BAM_DIR}/${TUMOR_BAM}"
    "${NORMAL_BAM_DIR}/${NORMAL_BAM}"
    "${REF_DIR}/${REF}"
    "${HIGH_CON_DIR}/${BASELINE_VCF_FILE_PATH}"
    "${HIGH_CON_DIR}/${BASELINE_BED_FILE_PATH}"
    # "${OUTPUT_DIR}"
)


for path in "${paths[@]}"; do
    if [ ! -e "$path" ]; then
        echo "[Error] not exist file: $path"
        exit 1
    fi
done 

for tumor_bam in "${TUMOR_BAM[@]}"; do
    echo "Processing tumor bam: $tumor_bam"

    OUTPUT_DIR="/big8_disk/liaoyoyo2001/ClairS_ssrs_output"
    mkdir -p ${OUTPUT_DIR}

    # platform: ont_r10_dorado_sup_5khz_ssrs
    docker run \
    -v ${TUMOR_BAM_DIR}:${TUMOR_BAM_DIR} \
    -v ${NORMAL_BAM_DIR}:${NORMAL_BAM_DIR} \
    -v ${REF_DIR}:${REF_DIR} \
    -v ${HIGH_CON_DIR}:${HIGH_CON_DIR} \
    -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
    hkubal/clairs:v0.4.1 \
    /opt/bin/run_clairs \
    --tumor_bam_fn ${TUMOR_BAM_DIR}/${tumor_bam} \
    --normal_bam_fn ${NORMAL_BAM_DIR}/${NORMAL_BAM} \
    --ref_fn ${REF_DIR}/${REF} \
    --threads 80 \
    --platform ont_r10_dorado_sup_5khz_ssrs \
    --enable_indel_calling \
    --output_dir ${OUTPUT_DIR} 
done

# docker run -it \
#   -v ${TUMOR_BAM_DIR}:${TUMOR_BAM_DIR} \
#   -v ${NORMAL_BAM_DIR}:${NORMAL_BAM_DIR} \
#   -v ${REF_DIR}:${REF_DIR} \
#   -v ${HIGH_CON_DIR}:${HIGH_CON_DIR} \
#   -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
#   hkubal/clairs:v0.4.0 \
#   python3 /opt/bin/clairs.py compare_vcf \
#      --truth_vcf_fn ${HIGH_CON_DIR}/${BASELINE_VCF_FILE_PATH} \
#      --input_vcf_fn ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
#      --bed_fn ${HIGH_CON_DIR}/${BASELINE_BED_FILE_PATH} \
#      --output_dir ${OUTPUT_DIR}/benchmark \
#      --input_filter_tag 'PASS' \