#!/bin/bash

# 用途
# 輸出 snv_tp snv_fp snv_fn  的 VCF 檔案
# 輸出 PASS 的 VCF 檔案，同時修改 GQ Type=Integer→Float
# 壓縮 VCF 並產生 .gz 檔
# 利用 tabix 產生 VCF 索引 (.tbi)
# 可選：若不需要保留未壓縮的 VCF 檔，可刪除之

### 設定區域 ###
# 預設參數 (可修改)
SNV_DIR="/big8_disk/liaoyoyo2001/data/longphase_somatic_result_ssrs/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter_HP3.out"
VCF_FILE="/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS_fixed.vcf"
OUTPUT_DIR="/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup"  # 預設輸出目錄為當前目錄

# 可修改的 SNV 檔案清單
declare -a SNV_FILES=("snv_tp.out" "snv_fp.out" "snv_fn.out")

### 允許使用者提供參數 ###
while [[ $# -gt 0 ]]; do
    case $1 in
        --snv-dir)
            SNV_DIR="$2"
            shift 2
            ;;
        --vcf-file)
            VCF_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --snv-files)
            IFS=',' read -r -a SNV_FILES <<< "$2"  # 允許以逗號分隔的字串作為 SNV 檔案清單
            shift 2
            ;;
        *)
            echo "未知參數: $1"
            echo "使用方式: $0 [--snv-dir 路徑] [--vcf-file 路徑] [--output-dir 路徑] [--snv-files 檔案1,檔案2,...]"
            exit 1
            ;;
    esac
done

### 確保輸出目錄存在 ###
mkdir -p "$OUTPUT_DIR"

# 決定解壓指令
if [[ "${VCF_FILE##*.}" == "gz" ]]; then
    DECOMP_CMD="zcat"
else
    DECOMP_CMD="cat"
fi

for snv in "${SNV_FILES[@]}"; do
    SNV_PATH="${SNV_DIR}/${snv}"
    OUTPUT_FILE="${OUTPUT_DIR}/filtered_${snv%.*}.vcf"

    echo "處理 ${SNV_PATH} → ${OUTPUT_FILE}"
    awk 'FNR==NR {
            if ($1 !~ /^#/) keys[$1"_"$2]=1
            next
        }
        /^#/ { print; next }
        { if (($1"_"$2) in keys) print }' \
        "$SNV_PATH" <($DECOMP_CMD "$VCF_FILE") > "$OUTPUT_FILE"

    bgzip -c "$OUTPUT_FILE" > "${OUTPUT_FILE}.gz"
    tabix -p vcf "${OUTPUT_FILE}.gz"

    echo "完成：${OUTPUT_FILE}.gz 及索引 ${OUTPUT_FILE}.gz.tbi"
done

# 只保留 chr19 的 snv_tp
bcftools view -r chr19 filtered_snv_tp.vcf.gz -Oz -o filtered_snv_tp_chr19.vcf.gz
tabix -p vcf filtered_snv_tp_chr19.vcf.gz
bcftools view -r chr19 filtered_snv_fp.vcf.gz -Oz -o filtered_snv_fp_chr19.vcf.gz
tabix -p vcf filtered_snv_fp_chr19.vcf.gz
bcftools view -r chr19 filtered_snv_fn.vcf.gz -Oz -o filtered_snv_fn_chr19.vcf.gz
tabix -p vcf filtered_snv_fn_chr19.vcf.gz