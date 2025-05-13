# 1. 抽 header + PASS
awk -F'\t' '/^#/ || $7=="PASS"' \
  /big8_disk/liaoyoyo2001/ClairS_ss_output/tmp/vcf_output/pileup_filter.vcf \
> /big8_disk/liaoyoyo2001/data/vcf/ClairS_ss/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS.vcf

# 2. 抽 header + PASS，同時修改 GQ Type=Integer→Float
awk -F'\t' '
  /^##FORMAT=<ID=GQ/ { sub("Type=Integer","Type=Float") } 
  /^#/ || $7=="PASS"
' \
  /big8_disk/liaoyoyo2001/ClairS_ss_output/tmp/vcf_output/pileup_filter.vcf \
> /big8_disk/liaoyoyo2001/data/vcf/ClairS_ss/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS_fixed.vcf

cd /big8_disk/liaoyoyo2001/data/vcf/ClairS_ss/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/
bcftools view -Oz HCC1395_methyl_PASS_fixed.vcf -o HCC1395_methyl_PASS_fixed.vcf.gz
bcftools index HCC1395_methyl_PASS_fixed.vcf.gz

./run_longphase_somatic_benchmark.sh > longphase_somatic_HCC1395_methyl_ss.log  2>&1 &

cd /big8_disk/liaoyoyo2001/longphase_somatic_output_ssrs/
samtools index HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam