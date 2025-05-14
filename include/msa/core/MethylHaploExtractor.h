#pragma once

#include <string>
#include <vector>
#include <htslib/sam.h>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief 甲基化與單倍型提取器，從BAM讀段中提取甲基化與單倍型信息
 */
class MethylHaploExtractor {
public:
    /**
     * @brief 建構函數
     * @param config 配置物件
     */
    MethylHaploExtractor(const msa::Config& config);
    
    /**
     * @brief 從讀段中提取甲基化與單倍型信息
     * @param read BAM讀段
     * @param target_variant 目標變異
     * @param bam_source_id BAM來源ID (Tumor/Normal)
     * @return std::vector<msa::MethylationSiteDetail> 提取的甲基化位點詳情
     */
    std::vector<msa::MethylationSiteDetail> extractFromRead(
        const bam1_t* read,
        const msa::VcfVariantInfo& target_variant,
        const std::string& bam_source_id
    );
    
private:
    const msa::Config& config_;  // 配置物件
    
    /**
     * @brief 從BAM讀段中提取甲基化信息
     * @param read BAM讀段
     * @param target_variant 目標變異
     * @param bam_source_id BAM來源ID
     * @param haplotype_tag 單倍型標籤
     * @param somatic_allele_type 等位基因類型
     * @param details 用於存儲提取的甲基化詳情
     */
    void extractMethylation(
        const bam1_t* read,
        const msa::VcfVariantInfo& target_variant,
        const std::string& bam_source_id,
        const std::string& haplotype_tag,
        const std::string& somatic_allele_type,
        const std::string& somatic_base_at_variant,
        const std::string& read_id,
        std::vector<msa::MethylationSiteDetail>& details
    );
    
    /**
     * @brief 從BAM讀段中提取單倍型標籤
     * @param read BAM讀段
     * @return std::string 單倍型標籤值
     */
    std::string extractHaplotypeTag(const bam1_t* read);
    
    /**
     * @brief 確定BAM讀段相對於變異的等位基因類型
     * @param read BAM讀段
     * @param target_variant 目標變異
     * @param somatic_base 用於存儲變異位點上的鹼基
     * @return std::string 等位基因類型 (ref/alt/unknown)
     */
    std::string determineAlleleType(
        const bam1_t* read,
        const msa::VcfVariantInfo& target_variant,
        std::string& somatic_base
    );
    
    /**
     * @brief 將甲基化水平分類為高、中、低或未知
     * @param meth_call 甲基化水平 (0.0-1.0)
     * @return std::string 甲基化狀態 (high/mid/low/unknown)
     */
    std::string classifyMethylationState(float meth_call);
    
    /**
     * @brief 將參考基因組位置轉換為讀段上的位置
     * @param read BAM讀段
     * @param ref_pos 參考基因組位置 (0-based)
     * @return int 讀段上的位置 (-1表示無法映射)
     */
    int refPosToReadPos(const bam1_t* read, int ref_pos);
    
    /**
     * @brief 獲取讀段上指定位置的鹼基
     * @param read BAM讀段
     * @param read_pos 讀段上的位置
     * @return char 鹼基 ('N'表示無效位置)
     */
    char getBaseAtReadPos(const bam1_t* read, int read_pos);
};

} // namespace msa::core