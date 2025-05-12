#pragma once

#include <string>
#include <vector>
#include <htslib/vcf.h>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief 變異加載器，從VCF檔案加載變異資訊
 */
class VariantLoader {
public:
    /**
     * @brief 建構函數
     */
    VariantLoader();
    
    /**
     * @brief 從VCF檔案加載變異資訊
     * @param vcfPaths VCF檔案路徑列表
     * @param bedPath BED檔案路徑，用於過濾變異（可為空）
     * @param config 配置物件
     * @return std::vector<msa::VcfVariantInfo> 變異資訊列表
     */
    std::vector<msa::VcfVariantInfo> loadVCFs(
        const std::vector<std::string>& vcfPaths,
        const std::string& bedPath,
        msa::Config& config
    );
    
private:
    /**
     * @brief 確定變異類型
     * @param vcf_record VCF記錄
     * @return std::string 變異類型（SNV, INDEL等）
     */
    std::string determineVariantType(const bcf1_t* vcf_record);
    
    /**
     * @brief 加載BED檔案以用於變異過濾
     * @param bedPath BED檔案路徑
     */
    void loadBedRegions(const std::string& bedPath);
    
    /**
     * @brief 檢查變異是否在BED區域內
     * @param chrom 染色體
     * @param pos 變異位置
     * @return bool 是否在BED區域內
     */
    bool isInBedRegions(const std::string& chrom, int pos);
    
    // BED區域存儲 
    // 簡單實現，使用<chrom, <start, end>>對的向量
    // 注意：實際實現可能需要更高效的數據結構，如區間樹
    struct BedRegion {
        std::string chrom;
        int start;  // 0-based
        int end;    // 0-based, exclusive
    };
    std::vector<BedRegion> bedRegions_;
    bool hasBedFile_ = false;
};

} // namespace msa::core