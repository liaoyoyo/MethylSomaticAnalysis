#pragma once

#include <string>
#include <htslib/sam.h>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief BAM檔案驗證器，檢查是否包含所需的標籤
 */
class BAMValidator {
public:
    /**
     * @brief 建構函數
     */
    BAMValidator();
    
    /**
     * @brief 驗證所有BAM檔案
     * @param config 配置物件，包含所需驗證的BAM檔案路徑
     * @return bool 是否成功 (檢查是否有致命性問題)
     */
    bool checkAllInputFiles(msa::Config& config);
    
    /**
     * @brief 驗證BAM檔案的甲基化標籤
     * @param bamPath BAM檔案路徑
     * @param fp BAM檔案指標
     * @param hdr BAM標頭
     * @param config 配置物件，用於記錄驗證結果
     * @return bool 是否找到甲基化標籤
     */
    bool checkMethylationTags(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, msa::Config& config);
    
    /**
     * @brief 驗證BAM檔案的單倍型標籤
     * @param bamPath BAM檔案路徑
     * @param fp BAM檔案指標
     * @param hdr BAM標頭
     * @param config 配置物件，用於記錄驗證結果
     * @return bool 是否找到單倍型標籤
     */
    bool checkHaplotypeTags(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, msa::Config& config);

private:
    /**
     * @brief 取樣檢查BAM讀段是否包含指定標籤
     * @param bamPath BAM檔案路徑
     * @param fp BAM檔案指標
     * @param hdr BAM標頭
     * @param tagName 標籤名稱
     * @param sampleSize 取樣大小
     * @param isTumorBam 是否為腫瘤樣本BAM
     * @return 找到標籤的讀段數量
     */
    int sampleCheckTag(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, 
                      const char* tagName, int sampleSize, bool isTumorBam);
};

} // namespace msa::core 