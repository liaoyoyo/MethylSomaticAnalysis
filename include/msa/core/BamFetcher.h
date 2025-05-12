#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <htslib/sam.h>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief BAM區域抓取器，用於從BAM檔案中取出指定區域的讀段
 */
class BamFetcher {
public:
    /**
     * @brief 建構函數
     * @param config 配置物件
     */
    BamFetcher(const msa::Config& config);
    
    /**
     * @brief 解構函數
     */
    ~BamFetcher();
    
    /**
     * @brief 打開腫瘤和正常BAM檔案
     * @return bool 是否成功打開
     */
    bool openBamFiles();
    
    /**
     * @brief 關閉BAM檔案
     */
    void closeBamFiles();
    
    /**
     * @brief 獲取變異周圍區域的讀段
     * @param variant 變異資訊
     * @param tumorReads 用於存儲腫瘤樣本的讀段
     * @param normalReads 用於存儲正常樣本的讀段
     * @return bool 是否成功獲取
     */
    bool fetchReadsAroundVariant(
        const msa::VcfVariantInfo& variant,
        std::vector<bam1_t*>& tumorReads,
        std::vector<bam1_t*>& normalReads
    );
    
    /**
     * @brief 獲取變異周圍區域的讀段（單樣本，用於 main.cpp）
     * @param variant 變異資訊
     * @param isTumor 是否為腫瘤樣本
     * @param window_size 搜索窗口大小，若為0則使用配置中的值
     * @return std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> 智能指針管理的讀段
     */
    std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> fetchReadsAroundVariant(
        const msa::VcfVariantInfo& variant,
        bool isTumor,
        int window_size = 0
    );
    
private:
    const msa::Config& config_;   // 配置物件
    
    // BAM檔案資源
    htsFile* tumor_fp_ = nullptr;
    bam_hdr_t* tumor_hdr_ = nullptr;
    hts_idx_t* tumor_idx_ = nullptr;
    
    htsFile* normal_fp_ = nullptr;
    bam_hdr_t* normal_hdr_ = nullptr;
    hts_idx_t* normal_idx_ = nullptr;
    
    /**
     * @brief 從BAM檔案抓取指定區域的讀段
     * @param fp BAM檔案指標
     * @param hdr BAM標頭
     * @param idx BAM索引
     * @param chrom 染色體名稱
     * @param start 起始位置 (0-based)
     * @param end 結束位置 (0-based)
     * @param reads 用於存儲抓取的讀段
     * @return bool 是否成功抓取
     */
    bool fetchReadsFromRegion(
        htsFile* fp,
        bam_hdr_t* hdr,
        hts_idx_t* idx,
        const std::string& chrom,
        int start,
        int end,
        std::vector<bam1_t*>& reads
    );
    
    /**
     * @brief 獲取BAM索引
     * @param bamPath BAM檔案路徑
     * @return hts_idx_t* BAM索引指標
     */
    hts_idx_t* getBamIndex(const std::string& bamPath);
    
    /**
     * @brief 檢查讀段是否符合條件（品質過濾等）
     * @param b BAM讀段
     * @return bool 是否符合條件
     */
    bool isReadValid(const bam1_t* b);
};

} // namespace msa::core 