#include "msa/core/BamFetcher.h"
#include "msa/utils/LogManager.h"
#include "msa/utils/MemoryPool.h"
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <functional>

// 使用正確的命名空間
using namespace msa::utils;
namespace fs = std::filesystem;

namespace msa::core {

/*
* 構造函數
*/
BamFetcher::BamFetcher(const msa::Config& config)
    : config_(config), tumor_fp_(nullptr), normal_fp_(nullptr),
      tumor_hdr_(nullptr), normal_hdr_(nullptr),
      tumor_idx_(nullptr), normal_idx_(nullptr) {
}

/*
* 析構函數
*/
BamFetcher::~BamFetcher() {
    closeBamFiles();
}

/*
* 開啟BAM檔案
* \return 是否成功開啟
*/
bool BamFetcher::openBamFiles() {
    // 開啟腫瘤BAM
    if (config_.tumor_bam != "-") {
        tumor_fp_ = sam_open(config_.tumor_bam.c_str(), "r");
        if (!tumor_fp_) {
            LOG_ERROR("BamFetcher", "無法開啟腫瘤BAM檔案: " + config_.tumor_bam);
            return false;
        }
        
        tumor_hdr_ = sam_hdr_read(tumor_fp_);
        if (!tumor_hdr_) {
            LOG_ERROR("BamFetcher", "無法讀取腫瘤BAM標頭: " + config_.tumor_bam);
            closeBamFiles();
            return false;
        }
        
        tumor_idx_ = getBamIndex(config_.tumor_bam);
        if (!tumor_idx_) {
            LOG_ERROR("BamFetcher", "無法載入腫瘤BAM索引: " + config_.tumor_bam);
            closeBamFiles();
            return false;
        }
    } else {
        LOG_WARN("BamFetcher", "腫瘤BAM指定為標準輸入，無法進行區域查詢");
    }
    
    // 開啟正常BAM
    if (config_.normal_bam != "-") {
        normal_fp_ = sam_open(config_.normal_bam.c_str(), "r");
        if (!normal_fp_) {
            LOG_ERROR("BamFetcher", "無法開啟正常BAM檔案: " + config_.normal_bam);
            closeBamFiles();
            return false;
        }
        
        normal_hdr_ = sam_hdr_read(normal_fp_);
        if (!normal_hdr_) {
            LOG_ERROR("BamFetcher", "無法讀取正常BAM標頭: " + config_.normal_bam);
            closeBamFiles();
            return false;
        }
        
        normal_idx_ = getBamIndex(config_.normal_bam);
        if (!normal_idx_) {
            LOG_ERROR("BamFetcher", "無法載入正常BAM索引: " + config_.normal_bam);
            closeBamFiles();
            return false;
        }
    } else {
        LOG_WARN("BamFetcher", "正常BAM指定為標準輸入，無法進行區域查詢");
    }
    
    return true;
}

/*
* 關閉BAM檔案
*/
void BamFetcher::closeBamFiles() {
    // 釋放腫瘤BAM資源
    if (tumor_idx_) {
        hts_idx_destroy(tumor_idx_);
        tumor_idx_ = nullptr;
    }
    
    if (tumor_hdr_) {
        bam_hdr_destroy(tumor_hdr_);
        tumor_hdr_ = nullptr;
    }
    
    if (tumor_fp_) {
        sam_close(tumor_fp_);
        tumor_fp_ = nullptr;
    }
    
    // 釋放正常BAM資源
    if (normal_idx_) {
        hts_idx_destroy(normal_idx_);
        normal_idx_ = nullptr;
    }
    
    if (normal_hdr_) {
        bam_hdr_destroy(normal_hdr_);
        normal_hdr_ = nullptr;
    }
    
    if (normal_fp_) {
        sam_close(normal_fp_);
        normal_fp_ = nullptr;
    }
}

/*
* 獲取變異周圍的讀段
* \param variant 變異資訊
* \param tumorReads 腫瘤讀段容器
* \param normalReads 正常讀段容器
* \return 是否成功獲取
*/
bool BamFetcher::fetchReadsAroundVariant(
    const msa::VcfVariantInfo& variant,
    std::vector<bam1_t*>& tumorReads,
    std::vector<bam1_t*>& normalReads) {
    
    // 清空輸入容器
    tumorReads.clear();
    normalReads.clear();
    
    // 計算查詢區域範圍 (0-based)
    int pos_0based = variant.pos - 1;  // 轉換為0-based位置
    int start = std::max(0, pos_0based - config_.window_size);
    int end = pos_0based + config_.window_size;
    
    LOG_DEBUG("BamFetcher", "使用窗口大小: " + std::to_string(config_.window_size) + 
             ", 區域: " + variant.chrom + ":" + std::to_string(start+1) + "-" + std::to_string(end+1));
    
    // 獲取腫瘤樣本的讀段
    if (tumor_fp_ && tumor_hdr_ && tumor_idx_) {
        if (!fetchReadsFromRegion(tumor_fp_, tumor_hdr_, tumor_idx_, 
                                variant.chrom, start, end, tumorReads)) {
            LOG_ERROR("BamFetcher", "無法取得腫瘤樣本區域: " + 
                              variant.chrom + ":" + std::to_string(start+1) + "-" + std::to_string(end+1));
            return false;
        }
    }
    
    // 獲取正常樣本的讀段
    if (normal_fp_ && normal_hdr_ && normal_idx_) {
        if (!fetchReadsFromRegion(normal_fp_, normal_hdr_, normal_idx_, 
                                variant.chrom, start, end, normalReads)) {
            LOG_ERROR("BamFetcher", "無法取得正常樣本區域: " + 
                              variant.chrom + ":" + std::to_string(start+1) + "-" + std::to_string(end+1));
            return false;
        }
    }
    
    // 記錄獲取的讀段數量
    std::ostringstream ss;
    ss << "變異 " << variant.chrom << ":" << variant.pos << " " << variant.ref << ">" << variant.alt 
       << " 區域取得: 腫瘤讀段=" << tumorReads.size() 
       << ", 正常讀段=" << normalReads.size();
    LOG_DEBUG("BamFetcher", ss.str());
    
    return true;
}

/*
* 從指定區域獲取讀段
* \param fp BAM檔案指針
* \param hdr BAM標頭指針
* \param idx BAM索引指針
* \param chrom 染色體名稱
* \param start 起始位置
* \param end 結束位置
* \param reads 讀段容器
* \return 是否成功獲取
*/
bool BamFetcher::fetchReadsFromRegion(
    htsFile* fp,
    bam_hdr_t* hdr,
    hts_idx_t* idx,
    const std::string& chrom,
    int start,
    int end,
    std::vector<bam1_t*>& reads) {
    
    // 獲取染色體ID
    int tid = bam_name2id(hdr, chrom.c_str());
    if (tid < 0) {
        LOG_ERROR("BamFetcher", "無法找到染色體: " + chrom);
        return false;
    }
    
    // 建立查詢迭代器
    hts_itr_t* iter = sam_itr_queryi(idx, tid, start, end);
    if (!iter) {
        LOG_ERROR("BamFetcher", "無法建立迭代器: " + chrom + ":" + 
                          std::to_string(start+1) + "-" + std::to_string(end+1));
        return false;
    }
    
    // 獲取記憶體池實例
    auto& memPool = MemoryPool::getInstance();
    
    // 讀取區域內的所有讀段
    bool max_depth_reached = false;
    int num_reads = 0;
    
    while (true) {
        // 檢查是否達到最大讀取深度
        if (num_reads >= config_.max_read_depth) {
            max_depth_reached = true;
            break;
        }
        
        // 從記憶體池獲取bam1_t物件
        bam1_t* b = memPool.getBam1();
        
        // 讀取下一條記錄
        int ret = sam_itr_next(fp, iter, b);
        if (ret < 0) {
            memPool.returnBam1(b);  // 返回記憶體池
            break;  // 讀取完畢或錯誤
        }
        
        // 檢查讀段是否有效
        if (isReadValid(b)) {
            reads.push_back(b);
            num_reads++;
        } else {
            memPool.returnBam1(b);  // 無效讀段返回記憶體池
        }
    }
    
    if (max_depth_reached) {
        LOG_WARN("BamFetcher", "區域 " + chrom + ":" + 
                          std::to_string(start+1) + "-" + std::to_string(end+1) + 
                          " 已達最大讀取深度: " + std::to_string(config_.max_read_depth));
    }
    
    // 釋放迭代器
    hts_itr_destroy(iter);
    
    return true;
}

/*
* 獲取BAM索引
* \param bamPath BAM檔案路徑
* \return BAM索引指針
*/
hts_idx_t* BamFetcher::getBamIndex(const std::string& bamPath) {
    hts_idx_t* idx = nullptr;
    
    // 如果是標準輸入，無法獲取索引
    if (bamPath != "-") {
        // 嘗試尋找 .bai 與 .bam.bai 兩種索引檔案
        std::string baiPath = bamPath + ".bai";
        std::string bamBaiPath = fs::path(bamPath).replace_extension(".bam.bai").string();
        
        idx = hts_idx_load(bamPath.c_str(), HTS_FMT_BAI);
        
        if (!idx) {
            LOG_ERROR("BamFetcher", "無法載入BAM索引: " + bamPath + 
                              " (嘗試了 " + baiPath + " 和 " + bamBaiPath + ")");
        }
    }
    
    return idx;
}

/*
* 檢查讀段是否有效
* \param b 讀段指針
* \return 是否有效
*/
bool BamFetcher::isReadValid(const bam1_t* b) {
    // 檢查是否為主要比對 (非輔助/次要比對)
    if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        return false;
    }
    
    // 檢查是否為有效比對 (非未比對)
    if (b->core.flag & BAM_FUNMAP) {
        return false;
    }
    
    // 檢查MAPQ (比對品質)
    if (b->core.qual < 10) {  // 基本MAPQ閾值
        return false;
    }
    
    return true;
}

/*
* 獲取變異周圍的讀段
* \param variant 變異資訊
* \param isTumor 是否為腫瘤樣本
* \param window_size 窗口大小
* \return 讀段容器
*/
std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> BamFetcher::fetchReadsAroundVariant(
    const msa::VcfVariantInfo& variant,
    bool isTumor,
    int window_size) {
    
    // 使用配置中的窗口大小（如果未提供）
    if (window_size <= 0) {
        window_size = config_.window_size;
    }
    
    // 計算查詢區域範圍 (0-based)
    int pos_0based = variant.pos - 1;  // 轉換為0-based位置
    int start = std::max(0, pos_0based - window_size);
    int end = pos_0based + window_size;
    
    // 獲取相應的BAM檔案資源
    htsFile* fp = isTumor ? tumor_fp_ : normal_fp_;
    bam_hdr_t* hdr = isTumor ? tumor_hdr_ : normal_hdr_;
    hts_idx_t* idx = isTumor ? tumor_idx_ : normal_idx_;
    
    // 創建結果容器
    std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> result;
    
    // 檢查資源是否有效
    if (!fp || !hdr || !idx) {
        LOG_ERROR("BamFetcher", "無效的BAM資源: " + std::string(isTumor ? "腫瘤" : "正常") + " 樣本");
        return result;
    }
    
    // 獲取染色體ID
    int tid = bam_name2id(hdr, variant.chrom.c_str());
    if (tid < 0) {
        LOG_ERROR("BamFetcher", "無法找到染色體: " + variant.chrom);
        return result;
    }
    
    // 建立查詢迭代器
    hts_itr_t* iter = sam_itr_queryi(idx, tid, start, end);
    if (!iter) {
        LOG_ERROR("BamFetcher", "無法建立迭代器: " + variant.chrom + ":" + 
                          std::to_string(start+1) + "-" + std::to_string(end+1));
        return result;
    }
    
    // 獲取記憶體池實例
    auto& memPool = MemoryPool::getInstance();
    
    // 定義自定義刪除器（用於智能指針）
    auto deleter = [&memPool](bam1_t* b) {
        memPool.returnBam1(b);
    };
    
    // 讀取區域內的所有讀段
    bool max_depth_reached = false;
    int num_reads = 0;
    
    while (true) {
        // 檢查是否達到最大讀取深度
        if (num_reads >= config_.max_read_depth) {
            max_depth_reached = true;
            break;
        }
        
        // 從記憶體池獲取bam1_t物件
        bam1_t* b = memPool.getBam1();
        
        // 讀取下一條記錄
        int ret = sam_itr_next(fp, iter, b);
        if (ret < 0) {
            memPool.returnBam1(b);  // 返回記憶體池
            break;  // 讀取完畢或錯誤
        }
        
        // 檢查讀段是否有效
        if (isReadValid(b)) {
            // 創建智能指針並添加到結果
            result.push_back(std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>(b, deleter));
            num_reads++;
        } else {
            memPool.returnBam1(b);  // 無效讀段返回記憶體池
        }
    }
    
    if (max_depth_reached) {
        LOG_WARN("BamFetcher", "區域 " + variant.chrom + ":" + 
                          std::to_string(start+1) + "-" + std::to_string(end+1) + 
                          " 已達最大讀取深度: " + std::to_string(config_.max_read_depth));
    }
    
    // 釋放迭代器
    hts_itr_destroy(iter);
    
    // 記錄獲取的讀段數量
    std::ostringstream ss;
    ss << "變異 " << variant.chrom << ":" << variant.pos << " " << variant.ref << ">" << variant.alt 
       << " 區域取得: " << (isTumor ? "腫瘤" : "正常") << "讀段=" << result.size();
    LOG_DEBUG("BamFetcher", ss.str());
    
    return result;
}

} // namespace msa::core 