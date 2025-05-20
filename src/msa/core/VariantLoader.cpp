#include "msa/core/VariantLoader.h"
#include "msa/core/ConfigParser.h"
#include "msa/utils/LogManager.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

// 使用正確的命名空間
using namespace msa::utils;
namespace fs = std::filesystem;

// 檢查 htslib 版本
#ifndef HTS_VERSION
#define HTSLIB_OLD_API
#endif

#ifdef HTSLIB_OLD_API
// 實現舊版本 htslib 缺少的函數
static inline int bcf_get_alleles_compat(const bcf_hdr_t *hdr, bcf1_t *line, const char ***alleles_dst, int *nals_dst) {
    int nals = line->n_allele;
    if (nals <= 0) return -1;
    
    const char **alleles = (const char**)malloc(nals * sizeof(const char *));
    if (!alleles) return -1;
    
    for (int i = 0; i < nals; i++) {
        alleles[i] = line->d.allele[i];
    }
    
    *alleles_dst = alleles;
    *nals_dst = nals;
    return nals;
}
#endif

// 聲明htslib函數，避免未定義的識別符錯誤
extern "C" {
#ifndef HTSLIB_OLD_API
    int bcf_get_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char ***alleles_dst, int *nals_dst);
#endif
}

namespace msa::core {

/*
* 構造函數
*/
VariantLoader::VariantLoader() : hasBedFile_(false) {}

/*
* 載入VCF檔案
*/
std::vector<msa::VcfVariantInfo> VariantLoader::loadVCFs(
    const std::vector<std::string>& vcfPaths,
    const std::string& bedPath,
    msa::Config& config) {
    
    // 載入BED檔案 (如果提供)
    if (!bedPath.empty()) {
        try {
            loadBedRegions(bedPath);
            hasBedFile_ = true;
            LOG_INFO("VariantLoader", "已載入BED檔案: " + bedPath + 
                               ", 共 " + std::to_string(bedRegions_.size()) + " 個區域");
        } catch (const std::exception& e) {
            LOG_ERROR("VariantLoader", "載入BED檔案失敗: " + std::string(e.what()));
            throw;
        }
    }
    
    std::vector<msa::VcfVariantInfo> variants;
    
    for (const auto& vcfPath : vcfPaths) {
        // 獲取VCF基本名稱作為識別碼
        std::string vcf_source_id = ConfigParser::getBasename(vcfPath);
        LOG_INFO("VariantLoader", "開始處理VCF檔案: " + vcfPath + 
                           " (source_id: " + vcf_source_id + ")");
        
        // 開啟VCF檔案
        htsFile* vcf_fp = bcf_open(vcfPath.c_str(), "r");
        if (!vcf_fp) {
            std::string err = "無法開啟VCF檔案: " + vcfPath;
            LOG_ERROR("VariantLoader", err);
            throw std::runtime_error(err);
        }
        
        // 讀取VCF標頭
        bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
        if (!vcf_hdr) {
            std::string err = "無法讀取VCF標頭: " + vcfPath;
            LOG_ERROR("VariantLoader", err);
            bcf_close(vcf_fp);
            throw std::runtime_error(err);
        }
        
        // 初始化VCF記錄
        bcf1_t* vcf_record = bcf_init();
        if (!vcf_record) {
            std::string err = "無法初始化VCF記錄結構";
            LOG_ERROR("VariantLoader", err);
            bcf_hdr_destroy(vcf_hdr);
            bcf_close(vcf_fp);
            throw std::runtime_error(err);
        }
        
        // 統計變數
        int total_variants = 0;
        int filtered_variants = 0;
        int variants_before = variants.size();
        
        // 迭代讀取VCF記錄
        while (bcf_read(vcf_fp, vcf_hdr, vcf_record) == 0) {
            total_variants++;
            
            // 跳過非PASS變異 (如果有FILTER欄位且不是PASS)
            if (bcf_has_filter(vcf_hdr, vcf_record, const_cast<char*>("PASS")) != 1) {
                filtered_variants++;
                continue;
            }
            
            // 獲取變異的基本資訊
            const char* chrom = bcf_hdr_id2name(vcf_hdr, vcf_record->rid);
            if (!chrom) {
                LOG_WARN("VariantLoader", "無法獲取染色體名稱，跳過變異");
                continue;
            }
            
            // 檢查是否在BED區域內
            if (hasBedFile_ && !isInBedRegions(chrom, vcf_record->pos)) {
                filtered_variants++;
                continue;
            }
            
            // 獲取REF和ALT等位基因
            char** alts = NULL;
            int n_alleles = 0;
            int ret;
                
            #ifdef HTSLIB_OLD_API
                ret = bcf_get_alleles_compat(vcf_hdr, vcf_record, (const char***)&alts, &n_alleles);
            #else
                ret = bcf_get_alleles(vcf_hdr, vcf_record, (const char***)&alts, &n_alleles);
            #endif
            
            if (ret <= 0 || n_alleles <= 1) {
                LOG_WARN("VariantLoader", "變異沒有ALT等位基因，跳過");
                continue;
            }
            
            // 對於每個ALT等位基因，建立一個變異記錄
            for (int i = 1; i < n_alleles; i++) {
                msa::VcfVariantInfo variant;
                variant.vcf_source_id = vcf_source_id;
                variant.chrom = chrom;
                variant.pos = vcf_record->pos + 1; // 轉換為1-based
                variant.ref = alts[0]; // REF等位基因
                variant.alt = alts[i]; // 目前的ALT等位基因
                
                // 判斷變異類型
                variant.variant_type = determineVariantType(vcf_record);
                
                // 檢查最小alt支持數 (如果配置指定)
                if (config.min_allele > 0) {
                    // 嘗試獲取AD標籤
                    int* dp = NULL;
                    int ndp = 0;
                    if (bcf_get_format_int32(vcf_hdr, vcf_record, "AD", &dp, &ndp) > 0) {
                        int alt_support = dp[i]; // 第i個ALT的支持數
                        
                        if (alt_support < config.min_allele) {
                            // 不符合最小支持數要求
                            filtered_variants++;
                            continue;
                        }
                    } else {
                        // 如果無法獲取ALT支持數，記錄警告但仍處理變異
                        LOG_WARN("VariantLoader", "無法獲取變異的ALT支持數，但仍保留此變異");
                    }
                    
                    if (dp) free(dp);
                }
                
                // 檢查變異品質 (如果有)
                if (vcf_record->qual > 0) {
                    variant.qual = vcf_record->qual;
                } else {
                    variant.qual = 0.0f;
                }
                
                variants.push_back(variant);
            }
            
            if (alts) free(alts);
        }
        
        // 清理資源
        bcf_destroy(vcf_record);
        bcf_hdr_destroy(vcf_hdr);
        bcf_close(vcf_fp);
        
        int variants_added = variants.size() - variants_before;
        LOG_INFO("VariantLoader", "已處理VCF檔案: " + vcfPath + 
                           ", 共 " + std::to_string(total_variants) + " 個變異, " + 
                           std::to_string(filtered_variants) + " 個被過濾, " + 
                           std::to_string(variants_added) + " 個保留");
    }
    
    // 對變異進行排序
    std::sort(variants.begin(), variants.end());
    
    return variants;
}

/*
* 確定變異類型
* \param vcf_record VCF記錄
* \return 變異類型
*/
std::string VariantLoader::determineVariantType(const bcf1_t* vcf_record) {
    // 從SVTYPE訊息欄位判斷是否為結構變異
    bcf_hdr_t* header = nullptr; // 這裡修改，header 需要從外部傳入，而不是從 vcf_record 獲取
    
    // 若要使用 SVTYPE，需要先獲取 header，這裡可能需要修改函數參數以接收 header
    // 暫時使用其他方法判斷變異類型
    
    // 檢查 REF 和 ALT 長度判斷變異類型
    int32_t* gt_arr = nullptr;
    int n_gt = 0;
    
    // 獲取 REF 和 ALT
    const char** alleles = nullptr;
    int n_allele = 0;
    
    // 注意：由於沒有 header，我們可能需要依靠 bcf_get_variant_types 函數
    // 或者直接使用 REF 和 ALT 長度比較來判斷變異類型
    
    if (vcf_record->n_allele == 2) {
        // 比較 REF 和 ALT 長度
        int refLen = strlen(vcf_record->d.allele[0]);
        int altLen = strlen(vcf_record->d.allele[1]);
        
        if (refLen == 1 && altLen == 1) {
            return "SNV";  // 單核苷酸變異
        } else if (refLen > altLen) {
            return "DEL";  // 刪除
        } else if (refLen < altLen) {
            return "INS";  // 插入
        } else {
            return "COMPLEX";  // 複雜變異（如多核苷酸替換）
        }
    } else if (vcf_record->n_allele > 2) {
        return "MULTI";  // 多個變異
    }
    
    return "UNKNOWN";  // 未知類型
}

/*
* 載入BED檔案
* \param bedPath BED檔案路徑
*/
void VariantLoader::loadBedRegions(const std::string& bedPath) {
    bedRegions_.clear();
    
    std::ifstream bedFile(bedPath);
    if (!bedFile.is_open()) {
        throw std::runtime_error("無法開啟BED檔案: " + bedPath);
    }
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(bedFile, line)) {
        line_num++;
        
        // 跳過注釋行與空行
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // 移除前導與尾隨空白
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        
        if (line.empty()) {
            continue;
        }
        
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        
        // 讀取染色體、起始位置和結束位置
        if (!(iss >> chrom >> start >> end)) {
            LOG_WARN("VariantLoader", "BED檔案格式錯誤，第 " + 
                               std::to_string(line_num) + " 行: " + line);
            continue;
        }
        
        // 檢查位置有效性
        if (start < 0 || end <= start) {
            continue;
        }
        
        // 創建並添加BED區域
        BedRegion region;
        region.chrom = chrom;
        region.start = start;
        region.end = end;
        
        bedRegions_.push_back(region);
    }
    
    bedFile.close();
}

/*
* 檢查變異是否在BED區域內
* \param chrom 染色體名稱
* \param pos 變異位置
* \return 是否在BED區域內
*/
bool VariantLoader::isInBedRegions(const std::string& chrom, int pos) {
    // 如果沒有BED檔案，則視為所有區域都在範圍內
    if (!hasBedFile_ || bedRegions_.empty()) {
        return true;
    }
    
    // 線性搜索所有BED區域
    // 注意：對於大量BED區域，應使用更高效的搜索結構（如區間樹）
    for (const auto& region : bedRegions_) {
        if (region.chrom == chrom && 
            pos >= region.start && pos < region.end) {
            return true;
        }
    }
    
    return false;
}

} // namespace msa::core