#include "msa/core/BAMValidator.h"
#include "msa/utils/LogManager.h"
#include "msa/utils/MemoryPool.h"
#include <sstream>
#include <stdexcept>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

// 使用正確的命名空間
using namespace msa::utils;

namespace msa::core {

BAMValidator::BAMValidator() = default;

bool BAMValidator::checkAllInputFiles(msa::Config& config) {
    bool success = true;
    
    // 檢查腫瘤BAM
    if (config.tumor_bam != "-") {
        htsFile* tumor_fp = sam_open(config.tumor_bam.c_str(), "r");
        if (!tumor_fp) {
            LOG_ERROR("BAMValidator", "無法開啟腫瘤BAM檔案: " + config.tumor_bam);
            return false;
        }
        
        bam_hdr_t* tumor_hdr = sam_hdr_read(tumor_fp);
        if (!tumor_hdr) {
            LOG_ERROR("BAMValidator", "無法讀取腫瘤BAM標頭: " + config.tumor_bam);
            hts_close(tumor_fp);
            return false;
        }
        
        bool tumor_has_meth = checkMethylationTags(config.tumor_bam, tumor_fp, tumor_hdr, config);
        if (!tumor_has_meth) {
            LOG_WARN("BAMValidator", "腫瘤BAM未檢測到甲基化標籤 (MM/ML)，這可能影響甲基化分析結果");
            // 非致命性警告，不要返回false
        }
        
        bool tumor_has_hp = checkHaplotypeTags(config.tumor_bam, tumor_fp, tumor_hdr, config);
        if (!tumor_has_hp) {
            LOG_WARN("BAMValidator", "腫瘤BAM未檢測到單倍型標籤 (HP/PS)，這可能影響單倍型分析結果");
            // 非致命性警告，不要返回false
        }
        
        bam_hdr_destroy(tumor_hdr);
        hts_close(tumor_fp);
    } else {
        LOG_WARN("BAMValidator", "腫瘤BAM指定為標準輸入，跳過標籤檢查");
    }
    
    // 檢查正常BAM
    if (config.normal_bam != "-") {
        htsFile* normal_fp = sam_open(config.normal_bam.c_str(), "r");
        if (!normal_fp) {
            LOG_ERROR("BAMValidator", "無法開啟正常BAM檔案: " + config.normal_bam);
            return false;
        }
        
        bam_hdr_t* normal_hdr = sam_hdr_read(normal_fp);
        if (!normal_hdr) {
            LOG_ERROR("BAMValidator", "無法讀取正常BAM標頭: " + config.normal_bam);
            hts_close(normal_fp);
            return false;
        }
        
        bool normal_has_meth = checkMethylationTags(config.normal_bam, normal_fp, normal_hdr, config);
        if (!normal_has_meth) {
            LOG_WARN("BAMValidator", "正常BAM未檢測到甲基化標籤 (MM/ML)，這可能影響甲基化分析結果");
            // 非致命性警告，不要返回false
        }
        
        bool normal_has_hp = checkHaplotypeTags(config.normal_bam, normal_fp, normal_hdr, config);
        if (!normal_has_hp) {
            LOG_WARN("BAMValidator", "正常BAM未檢測到單倍型標籤 (HP/PS)，這可能影響單倍型分析結果");
            // 非致命性警告，不要返回false
        }
        
        bam_hdr_destroy(normal_hdr);
        hts_close(normal_fp);
    } else {
        LOG_WARN("BAMValidator", "正常BAM指定為標準輸入，跳過標籤檢查");
    }
    
    return success;
}

bool BAMValidator::checkMethylationTags(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, msa::Config& config) {
       // 先重置文件指針
       if (bgzf_seek(fp->fp.bgzf, 0, SEEK_SET) < 0) {
           LOG_ERROR("BAMValidator", "無法重新定位至BAM檔案開頭: " + bamPath);
           return false;
       }
       
       // 重新讀取標頭以確保標頭指針也被重置
       bam_hdr_t* new_hdr = sam_hdr_read(fp);
       if (!new_hdr) {
           LOG_ERROR("BAMValidator", "無法重新讀取BAM標頭: " + bamPath);
           return false;
       }
       
       bam1_t* read = bam_init1();
       int found_mm = 0, found_ml = 0;
       int count = 0;
       
       while (sam_read1(fp, new_hdr, read) >= 0 && count < 100) {
           count++;
           
           // 檢查標籤
           uint8_t* mm_tag = bam_aux_get(read, "MM");
           uint8_t* ml_tag = bam_aux_get(read, "ML");
           
           if (mm_tag != nullptr) found_mm++;
           if (ml_tag != nullptr) found_ml++;
       }
       
       bam_destroy1(read);
       bam_hdr_destroy(new_hdr);  // 釋放新的標頭
       
       LOG_INFO("BAMValidator", "BAM檔案 " + bamPath + " 檢查結果: MM標籤發現於 " + 
                std::to_string(found_mm) + " 條讀段, ML標籤發現於 " + 
                std::to_string(found_ml) + " 條讀段 (取樣" + std::to_string(count) + "條)");
       
       bool has_methyl_tags = (found_mm > 0 && found_ml > 0);
       
       // 更新配置
       if (bamPath == config.tumor_bam) {
           config.tumor_has_methyl_tags = has_methyl_tags;
       } else if (bamPath == config.normal_bam) {
           config.normal_has_methyl_tags = has_methyl_tags;
       }
       
       return has_methyl_tags;
   }

bool BAMValidator::checkHaplotypeTags(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, msa::Config& config) {
    int hp_count = sampleCheckTag(bamPath, fp, hdr, "HP", 100, bamPath == config.tumor_bam);
    int ps_count = sampleCheckTag(bamPath, fp, hdr, "PS", 100, bamPath == config.tumor_bam);
    
    std::ostringstream ss;
    ss << "BAM檔案 " << bamPath << " 檢查結果: HP標籤發現於 " << hp_count 
       << " 條讀段, PS標籤發現於 " << ps_count << " 條讀段 (取樣100條)";
    LOG_INFO("BAMValidator", ss.str());
    
    bool has_hp_tags = (hp_count > 0);
    
    // 更新配置中的標籤狀態
    if (bamPath == config.tumor_bam) {
        config.tumor_has_hp_tags = has_hp_tags;
    } else if (bamPath == config.normal_bam) {
        config.normal_has_hp_tags = has_hp_tags;
    }
    
    return has_hp_tags;
}

int BAMValidator::sampleCheckTag(const std::string& bamPath, htsFile* fp, bam_hdr_t* hdr, 
                               const char* tagName, int sampleSize, bool isTumorBam) {
    // 重新定位到檔案開頭
    if (bgzf_seek(fp->fp.bgzf, 0, SEEK_SET) < 0) {
        std::ostringstream ss;
        ss << "無法重新定位至BAM檔案開頭: " << bamPath;
        LOG_ERROR("BAMValidator", ss.str());
        return 0;
    }
    
    // 重新讀取標頭以確保標頭指針也被重置
    bam_hdr_t* new_hdr = sam_hdr_read(fp);
    if (!new_hdr) {
        std::ostringstream ss;
        ss << "無法重新讀取BAM標頭: " << bamPath;
        LOG_ERROR("BAMValidator", ss.str());
        return 0;
    }
    
    bam1_t* read = bam_init1();
    int count = 0;
    int found = 0;
    
    // 讀取指定數量的記錄並檢查標籤
    while (sam_read1(fp, new_hdr, read) >= 0 && count < sampleSize) {
        // 隨機取樣
        if (rand() % 3 != 0 && count > 10) continue;
        
        count++;
        
        // 檢查是否有指定標籤
        uint8_t* tag_data = bam_aux_get(read, tagName);
        if (tag_data != nullptr) {
            found++;
        }
    }
    
    bam_destroy1(read);
    bam_hdr_destroy(new_hdr);  // 釋放新的標頭
    
    return found;
}

} // namespace msa::core 