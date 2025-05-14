#include "msa/core/MethylHaploExtractor.h"
#include "msa/utils/LogManager.h"
#include <sstream>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <cmath>

// 檢查 htslib 版本
#ifndef HTS_VERSION
#define HTSLIB_OLD_API
#endif

#ifdef HTSLIB_OLD_API
// 適配舊版 htslib，實現缺少的 bam_aux2array 函數
static inline void* bam_aux2array_compat(const uint8_t *aux, int *count) {
    if (aux[0] != 'B') return NULL;
    const char subtype = aux[1];
    if (subtype != 'c' && subtype != 'C' && subtype != 's' && subtype != 'S' && subtype != 'i' && subtype != 'I' && subtype != 'f') 
        return NULL;
    
    // 解析數組大小
    uint32_t nmemb = 0;
    memcpy(&nmemb, aux + 2, 4);
    *count = nmemb;
    
    // 分配記憶體並複製數據
    size_t elem_size;
    switch (subtype) {
        case 'c': case 'C': elem_size = 1; break;
        case 's': case 'S': elem_size = 2; break;
        case 'i': case 'I': case 'f': elem_size = 4; break;
        default: return NULL;
    }
    
    void* array = malloc(nmemb * elem_size);
    if (!array) return NULL;
    
    memcpy(array, aux + 6, nmemb * elem_size);
    return array;
}
#endif

// 使用正確的命名空間
using namespace msa::utils;

// 前置聲明htslib函數
extern "C" {
#ifndef HTSLIB_OLD_API
    void* bam_aux2array(const uint8_t *aux, int *count);
#endif
}

namespace msa::core {

MethylHaploExtractor::MethylHaploExtractor(const msa::Config& config)
    : config_(config) {
}

std::vector<msa::MethylationSiteDetail> MethylHaploExtractor::extractFromRead(
    const bam1_t* read,
    const msa::VcfVariantInfo& target_variant,
    const std::string& bam_source_id
) {
    std::vector<msa::MethylationSiteDetail> details;
    
    // 獲取讀段ID (QNAME)
    char* qname = bam_get_qname(read);
    std::string read_id(qname);
    
    // 獲取單倍型標籤
    std::string haplotype_tag = extractHaplotypeTag(read);
    
    // 確定對變異的支持 (ref/alt)
    std::string somatic_base;
    std::string somatic_allele_type = determineAlleleType(read, target_variant, somatic_base);
    
    // 提取甲基化信息
    extractMethylation(read, target_variant, bam_source_id, haplotype_tag, somatic_allele_type, details);
    
    // 為所有甲基化位點添加讀段ID
    for (auto& detail : details) {
        detail.read_id = read_id;
    }
    
    return details;
}

void MethylHaploExtractor::extractMethylation(
    const bam1_t* read,
    const msa::VcfVariantInfo& target_variant,
    const std::string& bam_source_id,
    const std::string& haplotype_tag,
    const std::string& somatic_allele_type,
    std::vector<msa::MethylationSiteDetail>& details
) {
    // 檢查MM和ML標籤是否存在
    uint8_t* mm_data = bam_aux_get(read, "MM");
    uint8_t* ml_data = bam_aux_get(read, "ML");
    
    // 如果MM或ML標籤不存在，則嘗試Mm和Ml標籤
    if (!mm_data || !ml_data) {
        mm_data = bam_aux_get(read, "Mm");
        ml_data = bam_aux_get(read, "Ml");
    }
    
    if (!mm_data || !ml_data) {
        // 沒有找到甲基化標籤，跳過此讀段
        return;
    }
    
    // 獲取MM標籤值
    const char* mm_str = bam_aux2Z(mm_data);
    if (!mm_str) {
        return;
    }
    
    // 獲取ML標籤值 (轉為概率數組)
    uint8_t* ml_array = NULL;
    int ml_len = 0;
    
    if (ml_data) {
        // 根據 htslib 版本選擇正確的函數
        #ifdef HTSLIB_OLD_API
            ml_array = (uint8_t*)bam_aux2array_compat(ml_data, &ml_len);
        #else
            ml_array = (uint8_t*)bam_aux2array(ml_data, &ml_len);
        #endif
    }
    
    if (!ml_array || ml_len == 0) {
        return;
    }
    
    // 獲取讀段的鏈方向
    bool is_reverse = bam_is_rev(read);
    char strand = is_reverse ? '-' : '+';
    
    // 解析MM字串
    // 格式: C+m,10,5,3;C+h,2,3;
    std::string mm_string(mm_str);
    std::istringstream mm_stream(mm_string);
    std::string token;
    
    // 分號分隔不同修飾類型
    while (std::getline(mm_stream, token, ';')) {
        if (token.empty()) continue;
        
        // 解析修飾類型
        size_t first_comma = token.find(',');
        if (first_comma == std::string::npos) continue;
        
        std::string mod_type = token.substr(0, first_comma);
        
        // 我們只關心C鹼基上的甲基化
        if (mod_type.size() < 3 || mod_type[0] != 'C') continue;
        
        // 檢查是否是C+m (5mC) 或 C+h (5hmC)
        bool is_valid_mod = (mod_type.substr(1, 2) == "+m" || mod_type.substr(1, 2) == "+h");
        if (!is_valid_mod) continue;
        
        std::string mod_pos_str = token.substr(first_comma + 1);
        std::istringstream mod_pos_stream(mod_pos_str);
        std::string pos_token;
        
        // 獲取第一個修飾位置
        std::getline(mod_pos_stream, pos_token, ',');
        int initial_offset = std::stoi(pos_token);
        
        // 計算參考基因組上的修飾位置
        uint32_t* cigar = bam_get_cigar(read);
        int n_cigar = read->core.n_cigar;
        
        // 起始參考位置 (0-based)
        int ref_pos = read->core.pos;
        int read_pos = 0;
        int current_c_count = 0;
        
        // 變異位置 (0-based)
        int variant_pos_0based = target_variant.pos - 1;
        
        // 遍歷CIGAR操作，找到所有C鹼基位置
        for (int i = 0; i < n_cigar; ++i) {
            uint32_t cigar_op = cigar[i];
            int op_len = bam_cigar_oplen(cigar_op);
            char op_char = bam_cigar_opchr(cigar_op);
            
            if (op_char == 'M' || op_char == '=' || op_char == 'X') {
                // 匹配區域，同時移動參考和讀段位置
                for (int j = 0; j < op_len; ++j) {
                    // 獲取讀段上當前位置的鹼基
                    char base = getBaseAtReadPos(read, read_pos + j);
                    
                    // 檢查是否為C鹼基
                    if (toupper(base) == 'C') {
                        current_c_count++;
                        
                        // 檢查是否是我們尋找的甲基化位置
                        if (current_c_count == initial_offset) {
                            // 找到了第一個甲基化位置
                            int methyl_ref_pos = ref_pos + j;
                            
                            // 檢查甲基化位置是否在變異位置的window_size範圍內
                            if (std::abs(methyl_ref_pos - variant_pos_0based) > config_.window_size) {
                                break;
                            }
                            
                            // 獲取ML值 (0-255的甲基化可能性)
                            float meth_prob = 0.0;
                            if (ml_array && ml_len > 0) {
                                meth_prob = ml_array[0] / 255.0; // 轉換為0-1範圍
                            }
                            
                            // 分類甲基化狀態
                            std::string meth_state = classifyMethylationState(meth_prob);
                            
                            // 建立甲基化位點詳情
                            msa::MethylationSiteDetail detail;
                            detail.chrom = target_variant.chrom;
                            detail.methyl_pos = methyl_ref_pos + 1; // 轉換為1-based
                            detail.somatic_pos = target_variant.pos;
                            detail.variant_type = target_variant.variant_type;
                            detail.vcf_source_id = target_variant.vcf_source_id;
                            detail.bam_source_id = bam_source_id;
                            detail.somatic_allele_type = somatic_allele_type;
                            detail.haplotype_tag = haplotype_tag;
                            detail.meth_call = meth_prob;
                            detail.meth_state = meth_state;
                            detail.strand = strand;
                            detail.somatic_base_at_variant = '?'; // 將在後面設置
                            
                            details.push_back(detail);
                            
                            // 處理剩餘的甲基化位置
                            int ml_index = 1;
                            while (std::getline(mod_pos_stream, pos_token, ',')) {
                                int next_offset = std::stoi(pos_token);
                                current_c_count += next_offset;
                                
                                // 計算參考基因組上的位置
                                int next_methyl_ref_pos = -1;
                                int next_read_pos = -1;
                                
                                // 使用refPosToReadPos的反向邏輯，計算參考位置
                                int c_count = 0;
                                for (int ci = 0; ci < n_cigar; ++ci) {
                                    uint32_t op = cigar[ci];
                                    int op_length = bam_cigar_oplen(op);
                                    char op_type = bam_cigar_opchr(op);
                                    
                                    if (op_type == 'M' || op_type == '=' || op_type == 'X') {
                                        for (int k = 0; k < op_length; ++k) {
                                            char next_base = getBaseAtReadPos(read, next_read_pos + k);
                                            if (toupper(next_base) == 'C') {
                                                c_count++;
                                                if (c_count == current_c_count) {
                                                    next_methyl_ref_pos = read->core.pos + k;
                                                    break;
                                                }
                                            }
                                        }
                                        if (next_methyl_ref_pos != -1) break;
                                        next_read_pos += op_length;
                                    } else if (op_type == 'I' || op_type == 'S') {
                                        // 插入或軟剪切，僅移動讀段位置
                                        next_read_pos += op_length;
                                    } else if (op_type == 'D' || op_type == 'N') {
                                        // 刪除或跳過，僅移動參考位置
                                    }
                                }
                                
                                if (next_methyl_ref_pos != -1) {
                                    // 檢查甲基化位置是否在變異位置的window_size範圍內
                                    if (std::abs(next_methyl_ref_pos - variant_pos_0based) > config_.window_size) {
                                        continue;
                                    }
                                    
                                    // 獲取ML值
                                    meth_prob = 0.0;
                                    if (ml_array && ml_len > ml_index) {
                                        meth_prob = ml_array[ml_index] / 255.0;
                                        ml_index++;
                                    }
                                    
                                    // 分類甲基化狀態
                                    meth_state = classifyMethylationState(meth_prob);
                                    
                                    // 建立甲基化位點詳情
                                    msa::MethylationSiteDetail next_detail;
                                    next_detail.chrom = target_variant.chrom;
                                    next_detail.methyl_pos = next_methyl_ref_pos + 1; // 轉換為1-based
                                    next_detail.somatic_pos = target_variant.pos;
                                    next_detail.variant_type = target_variant.variant_type;
                                    next_detail.vcf_source_id = target_variant.vcf_source_id;
                                    next_detail.bam_source_id = bam_source_id;
                                    next_detail.somatic_allele_type = somatic_allele_type;
                                    next_detail.haplotype_tag = haplotype_tag;
                                    next_detail.meth_call = meth_prob;
                                    next_detail.meth_state = meth_state;
                                    next_detail.strand = strand;
                                    next_detail.somatic_base_at_variant = '?'; // 將在後面設置
                                    
                                    details.push_back(next_detail);
                                }
                            }
                            
                            // 找到所有甲基化位置後，跳出循環
                            break;
                        }
                    }
                }
                
                ref_pos += op_len;
                read_pos += op_len;
            } else if (op_char == 'I' || op_char == 'S') {
                // 插入或軟剪切，僅移動讀段位置
                read_pos += op_len;
            } else if (op_char == 'D' || op_char == 'N') {
                // 刪除或跳過，僅移動參考位置
                ref_pos += op_len;
            }
        }
    }
    
    // 設置變異位點上的鹼基
    std::string somatic_base;
    determineAlleleType(read, target_variant, somatic_base);
    
    for (auto& detail : details) {
        detail.somatic_base_at_variant = somatic_base;
    }
}

std::string MethylHaploExtractor::extractHaplotypeTag(const bam1_t* read) {
    // 獲取HP標籤
    uint8_t* hp_data = bam_aux_get(read, "HP");
    if (hp_data) {
        int hp_value = bam_aux2i(hp_data);
        
        // 獲取PS標籤
        uint8_t* ps_data = bam_aux_get(read, "PS");
        if (ps_data) {
            int ps_value = bam_aux2i(ps_data);
            
            // 檢查是否為alt支持讀段且為HP=1或HP=2的突變單倍型
            if (hp_value == 1 || hp_value == 2) {
                return std::to_string(hp_value);
            }
        }
        
        // 如果沒有PS標籤但有HP標籤，仍然返回HP值
        return std::to_string(hp_value);
    }
    
    // 沒有HP標籤，返回"0"表示未知或其他
    return "0";
}

std::string MethylHaploExtractor::determineAlleleType(
    const bam1_t* read,
    const msa::VcfVariantInfo& target_variant,
    std::string& somatic_base
) {
    // 獲取變異位點的參考基因組位置 (0-based)
    int var_pos_0based = target_variant.pos - 1;
    
    // 將參考位置轉換為讀段位置
    int read_pos = refPosToReadPos(read, var_pos_0based);
    
    // 無法映射到讀段上
    if (read_pos < 0) {
        somatic_base = "?";
        return "unknown";
    }
    
    // 獲取讀段上變異位點的鹼基
    char base = getBaseAtReadPos(read, read_pos);
    somatic_base = std::string(1, base);
    
    // 比較讀段上的鹼基與參考和替代等位基因
    if (base == target_variant.ref[0]) {
        return "ref";
    } else if (base == target_variant.alt[0]) {
        return "alt";
    } else {
        return "unknown";
    }
}

std::string MethylHaploExtractor::classifyMethylationState(float meth_call) {
    if (meth_call >= config_.meth_high_threshold) {
        return "high";
    } else if (meth_call >= config_.meth_low_threshold) {
        return "mid";
    } else if (meth_call >= 0.0) {
        return "low";
    } else {
        return "unknown";
    }
}

int MethylHaploExtractor::refPosToReadPos(const bam1_t* read, int ref_pos) {
    // 讀段起始參考位置 (0-based)
    int start_ref_pos = read->core.pos;
    
    // 如果變異位點在讀段開始位置之前或結束位置之後，無法映射
    if (ref_pos < start_ref_pos) {
        return -1;
    }
    
    // 獲取CIGAR
    uint32_t* cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    
    int current_ref_pos = start_ref_pos;
    int current_read_pos = 0;
    
    // 遍歷CIGAR操作
    for (int i = 0; i < n_cigar; ++i) {
        uint32_t cigar_op = cigar[i];
        int op_len = bam_cigar_oplen(cigar_op);
        char op_char = bam_cigar_opchr(cigar_op);
        
        if (op_char == 'M' || op_char == '=' || op_char == 'X') {
            // 匹配區域，同時移動參考和讀段位置
            if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + op_len) {
                // 變異位點在當前匹配區域內
                return current_read_pos + (ref_pos - current_ref_pos);
            }
            
            current_ref_pos += op_len;
            current_read_pos += op_len;
        } else if (op_char == 'I' || op_char == 'S') {
            // 插入或軟剪切，僅移動讀段位置
            current_read_pos += op_len;
        } else if (op_char == 'D' || op_char == 'N') {
            // 刪除或跳過，僅移動參考位置
            if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + op_len) {
                // 變異位點在刪除區域內，無法映射
                return -1;
            }
            
            current_ref_pos += op_len;
        }
    }
    
    // 無法映射
    return -1;
}

char MethylHaploExtractor::getBaseAtReadPos(const bam1_t* read, int read_pos) {
    // 獲取讀段序列
    uint8_t* seq = bam_get_seq(read);
    
    // 檢查讀段位置是否有效
    if (read_pos < 0 || read_pos >= read->core.l_qseq) {
        return 'N';
    }
    
    // 獲取鹼基
    char base = seq_nt16_str[bam_seqi(seq, read_pos)];
    
    return base;
}

} // namespace msa::core 