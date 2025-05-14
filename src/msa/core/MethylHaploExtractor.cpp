#include "msa/core/MethylHaploExtractor.h"
#include "msa/utils/LogManager.h"
#include <sstream>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>


namespace msa::core {

// 輔助結構：記錄甲基化信息
struct MethylationRecord {
    int refPos;            // 1-based 參考座標
    double prob;           // 甲基化機率（0~1）
    int methyl_type;       // 甲基化類型: 1=高, 0=中, -1=低
    char strand;           // 鏈方向: '+', '-'
};

/**
 * @brief 從read到reference的位置映射表
 * @param aln BAM讀段
 * @return 讀段位置(0-based)到參考位置(1-based)的映射向量，-1表示不可映射
 */
static std::vector<int> buildReadToRefMap(const bam1_t *aln) {
    int readLength = aln->core.l_qseq;
    std::vector<int> readToRef(readLength, -1);
    int refPos = aln->core.pos + 1; // 轉為1-based
    uint32_t *cigar = bam_get_cigar(aln);
    int nCigar = aln->core.n_cigar;
    int readPos = 0; // 0-based讀段位置，從0開始
    
    // 記錄映射的起始位置，必要時處理非映射區域（如soft clip）
    std::string cigar_str = "";
    for (int i = 0; i < nCigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        char op_char = "MIDNSHP=X"[op];
        cigar_str += std::to_string(len) + op_char;
    }
    LOG_DEBUG("MethylHaploExtractor", "讀段 " + std::string(bam_get_qname(aln)) + 
              " CIGAR: " + cigar_str + ", 讀段長度: " + std::to_string(readLength) + 
              ", 起始參考位置: " + std::to_string(aln->core.pos));
    
    for (int i = 0; i < nCigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        char op_char = "MIDNSHP=X"[op];
        
        LOG_TRACE("MethylHaploExtractor", "處理CIGAR操作 " + std::string(1, op_char) + 
                 std::to_string(len) + " 讀段位置: " + std::to_string(readPos) + 
                 " 參考位置: " + std::to_string(refPos));
        
        switch (op) {
            case BAM_CMATCH:    // M - 序列匹配或不匹配
            case BAM_CEQUAL:    // = - 序列匹配
            case BAM_CDIFF:     // X - 序列不匹配
                // 匹配區域，同時更新讀段和參考座標
                for (int j = 0; j < len; j++) {
                    if (readPos < readLength) {
                        readToRef[readPos] = refPos;
                    }
                    readPos++;
                    refPos++;
                }
                break;
                
            case BAM_CINS:      // I - 插入到參考
                // 插入區域，只消耗讀段座標
                readPos += len;
                break;
                
            case BAM_CSOFT_CLIP: // S - 軟剪切（存在於SEQ中）
                // 軟剪切區域，只消耗讀段座標，但這些位置沒有對應的參考座標
                readPos += len;
                break;
                
            case BAM_CDEL:      // D - 從參考中刪除
            case BAM_CREF_SKIP: // N - 跳過參考（例如內含子）
                // 刪除或跳過區域，只消耗參考座標
                refPos += len;
                break;
                
            case BAM_CHARD_CLIP: // H - 硬剪切（不存在於SEQ中）
                // 硬剪切不影響座標計算，因為這些序列不存在於讀段SEQ中
                break;
                
            case BAM_CPAD:       // P - 填充（無意義，通常不使用）
                // 填充不影響座標計算
                break;
                
            default:
                // 其他情況忽略
                break;
        }
    }
    
    // 對於調試目的，輸出部分映射情況
    LOG_TRACE("MethylHaploExtractor", "讀段到參考的映射 (前10個位置):");
    for (int i = 0; i < std::min(10, readLength); i++) {
        LOG_TRACE("MethylHaploExtractor", "  讀段位置 " + std::to_string(i) + " -> 參考位置 " + std::to_string(readToRef[i]));
    }
    
    return readToRef;
}

/**
 * @brief 計算甲基化類型
 * @param prob 甲基化機率
 * @param meth_high_threshold 高甲基化閾值
 * @param meth_low_threshold 低甲基化閾值
 * @return 甲基化類型: 1=高, 0=中, -1=低
 */
static int calculateMethylType(double prob, float meth_high_threshold, float meth_low_threshold) {
    if (prob >= meth_high_threshold)
        return 1;
    else if (prob <= meth_low_threshold)
        return -1;
    else
        return 0;
}

/**
 * @brief 從BAM讀段中解析甲基化記錄
 * @param aln BAM讀段
 * @param config 配置物件
 * @return 甲基化記錄列表
 */
static std::vector<MethylationRecord> parseMethylationRecords(const bam1_t *aln, const msa::Config& config) {
    std::vector<MethylationRecord> records;
    
    // 獲取讀段ID以便日誌記錄
    std::string read_id(bam_get_qname(aln));
    LOG_DEBUG("MethylHaploExtractor", "開始解析讀段 " + read_id + " 的甲基化數據");
    
    // 檢查是否有MM和ML標籤
    uint8_t* mm_tag = bam_aux_get(aln, "MM");
    uint8_t* ml_tag = bam_aux_get(aln, "ML");
    
    // 如果MM或ML標籤不存在，則嘗試Mm和Ml標籤
    if (!mm_tag || !ml_tag) {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 沒有MM/ML標籤，嘗試Mm/Ml標籤");
        mm_tag = bam_aux_get(aln, "Mm");
        ml_tag = bam_aux_get(aln, "Ml");
    }
    
    // 檢查是否找到甲基化標籤
    if (mm_tag) {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 找到MM標籤: " + 
                 (mm_tag[0] == 'Z' ? std::string(bam_aux2Z(mm_tag)) : "非字符串類型"));
    } else {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 沒有找到MM或Mm標籤");
    }
    
    if (ml_tag) {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 找到ML標籤，類型: " + 
                 std::string(1, ml_tag[0]) + ", " + std::string(1, ml_tag[1]));
    } else {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 沒有找到ML或Ml標籤");
    }
    
    // 分配甲基化狀態物件
    hts_base_mod_state *modState = hts_base_mod_state_alloc();
    if (!modState) {
        LOG_ERROR("MethylHaploExtractor", "無法分配甲基化狀態物件");
        return records;
    }
    
    // 解析讀段中的甲基化信息
    int ret = bam_parse_basemod(aln, modState);
    if (ret < 0) {
        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 無甲基化標記或解析失敗，返回代碼: " + std::to_string(ret));
        hts_base_mod_state_free(modState);
        return records;
    }
    
    // 建立讀段到參考的映射
    auto readToRef = buildReadToRefMap(aln);
    int readLength = aln->core.l_qseq;
    
    // 獲取鏈方向
    char strand = ((aln->core.flag & BAM_FREVERSE) == 0) ? '+' : '-';
    
    // 用於存儲甲基化修飾
    hts_base_mod mods[10]; // 一次獲取最多10個修飾
    int readPos0;
    
    // 逐個獲取甲基化修飾
    // bam_next_basemod返回的readPos0是讀段上的0-based位置，表示甲基化修飾發生的位置
    // 這個位置需要通過readToRef映射到參考基因組上的位置
    int n = bam_next_basemod(aln, modState, mods, 10, &readPos0);
    LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 第一批甲基化修飾數量: " + std::to_string(n));
    
    int total_mods = 0;
    while (n > 0) {
        total_mods += n;
        for (int i = 0; i < n; i++) {
            // 我們只處理'C'鹼基上的'm'或'h'修飾，即5mC或5hmC
            if ((mods[i].modified_base == 'm' || mods[i].modified_base == 'h') && 
                (mods[i].canonical_base == 'C' || mods[i].canonical_base == 'c')) {
                
                // 確認修飾在讀段範圍內
                if (readPos0 >= 0 && readPos0 < readLength) {
                    // 獲取對應的參考座標
                    // readPos0是0-based讀段位置，我們從映射表中直接獲取對應的1-based參考位置
                    int refPos = readToRef[readPos0];
                    
                    if (refPos > 0) { // 有效參考座標 (1-based)
                        // 計算甲基化機率
                        double prob = static_cast<double>(mods[i].qual) / 255.0;
                        
                        // 計算甲基化類型
                        int methyl_type = calculateMethylType(prob, config.meth_high_threshold, config.meth_low_threshold);
                        
                        // 添加甲基化記錄
                        records.push_back({refPos, prob, methyl_type, strand});
                        
                        LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 發現甲基化位點: readPos=" + std::to_string(readPos0) + 
                                  ", refPos=" + std::to_string(refPos) + 
                                  ", prob=" + std::to_string(prob) + 
                                  ", type=" + std::to_string(methyl_type) + 
                                  ", strand=" + strand +
                                  ", base=" + std::string(1, mods[i].canonical_base) + 
                                  ", mod=" + std::string(1, mods[i].modified_base));
                    } else {
                        LOG_TRACE("MethylHaploExtractor", "讀段 " + read_id + " 在讀段位置 " + std::to_string(readPos0) + 
                                 " 的甲基化修飾無法映射到參考座標");
                    }
                }
            }
        }
        // 獲取下一批修飾
        n = bam_next_basemod(aln, modState, mods, 10, &readPos0);
    }
    
    // 釋放甲基化狀態物件
    hts_base_mod_state_free(modState);
    
    LOG_DEBUG("MethylHaploExtractor", "讀段 " + read_id + " 解析完成，總共 " + std::to_string(total_mods) + " 個甲基化修飾，提取出 " + std::to_string(records.size()) + " 個甲基化位點");
    return records;
}

MethylHaploExtractor::MethylHaploExtractor(const msa::Config& config)
    : config_(config) {
}

std::vector<msa::MethylationSiteDetail> MethylHaploExtractor::extractFromRead(
    const bam1_t* read,
    const msa::VcfVariantInfo& target_variant,
    const std::string& bam_source_id
) {
    std::vector<msa::MethylationSiteDetail> details;
    
    // 檢查是否為有效讀段
    if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
        LOG_TRACE("MethylHaploExtractor", "跳過無效讀段: " + std::string(bam_get_qname(read)));
        return details;
    }
    
    // 獲取讀段ID (QNAME)
    std::string read_id(bam_get_qname(read));
    
    // 獲取單倍型標籤
    std::string haplotype_tag = extractHaplotypeTag(read);
    
    // 確定對變異的支持 (ref/alt)
    std::string somatic_base;
    std::string somatic_allele_type = determineAlleleType(read, target_variant, somatic_base);
    
    // 如果無法確定等位基因類型，跳過
    if (somatic_allele_type == "unknown") {
        LOG_TRACE("MethylHaploExtractor", "無法確定讀段等位基因類型，跳過: " + read_id);
        return details;
    }
    
    // 獲取目標變異的0-based位置
    int target_pos_0based = target_variant.pos - 1;
    
    // 從讀段中提取甲基化記錄
    std::vector<MethylationRecord> methRecords = parseMethylationRecords(read, config_);
    
    // 如果沒有甲基化記錄，跳過
    if (methRecords.empty()) {
        LOG_TRACE("MethylHaploExtractor", "讀段無甲基化記錄，跳過: " + read_id);
        return details;
    }
    
    // 根據window_size篩選符合條件的甲基化記錄
    for (const auto& rec : methRecords) {
        int distance = std::abs(rec.refPos - target_variant.pos);
        
        // 只處理在window範圍內的甲基化位點
        if (distance <= config_.window_size) {
            // 創建甲基化位點詳情
            msa::MethylationSiteDetail detail;
            detail.chrom = target_variant.chrom;
            detail.methyl_pos = rec.refPos;
            detail.somatic_pos = target_variant.pos;
            detail.variant_type = target_variant.variant_type;
            detail.vcf_source_id = target_variant.vcf_source_id;
            detail.bam_source_id = bam_source_id;
            detail.somatic_allele_type = somatic_allele_type;
            detail.somatic_base_at_variant = somatic_base;
            detail.haplotype_tag = haplotype_tag;
            detail.meth_call = static_cast<float>(rec.prob);
            detail.strand = rec.strand;
            detail.read_id = read_id;
            
            // 設定甲基化狀態
            detail.meth_state = classifyMethylationState(detail.meth_call);
            
            // 將甲基化位點詳情添加到結果中
            details.push_back(detail);
        }
    }
    
    return details;
}

std::string MethylHaploExtractor::extractHaplotypeTag(const bam1_t* read) {
    // 獲取HP標籤
    uint8_t* hp_data = bam_aux_get(read, "HP");
    if (hp_data) {
        int hp_value = bam_aux2i(hp_data);
        
        // 檢查PS標籤是否存在（部分數據集可能沒有PS標籤）
        uint8_t* ps_data = bam_aux_get(read, "PS");
        
        // 有效的單倍型值為1或2
        if (hp_value == 1 || hp_value == 2) {
            return std::to_string(hp_value);
        } else if (hp_value == 1-1 || hp_value == 2-1){
            // 非標準單倍型值
            return std::to_string(hp_value);
        }else if (hp_value == 3){
            // 
            return std::to_string(hp_value);
        }
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
    } else if (meth_call > config_.meth_low_threshold) {
        return "mid";
    } else {
        return "low";
    }
}

int MethylHaploExtractor::refPosToReadPos(const bam1_t* read, int ref_pos) {
    // 讀段起始參考位置 (0-based)
    int start_ref_pos = read->core.pos;
    
    // 如果變異位點在讀段開始位置之前，無法映射
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
        uint32_t op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            // 匹配區域，同時移動參考和讀段位置
            if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len) {
                // 變異位點在當前匹配區域內
                return current_read_pos + (ref_pos - current_ref_pos);
            }
            
            current_ref_pos += len;
            current_read_pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            // 插入或軟剪切，僅移動讀段位置
            current_read_pos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            // 刪除或跳過，僅移動參考位置
            if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len) {
                // 變異位點在刪除區域內，無法映射
                return -1;
            }
            
            current_ref_pos += len;
        }
        // 硬剪切和填充不影響座標計算
    }
    
    // 無法映射（讀段可能不覆蓋變異位點）
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