#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <set>
#include <htslib/sam.h>

namespace msa {

/**
 * @brief 程式配置結構體
 */
struct Config {
    std::vector<std::string> vcf_files;  // VCF檔案列表
    std::string tumor_bam;               // 腫瘤樣本BAM檔案
    std::string normal_bam;              // 對照樣本BAM檔案
    std::string ref_file;                // 參考基因組檔案
    std::string bed_file;                // BED檔案（限制分析區域）
    std::string outdir;                  // 輸出目錄
    
    int window_size = 500;                // 甲基化窗口大小 (bp)
    float meth_high_threshold = 0.7f;     // 高甲基化閾值
    float meth_low_threshold = 0.3f;      // 低甲基化閾值
    float min_allele = 0.2f;              // 最小等位基因頻率
    int min_strand_reads = 3;             // 每條鏈上要求的最小讀段數
    int threads = 1;                      // 執行緒數
    bool gzip_output = true;              // 是否壓縮輸出
    int max_read_depth = 10000;           // 最大讀取深度
    int max_ram_gb = 32;                  // 最大RAM使用量(GB)
    std::string log_level = "TRACE";       // 日誌級別
    
    // BAM標籤檢測結果
    bool tumor_has_methyl_tags = false;   // 腫瘤BAM是否有甲基化標籤
    bool normal_has_methyl_tags = false;  // 正常BAM是否有甲基化標籤
    bool tumor_has_hp_tags = false;       // 腫瘤BAM是否有單倍型標籤
    bool normal_has_hp_tags = false;      // 正常BAM是否有單倍型標籤
    
    // 命令行選項
    bool help = false;                    // 顯示幫助信息
    bool version = false;                 // 顯示版本信息
};

/**
 * @brief VCF變異信息結構體
 */
struct VcfVariantInfo {
    std::string chrom;              // 染色體
    int pos;                        // 位置 (1-based)
    std::string ref;                // 參考基因組鹼基
    std::string alt;                // 變異鹼基
    std::string variant_type;       // 變異類型 (SNV, INS, DEL等)
    std::string vcf_source_id;      // VCF來源ID
    float allele_freq = 0.0f;       // 等位基因頻率
    float qual = 0.0f;              // 品質分數
    
    // 為了在容器中進行排序
    bool operator<(const VcfVariantInfo& other) const {
        if (chrom != other.chrom) return chrom < other.chrom;
        if (pos != other.pos) return pos < other.pos;
        return variant_type < other.variant_type;
    }
};

/**
 * @brief 甲基化位點詳細信息結構體
 */
struct MethylationSiteDetail {
    std::string chrom;                 // 染色體
    int methyl_pos = 0;                // 甲基化位點位置 (1-based)
    int somatic_pos = 0;               // 體細胞變異位置 (1-based)
    std::string variant_type;          // 變異類型
    std::string vcf_source_id;         // VCF來源ID
    std::string bam_source_id;         // BAM來源ID
    std::string somatic_allele_type;   // 體細胞等位基因類型 (ref/alt)
    std::string somatic_base_at_variant;  // 讀段在變異位置的鹼基
    std::string haplotype_tag;         // 單倍型標籤
    float meth_call = 0.0f;            // 甲基化程度 (0-1)
    std::string meth_state;            // 甲基化狀態 (high/mid/low)
    char strand = '.';                 // 鏈方向 (+/-)
    std::string read_id;               // 讀段ID
};

/**
 * @brief 體細胞變異甲基化摘要結構體
 */
struct SomaticVariantMethylationSummary {
    std::string chrom;                  // 染色體
    int somatic_pos = 0;                // 體細胞變異位置 (1-based)
    std::string variant_type;           // 變異類型
    std::string vcf_source_id;          // VCF來源ID
    std::string bam_source_id;          // BAM來源ID
    std::string somatic_allele_type;    // 體細胞等位基因類型 (ref/alt)
    std::string haplotype_tag;          // 單倍型標籤
    int supporting_read_count = 0;      // 支持該摘要的讀段數量
    int methyl_sites_count = 0;         // 甲基化位點數量
    float mean_methylation = 0.0f;      // 平均甲基化程度
    char strand = '.';                  // 主要鏈方向 (+/-/.)
};

/**
 * @brief 聚合單倍型統計結構體
 */
struct AggregatedHaplotypeStats {
    std::string haplotype_group;        // 單倍型組 (1/2/0等)
    std::string bam_source;             // BAM來源
    std::string variant_type_group;     // 變異類型組 (SNV/INDEL等)
    std::map<std::string, float> vcf_methylation_means;  // 每個VCF的平均甲基化
    float difference = 0.0f;            // VCF之間的甲基化差異
    float p_value = 1.0f;               // 統計顯著性p值
};

/**
 * @brief 全域摘要指標結構體
 */
struct GlobalSummaryMetrics {
    std::map<std::string, std::string> parameters;           // 分析參數
    std::map<std::string, std::string> numeric_metrics_str;  // 數值指標
};

/**
 * @brief 分析結果結構體
 */
struct AnalysisResults {
    std::vector<MethylationSiteDetail> level1_details;                 // Level 1: 原始甲基化詳情
    std::vector<SomaticVariantMethylationSummary> level2_summary;      // Level 2: 變異甲基化摘要
    std::vector<AggregatedHaplotypeStats> level3_stats;                // Level 3: 單倍型統計
    GlobalSummaryMetrics global_metrics;                               // 全域摘要指標
};

// 用於記憶體池的工作項目類型
enum class WorkItemType {
    None,       // 無任務
    Read,       // 讀段處理任務
    Exit        // 退出信號
};

// 工作項目結構體
struct WorkItem {
    WorkItemType type = WorkItemType::None;
    bam1_t* read = nullptr;
    VcfVariantInfo variant_info;
    std::string bam_source_id;
};

} // namespace msa 