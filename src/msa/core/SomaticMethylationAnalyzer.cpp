#include "msa/core/SomaticMethylationAnalyzer.h"
#include "msa/utils/LogManager.h"
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <numeric>
#include <iomanip>

// 使用正確的命名空間
using namespace msa::utils;

namespace msa::core {

SomaticMethylationAnalyzer::SomaticMethylationAnalyzer(const msa::Config& config)
    : config_(config) {
}

msa::AnalysisResults SomaticMethylationAnalyzer::analyze(const std::vector<msa::MethylationSiteDetail>& sites) {
    LOG_INFO("SomaticMethylationAnalyzer", "開始分析 " + std::to_string(sites.size()) + " 個甲基化位點");
    
    msa::AnalysisResults results;
    
    // 保存原始位點詳細資訊
    results.level1_details = sites;
    
    // 過濾雙股覆蓋不足的位點
    auto filtered_sites = filterSitesByStrandCoverage(sites);
    LOG_INFO("SomaticMethylationAnalyzer", "雙股覆蓋篩選後保留 " + std::to_string(filtered_sites.size()) + " 個位點");
    
    // 生成Level 2摘要統計
    results.level2_summary = generateLevel2Summary(filtered_sites);
    LOG_INFO("SomaticMethylationAnalyzer", "生成 " + std::to_string(results.level2_summary.size()) + " 個Level 2摘要記錄");
    
    // 生成Level 3聚合統計
    results.level3_stats = generateLevel3Statistics(results.level2_summary);
    LOG_INFO("SomaticMethylationAnalyzer", "生成 " + std::to_string(results.level3_stats.size()) + " 個Level 3聚合統計");
    
    // 計算全域摘要指標
    results.global_metrics = calculateGlobalMetrics(sites, results.level2_summary);
    
    return results;
}

std::vector<msa::MethylationSiteDetail> SomaticMethylationAnalyzer::filterSitesByStrandCoverage(
    const std::vector<msa::MethylationSiteDetail>& sites) {
    
    // 如果min_strand_reads為0，則不需要過濾
    if (config_.min_strand_reads <= 0) {
        return sites;
    }
    
    // 統計每個位點在正反鏈上的覆蓋
    std::map<std::string, std::map<char, int>> position_strand_counts;
    
    // 計算位點標識符的輔助函數
    auto getPositionKey = [](const msa::MethylationSiteDetail& site) {
        return site.chrom + ":" + std::to_string(site.methyl_pos) + ":" + site.bam_source_id;
    };
    
    // 統計每個位點在正反鏈的覆蓋數
    for (const auto& site : sites) {
        std::string pos_key = getPositionKey(site);
        position_strand_counts[pos_key][site.strand]++;
    }
    
    // 確定哪些位點符合雙股覆蓋要求
    std::set<std::string> valid_positions;
    for (const auto& [pos_key, strand_counts] : position_strand_counts) {
        // 檢查正鏈('+'）覆蓋
        bool plus_valid = strand_counts.count('+') && strand_counts.at('+') >= config_.min_strand_reads;
        // 檢查反鏈('-')覆蓋
        bool minus_valid = strand_counts.count('-') && strand_counts.at('-') >= config_.min_strand_reads;
        
        // 同時滿足正反鏈覆蓋要求
        if (plus_valid && minus_valid) {
            valid_positions.insert(pos_key);
        }
    }
    
    // 過濾出符合要求的位點
    std::vector<msa::MethylationSiteDetail> filtered_sites;
    for (const auto& site : sites) {
        std::string pos_key = getPositionKey(site);
        if (valid_positions.count(pos_key)) {
            filtered_sites.push_back(site);
        }
    }
    
    return filtered_sites;
}

std::vector<msa::SomaticVariantMethylationSummary> SomaticMethylationAnalyzer::generateLevel2Summary(
    const std::vector<msa::MethylationSiteDetail>& sites) {
    
    // 按分組鍵聚合數據
    std::map<std::string, std::vector<msa::MethylationSiteDetail>> groupedSites;
    std::map<std::string, std::set<std::string>> groupReadCounts;
    
    // 分組鍵格式: chrom:pos:variant_type:vcf_source_id:bam_source_id:somatic_allele_type:haplotype_tag
    for (const auto& site : sites) {
        std::ostringstream key;
        key << site.chrom << ":" << site.somatic_pos << ":" 
            << site.variant_type << ":" << site.vcf_source_id << ":" 
            << site.bam_source_id << ":" << site.somatic_allele_type << ":" 
            << site.haplotype_tag;
        
        std::string group_key = key.str();
        groupedSites[group_key].push_back(site);
        
        // 追蹤每個組中唯一的讀段
        groupReadCounts[group_key].insert(site.read_id);
    }
    
    // 生成Level 2摘要
    std::vector<msa::SomaticVariantMethylationSummary> summaries;
    
    for (const auto& [group_key, group_sites] : groupedSites) {
        // 跳過沒有足夠位點的組
        if (group_sites.empty()) {
            continue;
        }
        
        // 使用第一個站點獲取基本信息
        const auto& first_site = group_sites[0];
        
        // 創建摘要
        msa::SomaticVariantMethylationSummary summary;
        summary.chrom = first_site.chrom;
        summary.somatic_pos = first_site.somatic_pos;
        summary.variant_type = first_site.variant_type;
        summary.vcf_source_id = first_site.vcf_source_id;
        summary.bam_source_id = first_site.bam_source_id;
        summary.somatic_allele_type = first_site.somatic_allele_type;
        summary.haplotype_tag = first_site.haplotype_tag;
        
        // 計算支持該摘要的讀段數量
        summary.supporting_read_count = groupReadCounts[group_key].size();
        
        // 計算甲基化位點數量
        summary.methyl_sites_count = group_sites.size();
        
        // 計算平均甲基化水平
        double total_meth = 0.0;
        for (const auto& site : group_sites) {
            total_meth += site.meth_call;
        }
        summary.mean_methylation = group_sites.empty() ? 0.0f : static_cast<float>(total_meth / group_sites.size());
        
        // 確定主要鏈方向
        int plus_count = 0, minus_count = 0;
        for (const auto& site : group_sites) {
            if (site.strand == '+') plus_count++;
            else if (site.strand == '-') minus_count++;
        }
        
        if (plus_count > minus_count && plus_count > 0) {
            summary.strand = '+';
        } else if (minus_count > plus_count && minus_count > 0) {
            summary.strand = '-';
        } else {
            summary.strand = '.';  // 混合或未知
        }
        
        summaries.push_back(summary);
    }
    
    return summaries;
}

std::vector<msa::AggregatedHaplotypeStats> SomaticMethylationAnalyzer::generateLevel3Statistics(
    const std::vector<msa::SomaticVariantMethylationSummary>& level2Summary) {
    
    // 按單倍型組、樣本來源和變異類型分組
    std::map<std::string, std::map<std::string, std::vector<float>>> grouped_methylation;
    
    // 分組鍵格式: haplotype_group:bam_source:variant_type_group
    for (const auto& summary : level2Summary) {
        // 簡化變異類型分組（例如，合併所有INDELs）
        std::string variant_type_group = summary.variant_type;
        if (variant_type_group == "INS" || variant_type_group == "DEL") {
            variant_type_group = "INDEL";
        }
        
        // 創建分組鍵
        std::string group_key = summary.haplotype_tag + ":" + 
                              summary.bam_source_id + ":" + 
                              variant_type_group;
        
        // 將甲基化平均值添加到對應VCF的列表中
        grouped_methylation[group_key][summary.vcf_source_id].push_back(summary.mean_methylation);
    }
    
    // 生成Level 3統計
    std::vector<msa::AggregatedHaplotypeStats> stats;
    
    // 獲取唯一VCF來源
    std::set<std::string> vcf_sources;
    for (const auto& summary : level2Summary) {
        vcf_sources.insert(summary.vcf_source_id);
    }
    
    // 對每個分組生成統計
    for (const auto& [group_key, vcf_methyl_values] : grouped_methylation) {
        // 解析分組鍵
        std::string haplotype_group, bam_source, variant_type_group;
        {
            std::istringstream ss(group_key);
            std::getline(ss, haplotype_group, ':');
            std::getline(ss, bam_source, ':');
            std::getline(ss, variant_type_group, ':');
        }
        
        // 創建聚合統計
        msa::AggregatedHaplotypeStats stat;
        stat.haplotype_group = haplotype_group;
        stat.bam_source = bam_source;
        stat.variant_type_group = variant_type_group;
        
        // 計算每個VCF的平均甲基化
        for (const auto& vcf_source : vcf_sources) {
            if (vcf_methyl_values.count(vcf_source) && !vcf_methyl_values.at(vcf_source).empty()) {
                const auto& values = vcf_methyl_values.at(vcf_source);
                float mean = std::accumulate(values.begin(), values.end(), 0.0f) / values.size();
                stat.vcf_methylation_means[vcf_source] = mean;
            } else {
                // 此VCF在此分組中沒有數據
                stat.vcf_methylation_means[vcf_source] = -1.0f;  // 表示缺失值
            }
        }
        
        // 如果有多個VCF，計算它們之間的差異和p值
        if (vcf_sources.size() >= 2) {
            std::vector<std::string> vcf_list(vcf_sources.begin(), vcf_sources.end());
            const std::string& vcf1 = vcf_list[0];
            const std::string& vcf2 = vcf_list[1];
            
            // 檢查兩個VCF都有值
            if (stat.vcf_methylation_means[vcf1] >= 0 && stat.vcf_methylation_means[vcf2] >= 0) {
                stat.difference = stat.vcf_methylation_means[vcf1] - stat.vcf_methylation_means[vcf2];
                
                // 如果有原始數據，計算p值
                if (vcf_methyl_values.count(vcf1) && vcf_methyl_values.count(vcf2) &&
                    !vcf_methyl_values.at(vcf1).empty() && !vcf_methyl_values.at(vcf2).empty()) {
                    stat.p_value = calculatePValue(vcf_methyl_values.at(vcf1), vcf_methyl_values.at(vcf2));
                } else {
                    stat.p_value = 1.0f;  // 無法計算p值
                }
            } else {
                stat.difference = 0.0f;
                stat.p_value = 1.0f;
            }
        } else {
            stat.difference = 0.0f;
            stat.p_value = 1.0f;
        }
        
        stats.push_back(stat);
    }
    
    return stats;
}

msa::GlobalSummaryMetrics SomaticMethylationAnalyzer::calculateGlobalMetrics(
    const std::vector<msa::MethylationSiteDetail>& sites,
    const std::vector<msa::SomaticVariantMethylationSummary>& level2Summary) {
    
    msa::GlobalSummaryMetrics metrics;
    
    // 添加配置參數
    metrics.parameters["vcf_files"] = config_.vcf_files.empty() ? "None" : 
                                   config_.vcf_files.size() == 1 ? config_.vcf_files[0] :
                                   std::to_string(config_.vcf_files.size()) + " files";
    metrics.parameters["tumor_bam"] = config_.tumor_bam;
    metrics.parameters["normal_bam"] = config_.normal_bam;
    metrics.parameters["ref_file"] = config_.ref_file;
    metrics.parameters["bed_file"] = config_.bed_file.empty() ? "None" : config_.bed_file;
    metrics.parameters["window_size"] = std::to_string(config_.window_size);
    metrics.parameters["meth_high_threshold"] = std::to_string(config_.meth_high_threshold);
    metrics.parameters["meth_low_threshold"] = std::to_string(config_.meth_low_threshold);
    metrics.parameters["min_allele"] = std::to_string(config_.min_allele);
    metrics.parameters["min_strand_reads"] = std::to_string(config_.min_strand_reads);
    metrics.parameters["threads"] = std::to_string(config_.threads);
    
    // 收集統計指標
    std::map<std::string, std::map<std::string, int>> vcf_source_stats;  // [vcf_source][metric] = value
    std::map<std::string, std::map<std::string, double>> bam_source_meth_stats;  // [bam_source][metric] = value
    
    // 計算每個VCF源中的變異數
    std::set<std::string> variant_keys;
    for (const auto& site : sites) {
        std::string variant_key = site.vcf_source_id + ":" + site.chrom + ":" + 
                                std::to_string(site.somatic_pos) + ":" + site.variant_type;
        variant_keys.insert(variant_key);
        vcf_source_stats[site.vcf_source_id]["total_variants"]++;
    }
    
    // 計算處理的變異數
    for (const auto& key : variant_keys) {
        std::string vcf_source = key.substr(0, key.find(':'));
        vcf_source_stats[vcf_source]["processed_variants"]++;
    }
    
    // 計算甲基化位點數量和平均甲基化
    for (const auto& site : sites) {
        // 僅考慮高或中度甲基化的位點
        if (site.meth_state == "high" || site.meth_state == "mid") {
            bam_source_meth_stats[site.bam_source_id]["methyl_site_count"]++;
            bam_source_meth_stats[site.bam_source_id]["total_meth"] += site.meth_call;
        }
        bam_source_meth_stats[site.bam_source_id]["total_sites"]++;
    }
    
    // 將統計數據添加到指標
    for (const auto& [vcf_source, stats] : vcf_source_stats) {
        metrics.numeric_metrics_str[vcf_source + "_total_variants"] = std::to_string(stats.at("total_variants"));
        metrics.numeric_metrics_str[vcf_source + "_processed_variants"] = std::to_string(stats.at("processed_variants"));
    }
    
    for (const auto& [bam_source, stats] : bam_source_meth_stats) {
        std::ostringstream mean_meth_stream;
        double methyl_site_count = stats.at("methyl_site_count");
        double total_meth = stats.at("total_meth");
        double total_sites = stats.at("total_sites");
        double mean_meth = methyl_site_count > 0 ? total_meth / methyl_site_count : 0.0;
        
        // 格式化浮點數，保留四位小數
        mean_meth_stream << std::fixed << std::setprecision(4) << mean_meth;
        
        metrics.numeric_metrics_str[bam_source + "_methylated_site_count"] = std::to_string(static_cast<int>(methyl_site_count));
        metrics.numeric_metrics_str[bam_source + "_total_site_count"] = std::to_string(static_cast<int>(total_sites));
        metrics.numeric_metrics_str[bam_source + "_mean_methylation"] = mean_meth_stream.str();
    }
    
    return metrics;
}

float SomaticMethylationAnalyzer::calculatePValue(const std::vector<float>& group1, const std::vector<float>& group2) {
    // 簡單t檢驗實現
    // 注意：這是一個非常簡化的t檢驗實現，實際應用中應使用統計庫
    
    // 檢查樣本大小
    if (group1.size() < 2 || group2.size() < 2) {
        return 1.0f;  // 樣本太小無法計算
    }
    
    // 計算組1的均值和方差
    float mean1 = std::accumulate(group1.begin(), group1.end(), 0.0f) / group1.size();
    float var1 = 0.0f;
    for (float val : group1) {
        var1 += (val - mean1) * (val - mean1);
    }
    var1 /= (group1.size() - 1);
    
    // 計算組2的均值和方差
    float mean2 = std::accumulate(group2.begin(), group2.end(), 0.0f) / group2.size();
    float var2 = 0.0f;
    for (float val : group2) {
        var2 += (val - mean2) * (val - mean2);
    }
    var2 /= (group2.size() - 1);
    
    // 計算t統計量
    float t_stat = std::abs(mean1 - mean2) / 
                  std::sqrt((var1 / group1.size()) + (var2 / group2.size()));
    
    // 自由度 (approximation)
    float df = group1.size() + group2.size() - 2;
    
    // 簡化的p值計算 (使用一個近似公式)
    // 注意：這只是一個粗略的近似值，實際應用中應使用更準確的方法
    float p_value = 1.0f / (1.0f + t_stat * std::sqrt(df / 2.0f));
    
    return p_value;
}

} // namespace msa::core 