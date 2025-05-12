#pragma once

#include <vector>
#include <string>
#include <map>
#include <set>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief 體細胞甲基化分析器
 */
class SomaticMethylationAnalyzer {
public:
    /**
     * @brief 建構函數
     * @param config 配置物件
     */
    SomaticMethylationAnalyzer(const msa::Config& config);
    
    /**
     * @brief 分析甲基化位點數據
     * @param sites 甲基化位點詳情
     * @return msa::AnalysisResults 分析結果
     */
    msa::AnalysisResults analyze(const std::vector<msa::MethylationSiteDetail>& sites);
    
private:
    /**
     * @brief 根據雙股覆蓋條件過濾甲基化位點
     * @param sites 原始甲基化位點詳情
     * @return std::vector<msa::MethylationSiteDetail> 過濾後的位點
     */
    std::vector<msa::MethylationSiteDetail> filterSitesByStrandCoverage(
        const std::vector<msa::MethylationSiteDetail>& sites);
    
    /**
     * @brief 生成Level 2甲基化摘要
     * @param sites 甲基化位點詳情
     * @return std::vector<msa::SomaticVariantMethylationSummary> Level 2摘要
     */
    std::vector<msa::SomaticVariantMethylationSummary> generateLevel2Summary(
        const std::vector<msa::MethylationSiteDetail>& sites);
    
    /**
     * @brief 生成Level 3甲基化統計
     * @param level2Summary Level 2摘要
     * @return std::vector<msa::AggregatedHaplotypeStats> Level 3統計
     */
    std::vector<msa::AggregatedHaplotypeStats> generateLevel3Statistics(
        const std::vector<msa::SomaticVariantMethylationSummary>& level2Summary);
    
    /**
     * @brief 計算全域指標
     * @param sites 原始甲基化位點詳情
     * @param level2Summary Level 2摘要
     * @return msa::GlobalSummaryMetrics 全域指標
     */
    msa::GlobalSummaryMetrics calculateGlobalMetrics(
        const std::vector<msa::MethylationSiteDetail>& sites,
        const std::vector<msa::SomaticVariantMethylationSummary>& level2Summary);
    
    /**
     * @brief 計算統計p值
     * @param group1 第一組數據
     * @param group2 第二組數據
     * @return float p值
     */
    float calculatePValue(const std::vector<float>& group1, const std::vector<float>& group2);
    
    // 配置物件
    const msa::Config& config_;
};

} // namespace msa::core 