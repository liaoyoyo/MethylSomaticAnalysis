#pragma once

#include <string>
#include <vector>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief 結果匯出類，用於將分析結果輸出成TSV/JSON格式
 */
class ReportExporter {
public:
    /**
     * @brief 建構函數
     * @param config 配置物件
     */
    ReportExporter(const msa::Config& config);
    
    /**
     * @brief 匯出分析結果
     * @param results 分析結果
     * @param vcf_source_id VCF來源ID (用作輸出目錄名稱)
     * @return bool 匯出成功與否
     */
    bool exportResults(const msa::AnalysisResults& results, const std::string& vcf_source_id);
    
private:
    /**
     * @brief 匯出全域摘要指標
     * @param metrics 全域摘要指標
     * @param outputDir 輸出目錄
     * @return bool 匯出成功與否
     */
    bool exportGlobalSummary(const msa::GlobalSummaryMetrics& metrics, const std::string& outputDir);
    
    /**
     * @brief 匯出Level 1原始甲基化詳情
     * @param details 原始甲基化詳情
     * @param outputDir 輸出目錄
     * @return bool 匯出成功與否
     */
    bool exportLevel1Details(const std::vector<msa::MethylationSiteDetail>& details, const std::string& outputDir);
    
    /**
     * @brief 匯出Level 2變異甲基化摘要
     * @param summary 變異甲基化摘要
     * @param outputDir 輸出目錄
     * @return bool 匯出成功與否
     */
    bool exportLevel2Summary(const std::vector<msa::SomaticVariantMethylationSummary>& summary, const std::string& outputDir);
    
    /**
     * @brief 匯出Level 3單倍型統計
     * @param stats 單倍型統計
     * @param outputDir 輸出目錄
     * @return bool 匯出成功與否
     */
    bool exportLevel3Stats(const std::vector<msa::AggregatedHaplotypeStats>& stats, const std::string& outputDir);
    
    /**
     * @brief 創建輸出目錄
     * @param dirPath 目錄路徑
     * @return bool 創建成功與否
     */
    bool createDirectory(const std::string& dirPath);
    
    /**
     * @brief 壓縮檔案
     * @param inputPath 輸入檔案路徑
     * @param outputPath 輸出檔案路徑
     * @return bool 壓縮成功與否
     */
    bool compressFile(const std::string& inputPath, const std::string& outputPath);
    
    // 配置參數
    const msa::Config& config_;
};

} // namespace msa::core 