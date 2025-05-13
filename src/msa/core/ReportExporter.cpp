#include "msa/core/ReportExporter.h"
#include "msa/utils/LogManager.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <zlib.h>

// 使用正確的命名空間
using namespace msa::utils;
namespace fs = std::filesystem;

namespace msa::core {

ReportExporter::ReportExporter(const msa::Config& config)
    : config_(config) {
}

bool ReportExporter::exportResults(const msa::AnalysisResults& results, const std::string& vcf_source_id) {
    // 建立輸出目錄 (basedir/vcf_source_id)
    std::string outputDir = config_.outdir + "/" + vcf_source_id;
    
    if (!createDirectory(outputDir)) {
        LOG_ERROR("ReportExporter", "無法建立輸出目錄: " + outputDir);
        return false;
    }
    
    LOG_INFO("ReportExporter", "開始匯出VCF[" + vcf_source_id + "]的結果至目錄: " + outputDir);
    
    // 匯出全域摘要指標
    if (!exportGlobalSummary(results.global_metrics, outputDir)) {
        LOG_ERROR("ReportExporter", "匯出全域摘要指標失敗");
        return false;
    }
    
    // 匯出Level 1原始甲基化詳情
    if (!exportLevel1Details(results.level1_details, outputDir)) {
        LOG_ERROR("ReportExporter", "匯出Level 1原始甲基化詳情失敗");
        return false;
    }
    
    // 匯出Level 2變異甲基化摘要
    if (!exportLevel2Summary(results.level2_summary, outputDir)) {
        LOG_ERROR("ReportExporter", "匯出Level 2變異甲基化摘要失敗");
        return false;
    }
    
    // 匯出Level 3單倍型統計
    if (!exportLevel3Stats(results.level3_stats, outputDir)) {
        LOG_ERROR("ReportExporter", "匯出Level 3單倍型統計失敗");
        return false;
    }
    
    LOG_INFO("ReportExporter", "所有結果已成功匯出至目錄: " + outputDir);
    return true;
}

bool ReportExporter::exportGlobalSummary(const msa::GlobalSummaryMetrics& metrics, const std::string& outputDir) {
    std::string outputPath = outputDir + "/global_summary_metrics.tsv";
    std::ofstream outFile(outputPath);
    
    if (!outFile.is_open()) {
        LOG_ERROR("ReportExporter", "無法開啟輸出檔案: " + outputPath);
        return false;
    }
    
    // 寫入參數區段標題
    outFile << "# 參數\n";
    outFile << "parameter_name\tparameter_value\n";
    
    // 寫入參數
    for (const auto& [name, value] : metrics.parameters) {
        outFile << name << "\t" << value << "\n";
    }
    
    // 寫入統計數值區段標題
    outFile << "\n# 統計數值\n";
    outFile << "metric_name\tmetric_value\n";
    
    // 寫入統計數值
    for (const auto& [name, value] : metrics.numeric_metrics_str) {
        outFile << name << "\t" << value << "\n";
    }
    
    outFile.close();
    LOG_INFO("ReportExporter", "已匯出全域摘要指標: " + outputPath);
    return true;
}

bool ReportExporter::exportLevel1Details(const std::vector<msa::MethylationSiteDetail>& details, const std::string& outputDir) {
    std::string tempPath = outputDir + "/level1_raw_methylation_details.tsv";
    std::ofstream outFile(tempPath);
    
    if (!outFile.is_open()) {
        LOG_ERROR("ReportExporter", "無法開啟輸出檔案: " + tempPath);
        return false;
    }
    
    // 寫入標題列
    outFile << "chrom\tmethyl_pos\tsomatic_pos\tvariant_type\tvcf_source_id\tbam_source_id\t"
            << "somatic_allele_type\tsomatic_base_at_variant\thaplotype_tag\tmeth_call\t"
            << "meth_state\tstrand\tread_id\n";
    
    // 寫入數據
    for (const auto& detail : details) {
        outFile << detail.chrom << "\t"
                << detail.methyl_pos << "\t"
                << detail.somatic_pos << "\t"
                << detail.variant_type << "\t"
                << detail.vcf_source_id << "\t"
                << detail.bam_source_id << "\t"
                << detail.somatic_allele_type << "\t"
                << detail.somatic_base_at_variant << "\t"
                << detail.haplotype_tag << "\t"
                << std::fixed << std::setprecision(4) << detail.meth_call << "\t"
                << detail.meth_state << "\t"
                << detail.strand << "\t"
                << detail.read_id << "\n";
    }
    
    outFile.close();
    
    // 如果需要gzip壓縮
    if (config_.gzip_output) {
        std::string gzipPath = tempPath + ".gz";
        if (compressFile(tempPath, gzipPath)) {
            // 壓縮成功後刪除原檔案
            fs::remove(tempPath);
            LOG_INFO("ReportExporter", "已匯出Level 1原始甲基化詳情 (已壓縮): " + gzipPath);
        } else {
            LOG_ERROR("ReportExporter", "壓縮Level 1檔案失敗");
            return false;
        }
    } else {
        LOG_INFO("ReportExporter", "已匯出Level 1原始甲基化詳情: " + tempPath);
    }
    
    return true;
}

bool ReportExporter::exportLevel2Summary(const std::vector<msa::SomaticVariantMethylationSummary>& summaries, const std::string& outputDir) {
    std::string tempPath = outputDir + "/level2_somatic_variant_methylation_summary.tsv";
    std::ofstream outFile(tempPath);
    
    if (!outFile.is_open()) {
        LOG_ERROR("ReportExporter", "無法開啟輸出檔案: " + tempPath);
        return false;
    }
    
    // 寫入標題列
    outFile << "chrom\tsomatic_pos\tvariant_type\tvcf_source_id\tbam_source_id\t"
            << "somatic_allele_type\thaplotype_tag\tsupporting_read_count\t"
            << "methyl_sites_count\tmean_methylation\tstrand\n";
    
    // 寫入數據
    for (const auto& summary : summaries) {
        outFile << summary.chrom << "\t"
                << summary.somatic_pos << "\t"
                << summary.variant_type << "\t"
                << summary.vcf_source_id << "\t"
                << summary.bam_source_id << "\t"
                << summary.somatic_allele_type << "\t"
                << summary.haplotype_tag << "\t"
                << summary.supporting_read_count << "\t"
                << summary.methyl_sites_count << "\t"
                << std::fixed << std::setprecision(4) << summary.mean_methylation << "\t"
                << summary.strand << "\n";
    }
    
    outFile.close();
    
    // 如果需要gzip壓縮
    if (config_.gzip_output) {
        std::string gzipPath = tempPath + ".gz";
        if (compressFile(tempPath, gzipPath)) {
            // 壓縮成功後刪除原檔案
            fs::remove(tempPath);
            LOG_INFO("ReportExporter", "已匯出Level 2變異甲基化摘要 (已壓縮): " + gzipPath);
        } else {
            LOG_ERROR("ReportExporter", "壓縮Level 2檔案失敗");
            return false;
        }
    } else {
        LOG_INFO("ReportExporter", "已匯出Level 2變異甲基化摘要: " + tempPath);
    }
    
    return true;
}

bool ReportExporter::exportLevel3Stats(const std::vector<msa::AggregatedHaplotypeStats>& stats, const std::string& outputDir) {
    std::string outputPath = outputDir + "/level3_haplotype_group_statistics.tsv";
    std::ofstream outFile(outputPath);
    
    if (!outFile.is_open()) {
        LOG_ERROR("ReportExporter", "無法開啟輸出檔案: " + outputPath);
        return false;
    }
    
    // 獲取所有唯一的VCF來源ID
    std::set<std::string> vcf_sources;
    for (const auto& stat : stats) {
        for (const auto& [vcf_source, _] : stat.vcf_methylation_means) {
            vcf_sources.insert(vcf_source);
        }
    }
    
    // 寫入標題列
    outFile << "haplotype_group\tbam_source\tvariant_type_group";
    for (const auto& vcf_source : vcf_sources) {
        outFile << "\t" << vcf_source << "_mean_methylation";
    }
    if (vcf_sources.size() >= 2) {
        outFile << "\tdifference\tp_value";
    }
    outFile << "\n";
    
    // 寫入數據
    for (const auto& stat : stats) {
        outFile << stat.haplotype_group << "\t"
                << stat.bam_source << "\t"
                << stat.variant_type_group;
        
        for (const auto& vcf_source : vcf_sources) {
            outFile << "\t";
            if (stat.vcf_methylation_means.count(vcf_source) && 
                stat.vcf_methylation_means.at(vcf_source) >= 0) {
                outFile << std::fixed << std::setprecision(4) << stat.vcf_methylation_means.at(vcf_source);
            } else {
                outFile << "NA";
            }
        }
        
        if (vcf_sources.size() >= 2) {
            outFile << "\t" << std::fixed << std::setprecision(4) << stat.difference;
            outFile << "\t" << std::fixed << std::setprecision(6) << stat.p_value;
        }
        
        outFile << "\n";
    }
    
    outFile.close();
    LOG_INFO("ReportExporter", "已匯出Level 3單倍型統計: " + outputPath);
    return true;
}

bool ReportExporter::createDirectory(const std::string& dirPath) {
    try {
        fs::path path(dirPath);
        if (!fs::exists(path)) {
            return fs::create_directories(path);
        }
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("ReportExporter", "創建目錄時發生錯誤: " + std::string(e.what()));
        return false;
    }
}

bool ReportExporter::compressFile(const std::string& inputPath, const std::string& outputPath) {
    // 使用zlib進行gzip壓縮
    gzFile gzf = gzopen(outputPath.c_str(), "wb");
    if (!gzf) {
        LOG_ERROR("ReportExporter", "無法開啟gzip輸出檔案: " + outputPath);
        return false;
    }
    
    std::ifstream inFile(inputPath, std::ios::binary);
    if (!inFile) {
        gzclose(gzf);
        LOG_ERROR("ReportExporter", "無法開啟輸入檔案: " + inputPath);
        return false;
    }
    
    const int bufferSize = 8192;
    char buffer[bufferSize];
    
    while (inFile) {
        inFile.read(buffer, bufferSize);
        std::streamsize bytesRead = inFile.gcount();
        
        if (bytesRead > 0) {
            if (gzwrite(gzf, buffer, bytesRead) != bytesRead) {
                gzclose(gzf);
                inFile.close();
                LOG_ERROR("ReportExporter", "gzip寫入錯誤: " + outputPath);
                return false;
            }
        }
    }
    
    gzclose(gzf);
    inFile.close();
    return true;
}

} // namespace msa::core 