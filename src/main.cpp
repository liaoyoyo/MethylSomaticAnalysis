#include "msa/Types.h"
#include "msa/utils/LogManager.h"
#include "msa/utils/MemoryPool.h"
#include "msa/core/ConfigParser.h"
#include "msa/core/BAMValidator.h"
#include "msa/core/VariantLoader.h"
#include "msa/core/BamFetcher.h"
#include "msa/core/MethylHaploExtractor.h"
#include "msa/core/SomaticMethylationAnalyzer.h"
#include "msa/core/ReportExporter.h"

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include <filesystem>
#include <map>
#include <memory>
#include <functional>

// 定義版本信息
#define MSA_VERSION "1.0.0"

namespace fs = std::filesystem;
using namespace msa::utils;
using namespace msa::core;

// 宣告全域變數
std::vector<std::thread> worker_threads;

/**
 * @brief 顯示程式版本信息
 */
void showVersion() {
    std::cout << "MethylSomaticAnalysis v" << MSA_VERSION << std::endl;
    std::cout << "甲基化體細胞變異分析工具" << std::endl;
    std::cout << "構建日期: " << __DATE__ << " " << __TIME__ << std::endl;
}

/**
 * @brief 主程式入口點
 */
int main(int argc, char** argv) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // 顯示歡迎信息
    showVersion();
    std::cout << "---------------------------------------" << std::endl;
    
    // 初始化日誌管理器
    msa::utils::LogManager::getInstance().initialize(msa::utils::LogLevel::TRACE, "msa.log");
    LOG_INFO("Main", "初始化MethylSomaticAnalysis...");
    
    // 解析命令列參數
    ConfigParser config_parser;
    msa::Config config;
    
    try {
        config = config_parser.parse(argc, argv);
    } catch (const std::exception& e) {
        LOG_ERROR("Main", "解析參數錯誤: " + std::string(e.what()));
        std::cout << config_parser.getUsage() << std::endl;
        return 1;
    }
    
    // 檢查必要參數
    if (config.vcf_files.empty()) {
        LOG_ERROR("Main", "錯誤: 至少需要提供一個VCF檔案");
        std::cout << config_parser.getUsage() << std::endl;
        return 1;
    }
    
    if (config.tumor_bam.empty() && config.normal_bam.empty()) {
        LOG_ERROR("Main", "錯誤: 至少需要提供腫瘤或對照樣本BAM檔案");
        std::cout << config_parser.getUsage() << std::endl;
        return 1;
    }
    
    if (config.ref_file.empty()) {
        LOG_ERROR("Main", "錯誤: 必須提供參考基因組檔案");
        std::cout << config_parser.getUsage() << std::endl;
        return 1;
    }
    
    // 創建輸出目錄
    try {
        fs::path outdir_path(config.outdir);
        if (!fs::exists(outdir_path)) {
            LOG_INFO("Main", "創建輸出目錄: " + config.outdir);
            fs::create_directories(outdir_path);
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Main", "無法創建輸出目錄: " + std::string(e.what()));
        return 1;
    }
    
    // 驗證BAM檔案的甲基化和單倍型標籤
    BAMValidator validator;
    if (!validator.checkAllInputFiles(config)) {
        LOG_ERROR("Main", "BAM檔案驗證失敗，請檢查輸入文件是否有甲基化標籤");
        return 1;
    }
    
    LOG_INFO("Main", "BAM檔案驗證成功");
    LOG_INFO("Main", "腫瘤BAM甲基化標籤: " + std::string(config.tumor_has_methyl_tags ? "存在" : "不存在"));
    LOG_INFO("Main", "對照BAM甲基化標籤: " + std::string(config.normal_has_methyl_tags ? "存在" : "不存在"));
    LOG_INFO("Main", "腫瘤BAM單倍型標籤: " + std::string(config.tumor_has_hp_tags ? "存在" : "不存在"));
    LOG_INFO("Main", "對照BAM單倍型標籤: " + std::string(config.normal_has_hp_tags ? "存在" : "不存在"));
    
    // 設置使用的執行緒數
    if (config.threads <= 0) {
        config.threads = std::thread::hardware_concurrency();
        if (config.threads == 0) config.threads = 1;  // 防止硬體偵測失敗
    }
    LOG_INFO("Main", "使用 " + std::to_string(config.threads) + " 個執行緒");
    
    // 初始化記憶體池
    auto& mem_pool = MemoryPool::getInstance();
    mem_pool.initialize(100);  // 預分配100個bam1_t物件
    
    // 為每個VCF檔案執行分析
    LOG_INFO("Main", "共有 " + std::to_string(config.vcf_files.size()) + " 個VCF檔案需要處理");
    for (size_t vcf_idx = 0; vcf_idx < config.vcf_files.size(); ++vcf_idx) {
        const auto& vcf_file = config.vcf_files[vcf_idx];
        std::string vcf_base_name = ConfigParser::getBasename(vcf_file);
        LOG_INFO("Main", "開始處理VCF檔案 [" + std::to_string(vcf_idx+1) + "/" + 
                 std::to_string(config.vcf_files.size()) + "]: " + vcf_file);
        
        // 載入變異信息
        VariantLoader variant_loader;
        std::vector<msa::VcfVariantInfo> variants = variant_loader.loadVCFs({vcf_file}, config.bed_file, config);
        
        if (variants.empty()) {
            LOG_WARN("Main", "VCF檔案 " + vcf_file + " 未載入任何變異，跳過此檔案");
            continue;
        }
        
        LOG_INFO("Main", "已載入 " + std::to_string(variants.size()) + " 個變異");
        
        // 初始化BAM讀取器
        BamFetcher bam_fetcher(config);
        if (!bam_fetcher.openBamFiles()) {
            LOG_ERROR("Main", "無法開啟BAM檔案，跳過VCF檔案 " + vcf_file);
            continue;
        }
        
        // 初始化甲基化/單倍型提取器
        MethylHaploExtractor meth_extractor(config);
        
        // 初始化結果容器
        std::vector<msa::MethylationSiteDetail> all_methyl_sites;
        
        // 處理每個變異
        for (const auto& variant : variants) {
            LOG_INFO("Main", "處理變異: " + variant.chrom + ":" + std::to_string(variant.pos) + " " + variant.variant_type);
            
            // 提取腫瘤樣本中的甲基化位點
            if (!config.tumor_bam.empty() && config.tumor_has_methyl_tags) {
                std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> tumor_reads = 
                    bam_fetcher.fetchReadsAroundVariant(variant, true, config.window_size);
                LOG_INFO("Main", "腫瘤樣本有 " + std::to_string(tumor_reads.size()) + " 個讀段覆蓋此變異");
                
                for (const auto& read : tumor_reads) {
                    auto methyl_sites = meth_extractor.extractFromRead(read.get(), variant, "tumor");
                    all_methyl_sites.insert(all_methyl_sites.end(), methyl_sites.begin(), methyl_sites.end());
                }
            }
            
            // 提取對照樣本中的甲基化位點
            if (!config.normal_bam.empty() && config.normal_has_methyl_tags) {
                std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> normal_reads = 
                    bam_fetcher.fetchReadsAroundVariant(variant, false, config.window_size);
                LOG_INFO("Main", "對照樣本有 " + std::to_string(normal_reads.size()) + " 個讀段覆蓋此變異");
                
                for (const auto& read : normal_reads) {
                    auto methyl_sites = meth_extractor.extractFromRead(read.get(), variant, "normal");
                    all_methyl_sites.insert(all_methyl_sites.end(), methyl_sites.begin(), methyl_sites.end());
                }
            }
        }
        
        LOG_INFO("Main", "共提取 " + std::to_string(all_methyl_sites.size()) + " 個甲基化位點");
        
        // 關閉BAM檔案
        bam_fetcher.closeBamFiles();
        
        // 甲基化分析
        SomaticMethylationAnalyzer analyzer(config);
        msa::AnalysisResults results = analyzer.analyze(all_methyl_sites);
        
        LOG_INFO("Main", "分析完成，生成摘要報告");
        
        // 匯出結果
        ReportExporter exporter(config);
        if (!exporter.exportResults(results, vcf_base_name)) {
            LOG_ERROR("Main", "匯出結果失敗");
        } else {
            LOG_INFO("Main", "已成功匯出結果到 " + config.outdir + "/" + vcf_base_name);
        }
        
        LOG_INFO("Main", "VCF檔案 [" + std::to_string(vcf_idx+1) + "/" + 
                 std::to_string(config.vcf_files.size()) + "] 處理完成: " + vcf_file);
    }
    
    // 計算運行時間
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    
    LOG_INFO("Main", "所有VCF檔案分析完成! 總運行時間: " + std::to_string(duration) + " 秒");
    LOG_INFO("Main", "結果保存在目錄: " + config.outdir);
    
    return 0;
} 