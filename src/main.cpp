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

// 使用預處理器檢查是否編譯時啟用了OpenMP
#ifdef HAVE_OPENMP
#include <omp.h>  // 添加OpenMP頭文件
#endif

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
 * @brief 初始化OpenMP環境
 */
void initializeOpenMP(const msa::Config& config) {
#ifdef HAVE_OPENMP
    // 設置執行緒數量
    omp_set_num_threads(config.threads);
    
    // 禁用嵌套並行和動態執行緒
    omp_set_nested(0);
    omp_set_dynamic(0);
    
    LOG_INFO("Main", "OpenMP初始化完成，執行緒數: " + std::to_string(config.threads));
#else
    LOG_WARN("Main", "OpenMP未啟用，將使用單執行緒模式運行");
#endif
}

/**
 * @brief 並行處理單個變異，提取甲基化位點，並返回甲基化位點列表
 * \param variant 變異資訊
 * \param bam_fetcher BAM檔案讀取器
 * \param meth_extractor 甲基化/單倍型提取器
 * \param config 配置
 * \return 甲基化位點列表
 */
std::vector<msa::MethylationSiteDetail> processVariant(
    const msa::VcfVariantInfo& variant,
    BamFetcher& bam_fetcher,
    MethylHaploExtractor& meth_extractor,
    const msa::Config& config) {
    
    std::vector<msa::MethylationSiteDetail> variant_sites;
    
    LOG_DEBUG("Main", "處理變異: " + variant.chrom + ":" + std::to_string(variant.pos) + " " + variant.variant_type);
    
    // 提取腫瘤樣本中的甲基化位點
    if (!config.tumor_bam.empty() && config.tumor_has_methyl_tags) {
        std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> tumor_reads = 
            bam_fetcher.fetchReadsAroundVariant(variant, true, config.window_size);
        LOG_DEBUG("Main", "腫瘤樣本有 " + std::to_string(tumor_reads.size()) + " 個讀段覆蓋此變異");
        
        for (const auto& read : tumor_reads) {
            auto methyl_sites = meth_extractor.extractFromRead(read.get(), variant, "tumor");
            variant_sites.insert(variant_sites.end(), methyl_sites.begin(), methyl_sites.end());
        }
    }
    
    // 提取對照樣本中的甲基化位點
    if (!config.normal_bam.empty() && config.normal_has_methyl_tags) {
        std::vector<std::unique_ptr<bam1_t, std::function<void(bam1_t*)>>> normal_reads = 
            bam_fetcher.fetchReadsAroundVariant(variant, false, config.window_size);
        LOG_DEBUG("Main", "對照樣本有 " + std::to_string(normal_reads.size()) + " 個讀段覆蓋此變異");
        
        for (const auto& read : normal_reads) {
            auto methyl_sites = meth_extractor.extractFromRead(read.get(), variant, "normal");
            variant_sites.insert(variant_sites.end(), methyl_sites.begin(), methyl_sites.end());
        }
    }
    
    return variant_sites;
}

/**
 * @brief 並行處理變異的塊
 */
std::vector<msa::MethylationSiteDetail> processVariantsInParallel(
    const std::vector<msa::VcfVariantInfo>& variants,
    const msa::Config& config) {
    
    std::vector<msa::MethylationSiteDetail> all_methyl_sites;
    std::vector<std::vector<msa::MethylationSiteDetail>> thread_results(variants.size());
    
    // 使用OpenMP並行處理變異
#ifdef HAVE_OPENMP
    #pragma omp parallel
    {
        // 每個執行緒獨立開啟BAM檔案
        BamFetcher bam_fetcher(config);
        if (!bam_fetcher.openBamFiles()) {
            LOG_ERROR("Main", "執行緒無法開啟BAM檔案，跳過處理");
        } else {
            // 初始化甲基化/單倍型提取器
            MethylHaploExtractor meth_extractor(config);
            
            // 分配變異給各執行緒
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < variants.size(); ++i) {
                // 如果BAM開啟成功，則處理此變異
                thread_results[i] = processVariant(variants[i], bam_fetcher, meth_extractor, config);
            }
            
            // 關閉BAM檔案
            bam_fetcher.closeBamFiles();
        }
    }
#else
    // 非OpenMP模式 - 順序處理每個變異
    BamFetcher bam_fetcher(config);
    if (!bam_fetcher.openBamFiles()) {
        LOG_ERROR("Main", "無法開啟BAM檔案，跳過處理");
        return all_methyl_sites;
    }
    
    // 初始化甲基化/單倍型提取器
    MethylHaploExtractor meth_extractor(config);
    
    // 處理每個變異
    for (size_t i = 0; i < variants.size(); ++i) {
        thread_results[i] = processVariant(variants[i], bam_fetcher, meth_extractor, config);
    }
    
    // 關閉BAM檔案
    bam_fetcher.closeBamFiles();
#endif
    
    // 合併所有執行緒的結果
    for (const auto& results : thread_results) {
        all_methyl_sites.insert(all_methyl_sites.end(), results.begin(), results.end());
    }
    
    return all_methyl_sites;
}

/**
 * @brief 主程式入口點
 * \param argc 命令列參數數量
 * \param argv 命令列參數數組
 * \return 執行狀態碼
 */
int main(int argc, char** argv) {
    // 記錄程式開始執行時間
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // 顯示歡迎信息
    showVersion();
    std::cout << "---------------------------------------" << std::endl;
    
    // 解析命令列參數
    ConfigParser config_parser;
    msa::Config config;
    
    try {
        config = config_parser.parse(argc, argv);
    } catch (const std::exception& e) {
        std::cout << "解析參數錯誤: " << e.what() << std::endl;
        std::cout << config_parser.getUsage() << std::endl;
        return 1;
    }
    
    // 初始化日誌管理器
    msa::utils::LogManager::getInstance().initialize(msa::utils::LogLevel::INFO_Level, config.log_file);
    LOG_INFO("Main", "初始化MethylSomaticAnalysis...");
    
    // 使用設定的日誌級別重新初始化日誌管理器
    msa::utils::LogManager::getInstance().shutdown();
    msa::utils::LogManager::getInstance().initialize(
        msa::utils::LogManager::stringToLogLevel(config.log_level), 
        config.log_file);
    LOG_INFO("Main", "日誌系統以級別 " + config.log_level + " 初始化");
    
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
    
    // 初始化OpenMP環境
    initializeOpenMP(config);
    
    // 初始化記憶體池
    auto& mem_pool = MemoryPool::getInstance();
    mem_pool.initialize(100 * config.threads);  // 考慮執行緒數量預分配更多bam1_t物件
    
    // 為每個VCF檔案執行分析
    LOG_INFO("Main", "共有 " + std::to_string(config.vcf_files.size()) + " 個VCF檔案需要處理");
    
    // 使用OpenMP並行處理VCF檔案
#ifdef HAVE_OPENMP
    #pragma omp parallel for schedule(dynamic) if(config.vcf_files.size() > 1)
#endif
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
        
        // 並行處理變異
        std::vector<msa::MethylationSiteDetail> all_methyl_sites = processVariantsInParallel(variants, config);
        
        LOG_INFO("Main", "共提取 " + std::to_string(all_methyl_sites.size()) + " 個甲基化位點");
        
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