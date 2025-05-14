#include "msa/core/ConfigParser.h"
#include "msa/utils/LogManager.h"
#include "msa/external/cxxopts.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <thread>

using namespace msa::utils;

namespace fs = std::filesystem;
namespace msa::core {

ConfigParser::ConfigParser() = default;

msa::Config ConfigParser::parse(int argc, char** argv) {
    msa::Config config;
    
    // 使用cxxopts解析命令列參數
    cxxopts::Options options("MethylSomaticAnalysis", "Somatic variant methylation analysis tool");
    
    // 必要參數
    options.add_options()
        ("v,vcfs", "Somatic VCF檔案路徑 (必要，可提供多個)", cxxopts::value<std::vector<std::string>>())
        ("r,ref", "參考基因組路徑 (必要)", cxxopts::value<std::string>())
        ("t,tumor", "腫瘤BAM檔案路徑 (必要)", cxxopts::value<std::string>())
        ("n,normal", "正常BAM檔案路徑 (必要)", cxxopts::value<std::string>());
    
    // 可選參數
    options.add_options()
        ("w,window", "變異點擷取區域半徑(bp)", cxxopts::value<int>()->default_value("2000"))
        ("b,bed", "限定分析區域的BED檔案", cxxopts::value<std::string>())
        ("meth-high", "高甲基閾值 (0.01-1.0)", cxxopts::value<float>()->default_value("0.8"))
        ("meth-low", "低甲基閾值 (0.01-1.0)", cxxopts::value<float>()->default_value("0.2"))
        ("min-allele", "每個變異至少需有此數量Tumor BAM支持ALT讀數", cxxopts::value<int>()->default_value("0"))
        ("min-strand-reads", "每個CpG位點在正反鏈上各自至少需要的支持讀數", cxxopts::value<int>()->default_value("1"))
        ("log-level", "日誌級別 (trace/debug/info/warn/error/fatal)", cxxopts::value<std::string>()->default_value("info"))
        ("j,threads", "執行緒數", cxxopts::value<int>()->default_value("0"))
        ("o,outdir", "輸出總路徑", cxxopts::value<std::string>()->default_value("./results"))
        ("gzip-output", "是否gzip壓縮Level 1 & 2 TSV輸出", cxxopts::value<std::string>()->default_value("true"))
        ("max-read-depth", "最大讀取深度", cxxopts::value<int>()->default_value("10000"))
        ("max-ram-gb", "最大RAM使用量(GB)", cxxopts::value<int>()->default_value("32"))
        ("h,help", "顯示使用說明");

    // 設置需要參數值的選項
    options.positional_help("必要參數: --vcfs --ref --tumor --normal");
    options.show_positional_help();
    
    try {
        auto result = options.parse(argc, argv);
        
        // 檢查是否要顯示使用說明
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        
        // 解析必要參數
        if (result.count("vcfs")) {
            config.vcf_files = result["vcfs"].as<std::vector<std::string>>();
        } else {
            throw std::runtime_error("必須提供VCF檔案路徑 (--vcfs)");
        }
        
        if (result.count("ref")) {
            config.ref_file = result["ref"].as<std::string>();
        } else {
            throw std::runtime_error("必須提供參考基因組路徑 (--ref)");
        }
        
        if (result.count("tumor")) {
            config.tumor_bam = result["tumor"].as<std::string>();
        } else {
            throw std::runtime_error("必須提供腫瘤BAM檔案路徑 (--tumor)");
        }
        
        if (result.count("normal")) {
            config.normal_bam = result["normal"].as<std::string>();
        } else {
            throw std::runtime_error("必須提供正常BAM檔案路徑 (--normal)");
        }
        
        // 解析可選參數
        if (result.count("window")) {
            config.window_size = result["window"].as<int>();
        }
        
        if (result.count("bed")) {
            config.bed_file = result["bed"].as<std::string>();
        }
        
        if (result.count("meth-high")) {
            config.meth_high_threshold = result["meth-high"].as<float>();
        }
        
        if (result.count("meth-low")) {
            config.meth_low_threshold = result["meth-low"].as<float>();
        }
        
        if (result.count("min-allele")) {
            config.min_allele = result["min-allele"].as<int>();
        }
        
        if (result.count("min-strand-reads")) {
            config.min_strand_reads = result["min-strand-reads"].as<int>();
        }
        
        if (result.count("log-level")) {
            config.log_level = result["log-level"].as<std::string>();
        }
        
        if (result.count("threads")) {
            config.threads = result["threads"].as<int>();
            if (config.threads <= 0) {
                // 如果未指定或指定為0，則使用系統可用核心數
                config.threads = std::thread::hardware_concurrency();
                if (config.threads == 0) {
                    // 如果無法檢測到可用核心數，則預設為4
                    config.threads = 4;
                }
            }
        }
        
        if (result.count("outdir")) {
            config.outdir = result["outdir"].as<std::string>();
        }
        
        if (result.count("gzip-output")) {
            std::string gzip_value = result["gzip-output"].as<std::string>();
            // 轉換為小寫以便比較
            std::transform(gzip_value.begin(), gzip_value.end(), gzip_value.begin(),
                          [](unsigned char c){ return std::tolower(c); });
            // 判斷字串是否表示false
            config.gzip_output = !(gzip_value == "false" || gzip_value == "0" || gzip_value == "no" || gzip_value == "n" || gzip_value == "off");
            LOG_INFO("ConfigParser", "設定輸出壓縮: " + std::string(config.gzip_output ? "是" : "否") + " (原始值: " + gzip_value + ")");
        }
        
        if (result.count("max-read-depth")) {
            config.max_read_depth = result["max-read-depth"].as<int>();
        }
        
        if (result.count("max-ram-gb")) {
            config.max_ram_gb = result["max-ram-gb"].as<int>();
        }
        
        // 驗證配置是否合法
        validateConfig(config);
        
        // 檢查檔案存在性與索引
        checkVcfFiles(config.vcf_files);
        checkRefFile(config.ref_file);
        checkBamFile(config.tumor_bam, "Tumor BAM");
        checkBamFile(config.normal_bam, "Normal BAM");
        if (!config.bed_file.empty()) {
            checkFileExists(config.bed_file, "BED檔案", false);
        }
        
        // 創建輸出目錄
        try {
            fs::path outDirPath(config.outdir);
            if (!fs::exists(outDirPath)) {
                fs::create_directories(outDirPath);
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("無法創建輸出目錄: " + config.outdir + ", 錯誤: " + e.what());
        }
        
    } catch (const std::exception& e) {
        throw std::runtime_error("解析參數錯誤: " + std::string(e.what()) + "\n" + getUsage());
    }
    
    return config;
}

void ConfigParser::validateConfig(msa::Config& config) {
    // 檢查變異點擷取區域半徑
    if (config.window_size <= 0 || config.window_size > 100000) {
        throw std::runtime_error("變異點擷取區域半徑必須在1-100000範圍內");
    }
    
    // 檢查甲基閾值
    if (config.meth_high_threshold <= 0.01 || config.meth_high_threshold > 1.0) {
        throw std::runtime_error("高甲基閾值必須在0.01-1.0範圍內");
    }
    
    if (config.meth_low_threshold <= 0.01 || config.meth_low_threshold >= config.meth_high_threshold) {
        throw std::runtime_error("低甲基閾值必須在0.01-" + std::to_string(config.meth_high_threshold) + "範圍內");
    }
    
    // 檢查min-allele
    if (config.min_allele < 0) {
        throw std::runtime_error("min-allele必須大於等於0");
    }
    
    // 檢查min-strand-reads
    if (config.min_strand_reads < 0) {
        throw std::runtime_error("min-strand-reads必須大於等於0");
    }
    
    // 檢查max-read-depth
    if (config.max_read_depth < 100 || config.max_read_depth > 1000000) {
        throw std::runtime_error("max-read-depth必須在100-1000000範圍內");
    }
    
    // 檢查max-ram-gb
    if (config.max_ram_gb < 1 || config.max_ram_gb > 1024) {
        throw std::runtime_error("max-ram-gb必須在1-1024範圍內");
    }
    
    // 檢查日誌級別
    std::string level_lower = config.log_level;
    std::transform(level_lower.begin(), level_lower.end(), level_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (level_lower != "trace" && level_lower != "debug" && level_lower != "info" &&
        level_lower != "warn" && level_lower != "error" && level_lower != "fatal") {
        throw std::runtime_error("無效的日誌級別，必須是trace/debug/info/warn/error/fatal之一");
    }
}

void ConfigParser::checkFileExists(const std::string& filePath, const std::string& fileType, bool required) {
    // 處理標準輸入/輸出的特殊情況
    if (filePath == "-") {
        // 允許通過標準輸入/輸出進行串流處理
        return;
    }
    
    if (required && filePath.empty()) {
        throw std::runtime_error(fileType + "路徑不能為空");
    }
    
    if (!filePath.empty() && filePath != "-") {
        std::ifstream file(filePath);
        if (!file.good()) {
            if (required) {
                throw std::runtime_error(fileType + "不存在或無法讀取: " + filePath);
            } else {
                LOG_WARN("ConfigParser", fileType + "不存在或無法讀取: " + filePath);
            }
        }
    }
}

void ConfigParser::checkVcfFiles(const std::vector<std::string>& vcfPaths) {
    for (const auto& vcfPath : vcfPaths) {
        // 檢查VCF檔案
        checkFileExists(vcfPath, "VCF檔案");
        
        // 檢查VCF索引檔案
        std::string tbiPath = vcfPath + ".tbi";
        checkFileExists(tbiPath, "VCF索引檔案(tbi)");
    }
}

void ConfigParser::checkBamFile(const std::string& bamPath, const std::string& bamType, bool isRequired) {
    // 檢查BAM檔案
    checkFileExists(bamPath, bamType, isRequired);
    
    // 如果是標準輸入，不檢查索引
    if (bamPath == "-") {
        return;
    }
    
    // 檢查BAM索引檔案
    std::string baiPath = bamPath + ".bai";
    std::string bamBaiPath = bamPath.substr(0, bamPath.size() - 4) + ".bai"; // 假設bamPath以.bam結尾
    
    bool baiExists = false;
    std::ifstream baiFile(baiPath);
    if (baiFile.good()) {
        baiExists = true;
    } else {
        std::ifstream bamBaiFile(bamBaiPath);
        if (bamBaiFile.good()) {
            baiExists = true;
        }
    }
    
    if (!baiExists && isRequired) {
        throw std::runtime_error(bamType + "索引檔案(bai)不存在或無法讀取: 嘗試了" + baiPath + "和" + bamBaiPath);
    }
}

void ConfigParser::checkRefFile(const std::string& refPath) {
    // 檢查參考基因組檔案
    checkFileExists(refPath, "參考基因組檔案");
    
    // 檢查參考基因組索引檔案
    std::string faiPath = refPath + ".fai";
    checkFileExists(faiPath, "參考基因組索引檔案(fai)");
}

std::string ConfigParser::getUsage() const {
    return "用法: MethylSomaticAnalysis --vcfs <vcf_file1> [<vcf_file2> ...] --ref <ref_file> "
           "--tumor <tumor_bam> --normal <normal_bam> [選項]\n"
           "使用 --help 參數查看完整說明";
}

std::string ConfigParser::getBasename(const std::string& filepath) {
    // 獲取檔案名稱（移除路徑）
    fs::path path(filepath);
    std::string filename = path.filename().string();
    
    // 移除副檔名
    size_t lastDot = filename.find_last_of('.');
    if (lastDot != std::string::npos) {
        // 檢查是否為.vcf.gz這樣的雙重副檔名
        if (filename.size() > 7 && filename.substr(filename.size() - 7) == ".vcf.gz") {
            return filename.substr(0, filename.size() - 7);
        }
        // 檢查是否為.bam這樣的單一副檔名
        if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".bam") {
            return filename.substr(0, filename.size() - 4);
        }
        // 一般情況，移除單一副檔名
        return filename.substr(0, lastDot);
    }
    
    return filename;
}

} // namespace msa::core 