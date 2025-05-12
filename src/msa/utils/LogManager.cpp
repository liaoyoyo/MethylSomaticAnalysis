#include "msa/utils/LogManager.h"
#include <chrono>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <iostream>

namespace fs = std::filesystem;
namespace msa::utils {

// 單例實例獲取
LogManager& LogManager::getInstance() {
    static LogManager instance;
    return instance;
}

// 建構函數 - 預設初始化為 INFO 級別
LogManager::LogManager() : currentLevel_(LogLevel::INFO), initialized_(false) {
}

// 解構函數 - 關閉日誌檔案
LogManager::~LogManager() {
    shutdown();
}

// 初始化日誌系統
void LogManager::initialize(LogLevel logLevel, const std::string& logFile) {
    std::lock_guard<std::recursive_mutex> lock(logMutex_);
    
    // 設置日誌級別
    currentLevel_ = logLevel;
    
    // 如果指定了日誌檔案
    if (!logFile.empty()) {
        logFilePath_ = logFile;
        
        // 確保日誌目錄存在
        try {
            fs::path logPath(logFilePath_);
            fs::path logDir = logPath.parent_path();
            
            if (!logDir.empty() && !fs::exists(logDir)) {
                fs::create_directories(logDir);
            }
            
            // 開啟日誌檔案 (追加模式)
            logFileStream_.open(logFilePath_, std::ios::app);
            if (!logFileStream_.is_open()) {
                std::cerr << "[" << getCurrentTimeString() << "][ERROR][LogManager] "
                        << "無法開啟日誌檔案: " << logFilePath_ << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "[" << getCurrentTimeString() << "][ERROR][LogManager] "
                    << "建立日誌目錄錯誤: " << e.what() << std::endl;
        }
    }
    
    initialized_ = true;
    
    // 記錄初始化訊息
    log(LogLevel::INFO, "LogManager", "日誌系統已初始化，級別: " + logLevelToString(currentLevel_) +
        (logFilePath_.empty() ? "" : ", 檔案: " + logFilePath_));
}

// 關閉日誌系統
void LogManager::shutdown() {
    std::lock_guard<std::recursive_mutex> lock(logMutex_);
    if (initialized_) {
        if (logFileStream_.is_open()) {
            logFileStream_.close();
        }
        initialized_ = false;
    }
}

// 記錄日誌
void LogManager::log(LogLevel level, const std::string& module, const std::string& message) {
    // 檢查是否為可記錄的日誌級別
    if (level < currentLevel_) {
        return;
    }
    
    std::lock_guard<std::recursive_mutex> lock(logMutex_);
    
    // 格式化日誌訊息
    std::ostringstream logStream;
    logStream << "[" << getCurrentTimeString() << "]["
              << logLevelToString(level) << "]["
              << module << "] " << message;
    
    // 輸出到標準錯誤 (對於WARN以上級別)
    if (level >= LogLevel::WARN) {
        std::cerr << logStream.str() << std::endl;
    } else {
        // INFO和更低級別僅當初始化日誌系統時輸出到標準輸出
        if (initialized_) {
            std::cout << logStream.str() << std::endl;
        }
    }
    
    // 如果日誌檔案已開啟，輸出到檔案
    if (logFileStream_.is_open()) {
        logFileStream_ << logStream.str() << std::endl;
        logFileStream_.flush();  // 確保立即寫入檔案
    }
}

// 獲取當前時間的格式化字串
std::string LogManager::getCurrentTimeString() const {
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);
    auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    std::stringstream ss;
    ss << std::put_time(std::gmtime(&now_time_t), "%Y-%m-%dT%H:%M:%S");
    ss << '.' << std::setfill('0') << std::setw(3) << now_ms.count() << 'Z';
    
    return ss.str();
}

// 將字串轉換為LogLevel
LogLevel LogManager::stringToLogLevel(const std::string& levelStr) {
    std::string level_lower = levelStr;
    std::transform(level_lower.begin(), level_lower.end(), level_lower.begin(), 
        [](unsigned char c) { return std::tolower(c); });
    
    if (level_lower == "trace") return LogLevel::TRACE;
    if (level_lower == "debug") return LogLevel::DEBUG;
    if (level_lower == "info")  return LogLevel::INFO;
    if (level_lower == "warn")  return LogLevel::WARN;
    if (level_lower == "error") return LogLevel::ERROR;
    if (level_lower == "fatal") return LogLevel::FATAL;
    
    // 預設為 INFO
    return LogLevel::INFO;
}

// 將LogLevel轉換為字串
std::string LogManager::logLevelToString(LogLevel level) {
    switch(level) {
        case LogLevel::TRACE: return "TRACE";
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::INFO:  return "INFO";
        case LogLevel::WARN:  return "WARN";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::FATAL: return "FATAL";
        default:              return "UNKNOWN";
    }
}

} // namespace msa::utils 