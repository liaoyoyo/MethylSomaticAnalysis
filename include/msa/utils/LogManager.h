#pragma once

#include <string>
#include <mutex>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>

namespace msa::utils {

/**
 * @brief 日誌等級定義
 */
enum class LogLevel {
    TRACE_Level,  // 最詳細追蹤
    DEBUG_Level,  // 除錯資訊
    INFO_Level,   // 一般資訊
    WARN_Level,   // 警告
    ERROR_Level,  // 錯誤
    FATAL_Level   // 致命錯誤
};

/**
 * @brief 日誌管理類 - 使用單例模式
 */
class LogManager {
public:
    /**
     * @brief 獲取LogManager單例實例
     * @return LogManager& 單例引用
     */
    static LogManager& getInstance();

    /**
     * @brief 初始化日誌系統
     * @param logLevel 日誌級別
     * @param logFile 日誌檔案路徑，若為空則只輸出到標準錯誤
     */
    void initialize(LogLevel logLevel, const std::string& logFile = "");

    /**
     * @brief 記錄日誌
     * @param level 日誌級別
     * @param module 模組名稱
     * @param message 日誌訊息
     */
    void log(LogLevel level, const std::string& module, const std::string& message);

    /**
     * @brief 關閉日誌系統
     */
    void shutdown();

    /**
     * @brief 將字串轉換為LogLevel
     * @param levelStr 日誌級別字串
     * @return LogLevel 對應的LogLevel枚舉值
     */
    static LogLevel stringToLogLevel(const std::string& levelStr);

    /**
     * @brief 將LogLevel轉換為字串
     * @param level 日誌級別
     * @return std::string 對應的字串表示
     */
    static std::string logLevelToString(LogLevel level);

private:
    LogManager(); // 私有建構函數
    ~LogManager(); // 私有解構函數

    // 禁止複製與賦值
    LogManager(const LogManager&) = delete;
    LogManager& operator=(const LogManager&) = delete;

    LogLevel currentLevel_;          // 目前日誌級別
    std::string logFilePath_;        // 日誌檔案路徑
    std::ofstream logFileStream_;    // 日誌檔案串流
    std::recursive_mutex logMutex_;  // 用於多執行緒同步
    bool initialized_ = false;       // 初始化狀態

    /**
     * @brief 獲取當前時間字串，格式為YYYY-MM-DDTHH:MM:SSZ
     * @return std::string 時間字串
     */
    std::string getCurrentTimeString() const;
};

// 便利巨集，用於簡化日誌呼叫
#define LOG_TRACE(module, message) msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::TRACE_Level, module, message)
#define LOG_DEBUG(module, message) msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::DEBUG_Level, module, message)
#define LOG_INFO(module, message)  msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::INFO_Level, module, message)
#define LOG_WARN(module, message)  msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::WARN_Level, module, message)
#define LOG_ERROR(module, message) msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::ERROR_Level, module, message)
#define LOG_FATAL(module, message) msa::utils::LogManager::getInstance().log(msa::utils::LogLevel::FATAL_Level, module, message)

} // namespace msa::utils 