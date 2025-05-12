#pragma once

#include <string>
#include <vector>
#include "msa/Types.h"

namespace msa::core {

/**
 * @brief 命令列參數解析器
 */
class ConfigParser {
public:
    /**
     * @brief 建構函數
     */
    ConfigParser();

    /**
     * @brief 解析命令列參數
     * @param argc 參數數量
     * @param argv 參數陣列
     * @return msa::Config 解析後的配置物件
     * @throws std::runtime_error 若有參數錯誤或必要參數缺失
     */
    msa::Config parse(int argc, char** argv);

    /**
     * @brief 獲取使用說明
     * @return std::string 使用說明文字
     */
    std::string getUsage() const;

    /**
     * @brief 檢查配置合法性
     * @param config 待檢查的配置
     * @throws std::runtime_error 若有非法配置
     */
    void validateConfig(msa::Config& config);

    /**
     * @brief 從檔案路徑獲取基礎名稱 (移除目錄與副檔名)
     * @param filepath 檔案路徑
     * @return std::string 基礎名稱
     */
    static std::string getBasename(const std::string& filepath);

private:
    /**
     * @brief 確認文件存在與可讀取
     * @param filePath 檔案路徑
     * @param fileType 檔案類型描述（用於錯誤訊息）
     * @param required 是否為必要檔案，若為false則檔案不存在不會拋出例外
     * @throws std::runtime_error 若檔案不存在或不可讀取
     */
    void checkFileExists(const std::string& filePath, const std::string& fileType, bool required = true);

    /**
     * @brief 確認每個VCF檔案及其索引均存在
     * @param vcfPaths VCF檔案路徑列表
     * @throws std::runtime_error 若有VCF檔案或其索引不存在或不可讀取
     */
    void checkVcfFiles(const std::vector<std::string>& vcfPaths);

    /**
     * @brief 確認BAM檔案及其索引均存在
     * @param bamPath BAM檔案路徑
     * @param bamType BAM類型描述（用於錯誤訊息）
     * @param isRequired 是否為必要檔案
     * @throws std::runtime_error 若BAM檔案或其索引不存在或不可讀取
     */
    void checkBamFile(const std::string& bamPath, const std::string& bamType, bool isRequired = true);

    /**
     * @brief 確認參考基因組檔案及其索引均存在
     * @param refPath 參考基因組檔案路徑
     * @throws std::runtime_error 若參考基因組檔案或其索引不存在或不可讀取
     */
    void checkRefFile(const std::string& refPath);
};

} // namespace msa::core 