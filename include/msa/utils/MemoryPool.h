#pragma once

#include <memory>
#include <mutex>
#include <vector>
#include <queue>
#include <htslib/sam.h>
#include <thread>
#include <condition_variable>
#include <atomic>

namespace msa::utils {

/**
 * @brief 記憶體池，用於管理bam1_t物件的生命週期，減少頻繁分配與釋放
 */
class MemoryPool {
public:
    /**
     * @brief 獲取MemoryPool單例實例
     * @return MemoryPool& 單例引用
     */
    static MemoryPool& getInstance();

    /**
     * @brief 初始化記憶體池
     * @param initialCapacity 初始容量
     * @param maxCapacity 最大容量，0表示無限制
     */
    void initialize(size_t initialCapacity = 1000, size_t maxCapacity = 0);

    /**
     * @brief 從記憶體池中獲取一個bam1_t物件
     * @param waitIfEmpty 若為true且池為空時，等待而非動態分配
     * @return bam1_t* 已分配的bam1_t物件指標
     */
    bam1_t* getBam1(bool waitIfEmpty = false);

    /**
     * @brief 將bam1_t物件返回記憶體池
     * @param b bam1_t物件指標
     */
    void returnBam1(bam1_t* b);

    /**
     * @brief 釋放所有記憶體池中的物件
     */
    void releaseAll();

    /**
     * @brief 獲取目前池中物件數量
     * @return size_t 物件數量
     */
    size_t size() const;

    /**
     * @brief 獲取記憶體池當前統計資訊
     * @return std::string 統計信息字串
     */
    std::string getStats() const;

private:
    MemoryPool(); // 私有建構函數
    ~MemoryPool(); // 私有解構函數

    // 禁止複製與賦值
    MemoryPool(const MemoryPool&) = delete;
    MemoryPool& operator=(const MemoryPool&) = delete;

    // 自定義刪除器，釋放bam1_t物件
    struct Bam1tDeleter {
        void operator()(bam1_t* b) const {
            if (b) {
                bam_destroy1(b);
            }
        }
    };

    mutable std::mutex mutex_;                  // 互斥鎖保護共享資源
    std::condition_variable availableCondition_; // 用於阻塞等待
    std::queue<bam1_t*> availableBam1ts_;       // 可用的bam1_t物件隊列
    size_t maxCapacity_;                        // 最大容量
    std::atomic<size_t> totalAllocated_;        // 已分配物件總數
    std::atomic<size_t> currentlyInUse_;        // 當前使用中的物件數
    bool initialized_ = false;                  // 初始化狀態

    /**
     * @brief 建立一個新的bam1_t物件
     * @return bam1_t* 新建的bam1_t物件指標
     */
    bam1_t* createNewBam1();
};

} // namespace msa::utils 