#include "msa/utils/MemoryPool.h"
#include "msa/utils/LogManager.h"
#include <sstream>

namespace msa::utils {

// 單例實例獲取
MemoryPool& MemoryPool::getInstance() {
    static MemoryPool instance;
    return instance;
}

// 建構函數
MemoryPool::MemoryPool() : maxCapacity_(0), totalAllocated_(0), currentlyInUse_(0), initialized_(false) {
}

// 解構函數
MemoryPool::~MemoryPool() {
    releaseAll();
}

// 初始化記憶體池
void MemoryPool::initialize(size_t initialCapacity, size_t maxCapacity) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (initialized_) {
        LOG_WARN("MemoryPool", "記憶體池已初始化，忽略重複初始化");
        return;
    }
    
    maxCapacity_ = maxCapacity;
    
    // 預先分配指定數量的bam1_t物件
    for (size_t i = 0; i < initialCapacity; ++i) {
        availableBam1ts_.push(createNewBam1());
        totalAllocated_++;
    }
    
    initialized_ = true;
    
    std::ostringstream ss;
    ss << "記憶體池已初始化: 初始容量=" << initialCapacity 
       << ", 最大容量=" << (maxCapacity == 0 ? "無限制" : std::to_string(maxCapacity));
    LOG_INFO("MemoryPool", ss.str());
}

// 建立新的bam1_t物件
bam1_t* MemoryPool::createNewBam1() {
    bam1_t* b = bam_init1();
    if (!b) {
        throw std::bad_alloc();
    }
    return b;
}

// 獲取一個bam1_t物件
bam1_t* MemoryPool::getBam1(bool waitIfEmpty) {
    bam1_t* result = nullptr;
    
    {
        std::unique_lock<std::mutex> lock(mutex_);
        
        if (!initialized_) {
            // 如果尚未初始化，則自動初始化
            lock.unlock();
            initialize();
            lock.lock();
        }
        
        if (waitIfEmpty && availableBam1ts_.empty()) {
            // 如果設置為等待並且池中沒有可用物件，則阻塞直到有物件返回
            if (maxCapacity_ > 0 && totalAllocated_ >= maxCapacity_) {
                // 如果已達最大容量限制，等待直到有物件返回
                availableCondition_.wait(lock, [this]() {
                    return !availableBam1ts_.empty();
                });
            } else {
                // 如果未達最大容量，則建立新物件並返回
                result = createNewBam1();
                totalAllocated_++;
            }
        }
        
        if (!result) {  // 如果尚未創建新物件
            if (!availableBam1ts_.empty()) {
                // 從池中取出一個可用物件
                result = availableBam1ts_.front();
                availableBam1ts_.pop();
            } else {
                // 如果池為空，且未設置等待，則建立新物件
                result = createNewBam1();
                totalAllocated_++;
            }
        }
        
        // 更新使用計數
        currentlyInUse_++;
    }
    
    // 重置bam1_t物件的內容，確保沒有殘留數據
    if (result) {
        // 使用正確的htslib函數來重置，避免使用被棄用的函數
        if (result->data) {
            free(result->data);
            result->data = NULL;
        }
        result->l_data = 0;
        result->m_data = 0;
        result->core.pos = -1;
        result->core.mpos = -1;
        result->core.mtid = -1;
        result->core.tid = -1;
        result->core.flag = 0;
        result->core.n_cigar = 0;
        result->core.l_qname = 0;
        result->core.isize = 0;
        result->core.bin = 0;
        result->core.qual = 0;
    }
    
    return result;
}

// 返回一個bam1_t物件到池中
void MemoryPool::returnBam1(bam1_t* b) {
    if (!b) return;
    
    {
        std::lock_guard<std::mutex> lock(mutex_);
        
        // 將物件加入可用隊列
        availableBam1ts_.push(b);
        
        // 更新使用計數
        if (currentlyInUse_ > 0) {
            currentlyInUse_--;
        }
    }
    
    // 通知等待中的執行緒
    availableCondition_.notify_one();
}

// 釋放所有記憶體池中的物件
void MemoryPool::releaseAll() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // 釋放所有可用物件
    while (!availableBam1ts_.empty()) {
        bam1_t* b = availableBam1ts_.front();
        availableBam1ts_.pop();
        bam_destroy1(b);
    }
    
    // 輸出警告，如果仍有物件未返回
    if (currentlyInUse_ > 0) {
        std::ostringstream ss;
        ss << "記憶體池釋放時仍有 " << currentlyInUse_ << " 個物件未返回，可能導致記憶體洩漏";
        LOG_WARN("MemoryPool", ss.str());
    }
    
    // 重置計數
    totalAllocated_ = 0;
    currentlyInUse_ = 0;
    initialized_ = false;
    
    LOG_INFO("MemoryPool", "記憶體池已釋放所有物件");
}

// 獲取目前池中物件數量
size_t MemoryPool::size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return availableBam1ts_.size();
}

// 獲取統計資訊
std::string MemoryPool::getStats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::ostringstream ss;
    ss << "MemoryPool狀態: 總分配=" << totalAllocated_ 
       << ", 使用中=" << currentlyInUse_
       << ", 可用=" << availableBam1ts_.size()
       << ", 最大容量=" << (maxCapacity_ == 0 ? "無限制" : std::to_string(maxCapacity_));
    
    return ss.str();
}

} // namespace msa::utils 