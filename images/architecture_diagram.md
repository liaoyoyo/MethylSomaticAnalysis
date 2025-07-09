graph TB
    subgraph "主程式層"
        direction TB
        MAIN["main.cpp<br/>主程式入口"]
    end
    
    subgraph "核心模組層 (src/msa/core/)"
        direction TB
        CONFIG["ConfigParser<br/>配置解析"]
        VALIDATOR["BAMValidator<br/>BAM檔案驗證"]
        FETCHER["BamFetcher<br/>BAM讀段獲取"]
        LOADER["VariantLoader<br/>變異載入"]
        EXTRACTOR["MethylHaploExtractor<br/>甲基化和單倍型提取"]
        ANALYZER["SomaticMethylationAnalyzer<br/>體細胞甲基化分析"]
        EXPORTER["ReportExporter<br/>結果輸出"]
    end
    
    subgraph "工具模組層 (src/msa/utils/)"
        direction TB
        LOGGER["LogManager<br/>日誌管理"]
        MEMORY["MemoryPool<br/>記憶體池管理"]
    end
    
    subgraph "資料類型層 (include/msa/)"
        direction TB
        TYPES["Types.h<br/>資料結構定義"]
        STRUCTURES["MethylationSiteDetail<br/>SomaticVariantMethylationSummary<br/>AggregatedHaplotypeStats<br/>AnalysisResults"]
    end
    
    subgraph "外部依賴"
        direction LR
        HTSLIB["htslib<br/>BAM/VCF處理"]
        BOOST["Boost<br/>記憶體管理"]
        FMT["fmtlib<br/>格式化輸出"]
        OPENMP["OpenMP<br/>並行處理"]
    end
    
    subgraph "輸出層"
        direction TB
        LEVEL1["Level 1 輸出<br/>原始甲基化詳情"]
        LEVEL2["Level 2 輸出<br/>變異甲基化摘要"]
        LEVEL3["Level 3 輸出<br/>單倍型群組統計"]
    end
    
    %% 主程式依賴
    MAIN --> CONFIG
    MAIN --> LOGGER
    MAIN --> VALIDATOR
    MAIN --> MEMORY
    MAIN --> LOADER
    MAIN --> FETCHER
    MAIN --> EXTRACTOR
    MAIN --> ANALYZER
    MAIN --> EXPORTER
    
    %% 核心模組間依賴
    CONFIG --> TYPES
    VALIDATOR --> HTSLIB
    FETCHER --> HTSLIB
    FETCHER --> MEMORY
    LOADER --> HTSLIB
    EXTRACTOR --> FETCHER
    EXTRACTOR --> TYPES
    ANALYZER --> EXTRACTOR
    ANALYZER --> TYPES
    EXPORTER --> ANALYZER
    EXPORTER --> TYPES
    
    %% 工具模組依賴
    LOGGER --> FMT
    MEMORY --> BOOST
    
    %% 資料流
    TYPES --> STRUCTURES
    ANALYZER --> LEVEL1
    ANALYZER --> LEVEL2
    ANALYZER --> LEVEL3
    
    %% OpenMP 並行支援
    MAIN -.-> OPENMP
    ANALYZER -.-> OPENMP
    
    %% 樣式設定
    style MAIN fill:#e3f2fd
    style CONFIG fill:#fff3e0
    style VALIDATOR fill:#fff3e0
    style FETCHER fill:#fff3e0
    style LOADER fill:#fff3e0
    style EXTRACTOR fill:#f3e5f5
    style ANALYZER fill:#f3e5f5
    style EXPORTER fill:#e8f5e8
    style LOGGER fill:#fce4ec
    style MEMORY fill:#fce4ec
    style TYPES fill:#f1f8e9
    style STRUCTURES fill:#f1f8e9
    style LEVEL1 fill:#e0f2f1
    style LEVEL2 fill:#e0f2f1
    style LEVEL3 fill:#e0f2f1