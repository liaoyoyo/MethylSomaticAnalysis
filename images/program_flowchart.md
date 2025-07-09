flowchart TD
    A[["啟動 MSA 程式"]] --> B["解析命令列參數<br/>ConfigParser"]
    B --> C["初始化日誌系統<br/>LogManager"]
    C --> D["驗證輸入參數"]
    D --> E["創建輸出目錄"]
    E --> F["驗證 BAM 檔案<br/>BAMValidator<br/>檢查甲基化和單倍型標籤"]
    F --> G["設置 OpenMP<br/>執行緒數"]
    G --> H["初始化記憶體池<br/>MemoryPool"]
    H --> I["開始並行處理<br/>VCF 檔案"]
    
    I --> J["載入變異<br/>VariantLoader"]
    J --> K["並行處理變異<br/>processVariantsInParallel"]
    
    K --> L["針對每個變異<br/>processVariant"]
    L --> M["使用 BamFetcher<br/>獲取變異周圍讀段"]
    M --> N["使用 MethylHaploExtractor<br/>提取甲基化位點"]
    N --> O["收集甲基化資料"]
    
    O --> P["體細胞甲基化分析<br/>SomaticMethylationAnalyzer"]
    P --> Q["生成分析結果<br/>AnalysisResults"]
    
    Q --> R["結果輸出<br/>ReportExporter"]
    R --> S["輸出 Level 1<br/>原始甲基化詳情"]
    R --> T["輸出 Level 2<br/>變異甲基化摘要"]
    R --> U["輸出 Level 3<br/>單倍型群組統計"]
    
    S --> V[["分析完成"]]
    T --> V
    U --> V
    
    subgraph "並行處理區域"
        direction TB
        K
        L
        M
        N
        O
    end
    
    subgraph "輸出檔案"
        direction LR
        S
        T
        U
    end
    
    style A fill:#e1f5fe
    style V fill:#c8e6c9
    style F fill:#fff3e0
    style P fill:#f3e5f5
    style R fill:#e8f5e8