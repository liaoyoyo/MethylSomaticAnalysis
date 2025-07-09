# MethylSomaticAnalysis (MSA)

MethylSomaticAnalysis 是一個用於分析甲基化與體細胞變異關聯性的高效工具，基於 C++17 和 htslib 實現。

## 功能特點

- 整合分析體細胞變異與周圍甲基化位點的關聯性
- 支援單倍型（Haplotype）特異性甲基化分析
- 多層次分析結果：原始甲基化資料、變異摘要、單倍型聚合統計
- 針對高通量定序資料優化的記憶體管理和多執行緒處理
- 支援多個 VCF 檔案的批次處理

## 程式架構

MSA 採用模組化設計，各組件分工明確，具備良好的可擴展性和維護性：

![程式架構圖](images/architecture_diagram.md)

### 核心模組說明

- **主程式層**：負責整體流程控制和模組協調
- **核心模組層**：實現主要分析功能，包括資料載入、驗證、處理和輸出
- **工具模組層**：提供基礎設施支援，如記憶體管理和日誌系統
- **資料類型層**：定義標準化的資料結構，確保模組間資料一致性
- **外部依賴**：整合業界標準的高效能函式庫

## 系統需求

- Ubuntu 20.04+ 或相容 Linux 發行版
- C++17 相容的編譯器 (GCC 9.4.0+ 或 Clang 10+)
- htslib 1.17 或更高版本（處理 BAM/VCF 檔案）
- Boost 1.70.0+（用於記憶體管理）
- fmtlib 7.0.0+（用於日誌系統）
- CMake 3.16+ (建置系統)

## 安裝方式

### 1. 安裝相依套件

```bash
# 基本開發工具
sudo apt update
sudo apt install -y build-essential cmake git

# htslib 相依
sudo apt install -y libhts-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# Boost 相依
sudo apt install -y libboost-dev libboost-system-dev libboost-thread-dev

# fmtlib
sudo apt install -y libfmt-dev
```

### 2. 下載原始碼

```bash
git clone https://github.com/liaoyoyo/MethylSomaticAnalysis.git
cd MethylSomaticAnalysis
```

### 3. 使用 CMake 建置

```bash
# 創建並進入建置目錄
mkdir -p build && cd build

# 配置專案
cmake -DCMAKE_BUILD_TYPE=Release ..

# 編譯
make -j$(nproc)

# 安裝（可選）
sudo make install
```

## 程式執行流程

MSA 的完整執行流程如下圖所示，包含初始化、資料載入、並行處理和結果輸出四個主要階段：

![程式流程圖](images/program_flowchart.md)

### 流程說明

1. **初始化階段**：解析參數、設置日誌、驗證輸入檔案
2. **資料載入階段**：載入變異資訊、初始化記憶體池
3. **並行處理階段**：多執行緒處理每個變異，提取甲基化和單倍型資訊
4. **結果輸出階段**：生成三個層級的分析結果檔案

## 使用方法

### 基本用法

```bash
# 若已安裝至系統
msa --vcfs somatic.vcf.gz --ref reference.fa --tumor tumor.bam --normal normal.bam --outdir results

# 從建置目錄執行
./build/bin/msa --vcfs somatic.vcf.gz --ref reference.fa --tumor tumor.bam --normal normal.bam --outdir results
```

### 必要參數

| 參數 | 說明 |
|------|------|
| `--vcfs`, `-v` | Somatic VCF 檔案路徑（需有 .tbi 索引檔，可提供多個檔案） |
| `--ref`, `-r` | 參考基因組檔案路徑（需有 .fai 索引） |
| `--tumor`, `-t` | 腫瘤樣本 BAM 檔案路徑 |
| `--normal`, `-n` | 正常樣本 BAM 檔案路徑 |

### 常用選項

| 參數 | 預設值 | 說明 |
|------|--------|------|
| `--window`, `-w` | 2000 | 變異點擷取區域半徑 (bp) |
| `--bed`, `-b` | | 限定分析區域的 BED 檔案（可選） |
| `--meth-high` | 0.8 | 高甲基閾值 (0.01-1.0) |
| `--meth-low` | 0.2 | 低甲基閾值 (0.01-1.0) |
| `--min-allele`, `-a` | 0 | 每個變異至少需有此數量腫瘤 BAM 支持 ALT 讀數 |
| `--min-strand-reads` | 1 | 每個 CpG 位點在正反鏈上各自至少需要的支持讀數 |
| `--threads`, `-j` | [自動] | 使用的執行緒數 |
| `--outdir`, `-o` | ./results | 輸出目錄 |
| `--gzip-output` | true | 是否壓縮 Level 1/2 輸出 |
| `--log-level` | info | 日誌詳細程度 (trace/debug/info/warn/error/fatal) |

### 高級選項

| 參數 | 預設值 | 說明 |
|------|--------|------|
| `--max-read-depth` | 10000 | 每個區域最大讀取深度 |
| `--max-ram-gb` | 32 | 最大記憶體使用量 (GB) |

## 範例

### 單一 VCF 分析

```bash
MethylSomaticAnalysis \
  --vcfs /path/to/somatic_variants.vcf.gz \
  --ref /path/to/hg38.fa \
  --tumor /path/to/tumor.bam \
  --normal /path/to/normal.bam \
  --window 2000 \
  --min-strand-reads 2 \
  --threads 16 \
  --outdir ./results
```

### 多 VCF 比較分析

```bash
MethylSomaticAnalysis \
  --vcfs /path/to/tp.vcf.gz /path/to/fp.vcf.gz \
  --ref /path/to/hg38.fa \
  --tumor /path/to/tumor.bam \
  --normal /path/to/normal.bam \
  --bed /path/to/regions.bed \
  --meth-high 0.8 \
  --meth-low 0.2 \
  --min-strand-reads 2 \
  --log-level debug \
  --threads 16 \
  --outdir ./analysis_results
```

### 使用標準輸入串流處理

```bash
# 從標準輸入讀取 BAM
samtools view -h /path/to/tumor.bam chr1:1000-2000 | \
MethylSomaticAnalysis \
  --vcfs /path/to/chr1_variants.vcf.gz \
  --ref /path/to/hg38.fa \
  --tumor - \
  --normal /path/to/normal.bam \
  --outdir ./streamed_results
```

## 輸出檔案結構

每個輸入的 VCF 檔案會在輸出目錄中創建一個子目錄，包含以下檔案：

```
results/
└── {vcf_basename}/
    ├── global_summary_metrics.tsv          # 全域參數與甲基化摘要
    ├── level1_raw_methylation_details.tsv.gz  # 原始甲基化位點資料
    ├── level2_somatic_variant_methylation_summary.tsv.gz  # 每個變異周圍甲基化統計
    └── level3_haplotype_group_statistics.tsv  # 單倍型群組統計比較
```

### 輸出檔案說明

1. **global_summary_metrics.tsv**：包含執行參數和全域統計資訊
2. **level1_raw_methylation_details.tsv.gz**：每個甲基化位點的詳細資訊，包括染色體位置、甲基化比例、單倍型標籤等
3. **level2_somatic_variant_methylation_summary.tsv.gz**：針對每個體細胞變異位點，彙總周圍甲基化位點的統計數據
4. **level3_haplotype_group_statistics.tsv**：按單倍型和變異類型分組的聚合統計和比較

## 常見問題排解

### 編譯錯誤

**錯誤：找不到 htslib**
```bash
# 確認 htslib 安裝
pkg-config --cflags --libs hts

# 如果需要手動指定路徑
cmake -DCMAKE_BUILD_TYPE=Release -DHTS_ROOT=/path/to/hts ..
```

**錯誤：C++17 不支援**
```bash
# 檢查編譯器版本
g++ --version

# 更新 GCC
sudo apt install gcc-11 g++-11
export CC=gcc-11 CXX=g++-11
```

### 執行錯誤

**錯誤：BAM 檔案沒有甲基化標籤**
- 確認 BAM 檔案包含 MM 和 ML 標籤（ONT 資料）或 XM 標籤（Bismark 資料）
- 使用 `samtools view` 檢查 BAM 檔案格式

**錯誤：記憶體不足**
- 減少 `--max-ram-gb` 設定
- 減少 `--threads` 數量
- 使用 `--bed` 限制分析區域

**錯誤：索引檔案遺失**
```bash
# 建立 BAM 索引
samtools index input.bam

# 建立 VCF 索引
tabix -p vcf input.vcf.gz

# 建立參考基因組索引
samtools faidx reference.fa
```

## 效能調校

### 記憶體使用最佳化

1. **調整記憶體池大小**：使用 `--max-ram-gb` 根據系統記憶體調整
2. **限制分析區域**：使用 `--bed` 檔案只分析感興趣的區域
3. **調整讀取深度限制**：使用 `--max-read-depth` 避免高深度區域消耗過多記憶體

### 多執行緒效能

1. **最佳執行緒數**：通常設為 CPU 核心數的 0.8-1.0 倍
2. **I/O 密集任務**：對於網路儲存的 BAM 檔案，可能需要減少執行緒數
3. **記憶體頻寬**：高執行緒數可能受限於記憶體頻寬

### 建議的系統配置

| 資料大小 | RAM | CPU 核心 | 儲存 | 預估時間 |
|----------|-----|----------|------|----------|
| < 10GB BAM | 16GB | 8-16 | SSD | 1-2 小時 |
| 10-50GB BAM | 32GB | 16-32 | SSD | 2-6 小時 |
| > 50GB BAM | 64GB+ | 32+ | NVMe SSD | 6+ 小時 |

## 開發指南

### 程式碼結構

- `src/main.cpp`：主程式入口
- `src/msa/core/`：核心分析模組
- `src/msa/utils/`：工具函式和輔助類別
- `include/msa/`：標頭檔案和型別定義
- `tests/`：單元測試
- `scripts/`：建置和部署腳本

### 編碼標準

- 遵循 C++17 標準
- 使用 Google C++ Style Guide
- 所有公開 API 需要 Doxygen 註解
- 關鍵演算法需要單元測試覆蓋

### 貢獻指南

1. Fork 此儲存庫
2. 創建功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交變更 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 開啟 Pull Request

## 授權條款

此專案採用 MIT 授權條款。詳細資訊請參閱 [LICENSE](LICENSE) 檔案。

## 引用

如果您在研究中使用了 MethylSomaticAnalysis，請引用：

```
[待補充論文引用資訊]
```

## 聯絡資訊

- 作者：[您的姓名]
- Email：[您的Email]
- 專案首頁：https://github.com/liaoyoyo/MethylSomaticAnalysis
- 問題回報：https://github.com/liaoyoyo/MethylSomaticAnalysis/issues

## 更新日誌

### v1.0.0 (待發布)
- 初始版本發布
- 支援 ONT 和 Bismark 甲基化資料
- 多層次分析輸出
- 完整的並行處理支援