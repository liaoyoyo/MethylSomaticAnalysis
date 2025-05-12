# MethylSomaticAnalysis (MSA)

MethylSomaticAnalysis 是一個用於分析甲基化與體細胞變異關聯性的高效工具，基於 C++17 和 htslib 實現。

## 功能特點

- 整合分析體細胞變異與周圍甲基化位點的關聯性
- 支援單倍型（Haplotype）特異性甲基化分析
- 多層次分析結果：原始甲基化資料、變異摘要、單倍型聚合統計
- 針對高通量定序資料優化的記憶體管理和多執行緒處理
- 支援多個 VCF 檔案的批次處理

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

- **未定義的 htslib 函數**：如遇到 `undefined reference to 'bam_aux2array'` 或 `undefined reference to 'bcf_get_alleles'` 等錯誤，請確保安裝的 htslib 版本 ≥ 1.17。可使用以下命令檢查：
  ```bash
  pkg-config --modversion htslib
  ```
  若版本過低或使用系統庫路徑不正確，建議手動編譯安裝最新版 htslib：
  ```bash
  git clone https://github.com/samtools/htslib.git
  cd htslib
  git checkout 1.17 # 或更新版本
  autoreconf -i
  ./configure --prefix=/usr/local
  make
  sudo make install
  ```
  然後重新配置 CMake 時指定 htslib 路徑：
  ```bash
  cmake -DCMAKE_BUILD_TYPE=Release -DHTSLIB_ROOT=/usr/local ..
  ```

### BAM 檔案問題

- 確保 BAM 檔案已建立索引（.bai 檔案）
- 檢查 BAM 檔案是否包含甲基化標籤（MM/ML）和單倍型標籤（HP/PS）
- 若出現錯誤，可使用 `--log-level debug` 獲取更詳細的診斷信息

### 記憶體使用量調整

- 大型基因組分析可使用 `--max-read-depth` 限制每個區域的讀取深度
- 使用 `--max-ram-gb` 限制最大記憶體使用量

### 效能優化

- 增加 `--threads` 參數可提高處理速度
- 使用 `--bed` 參數限制分析區域，減少運算量
- 若僅關注特定變異類型，可先過濾 VCF 檔案

## 開發與貢獻

本專案遵循 C++17 標準和模組化設計原則。若要貢獻：

1. Fork 專案並克隆至本地
2. 建立功能分支：`git checkout -b feature/your-feature`
3. 提交變更：`git commit -m 'Add some feature'`
4. 推送到分支：`git push origin feature/your-feature`
5. 提交 Pull Request

## 授權

本專案採用 GNU GPLv3 授權條款發布。詳情請參閱 [LICENSE](LICENSE) 檔案。

## 引用

若在學術研究中使用本工具，請引用：

```
MethylSomaticAnalysis: A tool for analyzing the association between DNA methylation and somatic variants. 
(2025). GitHub repository, https://github.com/liaoyoyo/MethylSomaticAnalysis
```
