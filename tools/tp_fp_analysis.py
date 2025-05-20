#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def cohens_d(x, y):
    # compute Cohen's d
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx - 1) * x.std(ddof=1) ** 2 + (ny - 1) * y.std(ddof=1) ** 2) / dof)
    return (x.mean() - y.mean()) / pooled_std

def bootstrap_p(x, y, n_boot=1000):
    obs_diff = x.mean() - y.mean()
    pooled = np.concatenate((x, y))
    count = 0
    for _ in range(n_boot):
        bs_x = np.random.choice(pooled, size=len(x), replace=True)
        bs_y = np.random.choice(pooled, size=len(y), replace=True)
        if abs(bs_x.mean() - bs_y.mean()) >= abs(obs_diff):
            count += 1
    return count / n_boot

def plot_feature_comparison(tp_data, fp_data, feature, out_dir, group=None):
    plt.figure(figsize=(12, 6))
    
    # 準備數據
    tp_series = pd.Series(tp_data, name='TP')
    fp_series = pd.Series(fp_data, name='FP')
    
    # 箱型圖
    plt.subplot(1, 2, 1)
    data = pd.concat([tp_series, fp_series], axis=1)
    sns.boxplot(data=data)
    plt.title(f'Boxplot of {feature}')
    plt.ylabel(feature)
    
    # 小提琴圖
    plt.subplot(1, 2, 2)
    sns.violinplot(data=data)
    plt.title(f'Violin plot of {feature}')
    plt.ylabel(feature)
    
    plt.tight_layout()
    group_suffix = f"_{group}" if group else ""
    plt.savefig(os.path.join(out_dir, f'{feature}_comparison{group_suffix}.png'))
    plt.close()

def plot_feature_correlation(tp_df, fp_df, features, out_dir):
    # 確保所有特徵都存在
    common_features = [f for f in features if f in tp_df.columns and f in fp_df.columns]
    
    if not common_features:
        print("Warning: No common features found for correlation analysis")
        return
        
    # 計算相關係數矩陣
    tp_corr = tp_df[common_features].corr()
    fp_corr = fp_df[common_features].corr()
    
    # 繪製熱圖
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 2, 1)
    sns.heatmap(tp_corr, annot=True, cmap='coolwarm', center=0)
    plt.title('TP Feature Correlation')
    
    plt.subplot(1, 2, 2)
    sns.heatmap(fp_corr, annot=True, cmap='coolwarm', center=0)
    plt.title('FP Feature Correlation')
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'feature_correlation.png'))
    plt.close()

def main():
    p = argparse.ArgumentParser(description="TP vs FP feature comparison")
    p.add_argument('--tp', required=True)
    p.add_argument('--fp', required=True)
    p.add_argument('--feature-file', default=None,
                   help="CSV with list of features to test (one per row or column)")
    p.add_argument('--out-dir', required=True)
    p.add_argument('--group-by', default=None,
                   help="column name to stratify by; must exist in both TP and FP files")
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--n-boot', type=int, default=1000)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    # 讀取數據並進行預處理
    print("Loading data...")
    tp_df = pd.read_csv(args.tp, sep="\t")
    fp_df = pd.read_csv(args.fp, sep="\t")
    
    # 移除重複行
    tp_df = tp_df.drop_duplicates()
    fp_df = fp_df.drop_duplicates()
    
    print(f"TP data shape: {tp_df.shape}")
    print(f"FP data shape: {fp_df.shape}")

    # Check grouping
    if args.group_by:
        if args.group_by not in tp_df.columns or args.group_by not in fp_df.columns:
            print(f"Warning: group-by column '{args.group_by}' not found in both TP and FP. "
                  f"Proceeding without grouping.")
            args.group_by = None

    # Determine groups
    if args.group_by:
        groups = sorted(set(tp_df[args.group_by].dropna().unique()) |
                        set(fp_df[args.group_by].dropna().unique()))
    else:
        groups = [None]

    # Load features
    if args.feature_file and os.path.exists(args.feature_file):
        feats = pd.read_csv(args.feature_file, header=None).iloc[:,0].tolist()
    else:
        # auto-detect numeric columns excluding group-by
        feats = [c for c in tp_df.select_dtypes(include=[np.number]).columns
                 if c != args.group_by]
    
    print(f"Analyzing {len(feats)} features...")

    records = []
    for grp in groups:
        tp_sub = tp_df if grp is None else tp_df[tp_df[args.group_by] == grp]
        fp_sub = fp_df if grp is None else fp_df[fp_df[args.group_by] == grp]
        
        # 繪製特徵相關性圖
        if grp is None:
            plot_feature_correlation(tp_sub, fp_sub, feats, args.out_dir)
        
        for feat in feats:
            if feat not in tp_sub.columns or feat not in fp_sub.columns:
                continue
            x = tp_sub[feat].dropna().values
            y = fp_sub[feat].dropna().values
            if len(x) < 2 or len(y) < 2:
                continue
                
            # 繪製特徵比較圖
            plot_feature_comparison(x, y, feat, args.out_dir, grp)
            
            t_stat, p_val = ttest_ind(x, y, equal_var=False)
            d = cohens_d(x, y)
            p_boot = bootstrap_p(x, y, args.n_boot)
            records.append({'group': grp or 'ALL', 'feature': feat,
                            't_stat': t_stat, 'p_val': p_val,
                            'p_boot': p_boot, 'cohens_d': d})

    out_df = pd.DataFrame(records)
    # multiple testing correction
    out_df['p_adj'] = np.minimum(out_df['p_val'] * len(out_df), 1.0)
    out_df.to_csv(os.path.join(args.out_dir, 'tp_fp_feature_stats.csv'), index=False)
    print(f"Results saved to {args.out_dir}/tp_fp_feature_stats.csv")

if __name__ == '__main__':
    main()
