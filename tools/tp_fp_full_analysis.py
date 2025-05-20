#!/usr/bin/env python3
# tp_fp_full_analysis.py

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    roc_auc_score, accuracy_score,
    classification_report, confusion_matrix,
    roc_curve
)

def parse_args():
    p = argparse.ArgumentParser(description="TP/FP 全面统计与分类分析")
    p.add_argument('--tp-summary',  required=True, help="TP level2 summary TSV")
    p.add_argument('--fp-summary',  required=True, help="FP level2 summary TSV")
    p.add_argument('--tp-detail',   required=True, help="TP level1 detail TSV")
    p.add_argument('--fp-detail',   required=True, help="FP level1 detail TSV")
    p.add_argument('--out-dir',     required=True, help="结果输出目录")
    p.add_argument('--alpha', type=float, default=0.05, help="显著性水平")
    return p.parse_args()

def load_summary(path):
    df = pd.read_csv(path, sep='\t')
    # 如果已经包含 normal/tumor/delta，直接返回
    if set(['normal','tumor','delta']).issubset(df.columns):
        return df

    # 否则：用 bam_source_id 上的 mean_methylation pivot 成 normal/tumor
    if 'bam_source_id' in df.columns and 'mean_methylation' in df.columns:
        key = ['chrom','somatic_pos','variant_type',
               'vcf_source_id','somatic_allele_type','haplotype_tag']
        pivot = (
            df.pivot_table(
                index=key,
                columns='bam_source_id',
                values='mean_methylation',
                aggfunc='mean'
            )
            .reset_index()
        )
        if 'normal' not in pivot.columns or 'tumor' not in pivot.columns:
            raise KeyError(f"pivot 后缺少 normal/tumor 列，得到: {pivot.columns.tolist()}")
        pivot['delta'] = pivot['tumor'] - pivot['normal']
        return pivot

    raise KeyError(f"Summary 文件既无 normal/tumor/delta，也无法 pivot：{df.columns.tolist()}")

def compute_detail_features(detail_tsv):
    cols = ['chrom','somatic_pos','variant_type','vcf_source_id',
            'somatic_allele_type','haplotype_tag','methyl_pos',
            'meth_call','read_id']
    df = pd.read_csv(detail_tsv, sep='\t', usecols=cols, low_memory=False)
    grp = ['chrom','somatic_pos','variant_type',
           'vcf_source_id','somatic_allele_type','haplotype_tag']
    feats = (
        df
        .groupby(grp)
        .agg(
            supporting_read_count = ('read_id', 'nunique'),
            methyl_sites_count    = ('methyl_pos', 'nunique'),
            mean_methylation      = ('meth_call', 'mean')
        )
        .reset_index()
    )
    return feats

def merge_all(tp_sum, fp_sum, tp_feats, fp_feats):
    # 确保 summary 中一定有这三列
    for col in ('normal','tumor','delta'):
        if col not in tp_sum.columns:
            raise KeyError(f"TP summary 缺少必要列: {col}")
        if col not in fp_sum.columns:
            raise KeyError(f"FP summary 缺少必要列: {col}")

    # 只用 summary 和 detail 都有的列来 join
    key = [
        'chrom',
        'somatic_pos',
        'variant_type',
        'vcf_source_id',
        'somatic_allele_type',
        'haplotype_tag'
    ]
    tp = pd.merge(tp_sum, tp_feats, on=key, how='inner')
    fp = pd.merge(fp_sum, fp_feats, on=key, how='inner')

    tp['label'] = 'TP'
    fp['label'] = 'FP'
    df = pd.concat([tp, fp], ignore_index=True)

    # 最后再确认合并后所有特征都在
    required = [
        'supporting_read_count',
        'methyl_sites_count',
        'mean_methylation',
        'normal',
        'tumor',
        'delta'
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"合并后缺少列：{missing}\n可用列：{list(df.columns)}")
    return df



def descriptive_stats(df, features, out_dir):
    """输出每个特征的分组描述性统计并保存表格与直方图"""
    os.makedirs(out_dir, exist_ok=True)
    stats = []
    for feat in features:
        for grp in ['TP','FP']:
            arr = df[df.label==grp][feat].dropna()
            stats.append({
                'feature': feat,
                'group': grp,
                'count': len(arr),
                'mean': arr.mean(),
                'median': arr.median(),
                'std': arr.std(ddof=1),
                'min': arr.min(),
                'max': arr.max()
            })
            # 直方图
            plt.figure(figsize=(4,3))
            sns.histplot(arr, bins=30, kde=True)
            plt.title(f"{feat} distribution ({grp})")
            plt.xlabel(feat)
            plt.ylabel("Count")
            plt.tight_layout()
            plt.savefig(f"{out_dir}/{feat}_{grp}_hist.png", dpi=150)
            plt.close()
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(f"{out_dir}/descriptive_stats.csv", index=False)
    print(f"[+] 已输出描述性统计到 {out_dir}/descriptive_stats.csv")

def stat_tests(df, features, alpha, out_dir):
    """做 t-test + Cohen's d，并生成箱线图"""
    os.makedirs(out_dir, exist_ok=True)
    tests = []
    for feat in features:
        a = df[df.label=='TP'][feat].dropna()
        b = df[df.label=='FP'][feat].dropna()
        t_stat, p_val = ttest_ind(a, b, equal_var=False)
        pooled_sd = np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2) / 2)
        cohen_d = (a.mean() - b.mean()) / pooled_sd
        tests.append({
            'feature': feat,
            'TP_mean': a.mean(),
            'FP_mean': b.mean(),
            't_stat': t_stat,
            'p_val': p_val,
            'cohen_d': cohen_d,
            'significant': p_val < alpha
        })
        # 箱线图
        plt.figure(figsize=(4,3))
        sns.boxplot(x='label', y=feat, data=df, palette=['C0','C1'])
        plt.title(f"{feat} TP vs FP\np={p_val:.3g}, d={cohen_d:.2f}")
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{feat}_boxplot.png", dpi=150)
        plt.close()
    pd.DataFrame(tests).to_csv(f"{out_dir}/stat_tests.csv", index=False)
    print(f"[+] 已输出差异检验结果到 {out_dir}/stat_tests.csv")

def correlation_heatmaps(df, features, out_dir):
    """绘制 TP/FP 两组的特征相关性热图"""
    os.makedirs(out_dir, exist_ok=True)
    for grp in ['TP','FP']:
        sub = df[df.label==grp][features]
        corr = sub.corr()
        plt.figure(figsize=(5,4))
        sns.heatmap(corr, annot=True, fmt=".2f", cmap='coolwarm', vmin=-1, vmax=1)
        plt.title(f"{grp} feature correlation")
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{grp}_corr_heatmap.png", dpi=150)
        plt.close()
    print(f"[+] 已输出相关性热图到 {out_dir}/corr_heatmaps/")

def classification_eval(df, features, out_dir):
    """用逻辑回归评估特征区分度：ROC 曲线、特征权重"""
    os.makedirs(out_dir, exist_ok=True)
    X = df[features]
    imputer = SimpleImputer(strategy='median')
    X = pd.DataFrame(imputer.fit_transform(X), columns=features)
    y = (df.label == 'TP').astype(int)
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    X_tr, X_te, y_tr, y_te = train_test_split(
        Xs, y, test_size=0.3, stratify=y, random_state=42
    )
    clf = LogisticRegression(class_weight='balanced', max_iter=1000)
    clf.fit(X_tr, y_tr)
    y_pred = clf.predict(X_te)
    y_prob = clf.predict_proba(X_te)[:,1]

    # 打印
    print("=== Classification Report ===")
    print("Accuracy:", accuracy_score(y_te, y_pred))
    print("AUC:", roc_auc_score(y_te, y_prob))
    print(classification_report(y_te, y_pred))
    print("Confusion matrix:\n", confusion_matrix(y_te, y_pred))

    # ROC
    fpr, tpr, _ = roc_curve(y_te, y_prob)
    plt.figure(figsize=(4,4))
    plt.plot(fpr, tpr, label=f"AUC={roc_auc_score(y_te,y_prob):.2f}")
    plt.plot([0,1],[0,1],'--', color='grey')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_dir}/roc_curve.png", dpi=150)
    plt.close()

    # 特征权重
    coef = pd.Series(clf.coef_[0], index=features).sort_values(key=abs, ascending=False)
    coef.to_csv(f"{out_dir}/feature_coefficients.csv", header=['coef'])
    print(f"[+] 已输出 ROC 图和特征权重到 {out_dir}/")

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("[*] 加载 summary...")
    tp_sum = load_summary(args.tp_summary)
    fp_sum = load_summary(args.fp_summary)

    print("[*] 计算 detail features...")
    tp_feats = compute_detail_features(args.tp_detail)
    fp_feats = compute_detail_features(args.fp_detail)

    print("[*] 合并 TP/FP 数据集...")
    df = merge_all(tp_sum, fp_sum, tp_feats, fp_feats)

    # 要分析的特征列表
    features = [
        'supporting_read_count',
        'methyl_sites_count',
        'mean_methylation',
        'normal',
        'tumor',
        'delta'
    ]

    # 1. 描述性统计 & 直方图
    descriptive_stats(df, features, os.path.join(args.out_dir, 'descriptive'))

    # 2. 组间差异检验 & 箱线图
    stat_tests(df, features, args.alpha, os.path.join(args.out_dir, 'stat_tests'))

    # 3. 相关性热图
    correlation_heatmaps(df, features, os.path.join(args.out_dir, 'corr_heatmaps'))

    # 4. 分类效果评估
    classification_eval(df, features, os.path.join(args.out_dir, 'classification'))

    print("[✓] 全面分析完成！所有结果保存在", args.out_dir)

if __name__ == '__main__':
    main()
