#!/usr/bin/env python3
# threshold_grid_analysis.py

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    p = argparse.ArgumentParser(
        description="多参数阈值网格筛选：normal + delta + supporting_read_count_detail"
    )
    p.add_argument('--tp-summary',  required=True,
                   help="TP level2 summary TSV")
    p.add_argument('--fp-summary',  required=True,
                   help="FP level2 summary TSV")
    p.add_argument('--tp-detail',   required=True,
                   help="TP level1 raw methylation detail TSV")
    p.add_argument('--fp-detail',   required=True,
                   help="FP level1 raw methylation detail TSV")
    p.add_argument('--out-dir',     required=True,
                   help="结果输出目录")
    p.add_argument('--norm-ths', nargs='+', type=float,
                   default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
                   help="normal_mean 阈值列表")
    p.add_argument('--delta-ths', nargs='+', type=float,
                   default=[0.0, 0.05, 0.1, 0.15, 0.2],
                   help="delta 阈值列表")
    p.add_argument('--read-ths', nargs='+', type=int,
                   default=[0,2,5,10],
                   help="supporting_read_count_detail 阈值列表")
    return p.parse_args()

def load_summary(path):
    df = pd.read_csv(path, sep='\t')
    # 如果没有 normal/tumor/delta，就 pivot 出来
    if not {'normal','tumor','delta'}.issubset(df.columns):
        key = ['chrom','somatic_pos','variant_type',
               'vcf_source_id','somatic_allele_type','haplotype_tag']
        df = (
            df
            .pivot_table(
                index=key,
                columns='bam_source_id',
                values='mean_methylation',
                aggfunc='mean'
            )
            .reset_index()
            .rename_axis(None, axis=1)
        )
        df['delta'] = df['tumor'] - df['normal']
    return df

def compute_detail_features(detail_tsv):
    cols = ['chrom','somatic_pos','variant_type','vcf_source_id',
            'somatic_allele_type','haplotype_tag',
            'methyl_pos','meth_call','read_id']
    df = pd.read_csv(detail_tsv, sep='\t', usecols=cols, low_memory=False)
    grp = ['chrom','somatic_pos','variant_type',
           'vcf_source_id','somatic_allele_type','haplotype_tag']
    feats = (
        df
        .groupby(grp)
        .agg(
            read_count_detail   = ('read_id','nunique'),
            site_count_detail   = ('methyl_pos','nunique'),
            mean_meth_detail    = ('meth_call','mean')
        )
        .reset_index()
    )
    return feats

def merge_tp_fp(tp_sum, fp_sum, tp_det, fp_det):
    key = ['chrom','somatic_pos','variant_type',
           'vcf_source_id','somatic_allele_type','haplotype_tag']
    tp = pd.merge(tp_sum, tp_det, on=key, how='inner',
                  suffixes=('_summary','_detail'))
    fp = pd.merge(fp_sum, fp_det, on=key, how='inner',
                  suffixes=('_summary','_detail'))
    tp['label'] = 'TP'
    fp['label'] = 'FP'
    return pd.concat([tp, fp], ignore_index=True)

def grid_evaluate(df, norm_ths, delta_ths, read_ths, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    tpr_grids = {}
    fpr_grids = {}

    for r in read_ths:
        tpr = np.zeros((len(norm_ths), len(delta_ths)))
        fpr = np.zeros_like(tpr)
        for i, n_th in enumerate(norm_ths):
            for j, d_th in enumerate(delta_ths):
                sub = df[
                    (df['normal'] >= n_th) &
                    (df['delta'].abs() >= d_th) &
                    (df['read_count_detail'] >= r)
                ]
                tp_keep = (sub.label=='TP').sum()
                fp_keep = (sub.label=='FP').sum()
                # 原始总数
                tp_tot = (df.label=='TP').sum()
                fp_tot = (df.label=='FP').sum()
                tpr[i,j] = tp_keep / tp_tot
                fpr[i,j] = fp_keep / fp_tot

        tpr_grids[r] = tpr
        fpr_grids[r] = fpr

        # 绘制热图
        for name, grid in [('TPR', tpr), ('FPR', fpr)]:
            plt.figure(figsize=(5,4))
            sns.heatmap(
                grid,
                xticklabels=delta_ths,
                yticklabels=norm_ths,
                annot=True, fmt=".3f", cmap='viridis'
            )
            plt.xlabel('delta_th')
            plt.ylabel('normal_th')
            plt.title(f"{name} @ read≥{r}")
            plt.tight_layout()
            plt.savefig(f"{out_dir}/{name}_reads{r}.png", dpi=150)
            plt.close()

    # 将所有结果存成一个 CSV
    records = []
    for r in read_ths:
        for i, n_th in enumerate(norm_ths):
            for j, d_th in enumerate(delta_ths):
                records.append({
                    'read_th': r,
                    'normal_th': n_th,
                    'delta_th': d_th,
                    'TPR': tpr_grids[r][i,j],
                    'FPR': fpr_grids[r][i,j]
                })
    pd.DataFrame(records).to_csv(f"{out_dir}/threshold_grid_results.csv", index=False)
    print("已输出阈值网格结果到", out_dir)

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. 读 summary
    tp_sum = load_summary(args.tp_summary)
    fp_sum = load_summary(args.fp_summary)

    # 2. 读 detail & 计算特征
    tp_det = compute_detail_features(args.tp_detail)
    fp_det = compute_detail_features(args.fp_detail)

    # 3. 合并
    df = merge_tp_fp(tp_sum, fp_sum, tp_det, fp_det)

    # 4. 网格评估
    grid_evaluate(
        df,
        norm_ths  = args.norm_ths,
        delta_ths = args.delta_ths,
        read_ths  = args.read_ths,
        out_dir   = args.out_dir
    )

if __name__ == '__main__':
    main()
