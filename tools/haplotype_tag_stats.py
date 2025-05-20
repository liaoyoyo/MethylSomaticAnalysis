#!/usr/bin/env python3
# tools/haplotype_tag_stats.py

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute distribution of haplotype_tag & ASM modes in TP/FP and perform chi-square tests"
    )
    parser.add_argument('--tp-summary',    required=True, help="TP level2 summary TSV file")
    parser.add_argument('--fp-summary',    required=True, help="FP level2 summary TSV file")
    parser.add_argument('--out-prefix',    required=True, help="Output filename prefix (without extension)")
    parser.add_argument('--low-threshold',  type=float, default=0.2, help="Low methylation threshold (default: 0.2)")
    parser.add_argument('--high-threshold', type=float, default=0.8, help="High methylation threshold (default: 0.8)")
    return parser.parse_args()

def load_and_pivot(path):
    df = pd.read_csv(path, sep='\t')

    # If 'normal', 'tumor', 'delta' already exist, use them directly
    if set(['normal','tumor','delta']).issubset(df.columns):
        result = df.copy()
    else:
        keys = ['chrom','somatic_pos','variant_type',
                'vcf_source_id','somatic_allele_type','haplotype_tag']
        pivot = (
            df
            .pivot_table(
                index=keys,
                columns='bam_source_id',    # expects 'normal' / 'tumor'
                values='mean_methylation',
                aggfunc='mean'
            )
            .reset_index()
        )
        # Check for both columns
        for col in ('normal','tumor'):
            if col not in pivot.columns:
                raise KeyError(f"Missing column `{col}` after pivot; found: {pivot.columns.tolist()}")
        pivot['delta'] = pivot['tumor'] - pivot['normal']
        result = pivot

    # Ensure numeric types
    result['normal'] = pd.to_numeric(result['normal'], errors='coerce')
    result['tumor']  = pd.to_numeric(result['tumor'],  errors='coerce')

    # Print distribution summary
    name = os.path.basename(path)
    print(f"\n[{name}] normal distribution:\n{result['normal'].describe()}")
    print(f"[{name}] tumor  distribution:\n{result['tumor'].describe()}\n")

    return result

def asm_mode(row, low_thr, high_thr):
    n, t = row['normal'], row['tumor']
    # Handle missing values separately
    if pd.isna(n) and pd.isna(t):
        return 'miss→miss'
    elif pd.isna(n):
        # normal missing, classify tumor side
        if t > high_thr:
            return 'miss→high'
        elif t < low_thr:
            return 'miss→low'
        else:
            return 'miss→mid'
    elif pd.isna(t):
        # tumor missing, classify normal side
        if n > high_thr:
            return 'high→miss'
        elif n < low_thr:
            return 'low→miss'
        else:
            return 'mid→miss'
    # Both values present: determine ASM mode based on thresholds
    if n < low_thr and t > high_thr:
        return 'low→high'
    elif n > high_thr and t < low_thr:
        return 'high→low'
    elif n < low_thr and t < low_thr:
        return 'low→low'
    elif n > high_thr and t > high_thr:
        return 'high→high'
    elif low_thr <= n <= high_thr and low_thr <= t <= high_thr:
        return 'mid→mid'
    elif n < low_thr and low_thr <= t <= high_thr:
        return 'low→mid'
    elif n > high_thr and low_thr <= t <= high_thr:
        return 'high→mid'
    elif low_thr <= n <= high_thr and t < low_thr:
        return 'mid→low'
    elif low_thr <= n <= high_thr and t > high_thr:
        return 'mid→high'
    else:
        return 'other'

def plot_asm_mode_distribution(df, out_prefix):
    """Plot ASM mode distribution"""
    modes = df['asm_mode'].unique()
    labels = ['FP', 'TP']

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # 1. Bar chart
    mode_counts = df.groupby(['asm_mode', 'label']).size().unstack(fill_value=0).reindex(index=modes)
    mode_counts[labels].plot.bar(stacked=True, ax=axes[0,0])
    axes[0,0].set(title='ASM Mode Distribution (TP vs FP)', xlabel='ASM Mode', ylabel='Count')
    axes[0,0].tick_params(axis='x', rotation=45)

    # 2. Heatmap of proportions
    mode_props = mode_counts.div(mode_counts.sum(axis=1), axis=0)[labels]
    sns.heatmap(mode_props, annot=True, fmt='.2%', cmap='YlOrRd', ax=axes[0,1])
    axes[0,1].set(title='ASM Mode Proportions', xlabel='Label', ylabel='ASM Mode')

    # 3. Pie chart for TP
    tp_counts = df[df['label']=='TP']['asm_mode'].value_counts().reindex(index=modes).fillna(0)
    axes[1,0].pie(tp_counts, labels=tp_counts.index, autopct='%1.1f%%')
    axes[1,0].set_title('ASM Mode Breakdown in TP')

    # 4. Pie chart for FP
    fp_counts = df[df['label']=='FP']['asm_mode'].value_counts().reindex(index=modes).fillna(0)
    axes[1,1].pie(fp_counts, labels=fp_counts.index, autopct='%1.1f%%')
    axes[1,1].set_title('ASM Mode Breakdown in FP')

    plt.tight_layout()
    fig.savefig(f"{out_prefix}_asm_mode_distribution.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_haplotype_distribution(df, out_prefix):
    """Plot haplotype_tag distribution"""
    tags = df['haplotype_tag'].unique()
    labels = ['FP', 'TP']

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # 1. Bar chart
    haplo_counts = df.groupby(['haplotype_tag','label']).size().unstack(fill_value=0).reindex(index=tags)
    haplo_counts[labels].plot.bar(stacked=True, ax=axes[0,0])
    axes[0,0].set(title='Haplotype Distribution (TP vs FP)', xlabel='Haplotype Tag', ylabel='Count')
    axes[0,0].tick_params(axis='x', rotation=45)

    # 2. Heatmap
    haplo_props = haplo_counts.div(haplo_counts.sum(axis=1), axis=0)[labels]
    sns.heatmap(haplo_props, annot=True, fmt='.2%', cmap='YlOrRd', ax=axes[0,1])
    axes[0,1].set(title='Haplotype Proportions', xlabel='Label', ylabel='Haplotype Tag')

    # 3. Pie TP
    tp_counts = df[df['label']=='TP']['haplotype_tag'].value_counts().reindex(index=tags).fillna(0)
    axes[1,0].pie(tp_counts, labels=tp_counts.index, autopct='%1.1f%%')
    axes[1,0].set_title('Haplotype Breakdown in TP')

    # 4. Pie FP
    fp_counts = df[df['label']=='FP']['haplotype_tag'].value_counts().reindex(index=tags).fillna(0)
    axes[1,1].pie(fp_counts, labels=fp_counts.index, autopct='%1.1f%%')
    axes[1,1].set_title('Haplotype Breakdown in FP')

    plt.tight_layout()
    fig.savefig(f"{out_prefix}_haplotype_distribution.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_methylation_distribution(df, out_prefix):
    """Plot methylation level distributions"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # 1. Normal histogram
    sns.histplot(data=df, x='normal', hue='label', bins=50, kde=True, ax=axes[0,0])
    axes[0,0].set(title='Normal Sample Methylation Distribution', xlabel='Methylation Level', ylabel='Frequency')

    # 2. Tumor histogram
    sns.histplot(data=df, x='tumor', hue='label', bins=50, kde=True, ax=axes[0,1])
    axes[0,1].set(title='Tumor Sample Methylation Distribution', xlabel='Methylation Level', ylabel='Frequency')

    # 3. Delta histogram
    sns.histplot(data=df, x='delta', hue='label', bins=50, kde=True, ax=axes[1,0])
    axes[1,0].set(title='Delta Methylation (tumor - normal)', xlabel='Delta Methylation', ylabel='Frequency')

    # 4. Scatter plot
    sns.scatterplot(data=df, x='normal', y='tumor', hue='label', alpha=0.5, ax=axes[1,1])
    axes[1,1].set(title='Normal vs Tumor Methylation Levels', xlabel='Normal Level', ylabel='Tumor Level')

    plt.tight_layout()
    fig.savefig(f"{out_prefix}_methylation_distribution.png", dpi=300, bbox_inches='tight')
    plt.close(fig)

def main():
    args = parse_args()

    # Create output directory if needed
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    # 1. Load & pivot data
    tp = load_and_pivot(args.tp_summary)
    tp['label'] = 'TP'
    fp = load_and_pivot(args.fp_summary)
    fp['label'] = 'FP'

    # Keep relevant columns
    cols = ['chrom','somatic_pos','variant_type','vcf_source_id',
            'somatic_allele_type','haplotype_tag',
            'normal','tumor','delta','label']
    df = pd.concat([tp, fp], ignore_index=True)[cols]

    # Print data summary
    print("\nData Statistics:")
    print(f"Total rows: {len(df)}")
    print(f"Missing values in normal: {df['normal'].isna().sum()}")
    print(f"Missing values in tumor: {df['tumor'].isna().sum()}")
    print("\nValue ranges:")
    print(f"Normal range: {df['normal'].min()} to {df['normal'].max()}")
    print(f"Tumor range: {df['tumor'].min()} to {df['tumor'].max()}")

    # 2. Haplotype counts
    haplo_counts = (
        df
        .groupby(['haplotype_tag','label'])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    haplo_counts.to_csv(f"{args.out_prefix}_haplotype_counts.tsv", sep='\t', index=False)
    print(f"[+] Haplotype counts saved to: {args.out_prefix}_haplotype_counts.tsv")

    # Chi-square test for haplotype distribution
    chi2, p, dof, exp = chi2_contingency(haplo_counts[['FP','TP']].values)
    with open(f"{args.out_prefix}_haplotype_chi2.txt", 'w') as fo:
        fo.write(f"chi2 = {chi2:.4f}\n")
        fo.write(f"p-value = {p:.3g}\n")
        fo.write(f"dof = {dof}\n\n")
        fo.write("expected frequencies:\n")
        for row in exp:
            fo.write("  " + "\t".join(f"{x:.6g}" for x in row) + "\n")
    print(f"[+] Chi-square test results saved to: {args.out_prefix}_haplotype_chi2.txt")

    # 3. Compute ASM modes
    print("\nCalculating ASM modes...")
    df['asm_mode'] = df.apply(lambda r: asm_mode(r, args.low_threshold, args.high_threshold), axis=1)
    
    # Print ASM mode distribution
    print("\nASM mode counts:")
    print(df['asm_mode'].value_counts())
    print("\nASM mode counts by label:")
    print(df.groupby(['asm_mode', 'label']).size().unstack(fill_value=0))
    
    asm_counts = (
        df
        .groupby(['asm_mode','label'])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    asm_counts.to_csv(f"{args.out_prefix}_asm_mode_counts.tsv", sep='\t', index=False)
    print(f"[+] ASM mode counts saved to: {args.out_prefix}_asm_mode_counts.tsv")

    # Chi-square test for ASM modes
    chi2, p, dof, exp = chi2_contingency(asm_counts[['FP','TP']].values)
    with open(f"{args.out_prefix}_asm_mode_chi2.txt", 'w') as fo:
        fo.write(f"chi2 = {chi2:.4f}\n")
        fo.write(f"p-value = {p:.3g}\n")
        fo.write(f"dof = {dof}\n\n")
        fo.write("expected frequencies:\n")
        for row in exp:
            fo.write("  " + "\t".join(f"{x:.6g}" for x in row) + "\n")
    print(f"[+] ASM mode chi-square results saved to: {args.out_prefix}_asm_mode_chi2.txt")

    # 4. Generate plots
    print("\nGenerating plots...")
    plot_asm_mode_distribution(df, args.out_prefix)
    plot_haplotype_distribution(df, args.out_prefix)
    plot_methylation_distribution(df, args.out_prefix)
    print(f"[+] All plots generated")

if __name__ == '__main__':
    main()
