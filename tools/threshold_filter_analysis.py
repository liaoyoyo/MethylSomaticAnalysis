import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_thresholds(s):
    # parse comma-separated list of ints
    return [int(x) for x in s.split(',') if x.strip().isdigit()]


def main():
    parser = argparse.ArgumentParser(
        description='Analyze TP/FP counts across read/CpG thresholds and plot heatmaps.'
    )
    parser.add_argument('--tp', required=True,
                        help='Path to TP summary TSV')
    parser.add_argument('--fp', required=True,
                        help='Path to FP summary TSV')
    parser.add_argument('--min-reads', default='0,5,10,15,20,25,30,40,50',
                        help='Comma-separated min_reads thresholds')
    parser.add_argument('--min-cpg', default='0,1,2,3,4,5,10',
                        help='Comma-separated min_cpg thresholds')
    parser.add_argument('--out-dir', default='.', help='Directory to save results')
    args = parser.parse_args()

    # create output directory
    os.makedirs(args.out_dir, exist_ok=True)

    # parse thresholds
    reads_th = parse_thresholds(args.min_reads)
    cpg_th = parse_thresholds(args.min_cpg)

    # load data
    tp = pd.read_csv(args.tp, sep='\t')
    fp = pd.read_csv(args.fp, sep='\t')

    # prepare matrices
    tp_counts = np.zeros((len(cpg_th), len(reads_th)), dtype=int)
    fp_counts = np.zeros_like(tp_counts)

    # iterate thresholds
    for i, c in enumerate(cpg_th):
        for j, r in enumerate(reads_th):
            tp_filt = tp[(tp.supporting_read_count >= r) & (tp.methyl_sites_count >= c)]
            fp_filt = fp[(fp.supporting_read_count >= r) & (fp.methyl_sites_count >= c)]
            tp_counts[i, j] = len(tp_filt)
            fp_counts[i, j] = len(fp_filt)

    # save raw counts
    pd.DataFrame(tp_counts, index=cpg_th, columns=reads_th).to_csv(
        os.path.join(args.out_dir, 'tp_counts.csv'), sep='\t')
    pd.DataFrame(fp_counts, index=cpg_th, columns=reads_th).to_csv(
        os.path.join(args.out_dir, 'fp_counts.csv'), sep='\t')

    # compute retention rates (fraction of initial)
    tp_total = len(tp)
    fp_total = len(fp)
    tp_rate = tp_counts / tp_total
    fp_rate = fp_counts / fp_total

    # plot heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(tp_rate, ax=axes[0], xticklabels=reads_th, yticklabels=cpg_th,
                annot=True, fmt='.2f', cmap='viridis')
    axes[0].set_title('TP Retention Rate')
    axes[0].set_xlabel('min_reads')
    axes[0].set_ylabel('min_cpg')

    sns.heatmap(fp_rate, ax=axes[1], xticklabels=reads_th, yticklabels=cpg_th,
                annot=True, fmt='.2f', cmap='viridis')
    axes[1].set_title('FP Retention Rate')
    axes[1].set_xlabel('min_reads')
    axes[1].set_ylabel('min_cpg')

    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'threshold_heatmap.png'), dpi=300)
    plt.close()

    # also line plots for a fixed cpg or reads
    # e.g., for each reads threshold, total TP/FP
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(reads_th, tp_counts[0, :], marker='o', label='TP (cpg>=0)')
    ax.plot(reads_th, fp_counts[0, :], marker='o', label='FP (cpg>=0)')
    ax.set_xlabel('min_reads')
    ax.set_ylabel('count')
    ax.set_title('Counts vs min_reads (min_cpg=0)')
    ax.legend()
    plt.savefig(os.path.join(args.out_dir, 'counts_vs_reads.png'), dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(cpg_th, tp_counts[:, 0], marker='o', label='TP (reads>=0)')
    ax.plot(cpg_th, fp_counts[:, 0], marker='o', label='FP (reads>=0)')
    ax.set_xlabel('min_cpg')
    ax.set_ylabel('count')
    ax.set_title('Counts vs min_cpg (min_reads=0)')
    ax.legend()
    plt.savefig(os.path.join(args.out_dir, 'counts_vs_cpg.png'), dpi=300)
    plt.close()

    print(f'Results saved to {args.out_dir}')


if __name__ == '__main__':
    main()
