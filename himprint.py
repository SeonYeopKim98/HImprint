import sys
import os
import time
import argparse
import dask
import dask.multiprocessing
import numpy as np
import pandas as pd
from pomegranate.distributions import Categorical
from pomegranate.hmm import DenseHMM
from itertools import groupby
from datetime import datetime

dask.config.set(scheduler='processes')

VERSION = "1.0.0"

def get_args():
    parser = argparse.ArgumentParser(
        description="HImprint v{}: HMM-based Imprinted DMR Detection from Long-read data".format(VERSION),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:
  # From haplotype BED files:
  python himprint.py --hap1 sample.hap1.bed --hap2 sample.hap2.bed -o output.bed

  # From pre-calculated methylation difference file:
  python himprint.py -i sample.meth_diff.tsv -o output.bed

  # With custom parameters:
  python himprint.py --hap1 sample.hap1.bed --hap2 sample.hap2.bed \\
      -c 0.9 -l 10 --pe 0.6 --pt 0.95 -o output.bed
        """
    )

    input_group = parser.add_argument_group("Input options")
    input_group.add_argument(
        "-i", "--input",
        help="Pre-calculated methylation difference file (tab-separated, columns: chrom/start/end/meth_diff)"
    )
    input_group.add_argument(
        "--hap1",
        help="Haplotype 1 CpG methylation BED file (pb-CpG-tools format)"
    )
    input_group.add_argument(
        "--hap2",
        help="Haplotype 2 CpG methylation BED file (pb-CpG-tools format)"
    )
    input_group.add_argument(
        "-s", "--sample",
        help="Sample name (used in output filename and header; default: inferred from input filename)"
    )

    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument(
        "--min_cov", type=int, default=1,
        help="Minimum CpG coverage for both haplotypes (default: 1)"
    )
    filter_group.add_argument(
        "-c", "--cutoff", type=float, default=0.9,
        help="Methylation difference cutoff for imprinting classification\n"
             "  |meth_diff| >= cutoff → imprinted state (default: 0.9)"
    )

    hmm_group = parser.add_argument_group("HMM parameters")
    hmm_group.add_argument(
        "-l", "--length", type=int, default=10,
        help="Minimum number of CpGs in a detected DMR (default: 10)"
    )
    hmm_group.add_argument(
        "--pe", type=float, default=0.6,
        help="HMM emission probability for the dominant state\n"
             "  Higher → stricter CpG-level imprinting requirement (default: 0.6)"
    )
    hmm_group.add_argument(
        "--pt", type=float, default=0.95,
        help="HMM transition probability for staying in the same state\n"
             "  Higher → longer, more contiguous DMRs (default: 0.95)"
    )

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "-o", "--output",
        help="Output BED file path (default: auto-generated from sample name and parameters)"
    )
    output_group.add_argument(
        "-t", "--threads", type=int, default=16,
        help="Number of parallel threads for chromosome-level HMM (default: 16)"
    )

    args = parser.parse_args()

    if not args.input and (not args.hap1 or not args.hap2):
        parser.error("You must provide either -i/--input OR both --hap1 and --hap2.")

    return args


def calculate_methylation_diff(hap1_path, hap2_path, min_cov=1):
    print("[1/4] Reading haplotype BED files...")

    try:
        h1 = pd.read_csv(
            hap1_path, sep='\t', comment='#', header=None,
            usecols=[0, 1, 5, 8], names=['chrom', 'startp', 'cov1', 'meth1']
        )
        h2 = pd.read_csv(
            hap2_path, sep='\t', comment='#', header=None,
            usecols=[0, 1, 5, 8], names=['chrom', 'startp', 'cov2', 'meth2']
        )
    except Exception as e:
        print(f"  Error reading BED files: {e}")
        sys.exit(1)

    merged = pd.merge(h1, h2, on=['chrom', 'startp'])
    initial_len = len(merged)

    merged = merged[(merged['cov1'] >= min_cov) & (merged['cov2'] >= min_cov)]
    filtered_len = len(merged)
    print(f"  Removed {initial_len - filtered_len:,} low-coverage CpGs → {filtered_len:,} retained")

    merged['meth'] = (merged['meth1'] - merged['meth2']) / 100
    merged['endp'] = merged['startp'] + 1

    result_df = merged[['chrom', 'startp', 'endp', 'meth']].sort_values(['chrom', 'startp'])
    return result_df


def imprinted(diff, cutoff):
    if diff >= cutoff:
        return 1
    elif diff <= -cutoff:
        return 2
    else:
        return 0


def rle(in_seq, min_len):
    out_dict = {'char': [], 'startp': [], 'endp': []}
    n = 0
    for char, group in groupby(in_seq):
        char_list = list(group)
        count = len(char_list)
        start_p, end_p = n, n + count - 1
        if char != '0' and count >= min_len:
            out_dict['char'].append(char)
            out_dict['startp'].append(start_p)
            out_dict['endp'].append(end_p)
        n += count
    return pd.DataFrame(out_dict)


def run_hmm(seq_in, pe1, pt1):
    pe2 = (1 - pe1) / 2
    pt2 = (1 - pt1) / 2

    d1 = Categorical([[pe1, pe2, pe2]])
    d2 = Categorical([[pe2, pe1, pe2]])
    d3 = Categorical([[pe2, pe2, pe1]])

    model = DenseHMM()
    model.add_distributions([d1, d2, d3])

    model.add_edge(model.start, d1, 1/3)
    model.add_edge(model.start, d2, 1/3)
    model.add_edge(model.start, d3, 1/3)

    model.add_edge(d1, d1, pt1)
    model.add_edge(d1, d2, pt2)
    model.add_edge(d1, d3, pt2)

    model.add_edge(d2, d1, pt2)
    model.add_edge(d2, d2, pt1)
    model.add_edge(d2, d3, pt2)

    model.add_edge(d3, d1, pt2)
    model.add_edge(d3, d2, pt2)
    model.add_edge(d3, d3, pt1)

    X = np.array([[[['0', '1', '2'].index(char)] for char in seq_in]])
    y_hat = model.predict(X)
    return "".join([str(y.item()) for y in y_hat[0]])


@dask.delayed()
def run_hmm_chrom(df, chrom, params):
    output_list = []

    df_data = df[df['chrom'] == chrom].copy()
    df_data = df_data.sort_values('startp').reset_index(drop=True)

    imp_list = list(map(str, df_data['imprint'].tolist()))
    pred = run_hmm(imp_list, params['pe'], params['pt'])
    df_pred = rle(pred, params['length'])

    for _, row in df_pred.iterrows():
        s_idx, e_idx = int(row['startp']), int(row['endp'])

        if s_idx < len(df_data) and e_idx < len(df_data):
            genome_s = int(df_data['startp'].iloc[s_idx])
            genome_e = int(df_data['startp'].iloc[e_idx]) + 1

            dmr_slice = df_data.iloc[s_idx:e_idx + 1]
            n_cpg = len(dmr_slice)
            dmr_length = genome_e - genome_s

            state_int = int(row['char'])
            state_str = "hap1_methylated" if state_int == 1 else "hap2_methylated"

            output_list.append([
                chrom, genome_s, genome_e,
                state_str, n_cpg, dmr_length
            ])

    return output_list


def build_output_path(args):
    if args.output:
        return args.output

    if args.sample:
        sample = args.sample
    elif args.input:
        sample = os.path.basename(args.input).split(".")[0]
    else:
        sample = os.path.basename(args.hap1).split(".")[0]

    fname = (
        f"{sample}"
        f".himprint"
        f".c{args.cutoff}"
        f"_l{args.length}"
        f"_pe{args.pe}"
        f"_pt{args.pt}"
        f".bed"
    )
    return fname


def write_output(out_file, results, args, elapsed, total_input_cpgs):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sample_name = args.sample or (
        os.path.basename(args.input).split(".")[0] if args.input
        else os.path.basename(args.hap1).split(".")[0]
    )

    header_lines = [
        f"##HImprint v{VERSION}",
        f"##Date: {now}",
        f"##Sample: {sample_name}",
        f"##Parameters: cutoff={args.cutoff}, min_length={args.length}, "
        f"pe={args.pe}, pt={args.pt}, min_cov={args.min_cov}",
        f"##Input CpGs: {total_input_cpgs:,}",
        f"##Runtime: {elapsed:.2f}s",
        "#chrom\tstart\tend\tstate\tn_cpg\tdmr_length",
    ]

    count = 0
    with open(out_file, 'w') as f:
        for line in header_lines:
            f.write(line + "\n")
        for chrom_result in results[0]:
            if chrom_result:
                for row in chrom_result:
                    f.write("\t".join(str(v) for v in row) + "\n")
                    count += 1

    return count


def main():
    args = get_args()
    start_time = time.time()

    print(f"HImprint v{VERSION}")
    print("=" * 50)

    if args.input:
        print(f"[1/4] Reading input file: {args.input}")
        df = pd.read_csv(args.input, sep="\t", header=0, comment='#')
    else:
        df = calculate_methylation_diff(
            args.hap1, args.hap2,
            min_cov=args.min_cov
        )

    if len(df.columns) >= 4:
        df.columns = ['chrom', 'startp', 'endp', 'meth'] + list(df.columns[4:])

    total_input_cpgs = len(df)

    print(f"[2/4] Classifying CpGs (cutoff={args.cutoff})...")
    df['imprint'] = df['meth'].apply(lambda x: imprinted(x, args.cutoff))
    n_imp = (df['imprint'] != 0).sum()
    print(f"  Imprinted CpGs: {n_imp:,} / {total_input_cpgs:,} ({100*n_imp/total_input_cpgs:.1f}%)")

    chroms = list(df['chrom'].unique())
    print(f"[3/4] Running HMM on {len(chroms)} chromosomes ({args.threads} threads)...")

    params = {
        'cutoff': args.cutoff,
        'length': args.length,
        'pe': args.pe,
        'pt': args.pt,
    }

    futures = [run_hmm_chrom(df, chrom, params) for chrom in chroms]
    results = dask.compute(futures, num_workers=args.threads)

    out_file = build_output_path(args)

    end_time = time.time()
    elapsed = end_time - start_time

    print(f"[4/4] Writing output: {out_file}")
    count = write_output(out_file, results, args, elapsed, total_input_cpgs)

    print("=" * 50)
    print(f"Done. {count} DMRs detected → {out_file}")
    print(f"Runtime: {elapsed:.2f}s")


if __name__ == '__main__':
    main()
