# HImprint

> **Genome-wide imprinted DMR detection from haplotype-resolved long-read methylation, using a 3-state HMM**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python)](https://www.python.org/)
[![Version](https://img.shields.io/badge/Version-1.0.0-blue)]()
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey)]()

**HImprint** detects imprinted differentially methylated regions (DMRs) — regions where the two haplotypes carry opposite methylation — from phased long-read methylation data. It models the per-CpG methylation difference between haplotypes with a **directional 3-state hidden Markov model** (not-imprinted / haplotype 1–methylated / haplotype 2–methylated), so each detected region encodes *which* haplotype is methylated, not merely that the two differ. A run-length filter then keeps only contiguous segments, suppressing isolated single-CpG calls.

---

## Key Features

- **3-state HMM** — assigns each CpG to *not-imprinted*, *haplotype 1–methylated*, or *haplotype 2–methylated*, encoding the direction of the allelic difference.
- **Run-length filtering** — retains only runs of ≥ `-l` consecutive imprinted CpGs, removing short, isolated calls.
- **Annotation-free, genome-wide** — no ICR coordinate list or reference FASTA needed; it works from the per-CpG methylation difference alone.
- **Two input modes** — PacBio HiFi via pb-CpG-tools haplotype BEDs, or any platform via a pre-computed per-CpG difference file.
- **Chromosome-level parallelism** — Dask process-based scheduling, one worker per chromosome.
- **Per-CpG coverage filter** — drops CpGs below a minimum coverage on either haplotype before decoding (Mode 1 only).

---

## Installation

```bash
git clone https://github.com/<your-org>/HImprint.git
cd HImprint
pip install numpy pandas dask pomegranate
```

> Requires **Python ≥ 3.9** (tested on 3.10). `pomegranate` (≥ 1.0) provides the HMM and pulls in PyTorch as a dependency. A conda or virtual environment is recommended:

```bash
conda create -n himprint python=3.10
conda activate himprint
pip install numpy pandas dask pomegranate
```

---

## Usage

HImprint takes phased per-CpG methylation through one of two input modes.

### Mode 1 — From pb-CpG-tools haplotype BED files

Supply the two haplotype BED files produced by [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) (PacBio HiFi). HImprint computes the per-CpG difference internally.

```bash
python himprint.py \
  --hap1 sample.hap1.bed \      # Haplotype 1 BED (pb-CpG-tools format)
  --hap2 sample.hap2.bed \      # Haplotype 2 BED (pb-CpG-tools format)
  -s sample \                   # Sample name (optional; used in output name/header)
  --min_cov 1 \                 # Min coverage per CpG, both haplotypes (default: 1)
  -c 0.9 \                      # Per-CpG methylation-difference cutoff (default: 0.9)
  -l 10 \                       # Min consecutive CpGs per DMR (default: 10)
  --pe 0.6 \                    # HMM emission probability, dominant state (default: 0.6)
  --pt 0.95 \                   # HMM self-transition probability (default: 0.95)
  -t 16 \                       # Parallel workers / chromosomes (default: 16)
  -o sample.himprint.bed        # Output path (optional; auto-named if omitted)
```

> **Format note:** Mode 1 expects the **pb-CpG-tools** column layout (it reads columns 1, 2, 6, 9). Output from `modkit` or other bedMethyl producers places coverage and methylation in different columns and would be mis-parsed — for those, use **Mode 2** with a pre-computed difference.

### Mode 2 — From a pre-computed methylation-difference file

Platform-agnostic. Provide a per-CpG file containing the haplotype 1 − haplotype 2 difference (e.g. from ONT/`modkit` calls). This bypasses the pb-CpG-tools parsing.

```bash
python himprint.py \
  -i sample.meth_diff.tsv \     # Pre-computed difference (TSV with header)
  -c 0.9 -l 10 --pe 0.6 --pt 0.95 -t 16 \
  -o sample.himprint.bed
```

> `--min_cov` is applied only in Mode 1 (during difference computation); in Mode 2 the difference is taken as supplied.

---

## Parameters

| Argument | Short | Default | Description |
|---|---|---|---|
| `--input` | `-i` | — | Pre-computed per-CpG methylation-difference file (TSV, with header) |
| `--hap1` | — | — | Haplotype 1 BED (pb-CpG-tools format) |
| `--hap2` | — | — | Haplotype 2 BED (pb-CpG-tools format) |
| `--sample` | `-s` | inferred | Sample name for the output filename and header (default: from input filename) |
| `--min_cov` | — | `1` | Minimum coverage per CpG on **both** haplotypes (Mode 1 only) |
| `--cutoff` | `-c` | `0.9` | Per-CpG methylation-difference cutoff|
| `--length` | `-l` | `10` | Minimum number of consecutive CpGs to call a DMR |
| `--pe` | — | `0.6` | HMM emission probability for the dominant state (higher → stricter per-CpG requirement) |
| `--pt` | — | `0.95` | HMM self-transition probability (higher → longer, more contiguous DMRs) |
| `--output` | `-o` | auto | Output BED path |
| `--threads` | `-t` | `16` | Number of parallel workers (chromosome-level; Dask) |

Provide **either** `-i` **or** both `--hap1` and `--hap2`.

---

## Input Format

### Mode 1 — pb-CpG-tools haplotype BED

Tab-separated, no header. HImprint reads **columns 1, 2, 6, and 9**:

```
chr1    10468    10469    .    .    25    .    .    85.0
chr1    10470    10471    .    .    30    .    .    90.0
```

| Column | Field |
|---|---|
| 1 | Chromosome |
| 2 | Start position |
| 6 | Coverage |
| 9 | Methylation percentage (0–100) |

The two files are merged on chromosome + position, low-coverage CpGs are removed (`--min_cov`), and the difference is computed as `(meth_hap1 − meth_hap2) / 100`.

### Mode 2 — Pre-computed difference

Tab-separated **with a header line**; four columns. The difference is the haplotype 1 − haplotype 2 methylation, scaled to `−1.0` … `1.0`:

```
chrom    start    end      meth_diff
chr1     10468    10469    0.85
chr1     10470    10471    -0.72
```

| Column | Field |
|---|---|
| chrom | Chromosome |
| start | Start position |
| end | End position |
| meth_diff | (Hap1 − Hap2) / 100, range −1.0 … 1.0 |

---

## Output Format

A BED-like file with a metadata header (lines beginning `##`) followed by a column header and one row per detected DMR.

```
##HImprint v1.0.0
##Date: 2026-06-08 12:00:00
##Sample: sample
##Parameters: cutoff=0.9, min_length=10, pe=0.6, pt=0.95, min_cov=1
##Input CpGs: 28,415,002
##Runtime: 842.51s
#chrom  start      end        state             n_cpg   dmr_length
chr15   24954857   24956829   hap2_methylated   34      1972
chr7    50781028   50783615   hap1_methylated   41      2587
```

| Column | Field |
|---|---|
| chrom | Chromosome |
| start | DMR start (0-based, first CpG) |
| end | DMR end (last CpG + 1) |
| state | `hap1_methylated` or `hap2_methylated` |
| n_cpg | Number of CpGs in the DMR |
| dmr_length | DMR span in bp |

If `-o` is omitted, the file is named `{sample}.himprint.c{cutoff}_l{length}_pe{pe}_pt{pt}.bed`.

---

## Method Overview

```
<img width="1892" height="837" alt="Image" src="https://github.com/user-attachments/assets/d1722616-93fd-455b-83dc-836d099289e7" />
```

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| `numpy` | ≥ 1.24 | Numerical operations |
| `pandas` | ≥ 1.3 | Tabular I/O and merging |
| `dask` | ≥ 2022.0 | Chromosome-level parallelism |
| `pomegranate` | ≥ 1.0 | 3-state HMM (`DenseHMM`, `Categorical`); pulls in PyTorch |

---

## Citation

Manuscript in preparation. A release tag and Zenodo DOI will be added here on publication.

```
[Citation information coming soon]
```

---

## License

MIT — see [LICENSE](LICENSE).
