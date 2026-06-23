# HImprint

**HMM-based Imprinted DMR Detection from Long-read data**

HImprint detects allele-specific differentially methylated regions (DMRs) — the hallmark of genomic imprinting — from haplotype-phased long-read methylation calls. Per-CpG haplotype methylation differences are classified into three states (`hap1-methylated`, `hap2-methylated`, `unimprinted`) and stitched into DMRs with a 3-state HMM (Viterbi decoding) run per chromosome in parallel.

---

## Installation

HImprint requires **Python ≥ 3.9** on Linux.

### Option 1 — `pip` (from GitHub)

```bash
pip install "git+https://github.com/<org>/HImprint.git"
```

### Option 2 — `conda` + `pip`

```bash
git clone https://github.com/<org>/HImprint.git
cd HImprint
conda env create -f environment.yml
conda activate himprint
pip install .
```

> **Note on dependencies.** `pomegranate` (the HMM backend) depends on **PyTorch**, which `pip` will install automatically. This adds several hundred MB to the install footprint. A GPU is **not** required — HImprint runs on CPU.

---

## Usage

HImprint accepts two input modes.

### Mode A — Two haplotype BED files (pb-CpG-tools format)

```bash
himprint \
    --hap1 sample.hap1.bed \
    --hap2 sample.hap2.bed \
    -s sample \
    -o sample.imprinted.bed
```

### Mode B — Pre-computed methylation-difference TSV

The TSV must have a header and columns `chrom`, `start`, `end`, `meth_diff` (meth_diff in [-1, 1], hap1 − hap2).

```bash
himprint \
    -i sample.meth_diff.tsv \
    -s sample \
    -o sample.imprinted.bed
```

### With custom parameters

```bash
himprint \
    --hap1 sample.hap1.bed --hap2 sample.hap2.bed \
    -c 0.9 -l 10 --pe 0.6 --pt 0.95 \
    -t 8 \
    -o sample.imprinted.bed
```

---

## Default parameters

| Flag           | Default | Meaning                                                                 |
|----------------|---------|-------------------------------------------------------------------------|
| `--min_cov`    | `1`     | Minimum per-haplotype CpG coverage                                      |
| `-c, --cutoff` | `0.9`   | \|meth_diff\| ≥ cutoff → CpG is in an imprinted state                   |
| `-l, --length` | `10`    | Minimum number of CpGs in a reported DMR                                |
| `--pe`         | `0.6`   | HMM emission probability for the dominant state (higher = stricter)     |
| `--pt`         | `0.95`  | HMM transition probability for staying in the same state (higher = longer DMRs) |
| `-t, --threads`| `16`    | Parallel workers (one chromosome per worker)                            |

> **Threading.** HMM decoding is parallelized per chromosome via Dask. In our benchmarks throughput plateaus around **~8 threads**; larger values give diminishing returns.

---

## Output

A BED-like file with header comments:

```
##HImprint v1.0.0
##Date: ...
##Sample: ...
##Parameters: cutoff=0.9, min_length=10, pe=0.6, pt=0.95, min_cov=1
##Input CpGs: ...
##Runtime: ...
#chrom  start  end  state            n_cpg  dmr_length
chr1    11111  22222  hap1_methylated  42     11111
chr1    33333  44444  hap2_methylated  37     11111
...
```

---

## License

MIT — see [LICENSE](LICENSE).
