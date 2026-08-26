# genome_index

An implementation of genome indexing with a suffix array, read alignment, and mutation (SNP) calling.

---

## Overview

**genome_index** is a Python library that builds a suffix array from a raw reference genome FASTA file (hg38), uses it to align short reads back to the reference (exact and mismatch-tolerant), and calls point mutations (SNPs) from a set of aligned reads, filtering them out from random sequencing errors.

The algorithmic foundation is Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*, ch. 3 and 9:
- suffix array construction with **SA-IS** (Nong, Zhang & Chen, 2009), O(n)
- exact pattern matching via binary search over the suffix array
- mismatch-tolerant alignment via a **seed-and-extend** approach with a pigeonhole completeness guarantee
- mutation calling based on the **"bubble"** idea from assembly graphs: real variation is well-supported by multiple reads, a sequencing error is not

> A development/demonstration project — the mutation-calling thresholds (`min_coverage`, `min_variant_fraction`) are illustrative for demo-scale read sets, not clinically validated thresholds.

---

## Key Features

- Three suffix array construction backends: naive O(n² log n) (oracle for tests), pure-Python SA-IS O(n), and `divsufsort` (libdivsufsort) for whole-chromosome scale
- Mismatch-tolerant alignment on both strands (forward and reverse complement)
- Pigeonhole seed-and-extend: guaranteed to find every alignment within `max_mismatches`, while skipping seeds stuck in repeats/N-runs
- Pileup-based mutation calling with coverage and variant-allele-fraction thresholds, with correct base mapping from reverse-complement reads back to the forward reference
- Deduplication of identical reads during batch alignment
- A fuzz test that checks alignment output against a brute-force oracle across hundreds of randomized cases

---

## Tech Stack

- **Python** 3.13 (development/test version)
- **NumPy** — compact suffix array storage (`int32`/`int64` depending on sequence length)
- **pydivsufsort** (libdivsufsort) — fast suffix array backend at whole-chromosome scale
- **pytest** — unit and fuzz tests

---

## Project Structure

```text
genome_index/
├── genome_index/
│   ├── fasta_io.py            # parse a chromosome out of a (multi-)FASTA file
│   ├── suffix_array.py        # suffix array construction (SA-IS and divsufsort backend)
│   ├── alignment.py           # align a read/read set against the suffix array
│   └── mutation_calling.py    # pileup and SNP calling from aligned reads
├── tests/                     # unit and fuzz tests (pytest)
├── moj_primjer.py             # end-to-end example, from FASTA to mutations
└── requirements.txt
```

---

## Quickstart

### Prerequisites

- Python 3.13 (other 3.x versions likely work — the code doesn't use anything version-specific — but haven't been tested)
- a local copy of the hg38 reference genome (see below; not included in the repository)

### Install Dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` includes `numpy`, `pydivsufsort`, and `pytest`.

### Reference Genome

The repository does not include the reference genome (hg38 is several GB, so it's excluded via `.gitignore`). To run `moj_primjer.py` or the examples in `alignment.py`/`suffix_array.py`, you need locally:

```
hg38/hg38.u.fa
```

An hg38 assembly FASTA file (e.g. from UCSC, `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/`), unpacked. `parse_fasta_chromosome` reads a standard multi-FASTA format (sequence lines wrapped at a fixed width), so no extra preprocessing of the file is needed.

### Run Tests

```bash
pytest tests/ -v
```

Includes a comparison of SA-IS construction against a naive oracle, and a fuzz test that checks alignment against brute-force search across 560 randomized cases (exact, mutated, reverse-complement, and fully random reads).

### Run the Example

```bash
python moj_primjer.py
```

---

## How It Works

1. **Parsing** — `parse_fasta_chromosome` extracts a single chromosome (e.g. chrM) from the hg38 multi-FASTA file.
2. **Indexing** — `build_suffix_array_from_sequence` builds a suffix array for that sequence (SA-IS or divsufsort backend).
3. **Alignment** — `align_read`/`align_reads` search the suffix array via binary search for exact matches, or via seed-and-extend for matches with up to `max_mismatches` differences, on both strands.
4. **Mutation calling** — `call_mutations` builds a per-position pileup from all uniquely aligned reads and calls a mutation only where the alt base is sufficiently covered (`min_coverage`) and reproducibly supported (`min_variant_fraction`).

`moj_primjer.py` demonstrates this whole flow on chrM: aligning a single read with an introduced mutation back to its true position, then aligning 200 reads with planted mutations and random sequencing errors, and verifying that the planted mutations are found while the random errors are discarded.

---

## Engineering Highlights

- **SA-IS instead of DC3/skew** for linear-time suffix array construction — it only needs a single recursive step (sort LMS suffixes by induction, name them, recurse on collisions) instead of DC3's "split-then-merge" approach; this is also what `divsufsort` implements internally
- **Pigeonhole seed-and-extend**: splitting a read into `max_mismatches + 1` seeds guarantees at least one seed stays untouched by differences, so every alignment within the threshold is reliably found — confirmed by a fuzz test against a brute-force oracle (560 trials, 0 failures)
- Seeds with an excessive number of hits (repeat regions, N-runs — one all-N seed on chr21 had 6.6M hits during testing) are skipped rather than extended, trading recall for performance
- Multi-mapping reads (>1 alignment) are excluded from the pileup rather than inserted at every candidate position — trades recall for precision of called mutations
- The suffix array is stored as a NumPy `int32`/`int64` array (not a plain Python list) to save memory at whole-chromosome scale

---

## Future Improvements

- FM-index / BWT instead of a suffix array for a smaller memory footprint during search
- Support for indels (insertions/deletions); currently only point substitutions are called
- Parallelizing read alignment (batches are currently processed sequentially)
