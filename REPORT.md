# Comparison report: own alignment algorithm vs. STAR

**[English](#english)** · **[Hrvatski](#hrvatski)**

---

## English

**Summary.** The own read-alignment implementation (`genome_index/`) was
compared against STAR (v2.7.11b) at two scopes, both over the same real
dataset: all 14,988,127 single-end RNA-seq reads of *Drosophila
melanogaster* from `Fly/IFM48h_1.fastq`. At the scope of a single
chromosome (2L, 23,513,712 bases, via `run_pipeline.py`), strict
coordinate concordance with STAR is 97.5717%, and 99.9996% among the
reads the own aligner placed uniquely. At the whole-genome scope (all
1870 contigs merged into one index, via `star_comparison/scripts/`),
concordance is 96.42% strict, or 99.9999% among uniquely placed reads.
The two measurements agree within what their methodologies allow and
corroborate each other: the same set of positional disagreements appears
in both (7 of the 9 total disagreements are on chromosome 2L). The own
mutation-calling algorithm, run only at the single-chromosome scope,
called 10,445 point mutations that have not been independently
validated. Throughout this report, **measured values** are distinguished
from **derived conclusions**, which are labeled as such.

### Table of contents

1. [Comparison objective](#comparison-objective)
2. [Data](#data)
3. [Two comparison scopes](#two-comparison-scopes)
4. [Methodology](#methodology)
5. [Measurement conditions](#measurement-conditions)
6. [Commands used](#commands-used)
7. [Results on chromosome 2L](#results-on-chromosome-2l)
8. [Results on the whole genome](#results-on-the-whole-genome)
9. [Accuracy and reliability of the own algorithm](#accuracy-and-reliability-of-the-own-algorithm)
10. [Agreements and deviations relative to STAR](#agreements-and-deviations-relative-to-star)
11. [Possible causes of the observed differences](#possible-causes-of-the-observed-differences)
12. [Mutation-calling results](#mutation-calling-results)
13. [Strengths and limitations](#strengths-and-limitations)
14. [Assessment of practical applicability](#assessment-of-practical-applicability)
15. [Known issue: STAR on macOS (Apple Silicon)](#known-issue-star-on-macos-apple-silicon)
16. [Reproducing the pipeline](#reproducing-the-pipeline)
17. [Conclusion](#conclusion)

### Comparison objective

The own read-alignment implementation has no known ground truth it could
be checked against: the reads used are real RNA-seq data, not simulated,
so the true genomic origin of any given read is not known. Correctness is
therefore verified indirectly, by comparison with STAR, an independently
developed and widely used aligner.

The goal of the comparison is to answer the question: **does the own
aligner place reads on the same coordinates as STAR, and to what extent.**
STAR uses an uncompressed suffix array for its search (not an FM-index/
BWT), architecturally related to the approach used in the own
implementation, so the comparison is relevant to the question of
algorithmic correctness, not just speed.

STAR is used in this comparison as a **reference (baseline) aligner**, not
as a proven ground truth. The result is expressed as concordance with
STAR, not as absolute biological accuracy.

The comparison does not validate mutation calling. STAR does not call
variants, so the results in [that section](#mutation-calling-results) are
not independently confirmed.

### Data

| Item | Value |
|---|---|
| Reference genome | `Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa`, 146,261,227 B, 1870 records |
| Reads | `Fly/IFM48h_1.fastq`, 3,663,357,012 B, 14,988,127 single-end reads, 76 bases long |
| Reference aligner | STAR 2.7.11b |
| Tested aligner | own suffix-array aligner (`genome_index/`) |

The full set of reads was used at both comparison scopes, without
sampling.

### Two comparison scopes

The own aligner was tested in two separate setups, which cannot be
directly reduced to one number because they index different amounts of
the genome:

| | Chromosome 2L | Whole genome |
|---|---|---|
| Program | `run_pipeline.py` (root) | `star_comparison/scripts/` |
| Index | suffix array of chromosome 2L only (23,513,712 bases) | one suffix array over all 1870 contigs merged into a string (144 Mb), with 200bp N-spacers between contigs |
| Input to the comparison | all 14,988,127 reads, STAR filtered to chromosome 2L | all 14,988,127 reads, STAR filtered to the whole genome |
| STAR's comparison denominator | 1,755,830 (clean alignments on 2L) | 9,299,988 (clean alignments on any contig) |
| Purpose | matches the assignment's requirement literally (one chromosome, one suffix array) | verification at the scale of the whole genome, where repetitiveness and contig count matter more |

Both setups use the same `genome_index/` alignment code and the same
concordance criterion; only the indexed genome scope differs, and
consequently how many reads can even be covered by the comparison. Where
findings are consistent across the two scopes (e.g. [positional
disagreements](#agreements-and-deviations-relative-to-star)), this is
noted as additional confirmation.

### Methodology

#### Shared criteria: SAM fields and the CIGAR string

The standard SAM format (a header with `@HD`/`@SQ`, then one row per
alignment with 11 mandatory fields) is not described further here. What
matters for this comparison: `RNAME`+`POS` (the coordinate being
compared), `CIGAR` (the filtering criterion, below), and the optional
tags `NH` (number of loci the read was mapped to), `nM` (number of
mismatches in the alignment), and `AS` (STAR's local alignment score).
Example row:

```
ILLUMINA-075005_0053_FC:3:1:17967:18033#0  0  2L  3515966  255  76M  *  0  0  GAGAACTTCGCC...  hhhgehghhggg...  NH:i:1  HI:i:1  AS:i:74  nM:i:0
```

CIGAR describes an alignment as a sequence of `<number><operation>`
pairs, read left to right along the read:

| Op | Meaning | Consumes read bases | Consumes reference bases |
|---|---|---|---|
| `M` | aligned base (match or mismatch) | yes | yes |
| `I` | insertion in the read (bases absent from the reference) | yes | no |
| `D` | deletion (reference bases absent from the read) | no | yes |
| `N` | skipped reference region — in RNA-seq, an intron (splice) | no | yes |
| `S` | soft clipping — an unaligned end of the read, still present in SEQ | yes | no |
| `H` | hard clipping — a trimmed end, not present in SEQ either | no | no |

`M` **does not mean** the bases match, only that they are aligned without
a shift — mismatches are counted separately, in the `nM` tag. A read with
CIGAR `76M` can still carry an occasional point mutation or sequencing
error; that is exactly the scenario the own aligner supports
(mismatch-tolerant, without indels). The largest observed mismatch count
(`nM`) across the whole-genome run's SAM file is 10 on a 76bp read
(15,124,860 records checked).

The "aligned without a break" criterion, used at both comparison scopes:
the whole CIGAR must be exactly one `M` operation whose number matches
the read length — for 76bp reads, that is `76M`. Any `S`, `H`, `I`, `D`,
or `N` in the CIGAR automatically excludes the read from the comparison
(soft/hard clip, insertion, deletion, or a splice/intron gap). Examples
that get rejected:

```
...:17532:18022#0  16  X   2001069   255  75M1S       ...   -> 1 base soft-clipped
...:19344:18030#0  16  2R  13969542  255  61M69N15M   ...   -> read spans an intron (69bp splice)
```

Both scopes also exclude secondary (`FLAG & 0x100`), supplementary
(`FLAG & 0x800`), and unmapped (`FLAG & 0x4`) STAR alignments.

#### Note on the mismatch-count tag

The assignment mentions the SAM tag `NM:i:` for checking fully identical
alignments. Checking the entire SAM output (both scopes) found that STAR,
in this configuration, does not write that tag, but its own `nM:i:` tag
instead. The mismatch-count comparison therefore uses `nM:i:` throughout.

#### Accuracy criterion and tolerance

An alignment is considered concordant if the chromosome and start
coordinate match, with a tolerance of 0 bases. Zero tolerance was applied
without relaxation — it was checked (whole-genome scope, 9,299,988 reads)
that a tolerance of ±1 base gives an identical number of matches as
tolerance 0 ([details](#coordinate-tolerance-check)), which rules out a
systematic error in the coordinate conversion (0-indexed suffix array
position → 1-indexed SAM position).

At the chromosome-2L scope, the "same chromosome" condition is satisfied
by construction, since the own aligner's index contains only 2L —
misplacing a read from another chromosome would therefore show up as a
wrong coordinate within 2L, not as a wrong chromosome. At the
whole-genome scope, the chromosome is compared explicitly, since the
index covers every contig.

#### Procedure — chromosome 2L

The own pipeline (`run_pipeline.py`) builds a suffix array of chromosome
2L only and aligns all reads against it. STAR builds an index over the
whole genome and aligns the same reads, with SAM output. From STAR's SAM
output, alignments are taken that refer to chromosome 2L, satisfy the
[shared criterion](#shared-criteria-sam-fields-and-the-cigar-string), and
are not secondary/supplementary/unmapped. The comparison is run with
[`star_comparison/scripts/compare_chromosome.py`](star_comparison/scripts/compare_chromosome.py),
which does the SAM filtering and the comparison in a single pass.

#### Procedure — whole genome

1. **STAR alignment.** STAR aligns all 14,988,127 reads against the whole
   reference genome (`--outSAMtype SAM`), with full support for splicing
   and indels (default settings, no overrides).
2. **Filtering clean reads**
   ([`star_comparison/scripts/filter_clean_reads.py`](star_comparison/scripts/filter_clean_reads.py))
   by the [shared criterion](#shared-criteria-sam-fields-and-the-cigar-string).
3. **Alignment with the own aligner** against a single suffix array built
   over all 1870 contigs merged into one string — 200bp N-spacers between
   contigs prevent a read from falsely aligning across a contig boundary
   ([`star_comparison/scripts/build_my_index.py`](star_comparison/scripts/build_my_index.py),
   [`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py),
   `max_mismatches=3`).
4. **Coordinate comparison**
   ([`star_comparison/scripts/compare_coordinates.py`](star_comparison/scripts/compare_coordinates.py))
   between (RNAME, POS) from STAR and from the own aligner.

The comparison is deliberately asymmetric in scope: STAR handles every
read, including splicing and indels, while the own aligner only handles
the subset without those operations. That is a deliberate choice — it
measures the reliability of *positioning* on reads the own aligner should
be able to handle at all, not STAR's full RNA-seq functionality.

### Measurement conditions

| Item | Value |
|---|---|
| Machine | Apple M2 (Mac14,15), 8 cores, 8.6 GB RAM |
| OS | macOS 26.5 (build 25F71) |
| Compiler | Apple clang 15.0.0 |
| STAR | 2.7.11b, built from source with a macOS fix (see [Known issue](#known-issue-star-on-macos-apple-silicon)) |
| Python / numpy / pydivsufsort / pytest | 3.13.5 / 2.5.1 / 0.0.20 / 9.1.1 |
| Time/memory measurement | `/usr/bin/time -l` (BSD `time`, macOS); peak RSS is `maximum resident set size` in bytes |
| Date of the STAR run (whole-genome) | 2026-09-06, per `IFM48h_1.Log.out` |

Values refer to single runs on one machine and are not an average of
multiple repetitions. The runs for the two scopes were performed
separately, in different work sessions on the same machine.

### Commands used

**STAR — index build and alignment over the whole genome** (used for
both comparison scopes, since STAR always indexes the whole genome):

```bash
STAR --runMode genomeGenerate --runThreadN 6 \
  --genomeDir star_comparison/index \
  --genomeFastaFiles Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
  --genomeSAindexNbases 12 --genomeChrBinNbits 16

STAR --runMode alignReads --runThreadN 6 \
  --genomeDir star_comparison/index \
  --readFilesIn Fly/IFM48h_1.fastq \
  --outSAMtype SAM \
  --outFileNamePrefix star_comparison/results/star/IFM48h_1.
```

**Own pipeline — chromosome 2L:**

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_chromosome.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  results/chr2L/my_alignments_full.tsv \
  2L results/chr2L/comparison_chr2L.csv --tolerance 0
```

**Own pipeline — whole genome** (the full sequence with each step's
inputs/outputs is in
[Reproducing the pipeline](#reproducing-the-pipeline)):

```bash
python3 star_comparison/scripts/filter_clean_reads.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  star_comparison/results/filtered/clean_reads.tsv

python3 star_comparison/scripts/build_my_index.py

python3 star_comparison/scripts/run_my_aligner.py \
  star_comparison/results/my_aligner/clean_reads_input.fasta \
  star_comparison/results/my_aligner/my_alignments.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_coordinates.py \
  star_comparison/results/filtered/clean_reads.tsv \
  star_comparison/results/my_aligner/my_alignments.tsv \
  star_comparison/results/comparison/comparison.csv \
  --mismatches-out star_comparison/results/comparison/mismatches.csv \
  --tolerance 0
```

### Results on chromosome 2L

#### Alignment outcomes (2L)

The own aligner over all 14,988,127 reads, index of chromosome 2L,
`max_mismatches = 3`:

| Outcome | Number of reads | Share |
|---|---|---|
| `unique` | 2,321,574 | 15.49% |
| `multi` | 29,104 | 0.19% |
| `unmapped` | 12,637,449 | 84.32% |
| `skipped` | 0 | 0.00% |

The high share of unmapped reads is expected and **not a measure of
error**: the reads come from the whole genome, while the index covers
only chromosome 2L, whose length is about 16% of the assembly's total
length.

STAR over the whole genome, per its own `IFM48h_1.Log.final.out` report:

| Measure | Value |
|---|---|
| Uniquely mapped | 14,424,378 (96.24%) |
| Mapped to multiple loci | 197,559 (1.32%) |
| Unmapped | 342,852 (2.29%) |
| Number of splice junctions | 1,325,356 |
| Mismatch rate per base | 0.64% |

The splice-junction count confirms that the dataset contains a
significant share of reads that span exon boundaries, which the own
aligner cannot correctly align by construction.

#### Coordinate comparison (2L)

After filtering, STAR has 1,755,830 clean primary alignments on
chromosome 2L. That is the comparison's denominator.

| Measure | Count | Share |
|---|---|---|
| Both tools aligned uniquely | 1,713,201 | 97.58% |
| — same coordinate | 1,713,194 | 97.5717% |
| — different coordinate | 7 | 0.0004% |
| — same orientation | 1,713,201 | 100.00% |
| — different orientation | 0 | 0.00% |
| — same mismatch count (`nM`) | 1,710,110 | 99.82% |
| Own aligner: `multi` | 10,375 | 0.59% |
| Own aligner: `unmapped` | 32,254 | 1.84% |
| Missing from the own aligner's output | 0 | 0.00% |

Accuracy = reads with the same coordinate / STAR's clean alignments
= 1,713,194 / 1,755,830 = **97.5717%**.

Accuracy among reads the own aligner placed uniquely =
1,713,194 / 1,713,201 = **99.9996%**.

#### Runtime and memory (2L)

| Step | Tool | Threads | Reads | real | Peak RSS |
|---|---|---|---|---|---|
| Index build, whole genome | STAR | 6 | — | 53.49 s | 1.54 GB |
| Alignment, whole genome | STAR | 6 | 14,988,127 | 97.08 s | 1.98 GB |
| Suffix array build, 2L | own | 1 | — | 1.19 s | — |
| Alignment, 2L | own | 1 | 14,988,127 | 1529.34 s | 0.17 GB |
| Alignment and mutations, 2L | own | 1 | 14,988,127 | 1842.24 s | 2.90 GB |

The own aligner's alignment throughput was 9810 reads per second. The
own-aligner figures refer to building an index for one chromosome and
searching within it, while STAR's refer to the whole genome with six
threads. **The figures are therefore not directly comparable** and are
reported separately; a direct comparison at equal search scope is in
[Runtime (whole genome)](#runtime-whole-genome) below.

### Results on the whole genome

#### Filtering the SAM output

Distribution of CIGAR operations across all 14,621,937 primary
alignments, computed with an independent `awk` pass over the SAM (matches
the numbers reported by the Python filter,
[`star_comparison/scripts/filter_clean_reads.py`](star_comparison/scripts/filter_clean_reads.py)):

| Category | Reads | Share of primary |
|---|---|---|
| clean `<length>M` (input to the rest of the comparison) | 9,299,988 | 63.60% |
| contains `S` (soft clipping) | 4,236,843 | 28.98% |
| contains `N` (splice/intron) | 1,338,755 | 9.16% |
| contains `D` (deletion) | 51,549 | 0.35% |
| contains `I` (insertion) | 46,401 | 0.32% |
| contains `H` (hard clipping) | 0 | 0.00% |

Categories overlap — a single read can have both a splice and clipping —
so the shares do not sum to 100%.

#### Alignment outcomes (whole genome)

Over 9,299,988 clean reads, the own aligner
([`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py),
`max_mismatches=3`) classifies every read as uniquely aligned
(**unique**), **multi** (more than one alignment), or **unmapped**
(none):

| Outcome | Reads | Share |
|---|---|---|
| unique | 8,966,941 | 96.42% |
| multi-mapped | 179,897 | 1.93% |
| unmapped | 153,150 | 1.65% |

The share of "unique" outcomes is called the **unique placement rate**
(coverage). By itself it says nothing about accuracy — only how often the
aligner gives a single answer at all.

#### Coordinate concordance with STAR (whole genome)

For every clean read, (RNAME, POS) from the own aligner is compared
against STAR's (tolerance 0 bases):

| Measure | Numerator / denominator | Value | Definition |
|---|---|---|---|
| Unique placement rate (coverage) | 8,966,941 / 9,299,988 | 96.42% | share of clean reads for which the own aligner returns exactly one position |
| Concordance with STAR among unique calls | 8,966,932 / 8,966,941 | 99.9999% | of the reads the own aligner placed uniquely, the share that matches STAR |
| Strict coordinate concordance | 8,966,932 / 9,299,988 | 96.42% | matches relative to *all* clean reads — multi-mapped and unmapped reads count as non-matches |

These three measures are deliberately not collapsed into a single
"accuracy" number, since they answer different questions: the first
measures coverage, the second measures accuracy given that the aligner
answered, the third combines both. The gap between strict concordance
(96.42%) and the near-perfect concordance among unique calls (99.9999%)
shows that the shortfall comes almost entirely from lack of coverage, not
from incorrectly determined coordinates — of 8,966,941 unique calls, only
9 disagree with STAR
([details](#checking-the-nine-disagreements)).

#### Split by STAR's NH tag

STAR itself reports `NH:i:N` when a read is equally well aligned to N
locations and picks one as primary by its own internal rules. Of
9,299,988 clean reads, 109,411 (1.18%) have STAR `NH≥2`:

| Set (per STAR) | Reads | Own aligner: unique | — multi | — unmapped | Concordance with STAR |
|---|---|---|---|---|---|
| `NH:i:1` (STAR unique) | 9,190,577 | 8,962,784 (97.52%) | 78,951 (0.86%) | 148,842 (1.62%) | 97.52% strict, 99.9999% among unique calls |
| `NH≥2` (STAR multi-mapped) | 109,411 | 4,157 (3.80%) | 100,946 (92.26%) | 4,308 (3.94%) | 3.80% strict, 99.9759% among unique calls |

On the subset where STAR itself is confident (`NH:i:1`), the own aligner
achieves 97.52% strict concordance — higher than the overall average
(96.42%), because reads that are ambiguous by definition of the problem
are excluded. On the subset where STAR itself reported more than one
possible locus, the own aligner also most often (92.26%) returns
multi-mapping rather than a single position: the two independent tools
agree, in the large majority of cases, on which reads are ambiguous.

The only "unique" call within the `NH≥2` set that disagrees with STAR
(read `...3:14:17033:20323#0`) was also checked against STAR's secondary
locus: STAR reported 3L:28,103,368 (primary) and X:5,752,974 (secondary,
`NH:i:2`), while the own aligner reported a third position,
3L:18,780,015 — matching neither of the two.

#### Checking the nine disagreements

For all 9 reads where the own aligner settled on a unique position that
does not match STAR's, the actual number of mismatches (Hamming distance)
between the read and the reference genome was counted at both reported
positions
([`star_comparison/scripts/inspect_mismatches.py`](star_comparison/scripts/inspect_mismatches.py),
result in
[`star_comparison/results/comparison/mismatches_verified.csv`](star_comparison/results/comparison/mismatches_verified.csv)):

| Read | STAR position | STAR: mismatches | Own aligner's position | Own: mismatches |
|---|---|---|---|---|
| `...3:8:10065:3458#0` | 3L:3,901,348 | 5 | 3L:3,901,576 | 3 |
| `...3:14:17033:20323#0` | 3L:28,103,368 | 9 | 3L:18,780,015 | 3 |
| remaining 7 reads (identical pattern) | 2L:14,743,463 | 4 | 2L:14,743,493 | 3 |

These same 7 reads, at the same positions and with the same mismatch
count, also appear in the separate comparison limited to chromosome 2L
([Coordinate deviations](#coordinate-deviations)) — two independently
run comparisons give an identical finding.

At all 9 positions, the own aligner's position has a lower Hamming
distance than the one STAR reported (the count matches STAR's own `nM`
tag — `nM:i:5`, `nM:i:9`, `nM:i:4` for the reads shown — ruling out an
error in the counting itself). Eight of the 9 reads have STAR `NH:i:1`,
i.e. STAR itself considers them uniquely mapped.

**This is a measured fact** (lower Hamming distance = better-ranked
position by that criterion), **not proof** of the read's true biological
origin: the true coordinate of origin is not known, and a lower Hamming
distance can come from the true location, but also from a chance match to
a repeated sequence elsewhere in the genome.

Possible reasons for STAR choosing a different position include
differences in candidate generation, seed heuristics, and scoring
relative to the pigeonhole seed-and-extend used here — these are
hypotheses, not a confirmed cause. Since all 9 positions have CIGAR `76M`
on both sides of the comparison, the difference cannot specifically be
attributed to STAR's splice model.

The mismatch-counting method was also checked on 2000 random reads from
the set where the tools agree (a sanity check, not proof about the whole
set): 73% have 0 mismatches at the reported position, 94% at most 1.

#### Coordinate tolerance check

The comparison takes a `--tolerance` parameter (the largest position
difference still counted as a match). Run on all 9,299,988 clean reads
with tolerance 0 and tolerance ±1 base: both give an identical 8,966,932
matches. Widening the tolerance by one base changes no result, which
rules out a systematic off-by-one error in the coordinate conversion.

#### Runtime (whole genome)

| Phase | Tool | Threads | Reads | real | user | sys | Throughput |
|---|---|---|---|---|---|---|---|
| Genome index | STAR | 6 | — | 53.49 s | 189.83 s | 15.47 s | — |
| Genome index | own aligner | 1 | — | 9.06 s | 8.35 s | 0.48 s | — |
| Alignment | STAR | 6 | 14,988,127 | 97.08 s | 414.57 s | 52.20 s | ~154,400 reads/s |
| Alignment | own aligner | 1 | 9,299,988 | 1009.52 s (16m50s) | 968.57 s | 24.54 s | ~9,210 reads/s |

The thread count is confirmed from `--runThreadN 6` in the STAR
invocation and independently from the user/real time ratio: STAR's
alignment has user/real ≈ 4.27 (expected for 6 threads with some
non-parallel work), while the own aligner has user/real ≈ 0.96
(consistent with
[`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py)
using neither threading nor multiprocessing).

At the per-read level, STAR is ~16.8× faster during alignment (6.5
µs/read vs. 108.6 µs/read). This ratio **does not represent a normalized
comparison of the underlying algorithms** — it was measured in a
thread-asymmetric configuration (6 vs. 1), so part of the difference
comes from parallelism, not solely the algorithm. The own aligner is
also a simple, single-threaded Python program (the `divsufsort` C
backend is used only for suffix array construction; the per-read search
itself is a plain Python loop).

#### Peak memory (whole genome)

The value is `maximum resident set size` from `/usr/bin/time -l`. On
macOS it is expressed in bytes:

| Phase | Tool | Peak RSS (bytes) | Peak RSS (GB) |
|---|---|---|---|
| Genome index | STAR | 1,535,000,576 | 1.54 GB |
| Genome index | own aligner | 891,928,576 | 0.89 GB |
| Alignment | STAR | 1,978,810,368 | 1.98 GB |
| Alignment | own aligner | 1,562,148,864 | 1.56 GB |

On this machine and in this measurement, peak RSS for both tools is of
the same order of magnitude (1-2 GB), with the own aligner slightly more
economical during the alignment step. This is not a general claim about
either implementation's memory footprint — only a description of one
measurement on this genome (144 Mb): both tools spend memory dominantly
on the genome index.

### Accuracy and reliability of the own algorithm

**Measured.** Of the 1,713,201 reads both tools placed uniquely on 2L,
1,713,194 have an identical coordinate, and all 1,713,201 have the same
orientation. Not a single read from STAR's set is missing from the own
aligner's output.

**Measured.** Two independent runs of the own pipeline over 2L, with and
without mutation calling, produced identical alignment outcomes
(2,321,574 `unique`, 29,104 `multi`, 12,637,449 `unmapped`). The
algorithm is deterministic and does not depend on processing order.

**Measured.** An independent whole-genome comparison, over 9,299,988
clean reads, gave a strict accuracy of 96.42% and 99.9999% concordance
among uniquely placed reads. The seven positional disagreements found in
that comparison also appear in the comparison limited to 2L, at the same
positions and with the same mismatch counts — two independently run
comparisons give an identical finding.

**Measured.** Unit tests (16 tests) compare SA-IS against naive suffix
sorting on random strings, and alignment against an exhaustive search
procedure. Edge cases are covered as well: empty input, a read with no
alignment, reads on the negative strand, invalid input characters, a read
longer than the reference, and rejection of mutations at insufficient
coverage. The `moj_primjer.py` script additionally verifies the full
pipeline over data with a known ground truth.

**Derived conclusion.** Taken together, this supports the view that the
algorithm's core — binary search over the suffix array with
seed-and-extend under the pigeonhole guarantee — is correctly implemented
for the class of reads it handles (without indels, deletions, or
splicing).

**Limitation of this conclusion.** At the 2L scope, accuracy was
measured only on reads STAR cleanly aligned on that chromosome. The own
aligner placed 2,321,574 reads uniquely, of which 1,713,201 enter the
comparison; the remaining 608,373 **have not been verified by any
independent source**. At the whole-genome scope, the analogous unverified
gap is relatively smaller (the comparison covers 9,299,988 of all clean
reads), but it still excludes every read with splicing and indels, which
the own aligner does not handle by construction.

### Agreements and deviations relative to STAR

#### Coordinate deviations

Seven reads (at the 2L scope) were placed uniquely by both tools, but at
different positions. All the deviations are identical:

| Reads | Position per STAR (`nM`) | Position per own aligner (mismatches) | Difference |
|---|---|---|---|
| 7 | 2L:14,743,463 (4) | 2L:14,743,493 (3) | 30 bases |

These 7 are a subset of the total 9 positional disagreements found at the
whole-genome scope ([Checking the nine
disagreements](#checking-the-nine-disagreements)); the remaining 2
disagreements involve reads whose STAR alignment falls outside 2L (on
3L), so they do not appear in the comparison limited to 2L.

At all seven positions, the own aligner reports fewer mismatches than
STAR. This is measured and confirmed by STAR's own `nM` tag, but **does
not prove** that its position is biologically correct.

#### Mismatch-count deviations

Of the 1,713,201 reads both tools placed uniquely on 2L, 3,091 have a
different mismatch count. The cause was determined by directly checking
those reads' sequences in the input FASTQ file:

| Explanation | Reads | Share |
|---|---|---|
| The difference exactly matches the number of `N` bases in the read | 3,084 | 99.77% |
| Reads with a coordinate deviation | 7 | 0.23% |

STAR does not count an `N` base in the read as a mismatch, while the own
implementation compares character by character and counts it as one.
This is a difference in counting convention, not in read placement: the
coordinate is identical in all 3,084 cases.

The remaining seven reads have no `N` base at all, and the mismatch-count
difference is exactly −1, a consequence of the different position, not
of the convention.

### Possible causes of the observed differences

| Observed difference | Possible cause | Status |
|---|---|---|
| 3,084 reads (2L) with a different mismatch count | different convention for counting `N` bases | established by checking sequences |
| 7 reads (2L) at a different position | two copies of a repeated sequence 30 bases apart, where both tools find an acceptable alignment | hypothesis; supported by the identical difference across all seven reads and the low mismatch count at both positions |
| 32,254 reads (2L) the own aligner did not align, which STAR placed cleanly on 2L | seeds rejected for exceeding `max_seed_hits`, or the three-mismatch threshold exceeded | measured (see below) |
| 10,375 reads (2L) flagged as `multi` | repetitive regions within chromosome 2L | hypothesis |
| 153,150 reads (whole genome) the own aligner did not align | same mechanism as above, at whole-genome scale | hypothesis by analogy, not separately measured at whole-genome scale |
| Speed difference | single-threaded Python vs. a multi-threaded C++ program | partially measured; the shares due to parallelism vs. implementation language are not separated |

#### Why reads stay unmapped: measured (2L)

To distinguish two possible causes — a seed rejected for having too many
matches, versus a candidate that fails the mismatch threshold on
extension — a separate script was written
([`analyze_unmapped_chr2L.py`](analyze_unmapped_chr2L.py)) that, for
every `unmapped` read, calls the same internal seed-search functions as
`align_read` (`_split_into_seeds`, `_sa_search_bounds`), without modifying
`genome_index/alignment.py`. For every read, one of three states is
determined:

- **`no_exact_seed_match`** — none of the 4 seeds, on either of the two
  strands, has a single exact match in the sequence; no candidate
  position is generated at all.
- **`seed_rejected_frequent`** — at least one seed has exact matches, but
  all such matches were rejected for exceeding `max_seed_hits` (1000); no
  candidate is generated here either.
- **`extension_failed`** — at least one candidate was generated (a seed
  within `max_seed_hits`), but none satisfied the three-mismatch
  threshold when the whole read was compared.

The measurement was run at two scopes: over all 12,637,449 reads the own
aligner does not align on 2L, and separately over the subset of 32,254
that STAR nonetheless placed cleanly on 2L (so a genomic location is
known to exist for them within the indexed chromosome).

| Cause | All unmapped on 2L (12,637,449) | Subset STAR places cleanly (32,254) |
|---|---|---|
| `no_exact_seed_match` | 12,011,684 (95.05%) | 5,997 (18.59%) |
| `seed_rejected_frequent` | 0 (0.00%) | 0 (0.00%) |
| `extension_failed` | 625,765 (4.95%) | 26,257 (81.41%) |

`seed_rejected_frequent` measured as 0 at both scopes: no ~19-base seed
(a 76bp read split into 4 seeds for `max_mismatches = 3`) has more than
1000 exact matches on chromosome 2L, so `max_seed_hits` was not a
limiting factor in this measurement.

At the level of all 12.6 million unmapped reads, 95.05% have no exact
19-base match at all on 2L — expected, since most of those reads
originate from other chromosomes. Within the narrower, biologically
relevant subset that STAR places cleanly on 2L, the picture reverses:
81.41% fail because of the mismatch threshold on candidate extension.
This suggests that these reads are genomically located on 2L, but carry
more than three differences from the reference (real variants and/or
sequencing errors) — more than the own aligner, with its default
`max_mismatches = 3`, allows.

The full per-read result is in
`results/chr2L/unmapped_reason_breakdown.tsv` (gitignored due to size,
with the committed sample `unmapped_reason_breakdown.sample.tsv`). The
run took 1062.34 s. This measurement was not repeated at the
whole-genome scope.

### Mutation-calling results

Mutation calling was run only at the chromosome-2L scope, via
`run_pipeline.py`; the whole-genome pipeline in `star_comparison/` does
not call mutations.

Over 2,321,574 uniquely aligned reads, with `min_coverage = 3` and
`min_variant_fraction = 0.5`, **10,445 point mutations** were called in
264.24 s.

| Measure | Value |
|---|---|
| Coverage (min / median / max) | 3 / 6 / 9582 |
| VAF (min / median / max) | 0.5000 / 0.8468 / 1.0000 |
| Transitions | 5896 |
| Transversions | 4549 |
| Ti/Tv ratio | 1.296 |

Distribution by strength of evidence:

| Group | Count | Share |
|---|---|---|
| Coverage 3 | 2106 | 20.2% |
| Coverage below 5 | 3768 | 36.1% |
| Coverage 10 or more | 3951 | 37.8% |
| VAF equal to 1.0 | 4412 | 42.2% |
| Coverage 3 and VAF 1.0 | 1153 | 11.0% |

A Ti/Tv ratio above 1 is typical of real genomic variation. This is
**weak, set-wide corroboration** and says nothing about the correctness
of any individual call.

More than a third of the calls rest on coverage below five reads, and
1153 calls have coverage 3 with full agreement among the reads. At that
coverage, three reads sharing the same sequencing error satisfy the set
thresholds, so this group cannot be claimed to represent real variants.

**These results have not been independently validated.** STAR does not
call variants and cannot serve as their confirmation. A serious
assessment would require at least:

- comparing the positions against a known variant database for
  *Drosophila melanogaster*,
- calling variants with an established tool over the same alignments and
  comparing the sets,
- separating RNA editing from genomic variants, since this is RNA-seq
  data,
- testing the results' sensitivity to the coverage and VAF thresholds.

Until then, the called mutations should be treated as algorithm output,
not as a set of confirmed biological variants.

### Strengths and limitations

#### Strengths

- **Placement accuracy on the verified subset.** On reads the own
  aligner places uniquely, concordance with STAR is 99.9996% for
  chromosome 2L and 99.9999% in the whole-genome comparison.
- **Determinism.** Repeated runs give identical results.
- **Completeness within the given threshold.** The pigeonhole guarantee
  ensures the procedure does not miss alignments with at most
  `max_mismatches` mismatches, except for reads whose every seed was
  rejected for having too many matches (measured as 0 cases on 2L).
- **Modest memory use during alignment.** 0.17 GB (2L) to 1.56 GB (whole
  genome) to process 15 million reads, since reads are processed one at
  a time.

#### Limitations

- **No support for insertions, deletions, or splicing.** With 1,325,356
  splice junctions recorded by STAR, this limitation excludes a
  significant share of RNA-seq data (36.4% of primary STAR alignments
  contain at least one of these operations).
- **No alignment-confidence score.** There is no MAPQ equivalent, so
  alignments cannot be filtered by quality.
- **Speed.** At an equal search scope (whole genome), the own aligner was
  about 16.8 times slower than STAR, in a configuration of one thread
  versus six.
- **Memory use during mutation calling.** Peak RSS grows from 0.17 GB to
  2.90 GB (2L) because every read is kept in memory, including those that
  will not enter the pileup.
- **Unverified part of the output.** 608,373 unique alignments on 2L
  outside STAR's clean set have not been verified; the same holds,
  analogously, for the whole genome.
- **Measurements from a single run.** Time and memory are not an average
  of multiple repetitions, at either scope.
- **The macOS fix for STAR was checked on only one configuration**,
  untested on Linux
  ([Known issue](#known-issue-star-on-macos-apple-silicon)).

### Assessment of practical applicability

**For production RNA-seq analysis, the implementation is not
applicable.** The decisive reason is not accuracy, but the lack of
splicing support: STAR recorded 1,325,356 splice junctions on the same
dataset, and reads containing them cannot be correctly aligned by
ungapped alignment. On top of that, there is no alignment-confidence
score, without which results cannot be filtered in downstream analysis.

**For DNA data without indels, the implementation is applicable in
principle**, with the caveat that this scenario was not tested in this
work. The measured concordance concerns ungapped alignment sites, which
matches that class of problem, but no verification against DNA-seq data
was performed.

**As a correctness check of the implemented algorithms, the result is
successful, at both scopes.** The algorithms were run over the full real
dataset of 14,988,127 reads, both against one chromosome and against the
whole genome, not just against test examples, with concordance with STAR
that matches expectations for the class of reads handled in both cases.

Speed is a limitation of the implementation, not the algorithm. The
procedure is single-threaded and written in Python, with a C library
used only for suffix array construction, while the per-read search is a
plain Python loop.

### Known issue: STAR on macOS (Apple Silicon)

**Configuration where the issue was observed:** STAR 2.7.11b (Homebrew
`rna-star` build), macOS 26.5 (build 25F71), Apple Silicon (arm64), Apple
clang 15.0.0 (`arm64-apple-darwin25.5.0`), Apple's implementation of the
standard C++ library (libc++). On this configuration, `alignReads`
silently returns 0 aligned reads — the process exits with code 0, looking
like a successful run. Whether the same issue occurs on other macOS
versions, architectures, or STAR versions was not checked, so the finding
is not generalized beyond the configuration above, nor to every Homebrew
build on every Apple Silicon system.

**Cause.** `std::stringbuf::pubsetbuf()` is a standard C++ function, but
its specific behavior for `stringbuf` is left implementation-defined by
the C++ standard — there is no guarantee that the given `char*` becomes
the stream's actual internal buffer. STAR 2.7.11b relies on exactly that
happening. In the build tested, that worked with GNU's libstdc++ (a
behavior that particular implementation chooses, not something the
standard mandates), while Apple's libc++ treats the same call as a no-op
and never links the external buffer. Consequence: STAR's buffer for
reading reads and its buffer for writing SAM records stay empty
regardless of the streams' actual content, and the program still exits
with code 0 instead of an error.

**Fix.**
[`star_comparison/patches/star_2.7.11b_macos_libcxx_fix.patch`](star_comparison/patches/star_2.7.11b_macos_libcxx_fix.patch)
replaces the `pubsetbuf`-aliasing approach with an explicit
`std::stringbuf::str()` call — copying the content into the stringbuf
instead of relying on aliasing an external buffer (see also
[`star_comparison/scripts/build_star_macos.sh`](star_comparison/scripts/build_star_macos.sh)).
`stringbuf::str()` is standard, portable behavior, so the same functional
result is expected on Linux/libstdc++; this was not tested as part of
this work, so identical performance (time, memory copying) on that
platform is not claimed — only the expected functional equivalence.

### Reproducing the pipeline

The input data (`Fly/`) and intermediate results larger than a few
hundred KB are not in git (see `.gitignore`) — `.sample.*` versions are
committed (header + first ~2000 rows), and the full files are
regenerated locally.

**Chromosome 2L:**

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_chromosome.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  results/chr2L/my_alignments_full.tsv \
  2L results/chr2L/comparison_chr2L.csv --tolerance 0

python3 analyze_unmapped_chr2L.py
```

**Whole genome** (each command is run from the repository root, in
order):

```bash
star_comparison/scripts/build_star_macos.sh                     # build the patched STAR -> bin/STAR
star_comparison/scripts/run_star.sh                              # genome index + SAM     -> results/star/
star_comparison/scripts/filter_clean_reads.py  <sam> <out>       # -> results/filtered/clean_reads.tsv
star_comparison/scripts/make_fasta_input.py    <tsv> <out>       # -> results/my_aligner/clean_reads_input.fasta
star_comparison/scripts/build_my_index.py                        # -> my_index/
star_comparison/scripts/run_my_aligner.py      <fasta> <out>     # -> results/my_aligner/my_alignments.tsv
star_comparison/scripts/compare_coordinates.py <star> <my> <out> --mismatches-out <out2>
star_comparison/scripts/inspect_mismatches.py  <out2> <clean> <out3>
star_comparison/scripts/extract_nh.py          <sam> <out>       # -> results/filtered/clean_reads_nh.tsv
```

A short description of each step, and the exact command form with real
paths, is in [`star_comparison/README.md`](star_comparison/README.md),
section "Pipeline steps". Both pipelines are deterministic — rerunning
them on the same input data gives identical numbers.

### Conclusion

The own implementation of suffix array construction and read alignment
was run over the full real dataset, at two independent scopes — one
chromosome and the whole genome — and compared against STAR as a
reference solution.

At the chromosome-2L scope: of the 1,755,830 reads STAR cleanly aligned,
the own aligner placed 1,713,194 at an identical coordinate, a strict
accuracy of 97.5717%, 99.9996% among uniquely placed reads. At the
whole-genome scope: of 9,299,988 clean reads, strict accuracy 96.42%,
99.9999% among uniquely placed reads. The seven positional disagreements
from the 2L scope are identical to seven of the total nine disagreements
at the whole-genome scope — two independently run comparisons give the
same finding. The mismatch-count deviation is, in 99.77% of cases, a
consequence of the difference in the counting convention for `N` bases,
not a placement error.

From this follows the conclusion that the algorithm correctly places
reads of the class it is meant for, i.e. reads without gaps in the
alignment, at both tested scopes. The conclusion holds for the verified
subset at each scope; reads outside STAR's clean set (608,373 on 2L,
analogously on the whole genome) remain unverified, and the accuracy
claim cannot be extended to them without further measurement.

The mutation-calling algorithm, run at the single-chromosome scope,
returned 10,445 calls. These results have not been independently
validated and cannot be presented as confirmed biological variants,
especially since more than a third of the calls rest on coverage below
five reads.

A more complete assessment would require further tests: a comparison
over simulated reads with known origin, which would give an absolute
rather than relative accuracy measure; verification of the part of the
output not covered by STAR's clean set, at both scopes; and independent
validation of the called mutations against a known variant database.

---

## Hrvatski

**Sažetak.** Vlastita implementacija poravnanja readova (`genome_index/`)
uspoređena je sa STAR-om (v2.7.11b) u dva opsega, oba nad istim stvarnim
skupom podataka: svih 14.988.127 single-end RNA-seq reada *Drosophila
melanogaster* iz `Fly/IFM48h_1.fastq`. U opsegu jednog kromosoma (2L,
23.513.712 baza, preko `run_pipeline.py`) stroga podudarnost koordinata sa
STAR-om iznosi 97,5717 %, a među readovima koje je vlastiti aligner
jedinstveno smjestio 99,9996 %. U opsegu cijelog genoma (svih 1870
kontiga spojenih u jedan indeks, preko `star_comparison/scripts/`)
podudarnost iznosi 96,42 % strogo, odnosno 99,9999 % među jedinstveno
smještenim readovima. Oba mjerenja se slažu unutar granica koje njihove
metodologije dopuštaju i međusobno se potvrđuju: isti skup pozicijskih
neslaganja pojavljuje se u oba (7 od 9 ukupnih neslaganja nalazi se na
kromosomu 2L). Vlastiti algoritam pronalaženja mutacija, pokrenut samo u
opsegu jednog kromosoma, pozvao je 10.445 točkastih mutacija koje nisu
neovisno provjerene. U izvještaju se razlikuju **izmjerene vrijednosti** od
**izvedenih zaključaka**, koji su označeni kao takvi.

### Sadržaj

1. [Cilj usporedbe](#cilj-usporedbe)
2. [Podaci](#podaci)
3. [Dva opsega usporedbe](#dva-opsega-usporedbe)
4. [Metodologija](#metodologija)
5. [Uvjeti mjerenja](#uvjeti-mjerenja)
6. [Korištene naredbe](#korištene-naredbe)
7. [Rezultati na kromosomu 2L](#rezultati-na-kromosomu-2l)
8. [Rezultati na cijelom genomu](#rezultati-na-cijelom-genomu)
9. [Točnost i pouzdanost vlastitog algoritma](#točnost-i-pouzdanost-vlastitog-algoritma)
10. [Podudaranja i odstupanja u odnosu na STAR](#podudaranja-i-odstupanja-u-odnosu-na-star)
11. [Mogući uzroci uočenih razlika](#mogući-uzroci-uočenih-razlika)
12. [Rezultati pronalaženja mutacija](#rezultati-pronalaženja-mutacija)
13. [Prednosti i ograničenja](#prednosti-i-ograničenja)
14. [Procjena praktične primjenjivosti](#procjena-praktične-primjenjivosti)
15. [Poznati problem: STAR na macOS-u (Apple Silicon)](#poznati-problem-star-na-macos-u-apple-silicon)
16. [Reprodukcija pipelinea](#reprodukcija-pipelinea)
17. [Zaključak](#zaključak)

### Cilj usporedbe

Vlastita implementacija poravnanja readova nema poznato točno rješenje nad
kojim bi se mogla provjeriti: korišteni readovi su stvarni RNA-seq podatci,
a ne simulirani, pa stvarno genomsko podrijetlo pojedinog reada nije
poznato. Zbog toga se ispravnost provjerava neizravno, usporedbom s
neovisno razvijenim i široko korištenim alignerom STAR.

Cilj usporedbe je odgovoriti na pitanje: **smješta li vlastiti aligner
readove na iste koordinate kao STAR, i u kojoj mjeri.** STAR koristi
nekomprimirani suffix array za pretragu (ne FM-index/BWT), arhitekturalno
srodan pristupu korištenom u vlastitoj implementaciji, pa je usporedba
relevantna i za pitanje ispravnosti algoritma, ne samo brzine.

STAR se u ovoj usporedbi koristi kao **referentni (baseline) aligner**, ne
kao dokazana ground-truth vrijednost. Rezultat se izražava kao podudarnost
sa STAR-om (concordance), a ne kao apsolutna biološka točnost.

Usporedba ne provjerava pronalaženje mutacija. STAR ne poziva varijante, pa
rezultati iz [tog dijela](#rezultati-pronalaženja-mutacija) nisu neovisno
potvrđeni.

### Podaci

| Stavka | Vrijednost |
|---|---|
| Referentni genom | `Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa`, 146.261.227 B, 1870 zapisa |
| Readovi | `Fly/IFM48h_1.fastq`, 3.663.357.012 B, 14.988.127 single-end readova duljine 76 baza |
| Referentni aligner | STAR 2.7.11b |
| Testirani aligner | vlastiti suffix-array aligner (`genome_index/`) |

Korišten je cjelovit skup readova u oba opsega usporedbe, bez uzorkovanja.

### Dva opsega usporedbe

Vlastiti aligner je testiran u dva odvojena postava, koji se ne mogu
izravno svesti na jedan broj jer indeksiraju različitu količinu genoma:

| | Kromosom 2L | Cijeli genom |
|---|---|---|
| Program | `run_pipeline.py` (root) | `star_comparison/scripts/` |
| Indeks | suffix array samo kromosoma 2L (23.513.712 baza) | jedan suffix array svih 1870 kontiga spojenih u niz (144 Mb), s N-spacerima od 200bp između kontiga |
| Ulaz u usporedbu | svih 14.988.127 readova, STAR filtriran na kromosom 2L | svih 14.988.127 readova, STAR filtriran na cijeli genom |
| STAR-ov nazivnik usporedbe | 1.755.830 (čista poravnanja na 2L) | 9.299.988 (čista poravnanja na bilo kojem kontigu) |
| Namjena | odgovara doslovno zahtjevu rada (jedan kromosom, jedan suffix array) | provjera na skali cjelokupnog genoma, gdje repetitivnost i broj kontiga imaju veći utjecaj |

Oba postava koriste isti kod iz `genome_index/` za poravnanje i isti
kriterij podudarnosti; razlikuje se samo opseg indeksiranog genoma i
posljedično broj readova koji uopće mogu biti obuhvaćeni usporedbom. Ondje
gdje su nalazi konzistentni između dva opsega (npr. [pozicijska
neslaganja](#podudaranja-i-odstupanja-u-odnosu-na-star)), to je naznačeno
kao dodatna potvrda.

### Metodologija

#### Zajednički kriteriji: SAM polja i CIGAR string

Standardni SAM format (zaglavlje s `@HD`/`@SQ`, zatim jedan red po
alignmentu s 11 obaveznih polja) nije dalje opisan ovdje. Bitni su za ovu
usporedbu: `RNAME`+`POS` (koordinata koja se uspoređuje), `CIGAR`
(kriterij filtriranja, niže), i opcionalni tagovi `NH` (broj lokusa na
koje je read mapiran), `nM` (broj nepodudaranja u alignmentu) i `AS`
(STAR-ov lokalni alignment score). Primjer retka:

```
ILLUMINA-075005_0053_FC:3:1:17967:18033#0  0  2L  3515966  255  76M  *  0  0  GAGAACTTCGCC...  hhhgehghhggg...  NH:i:1  HI:i:1  AS:i:74  nM:i:0
```

CIGAR opisuje poravnanje kao niz parova `<broj><operacija>`, čitano slijeva
nadesno duž reada:

| Op | Značenje | Troši baze reada | Troši baze reference |
|---|---|---|---|
| `M` | poravnata baza (podudaranje ili nepodudaranje) | da | da |
| `I` | insercija u readu (baze kojih nema u referenci) | da | ne |
| `D` | delecija (baze reference kojih nema u readu) | ne | da |
| `N` | preskočena regija reference — kod RNA-seq to je intron (splice) | ne | da |
| `S` | soft clipping — kraj reada koji nije poravnat, ali je prisutan u SEQ | da | ne |
| `H` | hard clipping — odrezani kraj, nije ni prisutan u SEQ | ne | ne |

`M` **ne znači** da se baze podudaraju, nego samo da su poravnate bez
pomaka — nepodudaranja se broje odvojeno, u `nM` tagu. Read s CIGAR-om
`76M` može imati poneku point-mutaciju ili grešku sekvenciranja; to je
upravo scenarij koji vlastiti aligner podržava (mismatch-tolerant, bez
indela). Najveći opažen broj nepodudaranja (`nM`) u cijelom SAM-u iz
whole-genome runa je 10 na 76bp read (15.124.860 zapisa provjereno).

Kriterij "poravnan bez prekida", korišten u oba opsega usporedbe: cijeli
CIGAR mora biti točno jedna `M` operacija čiji broj odgovara duljini reada
— za readove od 76bp to je `76M`. Svaki `S`, `H`, `I`, `D` ili `N` u
CIGAR-u automatski isključuje read iz usporedbe (soft/hard clip, insercija,
delecija, ili splice/intron gap). Primjeri koji se odbacuju:

```
...:17532:18022#0  16  X   2001069   255  75M1S       ...   -> 1 baza soft-clipana
...:19344:18030#0  16  2R  13969542  255  61M69N15M   ...   -> read prelazi preko introna (69bp splice)
```

Uz to se u oba opsega isključuju sekundarna (`FLAG & 0x100`), supplementary
(`FLAG & 0x800`) i unmapped (`FLAG & 0x4`) STAR poravnanja.

#### Napomena o oznaci broja nepodudaranja

Zadatak spominje SAM oznaku `NM:i:` za provjeru potpuno identičnih
poravnanja. Provjerom cijelog SAM izlaza (oba opsega) utvrđeno je da STAR u
ovoj konfiguraciji tu oznaku ne zapisuje, nego vlastitu oznaku `nM:i:`.
Usporedba broja nepodudaranja stoga posvuda koristi `nM:i:`.

#### Kriterij točnosti i tolerancija

Poravnanje se smatra podudarnim ako se podudaraju kromosom i početna
koordinata, uz toleranciju od 0 baza. Tolerancija 0 primijenjena je bez
ublažavanja — provjereno je (opseg cijelog genoma, 9.299.988 readova) da
tolerancija ±1 baze daje identičan broj podudaranja kao tolerancija 0
([detalji](#provjera-tolerancije-koordinata)), čime je isključena sustavna
pogreška u pretvorbi koordinata (0-bazirana pozicija suffix arraya →
1-bazirana SAM pozicija).

U opsegu kromosoma 2L, uvjet "isti kromosom" zadovoljen je po konstrukciji,
jer indeks vlastitog alignera sadrži samo 2L — pogrešno smještanje reada s
drugog kromosoma stoga bi se očitovalo kao pogrešna koordinata unutar 2L,
ne kao pogrešan kromosom. U opsegu cijelog genoma kromosom se uspoređuje
eksplicitno, jer indeks obuhvaća sve kontige.

#### Postupak — kromosom 2L

Vlastiti pipeline (`run_pipeline.py`) gradi suffix array samo kromosoma 2L
i poravnava sve readove nad njim. STAR gradi indeks nad cijelim genomom i
poravnava iste readove, uz SAM izlaz. Iz STAR-ova SAM izlaza uzimaju se
poravnanja koja se odnose na kromosom 2L, zadovoljavaju
[zajednički kriterij](#zajednički-kriteriji-sam-polja-i-cigar-string) i
nisu sekundarna/supplementary/unmapped. Usporedba se izvodi skriptom
[`star_comparison/scripts/compare_chromosome.py`](star_comparison/scripts/compare_chromosome.py),
koja filtriranje SAM-a i usporedbu radi u jednom prolasku.

#### Postupak — cijeli genom

1. **STAR alignment.** STAR poravnava svih 14.988.127 reada nad cijelim
   referentnim genomom (`--outSAMtype SAM`), s punom podrškom za splicing i
   indele (zadane postavke, bez override-a).
2. **Filtriranje čistih reada**
   ([`star_comparison/scripts/filter_clean_reads.py`](star_comparison/scripts/filter_clean_reads.py))
   po [zajedničkom kriteriju](#zajednički-kriteriji-sam-polja-i-cigar-string).
3. **Poravnanje vlastitim alignerom** nad jednim suffix arrayem izgrađenim
   nad svim 1870 kontiga spojenim u jedan niz — N-spaceri od 200bp između
   kontiga sprječavaju lažno poravnanje preko granice kontiga
   ([`star_comparison/scripts/build_my_index.py`](star_comparison/scripts/build_my_index.py),
   [`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py),
   `max_mismatches=3`).
4. **Usporedba koordinata**
   ([`star_comparison/scripts/compare_coordinates.py`](star_comparison/scripts/compare_coordinates.py))
   između (RNAME, POS) iz STAR-a i vlastitog alignera.

Usporedba namjerno nije simetrična u opsegu: STAR obrađuje sve readove sa
splicingom i indelima, a vlastiti aligner samo podskup koji sam ne
podržava te operacije. To je svjesna odluka — mjeri se pouzdanost
*pozicioniranja* na readovima koje bi vlastiti aligner uopće trebao moći
obraditi, ne cjelokupna RNA-seq funkcionalnost STAR-a.

### Uvjeti mjerenja

| Stavka | Vrijednost |
|---|---|
| Stroj | Apple M2 (Mac14,15), 8 jezgara, 8,6 GB RAM |
| OS | macOS 26.5 (build 25F71) |
| Kompajler | Apple clang 15.0.0 |
| STAR | 2.7.11b, izgrađen iz izvornog koda uz popravak za macOS (vidi [Poznati problem](#poznati-problem-star-na-macos-u-apple-silicon)) |
| Python / numpy / pydivsufsort / pytest | 3.13.5 / 2.5.1 / 0.0.20 / 9.1.1 |
| Mjerenje vremena/memorije | `/usr/bin/time -l` (BSD `time`, macOS); peak RSS je `maximum resident set size` u bajtovima |
| Datum STAR runa (whole-genome) | 2026-09-06, prema `IFM48h_1.Log.out` |

Vrijednosti se odnose na pojedinačna izvođenja na jednom stroju i nisu
prosjek više ponavljanja. Runovi za oba opsega izvedeni su odvojeno, u
različitim sesijama rada na istom stroju.

### Korištene naredbe

**STAR — izgradnja indeksa i poravnanje nad cijelim genomom** (koristi se
za oba opsega usporedbe, jer STAR uvijek indeksira cijeli genom):

```bash
STAR --runMode genomeGenerate --runThreadN 6 \
  --genomeDir star_comparison/index \
  --genomeFastaFiles Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
  --genomeSAindexNbases 12 --genomeChrBinNbits 16

STAR --runMode alignReads --runThreadN 6 \
  --genomeDir star_comparison/index \
  --readFilesIn Fly/IFM48h_1.fastq \
  --outSAMtype SAM \
  --outFileNamePrefix star_comparison/results/star/IFM48h_1.
```

**Vlastiti pipeline — kromosom 2L:**

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_chromosome.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  results/chr2L/my_alignments_full.tsv \
  2L results/chr2L/comparison_chr2L.csv --tolerance 0
```

**Vlastiti pipeline — cijeli genom** (puni redoslijed s ulazima/izlazima
svakog koraka u [Reprodukciji pipelinea](#reprodukcija-pipelinea)):

```bash
python3 star_comparison/scripts/filter_clean_reads.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  star_comparison/results/filtered/clean_reads.tsv

python3 star_comparison/scripts/build_my_index.py

python3 star_comparison/scripts/run_my_aligner.py \
  star_comparison/results/my_aligner/clean_reads_input.fasta \
  star_comparison/results/my_aligner/my_alignments.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_coordinates.py \
  star_comparison/results/filtered/clean_reads.tsv \
  star_comparison/results/my_aligner/my_alignments.tsv \
  star_comparison/results/comparison/comparison.csv \
  --mismatches-out star_comparison/results/comparison/mismatches.csv \
  --tolerance 0
```

### Rezultati na kromosomu 2L

#### Ishodi poravnanja (2L)

Vlastiti aligner nad svih 14.988.127 readova, indeks kromosoma 2L,
`max_mismatches = 3`:

| Ishod | Broj readova | Udio |
|---|---|---|
| `unique` | 2.321.574 | 15,49 % |
| `multi` | 29.104 | 0,19 % |
| `unmapped` | 12.637.449 | 84,32 % |
| `skipped` | 0 | 0,00 % |

Visok udio nepovezanih readova očekivan je i **nije mjera pogreške**:
readovi potječu iz cijelog genoma, a indeks obuhvaća samo kromosom 2L, čija
duljina čini približno 16 % ukupne duljine sklopa.

STAR nad cijelim genomom, prema vlastitom izvještaju
`IFM48h_1.Log.final.out`:

| Mjera | Vrijednost |
|---|---|
| Jedinstveno poravnano | 14.424.378 (96,24 %) |
| Poravnano na više lokusa | 197.559 (1,32 %) |
| Nepovezano | 342.852 (2,29 %) |
| Broj splice spojeva | 1.325.356 |
| Stopa nepodudaranja po bazi | 0,64 % |

Broj splice spojeva potvrđuje da skup sadrži znatan udio readova koji
premošćuju granice egzona i koje vlastiti aligner po konstrukciji ne može
ispravno poravnati.

#### Usporedba koordinata (2L)

Nakon filtriranja, STAR ima 1.755.830 čistih primarnih poravnanja na
kromosomu 2L. To je nazivnik usporedbe.

| Mjera | Broj | Udio |
|---|---|---|
| Oba alata poravnala jedinstveno | 1.713.201 | 97,58 % |
| — ista koordinata | 1.713.194 | 97,5717 % |
| — različita koordinata | 7 | 0,0004 % |
| — ista orijentacija | 1.713.201 | 100,00 % |
| — različita orijentacija | 0 | 0,00 % |
| — isti broj nepodudaranja (`nM`) | 1.710.110 | 99,82 % |
| Vlastiti aligner: `multi` | 10.375 | 0,59 % |
| Vlastiti aligner: `unmapped` | 32.254 | 1,84 % |
| Bez zapisa u izlazu vlastitog alignera | 0 | 0,00 % |

Točnost = broj readova s istom koordinatom / broj čistih poravnanja STAR-a
= 1.713.194 / 1.755.830 = **97,5717 %**.

Točnost među readovima koje je vlastiti aligner jedinstveno smjestio =
1.713.194 / 1.713.201 = **99,9996 %**.

#### Vrijeme izvođenja i memorija (2L)

| Korak | Alat | Dretvi | Readova | real | Peak RSS |
|---|---|---|---|---|---|
| Izgradnja indeksa, cijeli genom | STAR | 6 | — | 53,49 s | 1,54 GB |
| Poravnanje, cijeli genom | STAR | 6 | 14.988.127 | 97,08 s | 1,98 GB |
| Izgradnja suffix arraya, 2L | vlastiti | 1 | — | 1,19 s | — |
| Poravnanje, 2L | vlastiti | 1 | 14.988.127 | 1529,34 s | 0,17 GB |
| Poravnanje i mutacije, 2L | vlastiti | 1 | 14.988.127 | 1842,24 s | 2,90 GB |

Propusnost poravnanja vlastitog alignera iznosila je 9810 readova u
sekundi. Vrijednosti za vlastiti aligner odnose se na izgradnju indeksa
jednog kromosoma i pretragu unutar njega, dok se STAR-ove odnose na cijeli
genom uz šest dretvi. **Brojke stoga nisu izravno usporedive** i navode se
odvojeno; izravna usporedba pri jednakom opsegu pretrage je u
[Vremenu izvršavanja](#vrijeme-izvršavanja-cijeli-genom) niže.

### Rezultati na cijelom genomu

#### Filtriranje SAM izlaza

Raspodjela CIGAR operacija po svih 14.621.937 primarnih alignmenata,
izračunata nezavisnim `awk` prolazom kroz SAM (poklapa se s brojkama koje
prijavljuje Python filter,
[`star_comparison/scripts/filter_clean_reads.py`](star_comparison/scripts/filter_clean_reads.py)):

| Kategorija | Broj reada | Udio primarnih |
|---|---|---|
| čisti `<duljina>M` (ulaz u ostatak usporedbe) | 9.299.988 | 63,60 % |
| sadrži `S` (soft clipping) | 4.236.843 | 28,98 % |
| sadrži `N` (splice/intron) | 1.338.755 | 9,16 % |
| sadrži `D` (delecija) | 51.549 | 0,35 % |
| sadrži `I` (insercija) | 46.401 | 0,32 % |
| sadrži `H` (hard clipping) | 0 | 0,00 % |

Kategorije se preklapaju — jedan read može istovremeno imati splice i
clipping — pa udjeli ne zbrajaju na 100 %.

#### Ishodi poravnanja (cijeli genom)

Nad 9.299.988 čistih reada, vlastiti aligner
([`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py),
`max_mismatches=3`) svaki read klasificira kao jedinstveno poravnat
(**unique**), **multi** (više poravnanja) ili **unmapped** (nijedno):

| Ishod | Broj reada | Udio |
|---|---|---|
| unique | 8.966.941 | 96,42 % |
| multi-mapped | 179.897 | 1,93 % |
| unmapped | 153.150 | 1,65 % |

Udio "unique" ishoda naziva se **stopom jedinstvenog pozicioniranja**
(coverage). Sam po sebi ne govori ništa o točnosti — samo o tome koliko
često aligner uopće daje jedan odgovor.

#### Podudarnost koordinata sa STAR-om (cijeli genom)

Za svaki čisti read uspoređuje se (RNAME, POS) vlastitog alignera sa
STAR-ovim (tolerancija 0 baza):

| Mjera | Brojnik / nazivnik | Vrijednost | Definicija |
|---|---|---|---|
| Stopa jedinstvenog pozicioniranja (coverage) | 8.966.941 / 9.299.988 | 96,42 % | udio čistih reada za koje vlastiti aligner vraća točno jednu poziciju |
| Podudarnost sa STAR-om među unique pozivima | 8.966.932 / 8.966.941 | 99,9999 % | od reada koje je vlastiti aligner jedinstveno smjestio, koliki udio se poklapa sa STAR-om |
| Stroga koordinatna podudarnost | 8.966.932 / 9.299.988 | 96,42 % | podudaranja u odnosu na *sve* čiste readove — multi-mapped i unmapped readovi računaju se kao nepodudarni |

Ove tri mjere namjerno se ne svode na jedan broj "točnosti", jer odgovaraju
na različita pitanja: prva mjeri pokrivenost, druga točnost uz
pretpostavku da je aligner dao odgovor, treća kombinira oboje. Razlika
između stroge podudarnosti (96,42 %) i gotovo savršene podudarnosti među
unique pozivima (99,9999 %) pokazuje da manjak dolazi gotovo isključivo iz
nedostatka pokrivenosti, ne iz pogrešno određenih koordinata — od
8.966.941 unique poziva, samo 9 se ne slaže sa STAR-om
([detalji](#provjera-devet-neslaganja)).

#### Podjela prema STAR-ovom NH tagu

STAR sam prijavljuje `NH:i:N` kad je read podjednako dobro poravnat na N
mjesta te jedno bira kao primarno prema vlastitim internim pravilima. Od
9.299.988 čistih reada, 109.411 (1,18 %) ima STAR `NH≥2`:

| Skup (prema STAR-u) | Broj reada | Vlastiti aligner: unique | — multi | — unmapped | Podudarnost sa STAR-om |
|---|---|---|---|---|---|
| `NH:i:1` (STAR jedinstven) | 9.190.577 | 8.962.784 (97,52 %) | 78.951 (0,86 %) | 148.842 (1,62 %) | 97,52 % strogo, 99,9999 % među unique pozivima |
| `NH≥2` (STAR multi-mapped) | 109.411 | 4.157 (3,80 %) | 100.946 (92,26 %) | 4.308 (3,94 %) | 3,80 % strogo, 99,9759 % među unique pozivima |

Na podskupu gdje je i STAR siguran (`NH:i:1`) vlastiti aligner postiže
97,52 % strogu podudarnost — više od ukupnog prosjeka (96,42 %), jer se
isključuju readovi koji su dvosmisleni po definiciji problema. Na podskupu
gdje je STAR sam prijavio više mogućih lokusa, vlastiti aligner također
najčešće (92,26 %) vraća multi-mapping umjesto jedne pozicije: dva neovisna
alata se u velikoj većini slučajeva slažu koji su readovi dvosmisleni.

Jedini "unique" poziv unutar `NH≥2` skupa koji se ne slaže sa STAR-om (read
`...3:14:17033:20323#0`) provjeren je i protiv STAR-ova sekundarnog
lokusa: STAR je prijavio 3L:28.103.368 (primarni) i X:5.752.974
(sekundarni, `NH:i:2`), a vlastiti aligner treću poziciju, 3L:18.780.015 —
ne poklapa se ni s jednom od te dvije.

#### Provjera devet neslaganja

Za svih 9 reada gdje se vlastiti aligner odlučio na jedinstvenu poziciju
koja se ne poklapa sa STAR-ovom, izbrojan je stvarni broj nepodudaranja
(Hammingova udaljenost) reada naspram referentnog genoma na obje
prijavljene pozicije
([`star_comparison/scripts/inspect_mismatches.py`](star_comparison/scripts/inspect_mismatches.py),
rezultat u
[`star_comparison/results/comparison/mismatches_verified.csv`](star_comparison/results/comparison/mismatches_verified.csv)):

| Read | STAR pozicija | STAR: nepodudaranja | Pozicija vlastitog alignera | Vlastiti: nepodudaranja |
|---|---|---|---|---|
| `...3:8:10065:3458#0` | 3L:3.901.348 | 5 | 3L:3.901.576 | 3 |
| `...3:14:17033:20323#0` | 3L:28.103.368 | 9 | 3L:18.780.015 | 3 |
| preostalih 7 reada (identičan obrazac) | 2L:14.743.463 | 4 | 2L:14.743.493 | 3 |

Ovih istih 7 readova, na istim pozicijama i s istim brojem nepodudaranja,
pojavljuje se i u zasebnoj usporedbi ograničenoj na kromosom 2L
([Odstupanja u koordinati](#odstupanja-u-koordinati)) — dvije neovisno
izvedene usporedbe daju identičan nalaz.

Na svih 9 pozicija, pozicija vlastitog alignera ima manju Hammingovu
udaljenost od pozicije koju je prijavio STAR (brojanje se poklapa sa
STAR-ovim vlastitim `nM` tagom — `nM:i:5`, `nM:i:9`, `nM:i:4` za prikazane
readove — pa je isključena greška u samom brojanju). Osam od 9 reada ima
STAR `NH:i:1`, tj. STAR ih sam smatra jedinstveno mapiranima.

**Ovo je izmjerena činjenica** (manja Hammingova udaljenost = bolje
rangirana pozicija prema tom kriteriju), **ne dokaz** o stvarnom biološkom
izvoru reada: prava koordinata podrijetla nije poznata, a niža Hammingova
udaljenost može doći od prave lokacije, ali i od slučajnog poklapanja s
ponovljenom sekvencom drugdje u genomu.

Mogući razlozi STAR-ova izbora druge pozicije uključuju razlike u
generiranju kandidata, seed-heuristikama i bodovanju u odnosu na
pigeonhole seed-and-extend korišten ovdje — to su hipoteze, ne potvrđen
uzrok. Budući da svih 9 pozicija ima CIGAR `76M` na obje strane usporedbe,
razlika se konkretno ne može pripisati STAR-ovom splice modelu.

Metoda brojanja nepodudaranja provjerena je i na 2000 nasumičnih reada iz
skupa na kojem se alati slažu (sanity check, ne dokaz o cijelom skupu):
73 % ima 0 nepodudaranja na prijavljenoj poziciji, 94 % najviše 1.

#### Provjera tolerancije koordinata

Usporedba prima parametar `--tolerance` (najveća dopuštena razlika
pozicija koja se još broji kao podudaranje). Pokrenuto na svih 9.299.988
čistih reada s tolerancijom 0 i s tolerancijom ±1 baza: oba daju identičnih
8.966.932 podudaranja. Proširenje tolerancije na susjednu bazu ne mijenja
nijedan rezultat, čime je isključena sustavna off-by-one greška u
pretvorbi koordinata.

#### Vrijeme izvršavanja (cijeli genom)

| Faza | Alat | Threadova | Reada | real | user | sys | Throughput |
|---|---|---|---|---|---|---|---|
| Genome index | STAR | 6 | — | 53,49 s | 189,83 s | 15,47 s | — |
| Genome index | vlastiti aligner | 1 | — | 9,06 s | 8,35 s | 0,48 s | — |
| Alignment | STAR | 6 | 14.988.127 | 97,08 s | 414,57 s | 52,20 s | ~154.400 reada/s |
| Alignment | vlastiti aligner | 1 | 9.299.988 | 1009,52 s (16m50s) | 968,57 s | 24,54 s | ~9.210 reada/s |

Broj threadova potvrđen je iz `--runThreadN 6` u pozivu STAR-a i neovisno
iz omjera user/real vremena: STAR-ov alignment ima user/real ≈ 4,27
(očekivano za 6 threadova uz dio neparalelnog rada), dok vlastiti aligner
ima user/real ≈ 0,96 (u skladu s time da
[`star_comparison/scripts/run_my_aligner.py`](star_comparison/scripts/run_my_aligner.py)
ne koristi ni threading ni multiprocessing).

Na razini pojedinačnog reada, STAR je ~16,8× brži tijekom alignmenta (6,5
µs/read naspram 108,6 µs/read). Ovaj omjer **ne predstavlja normaliziranu
usporedbu algoritama po jezgri** — mjeren je u threadovima asimetričnoj
konfiguraciji (6 naspram 1), pa dio razlike dolazi od paralelizma, ne
isključivo od algoritma. Vlastiti aligner je uz to jednostavan,
single-threaded Python (`divsufsort` C backend koristi se samo za
izgradnju suffix arraya, sama pretraga po readu je čista Python petlja).

#### Peak memory (cijeli genom)

Vrijednost je `maximum resident set size` iz `/usr/bin/time -l`. Na
macOS-u je izražena u bajtovima:

| Faza | Alat | Peak RSS (bajtova) | Peak RSS (GB) |
|---|---|---|---|
| Genome index | STAR | 1.535.000.576 | 1,54 GB |
| Genome index | vlastiti aligner | 891.928.576 | 0,89 GB |
| Alignment | STAR | 1.978.810.368 | 1,98 GB |
| Alignment | vlastiti aligner | 1.562.148.864 | 1,56 GB |

Na ovom stroju i u ovom mjerenju, peak RSS oba alata je istog reda
veličine (1–2 GB), a vlastiti aligner je u koraku poravnanja nešto
štedljiviji. Ovo nije opća tvrdnja o memory footprintu obje implementacije
— samo opis jednog mjerenja na ovom genomu (144 Mb): oba alata dominantno
troše memoriju na genomski indeks.

### Točnost i pouzdanost vlastitog algoritma

**Izmjereno.** Od 1.713.201 reada koje su oba alata smjestila jedinstveno
na 2L, 1.713.194 ima identičnu koordinatu, a svih 1.713.201 ima istu
orijentaciju. Nijedan read iz STAR-ova skupa nije izostao iz izlaza
vlastitog alignera.

**Izmjereno.** Dva neovisna izvođenja vlastitog pipelinea nad 2L, s
pozivanjem mutacija i bez njega, dala su identične ishode poravnanja
(2.321.574 `unique`, 29.104 `multi`, 12.637.449 `unmapped`). Algoritam je
determinističan i ne ovisi o redoslijedu obrade.

**Izmjereno.** Neovisna usporedba na razini cijelog genoma, nad 9.299.988
čistih readova, dala je strogu točnost od 96,42 % i podudarnost od
99,9999 % među jedinstveno smještenim readovima. Sedam pozicijskih
neslaganja pronađenih u toj usporedbi pojavljuje se i u usporedbi
ograničenoj na 2L, na istim pozicijama i s istim brojem nepodudaranja —
dvije neovisno izvedene usporedbe daju identičan nalaz.

**Izmjereno.** Jedinični testovi (16 testova) uspoređuju SA-IS s naivnim
sortiranjem sufiksa na nasumičnim nizovima, a poravnanje s iscrpnim
postupkom pretrage. Pokriveni su i rubni slučajevi: prazan ulaz, read bez
poravnanja, readovi s negativnog lanca, neispravni ulazni znakovi, read
dulji od reference te odbacivanje mutacija pri nedovoljnoj pokrivenosti.
Skripta `moj_primjer.py` dodatno provjerava cjelovit postupak nad podatcima
s poznatim rješenjem.

**Izvedeni zaključak.** Navedeno zajedno govori u prilog tome da je jezgra
algoritma — binarna pretraga nad suffix arrayem uz seed-and-extend s
pigeonhole jamstvom — implementirana ispravno za razred readova koje
obrađuje (bez indela, delecija i splicinga).

**Ograničenje ovog zaključka.** U opsegu 2L, točnost je izmjerena samo na
readovima koje je STAR čisto poravnao na tom kromosomu. Vlastiti aligner
jedinstveno je smjestio 2.321.574 reada, od kojih 1.713.201 ulazi u
usporedbu; preostalih 608.373 **nije provjereno nijednim neovisnim
izvorom**. U opsegu cijelog genoma analogna neprovjerena razlika je manja u
relativnom smislu (usporedba pokriva 9.299.988 od svih čistih reada), ali
i dalje isključuje sve readove sa splicingom i indelima, koje vlastiti
aligner po konstrukciji ne obrađuje.

### Podudaranja i odstupanja u odnosu na STAR

#### Odstupanja u koordinati

Sedam readova (u opsegu 2L) oba alata smjestila su jedinstveno, ali na
različite pozicije. Sva odstupanja identična su:

| Broj readova | Pozicija prema STAR-u (`nM`) | Pozicija prema vlastitom aligneru (nepodudaranja) | Razlika |
|---|---|---|---|
| 7 | 2L:14.743.463 (4) | 2L:14.743.493 (3) | 30 baza |

Ovih 7 je podskup od ukupno 9 pozicijskih neslaganja pronađenih u opsegu
cijelog genoma ([Provjera devet neslaganja](#provjera-devet-neslaganja));
preostala 2 neslaganja odnose se na readove čije se STAR-ovo poravnanje
nalazi izvan 2L (na 3L), pa se ne pojavljuju u usporedbi ograničenoj na
2L.

Na svih sedam pozicija vlastiti aligner prijavljuje manji broj
nepodudaranja od STAR-a. To je izmjereno i potvrđeno STAR-ovom vlastitom
oznakom `nM`, ali **ne dokazuje** da je njegova pozicija biološki ispravna.

#### Odstupanja u broju nepodudaranja

Od 1.713.201 zajednički jedinstveno smještenog reada na 2L, 3.091 ima
različit broj nepodudaranja. Uzrok je utvrđen izravnom provjerom sekvenci
tih readova u ulaznoj FASTQ datoteci:

| Objašnjenje | Broj readova | Udio |
|---|---|---|
| Razlika točno odgovara broju baza `N` u readu | 3.084 | 99,77 % |
| Readovi s odstupanjem u koordinati | 7 | 0,23 % |

STAR bazu `N` u readu ne broji kao nepodudaranje, dok je vlastita
implementacija uspoređuje znak po znak i broji kao nepodudaranje. Riječ je
o razlici u konvenciji brojanja, a ne o razlici u smještanju reada:
koordinata je u svih 3.084 slučaja identična.

Preostalih sedam readova nema nijednu bazu `N`, a razlika u broju
nepodudaranja iznosi točno −1, što je posljedica različite pozicije, a ne
konvencije.

### Mogući uzroci uočenih razlika

| Uočena razlika | Mogući uzrok | Status |
|---|---|---|
| 3.084 readova (2L) s različitim brojem nepodudaranja | različita konvencija brojanja baza `N` | utvrđeno provjerom sekvenci |
| 7 readova (2L) na različitoj poziciji | dvije kopije ponavljajuće sekvence udaljene 30 baza, na kojima oba alata nalaze prihvatljivo poravnanje | pretpostavka; podupire je jednaka razlika kod svih sedam readova i mali broj nepodudaranja na obje pozicije |
| 32.254 readova (2L) koje vlastiti aligner nije poravnao, a STAR ih je čisto smjestio na 2L | odbacivanje seedova s više od `max_seed_hits` podudaranja ili premašen prag od najviše tri nepodudaranja | izmjereno (vidi niže) |
| 10.375 readova (2L) označenih kao `multi` | ponavljajuće regije unutar kromosoma 2L | pretpostavka |
| 153.150 readova (cijeli genom) koje vlastiti aligner nije poravnao | isti mehanizam kao gore, na skali cijelog genoma | pretpostavka po analogiji, nije zasebno mjereno na cijelom genomu |
| Razlika u brzini | jednodretveni Python naspram višedretvenog C++ programa | djelomično izmjereno; udjeli paralelizma i implementacijskog jezika nisu razdvojeni |

#### Zašto readovi ostaju nepovezani: izmjereno (2L)

Za razlikovanje dva moguća uzroka — odbačen seed zbog prekomjernog broja
podudaranja naspram premašenog praga nepodudaranja pri proširenju —
napisana je zasebna skripta
([`analyze_unmapped_chr2L.py`](analyze_unmapped_chr2L.py)) koja za svaki
`unmapped` read poziva iste interne funkcije pretrage seedova kao
`align_read` (`_split_into_seeds`, `_sa_search_bounds`), bez izmjene
`genome_index/alignment.py`. Za svaki read utvrđuje se jedno od tri stanja:

- **`no_exact_seed_match`** — nijedan od 4 seeda, ni na jednom od dva
  lanca, nema nijedno egzaktno podudaranje u sekvenci; kandidatska pozicija
  se uopće ne generira.
- **`seed_rejected_frequent`** — barem jedan seed ima egzaktna
  podudaranja, ali su sva takva odbačena jer ih ima više od
  `max_seed_hits` (1000); kandidat se ni tu ne generira.
- **`extension_failed`** — barem jedan kandidat je generiran (seed unutar
  `max_seed_hits`), ali nijedan nije zadovoljio prag od najviše tri
  nepodudaranja pri usporedbi cijelog reada.

Mjerenje je provedeno u dva opsega: nad svih 12.637.449 readova koje
vlastiti aligner ne poravna na 2L, i zasebno nad podskupom od 32.254 koje
je STAR ipak čisto smjestio na 2L (pa se za njih zna da genomsko mjesto
postoji unutar indeksiranog kromosoma).

| Uzrok | Svi unmapped na 2L (12.637.449) | Podskup koji STAR čisto smjesti (32.254) |
|---|---|---|
| `no_exact_seed_match` | 12.011.684 (95,05 %) | 5.997 (18,59 %) |
| `seed_rejected_frequent` | 0 (0,00 %) | 0 (0,00 %) |
| `extension_failed` | 625.765 (4,95 %) | 26.257 (81,41 %) |

`seed_rejected_frequent` je izmjereno kao 0 u oba opsega: nijedan seed
duljine ~19 baza (76bp read podijeljen na 4 seeda za `max_mismatches = 3`)
nema više od 1000 egzaktnih podudaranja na kromosomu 2L, pa parametar
`max_seed_hits` u ovom mjerenju nije bio ograničavajući faktor.

Na razini svih 12,6 milijuna nepovezanih readova, 95,05 % nema nijedno
egzaktno 19-bazno podudaranje na 2L — očekivano, jer većina tih readova
potječe s drugih kromosoma. Unutar užeg, biološki relevantnog podskupa
koji STAR čisto smjesti na 2L, slika je obrnuta: 81,41 % otpada na
premašen prag nepodudaranja pri proširenju kandidata. To upućuje na to da
su ti readovi genomski smješteni na 2L, ali sadrže više od tri razlike
naspram reference (stvarne varijante i/ili greške sekvenciranja) — više
nego što vlastiti aligner sa zadanim `max_mismatches = 3` dopušta.

Puni rezultat po readu nalazi se u
`results/chr2L/unmapped_reason_breakdown.tsv` (gitignored zbog veličine, uz
commitani uzorak `unmapped_reason_breakdown.sample.tsv`). Izvođenje je
trajalo 1062,34 s. Ovo mjerenje nije ponovljeno za opseg cijelog genoma.

### Rezultati pronalaženja mutacija

Pozivanje mutacija izvedeno je samo u opsegu kromosoma 2L, preko
`run_pipeline.py`; whole-genome pipeline u `star_comparison/` ne poziva
mutacije.

Nad 2.321.574 jedinstveno poravnana reada, uz `min_coverage = 3` i
`min_variant_fraction = 0,5`, pozvano je **10.445 točkastih mutacija** u
264,24 s.

| Mjera | Vrijednost |
|---|---|
| Pokrivenost (najmanja / medijan / najveća) | 3 / 6 / 9582 |
| VAF (najmanji / medijan / najveći) | 0,5000 / 0,8468 / 1,0000 |
| Tranzicije | 5896 |
| Transverzije | 4549 |
| Omjer Ti/Tv | 1,296 |

Raspodjela prema jačini dokaza:

| Skupina | Broj | Udio |
|---|---|---|
| Pokrivenost 3 | 2106 | 20,2 % |
| Pokrivenost manja od 5 | 3768 | 36,1 % |
| Pokrivenost 10 ili više | 3951 | 37,8 % |
| VAF jednak 1,0 | 4412 | 42,2 % |
| Pokrivenost 3 i VAF 1,0 | 1153 | 11,0 % |

Omjer Ti/Tv veći od 1 uobičajen je za stvarnu genomsku varijaciju. To je
**slaba potvrda na razini cijelog skupa** i ne govori ništa o ispravnosti
pojedinačnog poziva.

Više od trećine poziva počiva na pokrivenosti manjoj od pet readova, a
1153 poziva ima pokrivenost 3 uz potpuno slaganje readova. Pri takvoj
pokrivenosti tri reada s istom greškom sekvenciranja zadovoljavaju
postavljene pragove, pa se za tu skupinu ne može tvrditi da predstavlja
stvarne varijante.

**Ovi rezultati nisu neovisno provjereni.** STAR ne poziva varijante i ne
može poslužiti kao njihova potvrda. Za ozbiljnu procjenu bilo bi potrebno:

- usporediti pozicije s poznatom bazom varijanti za *Drosophila
  melanogaster*,
- pozvati varijante utvrđenim alatom nad istim poravnanjima i usporediti
  skupove,
- odvojiti RNA editing od genomskih varijanti, budući da je riječ o
  RNA-seq podatcima,
- ispitati osjetljivost rezultata na pragove pokrivenosti i VAF-a.

Do tada pozvane mutacije treba smatrati izlazom algoritma, a ne skupom
potvrđenih bioloških varijanti.

### Prednosti i ograničenja

#### Prednosti

- **Točnost smještanja na provjerenom podskupu.** Na readovima koje
  vlastiti aligner jedinstveno smjesti, podudarnost sa STAR-om iznosi
  99,9996 % za kromosom 2L i 99,9999 % u usporedbi na razini cijelog genoma.
- **Determinizam.** Ponovljena izvođenja daju identične rezultate.
- **Potpunost unutar zadanog praga.** Pigeonhole jamstvo osigurava da
  postupak ne propušta poravnanja s najviše `max_mismatches` nepodudaranja,
  osim kod readova čiji su svi seedovi odbačeni zbog prevelikog broja
  podudaranja (izmjereno kao 0 slučajeva na 2L).
- **Skromna potrošnja memorije pri poravnanju.** 0,17 GB (2L) do 1,56 GB
  (cijeli genom) za obradu 15 milijuna readova, jer se readovi obrađuju
  pojedinačno.

#### Ograničenja

- **Nema podrške za insercije, delecije i splicing.** Pri 1.325.356 splice
  spojeva koje je STAR zabilježio, ovo ograničenje isključuje znatan dio
  RNA-seq podataka (36,4 % primarnih STAR poravnanja sadrži barem jednu od
  tih operacija).
- **Nema procjene pouzdanosti poravnanja.** Ne postoji ekvivalent MAPQ
  vrijednosti, pa se poravnanja ne mogu filtrirati po kvaliteti.
- **Brzina.** Pri jednakom opsegu pretrage (cijeli genom) vlastiti aligner
  bio je približno 16,8 puta sporiji od STAR-a, u konfiguraciji s jednom
  naspram šest dretvi.
- **Potrošnja memorije pri pozivanju mutacija.** Peak RSS raste s 0,17 GB
  na 2,90 GB (2L) jer se svi readovi zadržavaju u memoriji, uključujući one
  koji neće ući u pileup.
- **Neprovjeren dio izlaza.** 608.373 jedinstvenih poravnanja na 2L izvan
  STAR-ova čistog skupa nije provjereno; analogno vrijedi i za cijeli
  genom.
- **Mjerenja iz jednog izvođenja.** Vremena i memorija nisu prosjek više
  ponavljanja, u oba opsega.
- **Popravak STAR-a na macOS-u provjeren je samo na jednoj konfiguraciji**,
  netestiran na Linuxu
  ([Poznati problem](#poznati-problem-star-na-macos-u-apple-silicon)).

### Procjena praktične primjenjivosti

**Za produkcijsku analizu RNA-seq podataka implementacija nije
primjenjiva.** Odlučujući razlog nije točnost, nego nedostatak podrške za
splicing: STAR je na istom skupu zabilježio 1.325.356 splice spojeva, a
readovi koji ih sadrže ne mogu se ispravno poravnati neprekinutim
poravnanjem. Uz to nedostaje procjena pouzdanosti poravnanja, bez koje se
rezultati ne mogu filtrirati u daljnjoj analizi.

**Za analizu DNA podataka bez indela implementacija je načelno
upotrebljiva**, uz ogradu da takav scenarij nije testiran u ovom radu.
Izmjerena podudarnost odnosi se na neprekinuto poravnana mjesta, što
odgovara tom razredu problema, ali provjera nad DNA-seq podatcima nije
provedena.

**Kao provjera ispravnosti implementiranih algoritama rezultat je
uspješan, u oba opsega.** Algoritmi su izvedeni nad cjelovitim stvarnim
skupom od 14.988.127 readova, i nad jednim kromosomom i nad cijelim
genomom, a ne samo nad testnim primjerima, uz podudarnost sa STAR-om koja
odgovara očekivanjima za razred readova koji se obrađuje u oba slučaja.

Brzina je ograničenje implementacije, ne algoritma. Postupak je
jednodretven i izveden u Pythonu, uz C biblioteku samo za izgradnju suffix
arraya, dok je pretraga po readu Python petlja.

### Poznati problem: STAR na macOS-u (Apple Silicon)

**Konfiguracija na kojoj je problem opažen:** STAR 2.7.11b (Homebrew
`rna-star` build), macOS 26.5 (build 25F71), Apple Silicon (arm64), Apple
clang 15.0.0 (`arm64-apple-darwin25.5.0`), Appleova implementacija
standardne C++ biblioteke (libc++). Na ovoj konfiguraciji `alignReads`
tiho vraća 0 poravnatih reada — proces završava s exit kodom 0, izgleda
kao uspješan run. Nije provjereno javlja li se isti problem na drugim
verzijama macOS-a, arhitekturama ili verzijama STAR-a, pa se nalaz ne
generalizira izvan gore navedene konfiguracije niti na svaki Homebrew
build na svakom Apple Silicon sustavu.

**Uzrok.** `std::stringbuf::pubsetbuf()` je standardna C++ funkcija, ali
njeno konkretno ponašanje za `stringbuf` C++ standard ostavlja
implementation-defined — ne postoji jamstvo da će predani `char*` postati
stvarni interni buffer streama. STAR 2.7.11b se oslanja upravo na to da
hoće. U testiranom buildu to je funkcioniralo uz GNU-ov libstdc++
(ponašanje koje ta konkretna implementacija bira, ne nešto što standard
propisuje), dok Appleov libc++ isti poziv tretira kao no-op i nikad ne
povezuje vanjski buffer. Posljedica: STAR-ov buffer za čitanje readova i
buffer za pisanje SAM zapisa ostaju prazni bez obzira na stvarni sadržaj
streamova, a program ipak završava s exit kodom 0 umjesto greškom.

**Popravak.**
[`star_comparison/patches/star_2.7.11b_macos_libcxx_fix.patch`](star_comparison/patches/star_2.7.11b_macos_libcxx_fix.patch)
zamjenjuje `pubsetbuf`-aliasing pristup eksplicitnim
`std::stringbuf::str()` pozivom — kopira sadržaj u stringbuf umjesto da se
oslanja na aliasing vanjskog buffera (vidi i
[`star_comparison/scripts/build_star_macos.sh`](star_comparison/scripts/build_star_macos.sh)).
`stringbuf::str()` je standardno, prenosivo ponašanje, pa se na
Linuxu/libstdc++-u očekuje isti funkcionalni rezultat; to nije testirano u
sklopu ovog rada, pa se ne tvrdi identična izvedba (vrijeme, kopiranje
memorije) na toj platformi — samo očekivana funkcionalna ekvivalentnost.

### Reprodukcija pipelinea

Ulazni podatci (`Fly/`) i međurezultati veći od nekoliko stotina KB nisu u
gitu (vidi `.gitignore`) — commitane su `.sample.*` verzije (zaglavlje +
prvih ~2000 redaka), a puni fileovi regeneriraju se lokalno.

**Kromosom 2L:**

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3

python3 star_comparison/scripts/compare_chromosome.py \
  star_comparison/results/star/IFM48h_1.Aligned.out.sam \
  results/chr2L/my_alignments_full.tsv \
  2L results/chr2L/comparison_chr2L.csv --tolerance 0

python3 analyze_unmapped_chr2L.py
```

**Cijeli genom** (svaka naredba se pokreće iz roota repozitorija, redom):

```bash
star_comparison/scripts/build_star_macos.sh                     # build patchanog STAR-a -> bin/STAR
star_comparison/scripts/run_star.sh                              # genome index + SAM     -> results/star/
star_comparison/scripts/filter_clean_reads.py  <sam> <out>       # -> results/filtered/clean_reads.tsv
star_comparison/scripts/make_fasta_input.py    <tsv> <out>       # -> results/my_aligner/clean_reads_input.fasta
star_comparison/scripts/build_my_index.py                        # -> my_index/
star_comparison/scripts/run_my_aligner.py      <fasta> <out>     # -> results/my_aligner/my_alignments.tsv
star_comparison/scripts/compare_coordinates.py <star> <my> <out> --mismatches-out <out2>
star_comparison/scripts/inspect_mismatches.py  <out2> <clean> <out3>
star_comparison/scripts/extract_nh.py          <sam> <out>       # -> results/filtered/clean_reads_nh.tsv
```

Kratak opis svakog koraka i točan oblik naredbe s pravim putanjama je u
[`star_comparison/README.md`](star_comparison/README.md), sekcija "Koraci
pipelinea". Oba pipelinea su deterministička — ponovno pokretanje na istim
ulaznim podatcima daje identične brojke.

### Zaključak

Vlastita implementacija izgradnje suffix arraya i poravnanja readova
izvedena je nad cjelovitim stvarnim skupom podataka, u dva neovisna opsega
— jedan kromosom i cijeli genom — i uspoređena sa STAR-om kao referentnim
rješenjem.

U opsegu kromosoma 2L: od 1.755.830 readova koje je STAR čisto poravnao,
vlastiti aligner smjestio je 1.713.194 na identičnu koordinatu, stroga
točnost 97,5717 %, 99,9996 % među jedinstveno smještenima. U opsegu
cijelog genoma: od 9.299.988 čistih readova, stroga točnost 96,42 %,
99,9999 % među jedinstveno smještenima. Sedam pozicijskih neslaganja iz
opsega 2L identično je sedmorici od ukupno devet neslaganja u opsegu
cijelog genoma — dvije neovisno izvedene usporedbe daju isti nalaz. Odstupanje
u broju nepodudaranja u 99,77 % slučajeva posljedica je razlike u
konvenciji brojanja baza `N`, ne pogreške u smještanju.

Iz toga slijedi zaključak da algoritam ispravno smješta readove razreda za
koji je namijenjen, odnosno readove bez praznina u poravnanju, u oba
testirana opsega. Zaključak vrijedi za provjereni podskup u svakom opsegu;
readovi izvan STAR-ova čistog skupa (608.373 na 2L, analogno na cijelom
genomu) ostaju neprovjereni i za njih se tvrdnja o točnosti ne može
proširiti bez dodatnih mjerenja.

Algoritam pronalaženja mutacija, pokrenut u opsegu jednog kromosoma,
vratio je 10.445 poziva. Ti rezultati nisu neovisno provjereni i ne mogu se
predstaviti kao potvrđene biološke varijante, tim više što više od trećine
poziva počiva na pokrivenosti manjoj od pet readova.

Za potpuniju ocjenu bilo bi potrebno provesti dodatne testove: usporedbu
nad simuliranim readovima s poznatim podrijetlom, čime bi se dobila
apsolutna, a ne relativna mjera točnosti; provjeru dijela izlaza koji nije
pokriven STAR-ovim čistim skupom, u oba opsega; te neovisnu provjeru
pozvanih mutacija usporedbom s poznatom bazom varijanti.
