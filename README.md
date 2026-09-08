# 🧬 Genlyzee

Own implementation of three genomic data-processing algorithms: chromosome
**suffix array** construction, **read alignment** against that index, and
**point-mutation calling** from the resulting alignments. Written in
Python, without using any existing alignment tool.
[STAR](https://github.com/alexdobin/STAR) is used only as an external
reference solution to validate the results.

**[English](#english)** · **[Hrvatski](#hrvatski)**

---

## English

### Table of contents

1. [Purpose and scope](#purpose-and-scope)
2. [Architecture](#architecture)
3. [Prerequisites and installation](#prerequisites-and-installation)
4. [Repository layout](#repository-layout)
5. [Algorithm descriptions](#algorithm-descriptions)
   - [Suffix array construction](#suffix-array-construction)
   - [Read alignment](#read-alignment)
   - [Mutation calling](#mutation-calling)
6. [Programmatic interface](#programmatic-interface)
7. [Running it](#running-it)
8. [Input data](#input-data)
9. [Output data](#output-data)
10. [Example run](#example-run)
11. [Testing](#testing)
12. [STAR's role](#stars-role)
13. [Limitations](#limitations)

### Purpose and scope

The program takes a reference genome in FASTA format, a chromosome name,
and a reads file, and runs the full pipeline: it extracts the requested
chromosome, builds its suffix array, aligns every read against that index,
and calls point mutations from the resulting alignments.

The goal of the project is an own implementation of these algorithms and
verification of their correctness on real biological data. STAR is not
part of the implementation; it is used only as an independent reference
solution against which alignment coordinates are compared. Comparison
results are in [`REPORT.md`](REPORT.md).

### Architecture

All three stages are implemented in the `genome_index/` package. The
`run_pipeline.py` script wires them together into a single command-line
program and contains no algorithmic logic of its own.

```text
FASTA + chromosome name                    reads file
          │                                          │
          ▼                                          │
  parse_fasta_chromosome                             │
   (genome_index/fasta_io.py)                        │
          │  chromosome sequence                     │
          ▼                                          ▼
  build_suffix_array_from_sequence  ────────►  align_read
   (genome_index/suffix_array.py)        (genome_index/alignment.py)
          │                                          │
          │                                          │ alignments
          │                                          ▼
          └──────────────────────────────►  call_mutations
                                     (genome_index/mutation_calling.py)
                                                     │
                                                     ▼
                                    alignments TSV and mutations TSV
```

| Stage | Module | Input | Output |
|---|---|---|---|
| FASTA parsing | `fasta_io.py` | FASTA path, chromosome name | chromosome sequence |
| Suffix array construction | `suffix_array.py` | sequence | `numpy` array of suffix positions |
| Alignment | `alignment.py` | read, sequence, suffix array | list of alignments |
| Mutation calling | `mutation_calling.py` | sequence, reads, alignments | list of called mutations |

### Prerequisites and installation

| Component | Version used in this work | Role |
|---|---|---|
| Python | 3.13.5 | runtime |
| `numpy` | 2.5.1 | suffix array storage |
| `pydivsufsort` | 0.0.20 | `divsufsort` backend for suffix array construction |
| `pytest` | 9.1.1 | running the tests |

`requirements.txt` does not pin versions; the ones listed above are the
ones the results in [`REPORT.md`](REPORT.md) were obtained with.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The STAR comparison needs additional tools, described in
[`star_comparison/README.md`](star_comparison/README.md).

### Repository layout

```text
.
├── genome_index/                  own implementation
│   ├── fasta_io.py                extract a chromosome from a multi-FASTA file
│   ├── suffix_array.py            suffix array construction (naive, SA-IS, divsufsort)
│   ├── alignment.py               read alignment
│   └── mutation_calling.py        point-mutation calling
├── tests/                         unit and fuzz tests
│   ├── test_suffix_array.py
│   ├── test_alignment.py
│   ├── test_alignment_fuzz.py
│   └── test_mutation_calling.py
├── run_pipeline.py                CLI wiring all three stages together
├── moj_primjer.py                 end-to-end demo with a known ground truth
├── analyze_unmapped_chr2L.py      analysis of why reads stay unmapped, over run_pipeline.py's output
├── requirements.txt
├── results/chr2L/                 results of the run over chromosome 2L
├── REPORT.md                      STAR comparison report
└── star_comparison/               STAR comparison pipeline
    ├── scripts/                   pipeline scripts
    ├── patches/                   macOS source fix for STAR
    ├── results/                   comparison results
    └── README.md                  building and running STAR
```

The `Fly/` input-data directory, and the directories produced by building
STAR and its index, are not included in the repository due to size.

### Algorithm descriptions

#### Suffix array construction

A suffix array is an array of the starting positions of every suffix of a
sequence, sorted lexicographically. Suffixes sharing a prefix therefore
occupy one contiguous block of the array, which makes it possible to
search for a pattern with binary search instead of a linear scan of the
sequence.

Three procedures are implemented:

| Procedure | Function | Complexity | Purpose |
|---|---|---|---|
| Naive suffix sort | `build_suffix_array_naive` | O(n² log n) | reference solution for correctness checks in the tests |
| SA-IS | `_sa_is`, via `build_suffix_array_from_sequence(method="sais")` | O(n) | own implementation of the linear algorithm |
| `divsufsort` | `build_suffix_array_from_sequence(method="divsufsort")` | — | used for runs at real chromosome scale |

**SA-IS** (Nong, Zhang and Chen, 2009) builds the suffix array via induced
sorting. Procedure:

1. Every position is classified as S-type or L-type in a single
   right-to-left pass, since the type of position `i` depends only on the
   type of position `i+1`.
2. LMS positions are identified (an S-type position whose predecessor is
   L-type); consecutive LMS positions bound the LMS substrings.
3. LMS suffixes are placed at the end of their buckets, after which one
   L-pass and one S-pass of induced sorting correctly sort the LMS
   substrings relative to each other.
4. Every LMS substring is assigned a name according to that order.
5. If all names are distinct, the order of the LMS suffixes is known.
   Otherwise, the suffix array of the reduced name string — at most half
   the length of the input — is built recursively.
6. The correctly ordered LMS suffixes are placed at the end of the
   buckets again, and both induced-sorting passes are run once more,
   producing the final suffix array.

SA-IS was chosen over the DC3/skew algorithm because it needs only one
recursive step, instead of DC3's split-then-merge approach.

The pure-Python implementation is linear but too slow in practice for a
chromosome of 23.5 million bases, so real-data runs use `divsufsort` (the
`libdivsufsort` library via `pydivsufsort`). Measured construction of the
suffix array for chromosome 2L (23,513,712 bases) with that backend takes
1.19 s.

The result is a `numpy` array of type `int32`, or `int64` for sequences
longer than 2³¹ positions.

#### Read alignment

Alignment happens in two steps.

**Exact matching.** Binary search over the suffix array determines the
lower and upper bound of the block of suffixes that start with the given
pattern (`_sa_lower_bound`, `_sa_upper_bound`). Complexity is O(m log n),
where `m` is the pattern length and `n` the sequence length, plus O(k) to
read out `k` matches found.

**Mismatch-tolerant alignment.** A seed-and-extend approach based on the
pigeonhole principle is used:

1. The read is split into `max_mismatches + 1` seeds. If the read has at
   most `max_mismatches` mismatches, at least one seed is necessarily free
   of any difference, so it can be located by exact search. The procedure
   therefore never misses an alignment within the given threshold.
2. Exact matches are searched for every seed. The positions are shifted
   back by the seed's offset within the read, giving candidate start
   positions for the whole read.
3. Every candidate is checked by comparing the whole read against the
   corresponding reference window, base by base. Candidates whose mismatch
   count is within the threshold are returned as alignments.

A seed with more than `max_seed_hits` exact matches is skipped as
uninformative. Such seeds occur in repetitive regions and runs of `N`
bases, where checking every candidate would be infeasible, and
localization of the read is left to the remaining seeds.

The procedure runs on both strands: the read as given (`strand = "+"`) and
its reverse complement (`strand = "-"`). The position is reported on the
forward reference strand in both cases.

A read is classified by the number of alignments found: exactly one
(`unique`), more than one (`multi`), or none (`unmapped`).

#### Mutation calling

Only point substitutions are called. Insertions and deletions are not
supported, consistent with the alignment algorithm also not handling
them.

1. **Pileup construction.** For every read with exactly one alignment,
   its bases are mapped onto forward-strand reference coordinates. For an
   alignment on the negative strand, the reverse complement of the read is
   taken first. `N` bases in the read are skipped as uninformative. Reads
   with no alignment and multi-mapped reads are left out of the pileup,
   trading recall for precision.
2. **Position check.** For every position with coverage of at least
   `min_coverage`, the reference base is compared against the bases
   reported at that position by the aligned reads. Positions where the
   reference base is `N` are skipped.
3. **Mutation call.** Among the bases that differ from the reference, the
   most common one is taken as the alternate. If the fraction of reads
   supporting it (variant allele fraction, VAF) is at least
   `min_variant_fraction`, a mutation is called.

Both conditions are necessary: the fraction alone would call every
single-read sequencing error, whose VAF is trivially 1.0, while coverage
alone does not rule out several different errors accumulating by chance.

### Programmatic interface

| Function | Module | Signature |
|---|---|---|
| `parse_fasta_chromosome` | `fasta_io` | `(fasta_path, chromosome_name) -> str` |
| `build_suffix_array` | `suffix_array` | `(fasta_path, chromosome_name, method="sais") -> np.ndarray` |
| `build_suffix_array_from_sequence` | `suffix_array` | `(sequence, method="sais") -> np.ndarray` |
| `build_suffix_array_naive` | `suffix_array` | `(text) -> List[int]` |
| `find_exact_matches` | `alignment` | `(sequence, suffix_array, pattern) -> List[int]` |
| `align_read` | `alignment` | `(read, sequence, suffix_array, max_mismatches=2, max_seed_hits=1000) -> List[AlignmentResult]` |
| `align_reads` | `alignment` | `(reads, sequence, suffix_array, ...) -> Dict[str, List[AlignmentResult]]` |
| `align_reads_from_file` | `alignment` | `(path, sequence, suffix_array, ...) -> Dict[str, List[AlignmentResult]]` |
| `load_reads` | `alignment` | `(path) -> List[Tuple[str, str]]` |
| `reverse_complement` | `alignment` | `(sequence) -> str` |
| `call_mutations` | `mutation_calling` | `(reference_sequence, reads, alignment_results, min_coverage=3, min_variant_fraction=0.5) -> List[MutationCall]` |

`build_suffix_array` is the function that matches the assignment's
requirement: it takes a FASTA file and a chromosome name, and returns the
suffix array of that chromosome.

**`AlignmentResult`**

| Field | Description |
|---|---|
| `position` | 0-indexed offset on the forward reference strand |
| `mismatches` | number of mismatched bases between the read and the reference window |
| `strand` | `"+"` or `"-"` |

**`MutationCall`**

| Field | Description |
|---|---|
| `position` | 0-indexed offset, same coordinate space as `AlignmentResult` |
| `ref_base` | reference base |
| `alt_base` | alternate base |
| `supporting_reads` | number of reads supporting the alternate base |
| `coverage` | total number of reads at the position |
| `variant_allele_fraction` | `supporting_reads / coverage` |

Both structures are immutable (`frozen` dataclasses).

### Running it

```bash
python3 run_pipeline.py <fasta_path> <chromosome> <reads_path> \
  --alignments-out <path> [options]
```

**Positional arguments**

| Argument | Description |
|---|---|
| `fasta_path` | reference genome, multi-FASTA format |
| `chromosome` | name of the chromosome to extract from the FASTA file |
| `reads_path` | reads file, FASTQ or FASTA format |

**Options**

| Option | Default | Description |
|---|---|---|
| `--alignments-out` | required | output TSV with alignments |
| `--mutations-out` | not set | output TSV with mutations; without it, mutations are not called |
| `--max-mismatches` | 3 | maximum mismatches allowed per read |
| `--max-seed-hits` | 1000 | a seed with more matches than this is skipped as uninformative |
| `--min-coverage` | 3 | minimum coverage at a position required to call a mutation |
| `--min-variant-fraction` | 0.5 | minimum fraction of reads supporting the alternate base |
| `--sa-method` | `divsufsort` | suffix array construction backend (`divsufsort` or `sais`) |

The program prints its progress and per-stage measurements to standard
error.

### Input data

| File | Format | Note |
|---|---|---|
| reference genome | multi-FASTA | the chromosome name is matched against the first token of the header, case-insensitively |
| reads | FASTQ or FASTA | the format is detected from the first character of the first non-blank line |

Reads are read one at a time, so memory use does not grow with the input
file's size. The exception is mutation calling, which needs all reads at
once.

A read is skipped (`skipped`) if it is longer than the chromosome or
contains characters outside the `ACGTN` set.

### Output data

**Alignments file** (`--alignments-out`), one row per read:

| Column | Description |
|---|---|
| `read_id` | first token of the read header |
| `status` | `unique`, `multi`, `unmapped`, or `skipped` |
| `pos_1based` | alignment start coordinate; filled in only for `unique` |
| `strand` | `+` or `-`; filled in only for `unique` |
| `mismatches` | mismatch count; filled in only for `unique` |
| `num_hits` | number of alignments found |

**Mutations file** (`--mutations-out`), one row per called mutation:

| Column | Description |
|---|---|
| `chromosome` | chromosome name given at invocation |
| `pos_1based` | mutation position |
| `ref_base` | reference base |
| `alt_base` | alternate base |
| `supporting_reads` | number of reads supporting the alternate base |
| `coverage` | total number of reads at the position |
| `vaf` | `supporting_reads / coverage` |

The suffix array and the `AlignmentResult`/`MutationCall` structures use
0-indexed positions. Both output files convert them to 1-indexed, to
match the `POS` field in SAM format.

### Example run

The command that produced the results described in [`REPORT.md`](REPORT.md):

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3
```

Input files:

| File | Size | Content |
|---|---|---|
| `Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa` | 146,261,227 B | BDGP6 reference genome, 1870 records |
| `Fly/IFM48h_1.fastq` | 3,663,357,012 B | 14,988,127 single-end RNA-seq reads, 76 bases long |

Files produced in `results/chr2L/`:

| File | Content | In the repository |
|---|---|---|
| `my_alignments_full.tsv` | alignment of every read | no, with sample `my_alignments.sample.tsv` |
| `my_mutations.tsv` | 10,445 called mutations | yes |
| `comparison_chr2L.csv` | per-read comparison against STAR | no, with sample `comparison_chr2L.sample.csv` |
| `run_full_time_mem.log` | timing without mutation calling | yes |
| `run_full_with_mutations_time_mem.log` | timing for the full run | yes |

### Testing

```bash
pytest tests/ -v
```

16 tests, all passing.

| File | Test count | Content |
|---|---|---|
| `test_suffix_array.py` | 5 | known examples, SA-IS vs. the naive procedure on random strings, empty input, `divsufsort`/`sais` backend agreement |
| `test_alignment.py` | 6 | exact matching, alignment with mismatches, a nonexistent read, reverse complement, batch alignment with deduplication, input validation |
| `test_alignment_fuzz.py` | 1 | comparison against an exhaustive procedure on randomly generated cases |
| `test_mutation_calling.py` | 4 | a call with sufficient coverage and fraction, rejecting a single error, rejecting low coverage, reads on the negative strand |

Additionally, `moj_primjer.py` runs the full pipeline over chromosome chrM
from the hg38 assembly with a known ground truth: it aligns a read with an
inserted difference back to its known position, and, over a set of 200
reads with two planted mutations and five isolated errors, checks that
both mutations are called and neither error is. The script requires a
local `hg38/hg38.u.fa` file, which is not included in the repository.

```bash
python3 moj_primjer.py
```

### STAR's role

STAR is an external, independently developed aligner used in this project
solely as a reference solution. It has not been modified beyond a fix
needed to build it on macOS, and it plays no part in computing the own
implementation's results.

The comparison comes down to checking whether the own aligner places
reads on the same coordinates as STAR. Since the reads are not simulated,
the true biological origin of any given read is not known, so results are
expressed as concordance with STAR, not as absolute accuracy.

| Document | Content |
|---|---|
| [`REPORT.md`](REPORT.md) | full comparison report, both scopes: chromosome 2L and the whole genome |
| [`star_comparison/README.md`](star_comparison/README.md) | building and running STAR, comparison procedure |

### Limitations

| Limitation | Consequence |
|---|---|
| No support for insertions, deletions, or splicing | reads spanning exon boundaries cannot be aligned; comparison against STAR is only possible on reads without these operations |
| The index covers a single chromosome | reads from other chromosomes stay unmapped or may be placed incorrectly within the indexed chromosome |
| Alignment is single-threaded and written in Python | running it is substantially slower than STAR |
| Mutation calling keeps all reads in memory | memory use grows with the number of reads |
| No alignment-confidence score | there is no MAPQ equivalent, so alignments cannot be filtered by quality |
| Called mutations are not independently validated | the result is algorithm output, not a confirmed set of biological variants |

Measurements and a detailed analysis of these limitations are in
[`REPORT.md`](REPORT.md#strengths-and-limitations).

---

## Hrvatski

### Sadržaj

1. [Svrha i opseg](#svrha-i-opseg)
2. [Arhitektura](#arhitektura)
3. [Preduvjeti i instalacija](#preduvjeti-i-instalacija)
4. [Struktura repozitorija](#struktura-repozitorija)
5. [Opis algoritama](#opis-algoritama)
   - [Izgradnja suffix arraya](#izgradnja-suffix-arraya)
   - [Poravnanje readova](#poravnanje-readova)
   - [Pronalaženje mutacija](#pronalaženje-mutacija)
6. [Programsko sučelje](#programsko-sučelje)
7. [Pokretanje](#pokretanje)
8. [Ulazni podatci](#ulazni-podatci)
9. [Izlazni podatci](#izlazni-podatci)
10. [Primjer izvođenja](#primjer-izvođenja)
11. [Testiranje](#testiranje)
12. [Uloga STAR-a](#uloga-star-a)
13. [Ograničenja](#ograničenja)

### Svrha i opseg

Program prima referentni genom u FASTA formatu, naziv kromosoma i datoteku
s readovima, te izvodi cjelovit postupak: izdvaja traženi kromosom, gradi
njegov suffix array, poravnava svaki read nad tim indeksom i iz poravnanja
poziva točkaste mutacije.

Cilj projekta je vlastita implementacija navedenih algoritama i provjera
njihove ispravnosti na stvarnim bioloških podatcima. STAR nije dio
implementacije; koristi se samo kao neovisno referentno rješenje s kojim se
uspoređuju koordinate poravnanja. Rezultati usporedbe nalaze se u
[`REPORT.md`](REPORT.md).

### Arhitektura

Sve tri faze implementirane su u paketu `genome_index/`. Skripta
`run_pipeline.py` povezuje ih u jedan naredbeni program i ne sadrži vlastitu
algoritamsku logiku.

```text
FASTA + naziv kromosoma                    datoteka s readovima
          │                                          │
          ▼                                          │
  parse_fasta_chromosome                             │
   (genome_index/fasta_io.py)                        │
          │  sekvenca kromosoma                      │
          ▼                                          ▼
  build_suffix_array_from_sequence  ────────►  align_read
   (genome_index/suffix_array.py)        (genome_index/alignment.py)
          │                                          │
          │                                          │ poravnanja
          │                                          ▼
          └──────────────────────────────►  call_mutations
                                     (genome_index/mutation_calling.py)
                                                     │
                                                     ▼
                                    TSV s poravnanjima i TSV s mutacijama
```

| Faza | Modul | Ulaz | Izlaz |
|---|---|---|---|
| Parsiranje FASTA datoteke | `fasta_io.py` | putanja do FASTA datoteke, naziv kromosoma | sekvenca kromosoma |
| Izgradnja suffix arraya | `suffix_array.py` | sekvenca | `numpy` polje pozicija sufiksa |
| Poravnanje | `alignment.py` | read, sekvenca, suffix array | lista poravnanja |
| Pronalaženje mutacija | `mutation_calling.py` | sekvenca, readovi, poravnanja | lista pozvanih mutacija |

### Preduvjeti i instalacija

| Komponenta | Verzija korištena u radu | Uloga |
|---|---|---|
| Python | 3.13.5 | izvođenje |
| `numpy` | 2.5.1 | pohrana suffix arraya |
| `pydivsufsort` | 0.0.20 | backend `divsufsort` za izgradnju suffix arraya |
| `pytest` | 9.1.1 | pokretanje testova |

`requirements.txt` ne fiksira verzije; navedene su one na kojima su
dobiveni rezultati opisani u [`REPORT.md`](REPORT.md).

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Za usporedbu sa STAR-om potrebni su dodatni alati opisani u
[`star_comparison/README.md`](star_comparison/README.md).

### Struktura repozitorija

```text
.
├── genome_index/                  vlastita implementacija
│   ├── fasta_io.py                izdvajanje kromosoma iz multi-FASTA datoteke
│   ├── suffix_array.py            izgradnja suffix arraya (naivno, SA-IS, divsufsort)
│   ├── alignment.py               poravnanje readova
│   └── mutation_calling.py        pozivanje točkastih mutacija
├── tests/                         jedinični i fuzz testovi
│   ├── test_suffix_array.py
│   ├── test_alignment.py
│   ├── test_alignment_fuzz.py
│   └── test_mutation_calling.py
├── run_pipeline.py                naredbeni program koji povezuje sve tri faze
├── moj_primjer.py                 primjer cjelovitog postupka s poznatim rješenjem
├── analyze_unmapped_chr2L.py      analiza uzroka nepovezanih readova nad izlazom run_pipeline.py
├── requirements.txt
├── results/chr2L/                 rezultati izvođenja nad kromosomom 2L
├── REPORT.md                      izvještaj o usporedbi sa STAR-om
└── star_comparison/               pipeline za usporedbu sa STAR-om
    ├── scripts/                   skripte pipelinea
    ├── patches/                   popravak izvornog koda STAR-a za macOS
    ├── results/                   rezultati usporedbe
    └── README.md                  izgradnja i pokretanje STAR-a
```

Direktorij `Fly/` s ulaznim podatcima te direktoriji koji nastaju
izgradnjom STAR-a i indeksa nisu uključeni u repozitorij zbog veličine.

### Opis algoritama

#### Izgradnja suffix arraya

Suffix array je polje početnih pozicija svih sufiksa sekvence, poredanih
leksikografski. Sufiksi koji dijele isti prefiks time zauzimaju jedan
neprekinuti blok polja, što omogućuje pretraživanje uzorka binarnom
pretragom umjesto linearnim prolaskom kroz sekvencu.

Implementirana su tri postupka:

| Postupak | Funkcija | Složenost | Namjena |
|---|---|---|---|
| Naivno sortiranje sufiksa | `build_suffix_array_naive` | O(n² log n) | referentno rješenje za provjeru ispravnosti u testovima |
| SA-IS | `_sa_is`, preko `build_suffix_array_from_sequence(method="sais")` | O(n) | vlastita implementacija linearnog algoritma |
| `divsufsort` | `build_suffix_array_from_sequence(method="divsufsort")` | — | izvođenje na stvarnoj kromosomskoj skali |

**SA-IS** (Nong, Zhang i Chen, 2009) gradi suffix array induciranim
sortiranjem. Postupak:

1. Svaka pozicija klasificira se kao S-tip ili L-tip u jednom prolasku
   zdesna nalijevo, jer tip pozicije `i` ovisi samo o tipu pozicije `i+1`.
2. Određuju se LMS pozicije (S-tip pozicija čiji je prethodnik L-tip) koje
   omeđuju LMS podnizove.
3. LMS sufiksi se rasporede na kraj svojih kanti, nakon čega jedan L-prolaz
   i jedan S-prolaz induciranog sortiranja ispravno poredaju LMS podnizove
   međusobno.
4. Svakom LMS podnizu dodijeli se ime prema tom poretku.
5. Ako su sva imena različita, poredak LMS sufiksa je poznat. Inače se
   rekurzivno gradi suffix array reduciranog niza imena, koji je najviše
   upola kraći od ulaznog.
6. Ispravno poredani LMS sufiksi ponovno se rasporede na kraj kanti i
   izvedu se oba prolaza induciranog sortiranja, čime se dobiva konačni
   suffix array.

SA-IS je odabran umjesto DC3/skew algoritma jer zahtijeva samo jedan
rekurzivni korak umjesto pristupa s razdvajanjem i spajanjem.

Implementacija u čistom Pythonu je linearna, ali u praksi prespora za
kromosom veličine 23,5 milijuna baza, pa se za izvođenje nad stvarnim
podatcima koristi `divsufsort` (biblioteka `libdivsufsort` preko
`pydivsufsort`). Izmjerena izgradnja suffix arraya kromosoma 2L
(23.513.712 baza) tim backendom traje 1,19 s.

Rezultat je `numpy` polje tipa `int32`, odnosno `int64` za sekvence dulje
od 2³¹ pozicija.

#### Poravnanje readova

Poravnanje se izvodi u dva koraka.

**Egzaktno podudaranje.** Binarnom pretragom nad suffix arrayem određuju se
donja i gornja granica bloka sufiksa koji počinju zadanim uzorkom
(`_sa_lower_bound`, `_sa_upper_bound`). Složenost je O(m log n), gdje je `m`
duljina uzorka, a `n` duljina sekvence, uz dodatnih O(k) za očitavanje `k`
pronađenih pozicija.

**Poravnanje s dopuštenim nepodudaranjima.** Koristi se pristup
seed-and-extend zasnovan na pigeonhole principu:

1. Read se dijeli na `max_mismatches + 1` seedova. Ako read ima najviše
   `max_mismatches` nepodudaranja, barem jedan seed nužno je bez ijedne
   razlike, pa ga se može pronaći egzaktnom pretragom. Postupak time ne
   propušta poravnanja unutar zadanog praga.
2. Za svaki seed traže se egzaktna podudaranja. Pozicije se pomiču unatrag
   za pomak seeda unutar reada, čime se dobivaju kandidatske početne
   pozicije cijelog reada.
3. Svaki kandidat provjerava se usporedbom cijelog reada s odgovarajućim
   prozorom reference, baza po baza. Kandidati kod kojih je broj
   nepodudaranja unutar praga vraćaju se kao poravnanja.

Seed s više od `max_seed_hits` egzaktnih podudaranja preskače se kao
neinformativan. Takvi seedovi javljaju se u ponavljajućim regijama i
nizovima baza `N`, gdje bi provjera svih kandidata bila neizvediva, a
lokalizaciju reada preuzimaju preostali seedovi.

Postupak se izvodi nad oba lanca: nad readom kakav je zadan (`strand = "+"`)
i nad njegovim reverznim komplementom (`strand = "-"`). Pozicija se u oba
slučaja prijavljuje na forward lancu reference.

Read se klasificira prema broju pronađenih poravnanja: točno jedno
(`unique`), više njih (`multi`) ili nijedno (`unmapped`).

#### Pronalaženje mutacija

Pozivaju se isključivo točkaste supstitucije. Insercije i delecije nisu
podržane, u skladu s time da ih ni algoritam poravnanja ne obrađuje.

1. **Izgradnja pileupa.** Za svaki read s točno jednim poravnanjem baze se
   preslikavaju na koordinate forward lanca reference. Kod poravnanja na
   negativnom lancu prvo se uzima reverzni komplement reada. Baze `N` u
   readu se preskaču jer nisu informativne. Readovi bez poravnanja i
   multi-mapirani readovi izostavljaju se iz pileupa, čime se preciznost
   ostvaruje nauštrb odziva.
2. **Provjera pozicija.** Za svaku poziciju s pokrivenošću najmanje
   `min_coverage` uspoređuje se referentna baza s bazama koje na toj
   poziciji prijavljuju poravnati readovi. Pozicije na kojima je referentna
   baza `N` preskaču se.
3. **Poziv mutacije.** Među bazama koje se razlikuju od referentne uzima se
   najčešća kao alternativna. Ako je udio readova koji je podupiru
   (variant allele fraction, VAF) najmanje `min_variant_fraction`, poziva
   se mutacija.

Oba uvjeta su nužna: sam udio bi pozvao svaku pojedinačnu grešku
sekvenciranja, čiji je VAF trivijalno 1,0, dok sama pokrivenost ne
isključuje slučajno nakupljanje različitih grešaka.

### Programsko sučelje

| Funkcija | Modul | Potpis |
|---|---|---|
| `parse_fasta_chromosome` | `fasta_io` | `(fasta_path, chromosome_name) -> str` |
| `build_suffix_array` | `suffix_array` | `(fasta_path, chromosome_name, method="sais") -> np.ndarray` |
| `build_suffix_array_from_sequence` | `suffix_array` | `(sequence, method="sais") -> np.ndarray` |
| `build_suffix_array_naive` | `suffix_array` | `(text) -> List[int]` |
| `find_exact_matches` | `alignment` | `(sequence, suffix_array, pattern) -> List[int]` |
| `align_read` | `alignment` | `(read, sequence, suffix_array, max_mismatches=2, max_seed_hits=1000) -> List[AlignmentResult]` |
| `align_reads` | `alignment` | `(reads, sequence, suffix_array, ...) -> Dict[str, List[AlignmentResult]]` |
| `align_reads_from_file` | `alignment` | `(path, sequence, suffix_array, ...) -> Dict[str, List[AlignmentResult]]` |
| `load_reads` | `alignment` | `(path) -> List[Tuple[str, str]]` |
| `reverse_complement` | `alignment` | `(sequence) -> str` |
| `call_mutations` | `mutation_calling` | `(reference_sequence, reads, alignment_results, min_coverage=3, min_variant_fraction=0.5) -> List[MutationCall]` |

`build_suffix_array` je funkcija koja odgovara zahtjevu zadatka: prima
FASTA datoteku i naziv kromosoma, a vraća suffix array tog kromosoma.

**`AlignmentResult`**

| Polje | Opis |
|---|---|
| `position` | 0-indeksirani pomak na forward lancu reference |
| `mismatches` | broj nepodudarnih baza između reada i prozora reference |
| `strand` | `"+"` ili `"-"` |

**`MutationCall`**

| Polje | Opis |
|---|---|
| `position` | 0-indeksirani pomak, isti koordinatni prostor kao `AlignmentResult` |
| `ref_base` | baza u referenci |
| `alt_base` | alternativna baza |
| `supporting_reads` | broj readova koji podupiru alternativnu bazu |
| `coverage` | ukupan broj readova na poziciji |
| `variant_allele_fraction` | `supporting_reads / coverage` |

Obje strukture su nepromjenjive (`frozen` dataclass).

### Pokretanje

```bash
python3 run_pipeline.py <fasta_path> <chromosome> <reads_path> \
  --alignments-out <putanja> [opcije]
```

**Pozicijski argumenti**

| Argument | Opis |
|---|---|
| `fasta_path` | referentni genom u multi-FASTA formatu |
| `chromosome` | naziv kromosoma koji se izdvaja iz FASTA datoteke |
| `reads_path` | datoteka s readovima u FASTQ ili FASTA formatu |

**Opcije**

| Opcija | Zadana vrijednost | Opis |
|---|---|---|
| `--alignments-out` | obavezna | izlazna TSV datoteka s poravnanjima |
| `--mutations-out` | nije zadana | izlazna TSV datoteka s mutacijama; bez nje se mutacije ne pozivaju |
| `--max-mismatches` | 3 | najveći dopušteni broj nepodudaranja po readu |
| `--max-seed-hits` | 1000 | seed s više podudaranja od zadanog preskače se kao neinformativan |
| `--min-coverage` | 3 | najmanja pokrivenost pozicije potrebna za poziv mutacije |
| `--min-variant-fraction` | 0.5 | najmanji udio readova koji podupiru alternativnu bazu |
| `--sa-method` | `divsufsort` | backend za izgradnju suffix arraya (`divsufsort` ili `sais`) |

Program tijek izvođenja i mjerenja po fazama ispisuje na standardni izlaz
za pogreške.

### Ulazni podatci

| Datoteka | Format | Napomena |
|---|---|---|
| referentni genom | multi-FASTA | naziv kromosoma uspoređuje se s prvim tokenom zaglavlja, neovisno o veličini slova |
| readovi | FASTQ ili FASTA | format se prepoznaje iz prvog znaka prve neprazne linije |

Readovi se čitaju jedan po jedan, pa potrošnja memorije ne raste s
veličinom ulazne datoteke. Iznimka je pozivanje mutacija, koje zahtijeva
sve readove istovremeno.

Read se preskače (`skipped`) ako je dulji od kromosoma ili sadrži znakove
izvan skupa `ACGTN`.

### Izlazni podatci

**Datoteka s poravnanjima** (`--alignments-out`), jedan redak po readu:

| Stupac | Opis |
|---|---|
| `read_id` | prvi token zaglavlja reada |
| `status` | `unique`, `multi`, `unmapped` ili `skipped` |
| `pos_1based` | početna koordinata poravnanja; popunjeno samo za `unique` |
| `strand` | `+` ili `-`; popunjeno samo za `unique` |
| `mismatches` | broj nepodudaranja; popunjeno samo za `unique` |
| `num_hits` | broj pronađenih poravnanja |

**Datoteka s mutacijama** (`--mutations-out`), jedan redak po pozvanoj
mutaciji:

| Stupac | Opis |
|---|---|
| `chromosome` | naziv kromosoma zadan pri pokretanju |
| `pos_1based` | pozicija mutacije |
| `ref_base` | baza u referenci |
| `alt_base` | alternativna baza |
| `supporting_reads` | broj readova koji podupiru alternativnu bazu |
| `coverage` | ukupan broj readova na poziciji |
| `vaf` | `supporting_reads / coverage` |

Suffix array i strukture `AlignmentResult` i `MutationCall` koriste
0-indeksirane pozicije. Obje izlazne datoteke pretvaraju ih u
1-indeksirane, radi usklađenosti s poljem `POS` u SAM formatu.

### Primjer izvođenja

Naredba kojom su dobiveni rezultati opisani u [`REPORT.md`](REPORT.md):

```bash
python3 run_pipeline.py \
  Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa 2L Fly/IFM48h_1.fastq \
  --alignments-out results/chr2L/my_alignments_full.tsv \
  --mutations-out results/chr2L/my_mutations.tsv \
  --max-mismatches 3
```

Ulazne datoteke:

| Datoteka | Veličina | Sadržaj |
|---|---|---|
| `Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa` | 146.261.227 B | referentni genom BDGP6, 1870 zapisa |
| `Fly/IFM48h_1.fastq` | 3.663.357.012 B | 14.988.127 single-end RNA-seq readova duljine 76 baza |

Nastale datoteke u `results/chr2L/`:

| Datoteka | Sadržaj | U repozitoriju |
|---|---|---|
| `my_alignments_full.tsv` | poravnanje svih readova | ne, uz uzorak `my_alignments.sample.tsv` |
| `my_mutations.tsv` | 10.445 pozvanih mutacija | da |
| `comparison_chr2L.csv` | usporedba sa STAR-om po readu | ne, uz uzorak `comparison_chr2L.sample.csv` |
| `run_full_time_mem.log` | mjerenje izvođenja bez pozivanja mutacija | da |
| `run_full_with_mutations_time_mem.log` | mjerenje cjelovitog izvođenja | da |

### Testiranje

```bash
pytest tests/ -v
```

16 testova, svi prolaze.

| Datoteka | Broj testova | Sadržaj |
|---|---|---|
| `test_suffix_array.py` | 5 | poznati primjeri, usporedba SA-IS-a s naivnim postupkom na nasumičnim nizovima, prazan ulaz, podudarnost backenda `divsufsort` i `sais` |
| `test_alignment.py` | 6 | egzaktno podudaranje, poravnanje s nepodudaranjima, nepostojeći read, reverzni komplement, skupno poravnanje s uklanjanjem duplikata, provjera ulaznih vrijednosti |
| `test_alignment_fuzz.py` | 1 | usporedba s iscrpnim postupkom na nasumično generiranim slučajevima |
| `test_mutation_calling.py` | 4 | poziv uz dovoljnu pokrivenost i udio, odbacivanje pojedinačne greške, odbacivanje niske pokrivenosti, readovi s negativnog lanca |

Dodatno, `moj_primjer.py` izvodi cjelovit postupak nad kromosomom chrM iz
sklopa hg38 s poznatim rješenjem: poravnava read s umetnutom razlikom
natrag na poznatu poziciju te nad skupom od 200 readova s dvije zasađene
mutacije i pet pojedinačnih grešaka provjerava jesu li pozvane obje
mutacije i nijedna greška. Skripta zahtijeva lokalnu datoteku
`hg38/hg38.u.fa`, koja nije uključena u repozitorij.

```bash
python3 moj_primjer.py
```

### Uloga STAR-a

STAR je vanjski, neovisno razvijen aligner koji se u ovom projektu koristi
isključivo kao referentno rješenje. Nije mijenjan osim popravkom potrebnim
za izgradnju na macOS-u i ne sudjeluje u izračunu rezultata vlastite
implementacije.

Usporedba se svodi na provjeru smješta li vlastiti aligner readove na iste
koordinate kao STAR. Budući da readovi nisu simulirani, stvarno biološko
podrijetlo pojedinog reada nije poznato, pa se rezultati izražavaju kao
podudarnost sa STAR-om, a ne kao apsolutna točnost.

| Dokument | Sadržaj |
|---|---|
| [`REPORT.md`](REPORT.md) | puni izvještaj usporedbe, oba opsega: kromosom 2L i cijeli genom |
| [`star_comparison/README.md`](star_comparison/README.md) | izgradnja i pokretanje STAR-a, postupak usporedbe |

### Ograničenja

| Ograničenje | Posljedica |
|---|---|
| Nisu podržane insercije, delecije ni splicing | readovi koji premošćuju granice egzona ne mogu se poravnati; usporedba sa STAR-om moguća je samo na readovima bez tih operacija |
| Indeks obuhvaća jedan kromosom | readovi s ostalih kromosoma ostaju nepovezani ili se mogu pogrešno smjestiti unutar indeksiranog kromosoma |
| Poravnanje je jednodretveno i izvedeno u Pythonu | izvođenje je bitno sporije od STAR-a |
| Pozivanje mutacija drži sve readove u memoriji | potrošnja memorije raste s brojem readova |
| Nema procjene pouzdanosti poravnanja | ne postoji ekvivalent MAPQ vrijednosti, pa se poravnanja ne mogu filtrirati po kvaliteti |
| Pozvane mutacije nisu neovisno provjerene | rezultat je izlaz algoritma, a ne potvrđen skup bioloških varijanti |

Mjerenja i detaljna analiza ograničenja nalaze se u
[`REPORT.md`](REPORT.md#prednosti-i-ograničenja).
