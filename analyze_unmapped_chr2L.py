"""Zasebna analiza uzroka za readove koje run_pipeline.py klasificira kao 'unmapped' na 2L.

Ne mijenja genome_index/ - samo poziva njegove postojeće (interne) funkcije da razdvoji
dva moguća uzroka:
  A) nijedan seed (ni na jednom lancu) nije imao egzaktno podudaranje - "no_exact_seed_match"
  B) neki seed je imao podudaranja, ali sva takva su odbačena jer ih je bilo vise od
     max_seed_hits - "seed_rejected_frequent" (nijedan kandidat nikad nije generiran)
  C) barem jedan kandidat je generiran (seed unutar max_seed_hits), ali proširenje na
     cijeli read nije zadovoljilo max_mismatches (ili je ispao izvan granica sekvence)
     - "extension_failed"

Cita FASTQ i alignments TSV usporedno, redak po redak (isti redoslijed kao u run_pipeline.py),
i klasificira samo readove sa status=unmapped.
"""
from __future__ import annotations

import sys
import time

from genome_index.alignment import (
    _split_into_seeds,
    _sa_search_bounds,
    reverse_complement,
)
from genome_index.fasta_io import parse_fasta_chromosome
from genome_index.suffix_array import build_suffix_array_from_sequence

MAX_MISMATCHES = 3
MAX_SEED_HITS = 1000

FASTA_PATH = "Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
CHROMOSOME = "2L"
READS_PATH = "Fly/IFM48h_1.fastq"
ALIGNMENTS_PATH = "results/chr2L/my_alignments_full.tsv"
OUT_PATH = "results/chr2L/unmapped_reason_breakdown.tsv"


def iter_fastq(path):
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().rstrip()
            handle.readline()
            handle.readline()
            yield header[1:].split(None, 1)[0], seq.upper()


def classify(read: str, sequence: str, suffix_array) -> str:
    have_candidates = False
    any_seed_matched = False
    for strand_seq in (read, reverse_complement(read)):
        for offset, seed in _split_into_seeds(strand_seq, MAX_MISMATCHES):
            if not seed:
                continue
            first, last = _sa_search_bounds(sequence, suffix_array, seed)
            if first > last:
                continue
            any_seed_matched = True
            count = last - first + 1
            if count <= MAX_SEED_HITS:
                have_candidates = True
    if have_candidates:
        return "extension_failed"
    if any_seed_matched:
        return "seed_rejected_frequent"
    return "no_exact_seed_match"


def main() -> None:
    log = sys.stderr
    print(f"Izdvajam kromosom {CHROMOSOME!r} ...", file=log)
    sequence = parse_fasta_chromosome(FASTA_PATH, CHROMOSOME)
    print(f"  {len(sequence):,} baza", file=log)

    print("Gradim suffix array (divsufsort) ...", file=log)
    suffix_array = build_suffix_array_from_sequence(sequence, method="divsufsort")
    print("  gotovo", file=log)

    counts = {"no_exact_seed_match": 0, "seed_rejected_frequent": 0, "extension_failed": 0}
    total_unmapped = 0
    total_rows = 0

    fastq_iter = iter_fastq(READS_PATH)
    start = time.time()
    with open(ALIGNMENTS_PATH) as aln, open(OUT_PATH, "w") as out:
        header = aln.readline()
        assert header.startswith("read_id\tstatus"), header
        out.write("read_id\treason\n")

        for aln_line in aln:
            total_rows += 1
            read_id_aln, status = aln_line.split("\t", 2)[:2]
            fastq_id, fastq_seq = next(fastq_iter)
            assert fastq_id == read_id_aln, (fastq_id, read_id_aln, total_rows)

            if status == "unmapped":
                total_unmapped += 1
                reason = classify(fastq_seq, sequence, suffix_array)
                counts[reason] += 1
                out.write(f"{read_id_aln}\t{reason}\n")

            if total_rows % 1_000_000 == 0:
                elapsed = time.time() - start
                print(f"  ...{total_rows:,} redaka obradeno ({total_unmapped:,} unmapped do sad) "
                      f"za {elapsed:.0f}s", file=log)

    elapsed = time.time() - start
    print(f"\nGotovo za {elapsed:.2f}s. Ukupno redaka: {total_rows:,}, "
          f"unmapped: {total_unmapped:,}", file=log)
    for k, v in counts.items():
        pct = 100 * v / total_unmapped if total_unmapped else 0.0
        print(f"  {k}: {v:,} ({pct:.2f}%)", file=log)


if __name__ == "__main__":
    main()
