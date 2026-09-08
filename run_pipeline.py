"""Glavni program: FASTA + naziv kromosoma + datoteka s readovima -> poravnanja i mutacije

Izvodi sva tri vlastita algoritma nad jednim kromosomom:
  1. parse_fasta_chromosome  - izdvajanje tražene sekvence iz multi-FASTA datoteke
  2. build_suffix_array_from_sequence - izgradnja suffix arraya te sekvence
  3. align_read              - poravnanje svakog reada (seed-and-extend, oba lanca)
  4. call_mutations          - pozivanje SNP-ova iz pileupa jedinstveno poravnatih reada

Readovi se čitaju streamingom (FASTQ ili FASTA), pa memorija ne raste s veličinom ulazne
datoteke osim kad je uključeno pozivanje mutacija, koje po svojoj prirodi treba sve readove odjednom

Pozicije: suffix array radi s 0-baziranim indeksima, a SAM/izlaz ovog programa s 1-baziranima
"""
from __future__ import annotations

import argparse
import sys
import time
from typing import Dict, Iterator, List, TextIO, Tuple

from genome_index.alignment import AlignmentResult, align_read
from genome_index.fasta_io import parse_fasta_chromosome
from genome_index.mutation_calling import call_mutations
from genome_index.suffix_array import build_suffix_array_from_sequence


def iter_reads(path: str) -> Iterator[Tuple[str, str]]:
    """Čitaj (read_id, sekvenca) iz FASTQ ili FASTA datoteke, jedan po jedan
    Format se prepoznaje iz prvog znaka prve neprazne linije"""
    with open(path) as handle:
        first = handle.readline()
        while first and not first.strip():
            first = handle.readline()
        if not first:
            return

        if first.startswith("@"):  # FASTQ: 4 linije po readu
            while first:
                read_id = first[1:].split(None, 1)[0]
                sequence = handle.readline().rstrip()
                handle.readline()  # '+' linija
                handle.readline()  # kvalitete
                yield read_id, sequence.upper()
                first = handle.readline()
        elif first.startswith(">"):  # FASTA: zaglavlje + jedna ili više linija sekvence
            read_id = first[1:].split(None, 1)[0]
            chunks: List[str] = []
            for line in handle:
                if line.startswith(">"):
                    yield read_id, "".join(chunks).upper()
                    read_id = line[1:].split(None, 1)[0]
                    chunks = []
                else:
                    chunks.append(line.rstrip())
            yield read_id, "".join(chunks).upper()
        else:
            raise ValueError(f"Nepoznat format datoteke s readovima: {path!r} (očekivano '@' ili '>')")


def align_all(
    reads_path: str,
    sequence: str,
    suffix_array,
    out: TextIO,
    max_mismatches: int,
    max_seed_hits: int,
    keep_for_mutations: bool,
    log: TextIO,
) -> Tuple[Dict[str, int], List[str], Dict[str, List[AlignmentResult]]]:
    """Poravnaj sve readove i zapiši rezultate; vrati statistiku i (opcionalno) ulaz za pozivanje mutacija"""
    counts = {"total": 0, "unique": 0, "multi": 0, "unmapped": 0, "skipped": 0}
    reads_for_mutations: List[str] = []
    alignments_for_mutations: Dict[str, List[AlignmentResult]] = {}

    out.write("read_id\tstatus\tpos_1based\tstrand\tmismatches\tnum_hits\n")
    start = time.time()
    for read_id, sequence_read in iter_reads(reads_path):
        counts["total"] += 1

        # align_read odbija readove duže od reference i one s bazama izvan ACGTN
        if len(sequence_read) > len(sequence) or set(sequence_read) - set("ACGTN"):
            counts["skipped"] += 1
            out.write(f"{read_id}\tskipped\t\t\t\t0\n")
            continue

        hits = align_read(sequence_read, sequence, suffix_array,
                          max_mismatches=max_mismatches, max_seed_hits=max_seed_hits)

        if keep_for_mutations:
            reads_for_mutations.append(sequence_read)
            alignments_for_mutations[sequence_read] = hits

        if not hits:
            counts["unmapped"] += 1
            out.write(f"{read_id}\tunmapped\t\t\t\t0\n")
        elif len(hits) == 1:
            counts["unique"] += 1
            hit = hits[0]
            # 0-bazirana pozicija u suffix arrayu -> 1-bazirana, kao SAM POS
            out.write(f"{read_id}\tunique\t{hit.position + 1}\t{hit.strand}\t{hit.mismatches}\t1\n")
        else:
            counts["multi"] += 1
            out.write(f"{read_id}\tmulti\t\t\t\t{len(hits)}\n")

        if counts["total"] % 200_000 == 0:
            elapsed = time.time() - start
            print(f"  ...{counts['total']:,} reada za {elapsed:.0f}s "
                  f"({counts['total'] / elapsed:.0f} reada/s)", file=log)

    return counts, reads_for_mutations, alignments_for_mutations


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta_path", help="referentni genom u multi-FASTA formatu")
    parser.add_argument("chromosome", help="naziv kromosoma koji se izdvaja iz FASTA datoteke (npr. 2L)")
    parser.add_argument("reads_path", help="datoteka s readovima (FASTQ ili FASTA)")
    parser.add_argument("--alignments-out", required=True, help="izlazni TSV s poravnanjima")
    parser.add_argument("--mutations-out", help="izlazni TSV s mutacijama; bez ovoga se mutacije ne pozivaju")
    parser.add_argument("--max-mismatches", type=int, default=3,
                        help="najviše dopuštenih nepodudaranja po readu (zadano 3)")
    parser.add_argument("--max-seed-hits", type=int, default=1000,
                        help="seed s više pogodaka od ovoga se preskače kao neinformativan (zadano 1000)")
    parser.add_argument("--min-coverage", type=int, default=3,
                        help="najmanja pokrivenost pozicije za pozivanje mutacije (zadano 3)")
    parser.add_argument("--min-variant-fraction", type=float, default=0.5,
                        help="najmanji udio reada koji podupiru alternativnu bazu (zadano 0.5)")
    parser.add_argument("--sa-method", default="divsufsort", choices=["divsufsort", "sais"],
                        help="backend za izgradnju suffix arraya (zadano divsufsort)")
    args = parser.parse_args()

    log = sys.stderr

    print(f"[1/4] Izdvajam kromosom {args.chromosome!r} iz {args.fasta_path} ...", file=log)
    start = time.time()
    sequence = parse_fasta_chromosome(args.fasta_path, args.chromosome)
    parse_seconds = time.time() - start
    print(f"      {len(sequence):,} baza, {parse_seconds:.2f}s", file=log)

    print(f"[2/4] Gradim suffix array ({args.sa_method}) ...", file=log)
    start = time.time()
    suffix_array = build_suffix_array_from_sequence(sequence, method=args.sa_method)
    index_seconds = time.time() - start
    print(f"      gotovo za {index_seconds:.2f}s", file=log)

    print(f"[3/4] Poravnavam readove iz {args.reads_path} (max_mismatches={args.max_mismatches}) ...", file=log)
    start = time.time()
    with open(args.alignments_out, "w") as out:
        counts, reads, alignments = align_all(
            args.reads_path, sequence, suffix_array, out,
            args.max_mismatches, args.max_seed_hits,
            keep_for_mutations=bool(args.mutations_out), log=log,
        )
    align_seconds = time.time() - start

    total = counts["total"]
    print(f"      {total:,} reada za {align_seconds:.2f}s "
          f"({total / align_seconds:.0f} reada/s)", file=log)
    for status in ("unique", "multi", "unmapped", "skipped"):
        print(f"      {status}: {counts[status]:,} ({100 * counts[status] / total:.2f}%)", file=log)

    mutation_seconds = 0.0
    num_mutations = 0
    if args.mutations_out:
        print(f"[4/4] Pozivam mutacije (min_coverage={args.min_coverage}, "
              f"min_variant_fraction={args.min_variant_fraction}) ...", file=log)
        start = time.time()
        calls = call_mutations(sequence, reads, alignments,
                               min_coverage=args.min_coverage,
                               min_variant_fraction=args.min_variant_fraction)
        mutation_seconds = time.time() - start
        num_mutations = len(calls)
        with open(args.mutations_out, "w") as out:
            out.write("chromosome\tpos_1based\tref_base\talt_base\tsupporting_reads\tcoverage\tvaf\n")
            for call in calls:
                out.write(f"{args.chromosome}\t{call.position + 1}\t{call.ref_base}\t{call.alt_base}\t"
                          f"{call.supporting_reads}\t{call.coverage}\t"
                          f"{call.variant_allele_fraction:.4f}\n")
        print(f"      {num_mutations:,} mutacija za {mutation_seconds:.2f}s", file=log)
    else:
        print("[4/4] Pozivanje mutacija preskočeno (nema --mutations-out)", file=log)

    print(f"\nUkupno: {parse_seconds + index_seconds + align_seconds + mutation_seconds:.2f}s "
          f"(FASTA {parse_seconds:.2f}s + indeks {index_seconds:.2f}s + "
          f"poravnanje {align_seconds:.2f}s + mutacije {mutation_seconds:.2f}s)", file=log)
