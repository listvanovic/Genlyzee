"""Poravnanje čistih readova vlastitim alignerom nad prethodno izgrađenim indeksom (build_my_index.py)

Za svaki read: točno jedno poravnanje unutar max_mismatches -> "unique" (pozicija se bilježi),
više njih -> "multi", nijedno -> "unmapped"
Ista podjela kao STAR-ova unique/multi-mapped/unmapped, da se izlazi poklope 1:1 pri usporedbi koordinata
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from genome_index.alignment import align_read, load_reads
from genome_concat import global_pos_to_contig

INDEX_DIR = "star_comparison/my_index"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", help="ulazni FASTA s readovima")
    parser.add_argument("out_path", help="izlazni TSV s pozicijama")
    parser.add_argument("--max-mismatches", type=int, default=3)
    parser.add_argument("--max-seed-hits", type=int, default=1000)
    parser.add_argument("--index-dir", default=INDEX_DIR)
    args = parser.parse_args()

    print("Učitavam indeks ...", file=sys.stderr)
    start = time.time()
    suffix_array = np.load(f"{args.index_dir}/suffix_array.npy")
    with open(f"{args.index_dir}/contigs.json") as handle:
        contigs = json.load(handle)
    names, offsets, lengths = contigs["names"], contigs["offsets"], contigs["lengths"]
    with open(f"{args.index_dir}/genome.txt") as handle:
        genome = handle.read()
    print(f"  učitano za {time.time() - start:.1f}s ({len(genome):,} baza)", file=sys.stderr)

    reads = load_reads(args.fasta_path)
    print(f"Poravnavam {len(reads):,} reada (max_mismatches={args.max_mismatches}) ...", file=sys.stderr)

    unique = multi = unmapped = 0
    start = time.time()
    with open(args.out_path, "w") as out:
        out.write("read_id\tmy_status\tmy_rname\tmy_pos\tmy_mismatches\tmy_num_hits\n")
        for i, (read_id, sequence) in enumerate(reads):
            hits = align_read(sequence, genome, suffix_array,
                              max_mismatches=args.max_mismatches,
                              max_seed_hits=args.max_seed_hits)
            if not hits:
                unmapped += 1
                out.write(f"{read_id}\tunmapped\t\t\t\t0\n")
            elif len(hits) == 1:
                unique += 1
                best = hits[0]
                mapped = global_pos_to_contig(best.position, names, offsets, lengths)
                if mapped is None:
                    out.write(f"{read_id}\tunique_in_spacer\t\t\t{best.mismatches}\t1\n")
                else:
                    rname, pos = mapped
                    out.write(f"{read_id}\tunique\t{rname}\t{pos}\t{best.mismatches}\t1\n")
            else:
                multi += 1
                out.write(f"{read_id}\tmulti\t\t\t\t{len(hits)}\n")

            if (i + 1) % 500_000 == 0:
                elapsed = time.time() - start
                print(f"  ...{i + 1:,} reada za {elapsed:.0f}s ({(i + 1) / elapsed:.0f} reada/s)",
                      file=sys.stderr)

    elapsed = time.time() - start
    total = len(reads)
    print(f"\nGotovo za {elapsed:.1f}s ({total / elapsed:.0f} reada/s)", file=sys.stderr)
    print(f"Unique:   {unique:,} ({100 * unique / total:.2f}%)", file=sys.stderr)
    print(f"Multi:    {multi:,} ({100 * multi / total:.2f}%)", file=sys.stderr)
    print(f"Unmapped: {unmapped:,} ({100 * unmapped / total:.2f}%)", file=sys.stderr)
