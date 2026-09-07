"""Pretvorba filtriranih čistih readova (clean_reads.tsv) u FASTA, format koji vlastiti aligner očekuje
preko genome_index.alignment.load_reads()"""
from __future__ import annotations

import argparse
import sys


def write_fasta(tsv_path: str, fasta_path: str) -> int:
    count = 0
    with open(tsv_path) as tsv, open(fasta_path, "w") as fasta:
        next(tsv)  # zaglavlje
        for line in tsv:
            read_id, sequence, _rname, _pos = line.rstrip("\n").split("\t")
            fasta.write(f">{read_id}\n{sequence}\n")
            count += 1
    return count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv_path", help="ulazni clean_reads.tsv")
    parser.add_argument("fasta_path", help="izlazni FASTA file")
    args = parser.parse_args()

    count = write_fasta(args.tsv_path, args.fasta_path)
    print(f"Zapisano {count:,} reada u {args.fasta_path}", file=sys.stderr)
