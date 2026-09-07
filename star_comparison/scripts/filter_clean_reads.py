"""Izdvajanje "čistih" readova iz STAR-ovog SAM-a: onih čiji je CIGAR točno "<duljina reada>M",
dakle bez insercija, delecija, soft/hard clippinga i splicanja

Gledaju se samo primarni alignmenti (FLAG & 0x100 == 0 i FLAG & 0x800 == 0), pa multi-mapped ili
chimerični read daje najviše jedan redak
"""
from __future__ import annotations

import argparse
import re
import sys
from typing import TextIO

CIGAR_FULL_MATCH = re.compile(r"^(\d+)M$")

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def is_clean_cigar(cigar: str, read_length: int) -> bool:
    match = CIGAR_FULL_MATCH.match(cigar)
    return match is not None and int(match.group(1)) == read_length


def filter_sam(sam_path: str, out_path: str, log: TextIO = sys.stderr) -> None:
    total_alignments = 0
    primary_alignments = 0
    clean_reads = 0

    with open(sam_path) as sam, open(out_path, "w") as out:
        out.write("read_id\tsequence\trname\tpos\n")
        for line in sam:
            if line.startswith("@"):
                continue
            total_alignments += 1

            fields = line.rstrip("\n").split("\t")
            qname, flag, rname, pos, _mapq, cigar = fields[0:6]
            sequence = fields[9]
            flag = int(flag)

            if flag & (FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_SUPPLEMENTARY):
                continue
            primary_alignments += 1

            if not is_clean_cigar(cigar, len(sequence)):
                continue
            clean_reads += 1
            out.write(f"{qname}\t{sequence}\t{rname}\t{pos}\n")

            if total_alignments % 2_000_000 == 0:
                print(f"  ...obrađeno {total_alignments:,} SAM linija", file=log)

    print(f"Ukupno SAM alignment linija: {total_alignments:,}", file=log)
    print(f"Primarnih alignmenata:       {primary_alignments:,}", file=log)
    print(f"Čistih (full-match) reada:   {clean_reads:,}", file=log)
    print(f"Udio čistih u primarnima:    {100 * clean_reads / primary_alignments:.2f}%", file=log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sam_path", help="ulazni STAR .sam file")
    parser.add_argument("out_path", help="izlazni TSV: read_id, sequence, rname, pos")
    args = parser.parse_args()

    filter_sam(args.sam_path, args.out_path)
