"""Izdvajanje STAR-ovog NH taga (broj lokusa na koje je read mapiran) za svaki čisti read

Koristi identičan filter kao filter_clean_reads.py (primarni, CIGAR=<duljina>M), pa je izlaz
u istom poretku redaka kao clean_reads.tsv i može se uz njega zipati bez hash joina

Treba za provjeru je li STAR-ova referentna pozicija bila jedina koju je STAR razmatrao (NH:i:1)
ili je STAR i sam bio neodlučan (NH>=2), jer u drugom slučaju stroga usporedba s jednim STAR POS-om
može biti nepravedna prema drugom alignerau
"""
from __future__ import annotations

import argparse
import re
import sys

CIGAR_FULL_MATCH = re.compile(r"^(\d+)M$")

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def extract_nh(fields: list[str]) -> int:
    for tag in fields[11:]:
        if tag.startswith("NH:i:"):
            return int(tag[5:])
    raise ValueError(f"nema NH taga u zapisu za {fields[0]!r}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sam_path")
    parser.add_argument("out_path")
    args = parser.parse_args()

    n = 0
    with open(args.sam_path) as sam, open(args.out_path, "w") as out:
        out.write("read_id\tnh\n")
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & (FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_SUPPLEMENTARY):
                continue
            cigar = fields[5]
            sequence = fields[9]
            match = CIGAR_FULL_MATCH.match(cigar)
            if not (match and int(match.group(1)) == len(sequence)):
                continue
            out.write(f"{fields[0]}\t{extract_nh(fields)}\n")
            n += 1
            if n % 2_000_000 == 0:
                print(f"  ...obrađeno {n:,} čistih reada", file=sys.stderr)

    print(f"Zapisano NH za {n:,} čistih reada", file=sys.stderr)
