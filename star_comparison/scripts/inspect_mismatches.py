"""Provjera tko je u pravu kod neslaganja: za svaki read koji su oba alignera jedinstveno smjestila,
ali na različite pozicije, broji stvarna nepodudaranja reada naspram genoma na objema pozicijama

Ovo je jedina provjera koja ne uzima STAR-ovu poziciju zdravo za gotovo — bez nje se ne može reći
je li razlika greška vlastitog alignera ili je STAR odabrao lošije poravnanje
"""
from __future__ import annotations

import argparse
import csv
import sys

sys.path.insert(0, ".")

from genome_index.alignment import reverse_complement
from genome_concat import load_genome_concat

FASTA_PATH = "Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"


def count_mismatches(genome: str, contig_start: int, pos: int, read: str) -> int | None:
    """Nepodudaranja reada naspram genoma na 1-baziranoj poziciji `pos` unutar kontiga
    Uzima bolju od dvije orijentacije, jer SAM sekvenca je već okrenuta prema forward lancu"""
    start = contig_start + pos - 1
    window = genome[start:start + len(read)]
    if len(window) != len(read):
        return None
    forward = sum(1 for a, b in zip(read, window) if a != b)
    reverse = sum(1 for a, b in zip(reverse_complement(read), window) if a != b)
    return min(forward, reverse)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mismatches_path", help="CSV s neslaganjima (iz compare_coordinates.py)")
    parser.add_argument("clean_reads_path", help="clean_reads.tsv, radi sekvenci readova")
    parser.add_argument("out_path", help="izlazna CSV s brojem nepodudaranja i presudom")
    parser.add_argument("--fasta", default=FASTA_PATH)
    args = parser.parse_args()

    with open(args.mismatches_path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        disagreements = list(reader)

    if not disagreements:
        print("Nema neslaganja za provjeru.")
        sys.exit()

    wanted = {row[0] for row in disagreements}
    sequences = {}
    with open(args.clean_reads_path) as handle:
        next(handle)
        for line in handle:
            read_id, sequence, _rname, _pos = line.rstrip("\n").split("\t")
            if read_id in wanted:
                sequences[read_id] = sequence
                if len(sequences) == len(wanted):
                    break

    print("Učitavam genom ...", file=sys.stderr)
    names, offsets, _lengths, genome = load_genome_concat(args.fasta)
    contig_start = dict(zip(names, offsets))

    star_better = my_better = tied = 0
    with open(args.out_path, "w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["read_id", "star_rname", "star_pos", "star_mismatches",
                         "my_rname", "my_pos", "my_mismatches", "bolja_pozicija"])

        for read_id, star_rname, star_pos, my_rname, my_pos, _delta in disagreements:
            read = sequences[read_id]
            star_mm = count_mismatches(genome, contig_start[star_rname], int(star_pos), read)
            my_mm = count_mismatches(genome, contig_start[my_rname], int(my_pos), read)

            if star_mm < my_mm:
                verdict = "STAR"
                star_better += 1
            elif my_mm < star_mm:
                verdict = "vlastiti"
                my_better += 1
            else:
                verdict = "izjednačeno"
                tied += 1

            writer.writerow([read_id, star_rname, star_pos, star_mm,
                             my_rname, my_pos, my_mm, verdict])

    total = len(disagreements)
    print(f"Provjereno neslaganja: {total}")
    print(f"  STAR ima manje nepodudaranja:            {star_better} ({100 * star_better / total:.1f}%)")
    print(f"  vlastiti aligner ima manje nepodudaranja: {my_better} ({100 * my_better / total:.1f}%)")
    print(f"  izjednačeno:                             {tied} ({100 * tied / total:.1f}%)")
