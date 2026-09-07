"""Usporedba koordinata koje je za svaki čisti read vratio vlastiti aligner s onima iz STAR-a

Read se broji kao pogodak ako je RNAME isti i |my_pos - star_pos| <= --tolerance (zadano 0, dakle identična pozicija)

Prijavljuju se dvije točnosti jer odgovaraju na različita pitanja:
  - stroga     = pogoci / svi čisti readovi; read koji vlastiti aligner nije jedinstveno smjestio
                 (multi ili unmapped) računa se kao promašaj, jednako kao i kriva pozicija
  - unique-only = pogoci / readovi za koje se vlastiti aligner odlučio na jednu poziciju;
                 mjeri koliko je točan kad se odluči, odvojeno od toga koliko se često odlučuje

clean_reads.tsv i my_alignments.tsv su u istom poretku redaka (FASTA je nastao iz clean_reads.tsv bez presortiranja,
a aligner ga čita sekvencijalno), pa se ide red po red umjesto hash joina, uz provjeru read_id-a kao osigurač
"""
from __future__ import annotations

import argparse
import csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("star_path", help="clean_reads.tsv (STAR-ove pozicije)")
    parser.add_argument("my_path", help="my_alignments.tsv (izlaz vlastitog alignera)")
    parser.add_argument("out_path", help="izlazna CSV tablica usporedbe")
    parser.add_argument("--mismatches-out",
                        help="zasebna CSV s readovima koje su oba alignera jedinstveno smjestila, "
                             "ali na različite pozicije")
    parser.add_argument("--tolerance", type=int, default=0,
                        help="najveća dopuštena razlika pozicija koja se još broji kao pogodak")
    args = parser.parse_args()

    total = matches = unique_calls = unique_matches = 0
    status_counts: dict[str, int] = {}
    disagreements = []

    # lineterminator="\n": bez toga csv.writer piše CRLF, na čemu se lome awk/grep/cut
    with open(args.star_path) as star_file, open(args.my_path) as my_file, \
            open(args.out_path, "w", newline="") as out:
        next(star_file)
        next(my_file)
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["read_id", "star_rname", "star_pos", "my_status",
                         "my_rname", "my_pos", "match"])

        for star_line, my_line in zip(star_file, my_file):
            read_id, _sequence, star_rname, star_pos = star_line.rstrip("\n").split("\t")
            my_read_id, status, my_rname, my_pos, _mismatches, _hits = my_line.rstrip("\n").split("\t")
            if read_id != my_read_id:
                raise ValueError(
                    f"Redoslijed redaka se razišao: {read_id!r} u {args.star_path!r} "
                    f"naspram {my_read_id!r} u {args.my_path!r}"
                )

            total += 1
            status_counts[status] = status_counts.get(status, 0) + 1

            is_match = False
            if status == "unique":
                unique_calls += 1
                if my_rname == star_rname and abs(int(my_pos) - int(star_pos)) <= args.tolerance:
                    is_match = True
                    matches += 1
                    unique_matches += 1
                else:
                    disagreements.append([read_id, star_rname, star_pos, my_rname, my_pos])

            writer.writerow([read_id, star_rname, star_pos, status, my_rname, my_pos,
                             "yes" if is_match else "no"])

    if args.mismatches_out:
        with open(args.mismatches_out, "w", newline="") as out:
            writer = csv.writer(out, lineterminator="\n")
            writer.writerow(["read_id", "star_rname", "star_pos", "my_rname", "my_pos", "delta"])
            for read_id, star_rname, star_pos, my_rname, my_pos in disagreements:
                delta = int(my_pos) - int(star_pos) if my_rname == star_rname else ""
                writer.writerow([read_id, star_rname, star_pos, my_rname, my_pos, delta])

    print(f"Tolerancija: {args.tolerance} baza")
    print(f"Ukupno čistih reada (STAR referenca): {total:,}")
    print(f"Pogodaka: {matches:,}")
    print(f"STROGA točnost (pogoci / svi čisti readovi): {100 * matches / total:.2f}%")
    print()
    print("Raspodjela ishoda vlastitog alignera:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count:,} ({100 * count / total:.2f}%)")
    print()
    print(f"Od {unique_calls:,} jedinstvenih poravnanja, {unique_matches:,} se poklapa sa STAR-om "
          f"({100 * unique_matches / unique_calls:.2f}%) [unique-only točnost]")
