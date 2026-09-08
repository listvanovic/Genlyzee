"""Usporedba poravnanja vlastitog alignera (run_pipeline.py) sa STAR-ovim, za jedan kromosom

Iz STAR-ovog SAM-a uzima samo primarna poravnanja na zadanom kromosomu (FLAG & 0x100 == 0,
FLAG & 0x800 == 0, FLAG & 0x4 == 0) čiji je CIGAR točno <duljina reada>M, i uspoređuje ih s
izlazom vlastitog alignera po: koordinati (POS), orijentaciji (lanac) i broju nepodudaranja

STAR skraćuje mate-sufiks iz FASTQ zaglavlja ("...#0/1" -> "...#0"), pa se ID normalizira
prije spajanja
"""
from __future__ import annotations

import argparse
import csv
import re
from typing import Dict, Tuple

CIGAR_FULL_MATCH = re.compile(r"^(\d+)M$")

FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def normalize_read_id(read_id: str) -> str:
    """Makni mate-sufiks ('/1', '/2') koji STAR ne zadržava u QNAME polju"""
    return read_id[:-2] if read_id.endswith(("/1", "/2")) else read_id


def load_star_alignments(sam_path: str, chromosome: str) -> Dict[str, Tuple[int, str, int]]:
    """Vrati {read_id: (pos_1based, strand, nM)} za čiste primarne alignmente na `chromosome`"""
    star: Dict[str, Tuple[int, str, int]] = {}
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != chromosome:
                continue
            flag = int(fields[1])
            if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED):
                continue
            match = CIGAR_FULL_MATCH.match(fields[5])
            if not (match and int(match.group(1)) == len(fields[9])):
                continue

            mismatches = -1
            for tag in fields[11:]:
                if tag.startswith("nM:i:"):
                    mismatches = int(tag[5:])
                    break
            strand = "-" if flag & FLAG_REVERSE else "+"
            star[normalize_read_id(fields[0])] = (int(fields[3]), strand, mismatches)
    return star


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sam_path", help="STAR .sam")
    parser.add_argument("my_path", help="TSV vlastitog alignera (run_pipeline.py --alignments-out)")
    parser.add_argument("chromosome", help="kromosom na kojem se uspoređuje (npr. 2L)")
    parser.add_argument("out_path", help="izlazna CSV s usporedbom po readu")
    parser.add_argument("--tolerance", type=int, default=0,
                        help="dopuštena razlika koordinata koja se još broji kao pogodak (zadano 0)")
    args = parser.parse_args()

    star = load_star_alignments(args.sam_path, args.chromosome)
    print(f"STAR: {len(star):,} čistih primarnih poravnanja na kromosomu {args.chromosome}")

    stats = {
        "star_total": len(star),
        "both_aligned": 0,
        "same_pos": 0,
        "diff_pos": 0,
        "same_strand": 0,
        "diff_strand": 0,
        "same_mismatches": 0,
        "my_multi": 0,
        "my_unmapped": 0,
        "my_missing": 0,
    }

    seen = set()
    with open(args.my_path) as my_file, open(args.out_path, "w", newline="") as out:
        next(my_file)
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["read_id", "star_pos", "star_strand", "star_nM",
                         "my_status", "my_pos", "my_strand", "my_mismatches", "match"])

        for line in my_file:
            read_id, status, pos, strand, mismatches, _hits = line.rstrip("\n").split("\t")
            key = normalize_read_id(read_id)
            if key not in star:
                continue  # STAR ovaj read nije čisto poravnao na ovaj kromosom
            seen.add(key)
            star_pos, star_strand, star_nm = star[key]

            is_match = False
            if status == "unique":
                stats["both_aligned"] += 1
                my_pos = int(pos)
                if abs(my_pos - star_pos) <= args.tolerance:
                    stats["same_pos"] += 1
                    is_match = True
                else:
                    stats["diff_pos"] += 1
                if strand == star_strand:
                    stats["same_strand"] += 1
                else:
                    stats["diff_strand"] += 1
                if int(mismatches) == star_nm:
                    stats["same_mismatches"] += 1
            elif status == "multi":
                stats["my_multi"] += 1
            else:
                stats["my_unmapped"] += 1

            writer.writerow([key, star_pos, star_strand, star_nm,
                             status, pos, strand, mismatches, "yes" if is_match else "no"])

    stats["my_missing"] = stats["star_total"] - len(seen)

    total = stats["star_total"]
    print(f"\nUsporedba na kromosomu {args.chromosome} (tolerancija {args.tolerance} baza)")
    print(f"  STAR čistih poravnanja (nazivnik):      {total:,}")
    print(f"  oba programa poravnala jedinstveno:     {stats['both_aligned']:,}")
    print(f"  ista koordinata:                        {stats['same_pos']:,}")
    print(f"  različita koordinata:                   {stats['diff_pos']:,}")
    print(f"  ista orijentacija:                      {stats['same_strand']:,}")
    print(f"  različita orijentacija:                 {stats['diff_strand']:,}")
    print(f"  isti broj mismatcha:                    {stats['same_mismatches']:,}")
    print(f"  vlastiti: multi-mapped:                 {stats['my_multi']:,}")
    print(f"  vlastiti: nije poravnao:                {stats['my_unmapped']:,}")
    print(f"  nema zapisa u izlazu vlastitog:         {stats['my_missing']:,}")
    print(f"\n  Točnost = ista koordinata / STAR čistih poravnanja = "
          f"{stats['same_pos']:,} / {total:,} = {100 * stats['same_pos'] / total:.4f}%")
    if stats["both_aligned"]:
        print(f"  Točnost među jedinstveno poravnatima = {stats['same_pos']:,} / "
              f"{stats['both_aligned']:,} = {100 * stats['same_pos'] / stats['both_aligned']:.4f}%")
