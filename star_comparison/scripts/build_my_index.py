"""Izgradnja suffix arraya vlastitog alignera nad cijelim genomom i spremanje na disk,
da ga svako poravnanje ne mora graditi ispočetka

Pandan STAR-ovom genomeGenerate koraku; vrijeme i memorija se mjere odvojeno od samog poravnanja (run_my_aligner.py),
jednako kao što je i STAR pipeline mjerio genomeGenerate odvojeno od alignReads
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from genome_index.suffix_array import build_suffix_array_from_sequence
from genome_concat import load_genome_concat

FASTA_PATH = "Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
INDEX_DIR = "star_comparison/my_index"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default=FASTA_PATH, help="referentni genom (multi-FASTA)")
    parser.add_argument("--index-dir", default=INDEX_DIR, help="direktorij u koji se sprema indeks")
    args = parser.parse_args()

    print("Učitavam genom i spajam kontige ...", file=sys.stderr)
    start = time.time()
    names, offsets, lengths, genome = load_genome_concat(args.fasta)
    print(f"  {len(names)} kontiga, {len(genome):,} baza (uklj. spacere), {time.time() - start:.1f}s",
          file=sys.stderr)

    print("Gradim suffix array (divsufsort) ...", file=sys.stderr)
    start = time.time()
    suffix_array = build_suffix_array_from_sequence(genome, method="divsufsort")
    print(f"  gotovo za {time.time() - start:.1f}s", file=sys.stderr)

    np.save(f"{args.index_dir}/suffix_array.npy", suffix_array)
    with open(f"{args.index_dir}/contigs.json", "w") as handle:
        json.dump({"names": names, "offsets": offsets, "lengths": lengths}, handle)
    with open(f"{args.index_dir}/genome.txt", "w") as handle:
        handle.write(genome)

    print(f"Indeks zapisan u {args.index_dir}/", file=sys.stderr)
