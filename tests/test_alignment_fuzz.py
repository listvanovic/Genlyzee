""" Randomizirani fuzz test za genome_index.alignment

Uspoređuje seed-and-extend izlaz od `align_read` s brute-force oracleom (iscrpno pretraživanje svakog prozora, oba lanca) kroz nasumično generirane
sekvence i readove (egzaktni podstringovi, mutirani forward/reverse-complement readovi, i čisto nasumični readovi). Fiksni seed radi
reproducibilnosti -- ovo je ono što je tijekom razvoja potvrdilo garanciju potpunosti pigeonhole principa (560 pokušaja, 0 neuspjeha) i sada čuva
od regresija
"""
from __future__ import annotations

import os
import random
import sys
from typing import List, Set, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome_index.alignment import align_read, reverse_complement
from genome_index.suffix_array import build_suffix_array_from_sequence

FUZZ_SEED = 20240115  # proizvoljan, fiksiran radi reproducibilnosti
NUM_TRIALS = 560
BASES = "ACGT"


def _brute_force_align(
    read: str, sequence: str, max_mismatches: int
) -> Set[Tuple[int, int, str]]:
    """Oracle: provjeri svaki prozor na oba lanca izravnom usporedbom, bez ikakvog suffix arraya. Vraća torke (position, mismatches, strand)
    (odražava AlignmentResult) kao skup, za usporedbu neovisnu o redoslijedu"""
    results: Set[Tuple[int, int, str]] = set()
    read_len = len(read)
    rev_read = reverse_complement(read)
    for start in range(len(sequence) - read_len + 1):
        window = sequence[start:start + read_len]
        fwd_mismatches = sum(1 for a, b in zip(read, window) if a != b)
        if fwd_mismatches <= max_mismatches:
            results.add((start, fwd_mismatches, "+"))
        rev_mismatches = sum(1 for a, b in zip(rev_read, window) if a != b)
        if rev_mismatches <= max_mismatches:
            results.add((start, rev_mismatches, "-"))
    return results


def _random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def _mutate(bases: List[str], rng: random.Random, num_mutations: int) -> List[str]:
    mutated = bases.copy()
    for site in rng.sample(range(len(mutated)), num_mutations):
        original = mutated[site]
        mutated[site] = rng.choice([b for b in BASES if b != original])
    return mutated


def test_fuzz_matches_brute_force_oracle() -> None:
    # rezultat od align_read mora se točno poklapati s onim od brute-force oraclea, kroz NUM_TRIALS randomiziranih pokušaja koji pokrivaju
    # egzaktne/mutirane/RC/nasumične readove
    rng = random.Random(FUZZ_SEED)

    for trial in range(NUM_TRIALS):
        seq_len = rng.randint(200, 600)
        sequence = _random_sequence(rng, seq_len)
        sa = build_suffix_array_from_sequence(sequence, method="sais")

        read_len = rng.randint(10, 40)
        max_mismatches = rng.randint(0, 3)
        kind = rng.choice(["exact", "forward_mutated", "reverse_complement", "random"])

        true_start = rng.randint(0, seq_len - read_len)
        window = list(sequence[true_start:true_start + read_len])

        if kind == "exact":
            read = "".join(window)
        elif kind == "forward_mutated":
            num_mutations = rng.randint(0, max_mismatches)
            read = "".join(_mutate(window, rng, num_mutations))
        elif kind == "reverse_complement":
            num_mutations = rng.randint(0, max_mismatches)
            mutated = _mutate(window, rng, num_mutations)
            read = reverse_complement("".join(mutated))
        else:  # "random": neovisne nasumične baze, ne izrezane iz sequence
            read = _random_sequence(rng, read_len)

        actual = {
            (r.position, r.mismatches, r.strand)
            for r in align_read(read, sequence, sa, max_mismatches=max_mismatches)
        }
        expected = _brute_force_align(read, sequence, max_mismatches)

        assert actual == expected, (
            f"trial {trial} (kind={kind}, seq_len={seq_len}, read_len={read_len}, "
            f"max_mismatches={max_mismatches}, read={read!r}): "
            f"align_read output {actual} != brute-force oracle {expected}"
        )


if __name__ == "__main__":
    test_fuzz_matches_brute_force_oracle()
    print(f"Fuzz test passed: {NUM_TRIALS} trials, seed={FUZZ_SEED}.")
