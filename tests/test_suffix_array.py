# Jedinični testovi za genome_index.suffix_array
from __future__ import annotations

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome_index.suffix_array import (
    build_suffix_array_from_sequence,
    build_suffix_array_naive,
)


def test_banana() -> None:
    """Primjer iz knjige: sufiksi od "banana$", sortirani leksikografski ('$' je najmanji), su:

        $        -> 6
        a$       -> 5
        ana$     -> 3
        anana$   -> 1
        banana$  -> 0
        na$      -> 4
        nana$    -> 2

    pa je SUFFIXARRAY("banana$") = [6, 5, 3, 1, 0, 4, 2].

    `build_suffix_array_from_sequence` prima string "banana" bez sentinela i sam upravlja svojim internim sentinelom, pa njen rezultat izostavlja
    poziciju samog sentinela (6) ali se inače slaže:
    SUFFIXARRAY("banana") = [5, 3, 1, 0, 4, 2].
    """
    assert build_suffix_array_naive("banana$") == [6, 5, 3, 1, 0, 4, 2]
    result = build_suffix_array_from_sequence("banana", method="sais")
    assert result.tolist() == [5, 3, 1, 0, 4, 2]


def test_panamabananas_from_chapter_9() -> None:
    """Slika 9.7 iz Compeau & Pevzner, Poglavlje 9:
    SUFFIXARRAY("panamabananas$") = [13, 5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12]

    Kao i u `test_banana`, `build_suffix_array_from_sequence` dobiva "panamabananas" bez sentinela, pa njen rezultat izostavlja poziciju samog
    sentinela (13).
    """
    expected_with_sentinel = [13, 5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12]
    assert build_suffix_array_naive("panamabananas$") == expected_with_sentinel

    expected = [p for p in expected_with_sentinel if p != 13]
    result = build_suffix_array_from_sequence("panamabananas", method="sais")
    assert result.tolist() == expected


def test_sais_matches_naive_on_random_strings() -> None:
    """Fuzz-testiraj SA-IS izgradnju u linearnom vremenu naspram naivne, izravno sortirane izgradnje na puno malih nasumičnih DNA-sličnih stringova,
    uključujući rubne slučajeve poput visoko repetitivnih sekvenci."""
    random.seed(0)
    alphabets = ["ACGT", "AB", "ACGTN"]
    lengths = list(range(1, 20)) + [50, 200]

    for letters in alphabets:
        for length in lengths:
            for _ in range(10):
                text = "".join(random.choice(letters) for _ in range(length))
                # izbaci poziciju samog sentinela (indeks `length`) da se
                # usporedi s build_suffix_array_from_sequence, koja je
                # već izbacuje
                expected = [
                    p for p in build_suffix_array_naive(text + "$") if p != length
                ]
                got = build_suffix_array_from_sequence(text, method="sais")
                assert got.tolist() == expected, (text, expected, got.tolist())

    # visoko repetitivni stringovi testiraju rekurzivni slučaj SA-IS-a
    for length in [1, 2, 3, 10, 100, 500]:
        text = "A" * length
        expected = [p for p in build_suffix_array_naive(text + "$") if p != length]
        got = build_suffix_array_from_sequence(text, method="sais")
        assert got.tolist() == expected


def test_empty_sequence() -> None:
    result = build_suffix_array_from_sequence("", method="sais")
    assert result.tolist() == []


def test_divsufsort_backend_matches_sais_if_available() -> None:
    try:
        import pydivsufsort  # noqa: F401
    except ImportError:
        return  # opcionalna ovisnost nije instalirana; preskoči

    random.seed(1)
    for _ in range(20):
        length = random.randint(1, 300)
        text = "".join(random.choice("ACGT") for _ in range(length))
        sais_result = build_suffix_array_from_sequence(text, method="sais")
        divsufsort_result = build_suffix_array_from_sequence(text, method="divsufsort")
        assert sais_result.tolist() == divsufsort_result.tolist()


if __name__ == "__main__":
    test_banana()
    test_panamabananas_from_chapter_9()
    test_sais_matches_naive_on_random_strings()
    test_empty_sequence()
    test_divsufsort_backend_matches_sais_if_available()
    print("all suffix_array tests OK")
