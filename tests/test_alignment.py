# Jedinični testovi za genome_index.alignment
from __future__ import annotations

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome_index.alignment import align_read, align_reads, reverse_complement
from genome_index.suffix_array import build_suffix_array_from_sequence

SEQUENCE = (
    "ACGTTGCATGCATGCACGTAGCTAGCTGATCGATCGTAGCTAGCATCGATCGATGCATGC"
    "TAGCTAGCTAGCATCGATGCATCGTAGCTAGCATCGATCGATCGTAGCATCGATCGATCG"
)


def _sa(sequence: str = SEQUENCE):
    return build_suffix_array_from_sequence(sequence, method="sais")


def test_exact_match() -> None:
    """Read koji je egzaktan podstring sekvence mora se poravnati na ispravnoj poziciji s nula razlika, na '+' lancu"""
    true_pos = 10
    read = SEQUENCE[true_pos:true_pos + 25]
    sa = _sa()
    results = align_read(read, SEQUENCE, sa, max_mismatches=2)

    assert len(results) >= 1
    exact_hits = [r for r in results if r.position == true_pos]
    assert len(exact_hits) == 1
    assert exact_hits[0].mismatches == 0
    assert exact_hits[0].strand == "+"


def test_read_with_introduced_mismatches_still_aligns() -> None:
    """Read izrezan iz sekvence s uvedenih 1-2 točkastih mutacija i dalje se mora poravnati na svoju pravu poziciju, s ispravnim brojem razlika,
    sve dok ga max_mismatches pokriva"""
    true_pos = 40
    read_len = 30
    bases = list(SEQUENCE[true_pos:true_pos + read_len])

    for num_mutations in (1, 2):
        mutated = bases.copy()
        mutation_sites = random.Random(num_mutations).sample(range(read_len), num_mutations)
        for site in mutation_sites:
            original = mutated[site]
            mutated[site] = next(b for b in "ACGT" if b != original)
        read = "".join(mutated)

        sa = _sa()
        results = align_read(read, SEQUENCE, sa, max_mismatches=2)
        hits_at_true_pos = [r for r in results if r.position == true_pos and r.strand == "+"]

        assert len(hits_at_true_pos) == 1, (num_mutations, results)
        assert hits_at_true_pos[0].mismatches == num_mutations


def test_nonexistent_read_returns_no_alignment() -> None:
    """Read koji se nigdje ne pojavljuje u sekvenci (čak i uz dopušteni max_mismatches) mora dati praznu listu rezultata"""
    # SEQUENCE nema niz od 30+ 'G' slova, pa niz od 30 'G' ne može
    # odgovarati nijednom prozoru unutar max_mismatches=2.
    read = "G" * 30
    sa = _sa()
    results = align_read(read, SEQUENCE, sa, max_mismatches=2)
    assert results == []


def test_reverse_complement_alignment() -> None:
    """Read jednak reverse complementu podstringa sekvence mora se poravnati na forward-lanac početnu poziciju tog podstringa, označen
    s strand='-'"""
    true_pos = 15
    read_len = 20
    forward_substring = SEQUENCE[true_pos:true_pos + read_len]
    read = reverse_complement(forward_substring)

    sa = _sa()
    results = align_read(read, SEQUENCE, sa, max_mismatches=2)

    reverse_hits = [r for r in results if r.strand == "-"]
    assert len(reverse_hits) == 1
    assert reverse_hits[0].position == true_pos
    assert reverse_hits[0].mismatches == 0


def test_align_reads_batch_and_dedup() -> None:
    """align_reads treba ključati rezultate po stringu reada i računati svaki različit read samo jednom, čak i ako se pojavljuje više puta u ulazu"""
    true_pos = 5
    read = SEQUENCE[true_pos:true_pos + 15]
    sa = _sa()

    batch_results = align_reads([read, read, "G" * 30],
                                 SEQUENCE, sa, max_mismatches=1)

    assert read in batch_results
    assert len(batch_results) == 2  # duplicirani read sažet u jedan ključ
    assert any(r.position == true_pos for r in batch_results[read])


def test_input_validation() -> None:
    sa = _sa()
    import pytest

    with pytest.raises(ValueError):
        align_read("", SEQUENCE, sa, max_mismatches=1)
    with pytest.raises(ValueError):
        align_read("ACGTN" * 100, SEQUENCE, sa, max_mismatches=1)  # duži od SEQUENCE
    with pytest.raises(ValueError):
        align_read("ACGXT", SEQUENCE, sa, max_mismatches=1)  # nevažeći znak
    with pytest.raises(ValueError):
        align_read("ACG", SEQUENCE, sa, max_mismatches=5)  # max_mismatches >= len(read)


if __name__ == "__main__":
    tests = [
        test_exact_match,
        test_read_with_introduced_mismatches_still_aligns,
        test_nonexistent_read_returns_no_alignment,
        test_reverse_complement_alignment,
        test_align_reads_batch_and_dedup,
    ]
    for test in tests:
        test()
        print(f"{test.__name__}: OK")

    # test_input_validation koristi pytest.raises; pokreni ga samo ako je
    # pytest importabilan, inače ga preskoči u ovom samostalnom pokretaču.
    try:
        import pytest  # noqa: F401
        test_input_validation()
        print("test_input_validation: OK")
    except ImportError:
        print("test_input_validation: SKIPPED (pytest not installed)")

    print("All tests passed.")
