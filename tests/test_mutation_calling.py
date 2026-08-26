""" Jedinični testovi za genome_index.mutation_calling

`AlignmentResult`-i se grade izravno tako da svaki test izolira jedan dio pileup/calling logike
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome_index.alignment import AlignmentResult, reverse_complement
from genome_index.mutation_calling import call_mutations

REFERENCE = "ACGTACGTACGTACGTACGTACGTACGTAC"  # 31 baza


def test_mutation_called_with_sufficient_coverage_and_variant_fraction() -> None:
    """Pozicija na kojoj se svaki od nekoliko neovisnih readova slaže oko iste alternativne baze mora biti pozvana, s ispravnim ref/alt bazama,
    brojem podržavajućih readova, coverageom i variant allele fraction"""
    start = 5
    true_pos = 10  # unutar [start, start+10), offset 5 u readu
    window = list(REFERENCE[start:start + 10])
    alt_base = next(b for b in "ACGT" if b != window[true_pos - start])
    window[true_pos - start] = alt_base
    mutated_read = "".join(window)

    # 5 neovisnih fizičkih readova koji svi slučajno nose ovu mutaciju
    # (identična sekvenca, pa se sažmu u jedan ključ u alignment_results
    # ali se svaki i dalje mora zasebno brojati u coverage)
    reads = [mutated_read] * 5
    alignment_results = {mutated_read: [AlignmentResult(start, 1, "+")]}

    calls = call_mutations(REFERENCE, reads, alignment_results)

    assert len(calls) == 1
    call = calls[0]
    assert call.position == true_pos
    assert call.ref_base == REFERENCE[true_pos]
    assert call.alt_base == alt_base
    assert call.supporting_reads == 5
    assert call.coverage == 5
    assert call.variant_allele_fraction == 1.0


def test_single_read_error_not_called() -> None:
    """ Pozicija s dovoljnim ukupnim coverageom, ali gdje se samo jedan read ne slaže s referencom (moguća greška sekvenciranja), a ostali se slažu
    s njom, NE smije biti pozvana: udio varijante za tu usamljenu alternativnu bazu ostaje daleko ispod zadanog praga 0.5.
    """
    true_pos = 10

    # 3 reda, svaki drukčiji prozor od 10 baza koji i dalje pokriva true_pos,
    # svi se točno poklapaju s referencom tamo (0 razlika).
    correct_starts = [2, 5, 8]
    correct_reads = [REFERENCE[s:s + 10] for s in correct_starts]
    alignment_results = {
        read: [AlignmentResult(s, 0, "+")] for s, read in zip(correct_starts, correct_reads)
    }

    # 4. read, iz još jednog prozora, sa zasađenom greškom jedne baze
    # točno na true_pos.
    error_start = 3
    offset = true_pos - error_start
    error_bases = list(REFERENCE[error_start:error_start + 10])
    error_bases[offset] = next(b for b in "ACGT" if b != error_bases[offset])
    erroneous_read = "".join(error_bases)
    alignment_results[erroneous_read] = [AlignmentResult(error_start, 1, "+")]

    # coverage na true_pos = 4 (3 ispravna + 1 pogrešan), ali samo
    # 1/4 = 0.25 < 0.5 podržava alternativnu bazu greške.
    reads = correct_reads + [erroneous_read]

    calls = call_mutations(REFERENCE, reads, alignment_results)

    assert calls == []


def test_low_coverage_not_called_even_with_full_agreement() -> None:
    """Pozicija na kojoj se svi readovi slažu oko alternativne baze, ali je ukupan broj readova koji je pokrivaju ispod min_coverage, NE smije
    biti pozvana - samo slaganje nije dovoljan dokaz uz vrlo nisku dubinu"""
    start = 0
    true_pos = 4
    window = list(REFERENCE[start:start + 8])
    alt_base = next(b for b in "ACGT" if b != window[true_pos - start])
    window[true_pos - start] = alt_base
    mutated_read = "".join(window)

    # Samo 2 reda, oba se slažu oko alt_base (VAF=1.0), ali je zadani
    # min_coverage 3
    reads = [mutated_read, mutated_read]
    alignment_results = {mutated_read: [AlignmentResult(start, 1, "+")]}

    calls = call_mutations(REFERENCE, reads, alignment_results, min_coverage=3)

    assert calls == []


def test_reverse_strand_reads_mapped_back_to_forward_reference_base() -> None:
    """ Readovima poravnatim na '-' lancu baze se moraju ispravno mapirati natrag na forward-lanac referencu: sam string reada je reverse
    complement forward-lanac prozora s kojim se poklapa, pa ga pileup mora "odmotati" (un-reverse-complement) prije usporedbe s
    `reference_sequence`, i za poziciju i za identitet baze
    """
    start = 12
    true_pos = 18
    offset = true_pos - start
    forward_window = list(REFERENCE[start:start + 10])
    ref_base = forward_window[offset]
    alt_base = next(b for b in "ACGT" if b != ref_base)
    forward_window[offset] = alt_base

    # Fizički read, kakav je sekvenciran s suprotnog lanca, je reverse
    # complement (mutiranog) forward-lanac prozora
    read = reverse_complement("".join(forward_window))

    reads = [read] * 4
    alignment_results = {read: [AlignmentResult(start, 1, "-")]}

    calls = call_mutations(REFERENCE, reads, alignment_results)

    assert len(calls) == 1
    call = calls[0]
    assert call.position == true_pos
    assert call.ref_base == ref_base
    assert call.alt_base == alt_base  # identitet forward-lanca, ne njegov komplement
    assert call.supporting_reads == 4
    assert call.coverage == 4


if __name__ == "__main__":
    for test in (
        test_mutation_called_with_sufficient_coverage_and_variant_fraction,
        test_single_read_error_not_called,
        test_low_coverage_not_called_even_with_full_agreement,
        test_reverse_strand_reads_mapped_back_to_forward_reference_base,
    ):
        test()
        print(test.__name__, "OK")
    print("done.")
