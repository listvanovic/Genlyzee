"""
Pozivanje mutacija (SNP-ova) iz poravnatog reada

Compeau & Pevzner pogl. 9 tvrdi da mismatchevi u poravnatim readovima ukazuju na SNP-ove, ali stanu na *pronalaženju* podudaranja
Ne razdvajaju pravi SNP od greške sekvenciranja
Ideja je iz "bubbles" koncepta assembly grafova iz pogl. 3 (greške se pojavljuju uz nisku pokrivenost, prava varijacija je dobro poduprta):
izgradi pileup po poziciji i pozovi mutaciju samo tamo gdje je alt baza reproducibilno podržana od strane više readova

`call_mutations` zahtijeva coverage >= min_coverage I udio većinske alt baze >= min_variant_fraction
Oboje bitno: sam udio bi pozvao svaku grešku s jednim readom (njen vlastiti VAF je trivijalno 1.0)
Sam coverage ne isključuje da se nekoliko različitih grešaka nagomila slučajno
Zadane vrijednosti (3, 0.5) su ilustrativne za readove demo-skale, nisu pragovi

Rukovanje lancem: `position` u AlignmentResultu na "-" lancu je mjesto gdje se reverse_complement(read), a ne sam `read`, poklapa s forward
referencom (alignment.py)
`_mapped_bases` odmotava komplement počevši od `position`

Multi-mapping readovi (dvosmislenog podrijetla, repeat) isključeni su iz pileupa umjesto da se ubace na svaku kandidatsku poziciju - ideja kao s
 `max_seed_hits` u alignment.py (žrtvuje se recall za preciznost)
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from genome_index.alignment import AlignmentResult, reverse_complement

_VALID_BASES = frozenset("ACGTN")


@dataclass(frozen=True)
class MutationCall:
    """Jedna pozvana mutacija (točkasta supstitucija) u odnosu na referentnu sekvencu

    position: 0-indeksirani offset, isti koordinatni prostor kao AlignmentResult.position
    ref_base / alt_base: referentna baza na `position` / dogovorena alternativa
    supporting_reads: readovi u pileupu na `position` koji prijavljuju `alt_base`
    coverage: ukupno readova u pileupu na `position` (bilo koja baza)
    variant_allele_fraction: supporting_reads / coverage
    """
    position: int
    ref_base: str
    alt_base: str
    supporting_reads: int
    coverage: int
    variant_allele_fraction: float


def _validate_inputs(
    reference_sequence: str,
    reads: Sequence[str],
    alignment_results: Dict[str, List[AlignmentResult]],
    min_coverage: int,
    min_variant_fraction: float,
) -> None:
    # Provjeri ulaze prije izgradnje pileupa
    if not isinstance(reference_sequence, str) or len(reference_sequence) == 0:
        raise ValueError("reference_sequence must be a non-empty string")

    invalid_chars = set(reference_sequence) - _VALID_BASES
    if invalid_chars:
        raise ValueError(
            f"reference_sequence contains invalid characters "
            f"{sorted(invalid_chars)}; only A/C/G/T/N are allowed"
        )
    if min_coverage < 1:
        raise ValueError(f"min_coverage must be >= 1, got {min_coverage}")
    if not (0.0 < min_variant_fraction <= 1.0):
        raise ValueError(
            f"min_variant_fraction must be in (0.0, 1.0], got {min_variant_fraction}"
        )
    missing = [read for read in reads if read not in alignment_results]
    if missing:
        raise ValueError(
            f"{len(missing)} read(s) in `reads` have no entry in "
            f"`alignment_results` (e.g. {missing[0]!r}); `alignment_results` "
            f"must be the output of aligning exactly this read set (e.g. via "
            f"`align_reads(reads, ...)`)"
        )


def _mapped_bases(
    read: str, alignment: AlignmentResult, reference_length: int
) -> List[Tuple[int, str]]:
    """Mapiraj baze reada na koordinate forward-lanca reference
    'N' baze se izbacuju - nisu informativan poziv
    """
    read = read.upper()
    mapped_read = reverse_complement(read) if alignment.strand == "-" else read
    start = alignment.position
    end = start + len(mapped_read)
    if start < 0 or end > reference_length:
        raise ValueError(
            f"alignment position {start} (length {len(mapped_read)}) falls "
            f"outside the reference sequence (length {reference_length}); "
            f"`alignment_results` must have been computed against this same "
            f"`reference_sequence`"
        )
    return [
        (start + offset, base)
        for offset, base in enumerate(mapped_read)
        if base != "N"
    ]


def _build_pileup(
    reference_sequence: str,
    reads: Sequence[str],
    alignment_results: Dict[str, List[AlignmentResult]],
) -> Dict[int, List[str]]:
    #Izgradi pileup - za svaku poziciju u referenci, baze koje tamo prijavljuje svaki jedinstveno poravnati read
    pileup: Dict[int, List[str]] = defaultdict(list)
    for read in reads:
        alignments = alignment_results[read]
        if len(alignments) != 1:
            continue  # 0 = nepoklopljen, >1 = multi-mapping; oba se isključuju
        for position, base in _mapped_bases(read, alignments[0], len(reference_sequence)):
            pileup[position].append(base)
    return pileup


def call_mutations(
    reference_sequence: str,
    reads: List[str],
    alignment_results: Dict[str, List[AlignmentResult]],
    min_coverage: int = 3,
    min_variant_fraction: float = 0.5,
) -> List[MutationCall]:
    """Pozovi točkaste mutacije (SNP-ove) iz `reads` i `alignment_results`

    `reads` može sadržavati duplikate (svaki se posebno broji u coverage)
    `alignment_results` je `{read: [AlignmentResult, ...]}` npr. iz `alignment.align_reads`, sa zapisom za svaki read
    """
    _validate_inputs(reference_sequence, reads, alignment_results, min_coverage, min_variant_fraction)

    pileup = _build_pileup(reference_sequence, reads, alignment_results)

    calls: List[MutationCall] = []
    for position in sorted(pileup):
        bases = pileup[position]
        coverage = len(bases)
        if coverage < min_coverage:
            continue

        ref_base = reference_sequence[position]
        if ref_base == "N":
            continue  # nerazriješena referentna baza: nema s čim usporediti

        variant_bases = [base for base in bases if base != ref_base]
        if not variant_bases:
            continue

        alt_base, supporting_reads = Counter(variant_bases).most_common(1)[0]
        variant_allele_fraction = supporting_reads / coverage
        if variant_allele_fraction >= min_variant_fraction:
            calls.append(
                MutationCall(
                    position=position,
                    ref_base=ref_base,
                    alt_base=alt_base,
                    supporting_reads=supporting_reads,
                    coverage=coverage,
                    variant_allele_fraction=variant_allele_fraction,
                )
            )
    return calls


# Primjer - suffix array -> poravnanje -> pozivanje mutacija na pravom isječku chrM, sa zasađenim mutacijama i greškama sekvenciranja na 
# pojedinačnim readovima

if __name__ == "__main__":
    import random

    from genome_index.alignment import align_reads
    from genome_index.fasta_io import parse_fasta_chromosome
    from genome_index.suffix_array import build_suffix_array_from_sequence

    fasta_path = "hg38/hg38.u.fa"
    rng = random.Random(0)

    print("Loading chrM and slicing out a 4,000 bp reference window ...")
    chrm_sequence = parse_fasta_chromosome(fasta_path, "chrM")
    window_start = 3000
    reference_sequence = chrm_sequence[window_start:window_start + 4000]

    # Zasadi dvije prave mutacije - readovi koji pokrivaju te pozicije dobiju istu alt bazu
    true_mutation_positions = [800, 2500]
    true_mutations = {}
    for pos in true_mutation_positions:
        ref_base = reference_sequence[pos]
        alt_base = next(b for b in "ACGT" if b != ref_base)
        true_mutations[pos] = (ref_base, alt_base)

    read_length = 100
    num_reads = 200
    reads: List[str] = []
    for _ in range(num_reads):
        start = rng.randint(0, len(reference_sequence) - read_length)
        bases = list(reference_sequence[start:start + read_length])

        for pos, (ref_base, alt_base) in true_mutations.items():
            if start <= pos < start + read_length:
                bases[pos - start] = alt_base

        reads.append("".join(bases))

    # min grešaka na pojedinačnim readovima- ne smiju se pozvati (nemaju podršku više readova)
    single_read_error_indices = rng.sample(range(num_reads), 5)
    for i in single_read_error_indices:
        bases = list(reads[i])
        error_offset = rng.randrange(read_length)
        original = bases[error_offset]
        bases[error_offset] = rng.choice([b for b in "ACGT" if b != original])
        reads[i] = "".join(bases)

    print(f"Building suffix array for the {len(reference_sequence)} bp reference window ...")
    suffix_array = build_suffix_array_from_sequence(reference_sequence, method="sais")

    print(f"Aligning {len(reads)} reads (length {read_length}, max_mismatches=2) ...")
    alignment_results = align_reads(reads, reference_sequence, suffix_array, max_mismatches=2)

    print("Calling mutations (min_coverage=3, min_variant_fraction=0.5) ...\n")
    calls = call_mutations(reference_sequence, reads, alignment_results)

    called_positions = {call.position for call in calls}
    print(f"{len(calls)} mutation(s) called:")
    for call in calls:
        print(
            f"  position={call.position:>5}  {call.ref_base}->{call.alt_base}  "
            f"supporting_reads={call.supporting_reads:>3}/{call.coverage:<3}  "
            f"VAF={call.variant_allele_fraction:.2f}"
        )

    print("\nValidation against ground truth:")
    for pos, (ref_base, alt_base) in true_mutations.items():
        status = "OK, called" if pos in called_positions else "MISSED"
        print(f"  planted true mutation at {pos} ({ref_base}->{alt_base}): [{status}]")

    spurious = called_positions - set(true_mutations)
    status = "OK, none called" if not spurious else f"UNEXPECTED calls: {spurious}"
    print(f"  single-read sequencing errors correctly filtered out: [{status}]")
