"""
Poravnanje reada/fragmenta nad suffix arrayjem

Egzaktna podudaranja: binarna pretraga nad suffix arrayjem (PatternMatchingWithSuffixArray, Compeau & Pevzner pogl. 9) - sufiksi koji dijele prefiks
sortiraju se u jedan kontinuirani blok, pronađen u O(m log n) preko `_sa_lower_bound`/`_sa_upper_bound`

Približna podudaranja: seed-and-extend preko pigeonhole teorema iz istog poglavlja -- dijeljenje reada na d+1 seedova garantira da barem jedan ostane
netaknut s <= d razlika. Odabrano umjesto BWT/branching alternative jer smo izgradili samo suffix array, ne FM-index

Pokušavaju se oba lanca; pogoci na reverse-complementu prijavljuju se na poziciji forward lanca (AlignmentResult)
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Sequence, Set, Tuple

from genome_index.fasta_io import parse_fasta_chromosome
from genome_index.suffix_array import build_suffix_array

_VALID_BASES = frozenset("ACGTN")
_COMPLEMENT_TABLE = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class AlignmentResult:
    """
    position: 0-indeksirani offset u forward referenci (bez obzira na lanac)
    mismatches: broj neusklađenih baza između (eventualno RC-anog) reada i prozora
    strand: "+" ako se poklapa `read` kakav je zadan, "-" ako se poklapa njegov reverse complement
    """
    position: int
    mismatches: int
    strand: str


def reverse_complement(sequence: str) -> str:
    # Vrati reverse complement DNA stringa (A<->T, C<->G)
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


def _validate_read(read: str, sequence: str, max_mismatches: int) -> str:
    if not isinstance(read, str) or len(read) == 0:
        raise ValueError("read must be a non-empty string")

    read = read.upper()
    invalid_chars = set(read) - _VALID_BASES
    if invalid_chars:
        raise ValueError(
            f"read contains invalid characters {sorted(invalid_chars)}; "
            f"only A/C/G/T/N are allowed"
        )
    if len(read) > len(sequence):
        raise ValueError(
            f"read length ({len(read)}) exceeds sequence length "
            f"({len(sequence)}); a read cannot align to a shorter sequence"
        )
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be >= 0")
    if max_mismatches >= len(read):
        raise ValueError(
            f"max_mismatches ({max_mismatches}) must be smaller than the "
            f"read length ({len(read)})"
        )
    return read


#Egzaktno podudaranje: PatternMatchingWithSuffixArray (Compeau & Pevzner, pogl. 9)

def _suffix_prefix(sequence: str, suffix_array: Sequence[int], index: int, length: int) -> str:
    # Prvih `length` znakova sufiksa na suffix_array[index] (manje ako je taj sufiks kraći, tj. blizu kraja `sequence`)
    start = int(suffix_array[index])
    return sequence[start:start + length]


def _sa_lower_bound(sequence: str, suffix_array: Sequence[int], pattern: str) -> int:
    #Prvi indeks i takav da je sufiks od suffix_array[i] >= pattern
    lo, hi = 0, len(suffix_array)
    m = len(pattern)
    while lo < hi:
        mid = (lo + hi) // 2
        if _suffix_prefix(sequence, suffix_array, mid, m) < pattern:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _sa_upper_bound(sequence: str, suffix_array: Sequence[int], pattern: str) -> int:
    """Prvi indeks i takav da sufiks od suffix_array[i] više ne počinje s `pattern` (jedno mjesto iza zadnjeg podudarajućeg sufiksa)
    Skraćena jednakost se ovdje i dalje broji kao podudaranje, jer stvarni sufiks počinje s pattern bez obzira što slijedi"""
    lo, hi = 0, len(suffix_array)
    m = len(pattern)
    while lo < hi:
        mid = (lo + hi) // 2
        if _suffix_prefix(sequence, suffix_array, mid, m) <= pattern:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _sa_search_bounds(sequence: str, suffix_array: Sequence[int], pattern: str) -> Tuple[int, int]:
    """(first, last) uključivi raspon u suffix_array sufiksa koji počinju s pattern; first > last znači da nema podudaranja
    Izdvojeno kako bi pozivatelji mogli provjeriti broj podudaranja bez materijalizacije svih (vidi `_candidate_start_positions`)"""
    first = _sa_lower_bound(sequence, suffix_array, pattern)
    last = _sa_upper_bound(sequence, suffix_array, pattern) - 1
    return first, last


def find_exact_matches(sequence: str, suffix_array: Sequence[int], pattern: str) -> List[int]:
    """Sve početne pozicije `pattern` u `sequence`, putem binarne pretrage nad `suffix_array`
     O(len(pattern) * log(len(sequence))) plus O(k) za očitavanje k podudaranja"""
    if len(pattern) == 0:
        return []
    first, last = _sa_search_bounds(sequence, suffix_array, pattern)
    if first > last:
        return []
    return [int(p) for p in suffix_array[first:last + 1]]


# Približno podudaranje: seed-and-extend (pogl. 9, "Mismatch-Tolerant Read Mapping")

def _split_into_seeds(read: str, max_mismatches: int) -> List[Tuple[int, str]]:
    """Podijeli `read` na max_mismatches+1 seedova duljine k = len(read) // num_seeds, zadnji seed uzima ostatak
    Vraća parove (offset_in_read, seed_sequence)"""
    num_seeds = max_mismatches + 1
    k = len(read) // num_seeds
    seeds = [(i * k, read[i * k:(i + 1) * k]) for i in range(num_seeds - 1)]
    last_offset = (num_seeds - 1) * k
    seeds.append((last_offset, read[last_offset:]))
    return seeds


def _candidate_start_positions(
    sequence: str,
    suffix_array: Sequence[int],
    read: str,
    max_mismatches: int,
    max_seed_hits: int,
) -> Set[int]:
    """Detekcija seedova:
    Nađi egzaktna podudaranja svakog seeda, pomakni ih unatrag za njegov offset da dobiješ kandidatske početne pozicije za cijeli read

    Seed unutar dugog repeata ili N-runa može imati milijune pogodaka (seed od samih N-ova na chr21 je tijekom testiranja imao 6.6M pogodaka)
    preskoči svaki seed s više od `max_seed_hits` pogodaka umjesto da ga proširuješ, i pusti ostale seedove da lokaliziraju read
    Read kojemu je svaki seed prekomjerno zastupljen jednostavno se neće naći
    """
    candidates: Set[int] = set()
    for offset, seed in _split_into_seeds(read, max_mismatches):
        if not seed:
            continue
        first, last = _sa_search_bounds(sequence, suffix_array, seed)
        if first > last or (last - first + 1) > max_seed_hits:
            continue
        for hit in suffix_array[first:last + 1]:
            candidates.add(int(hit) - offset)
    return candidates


def _count_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _align_one_strand(
    sequence: str,
    suffix_array: Sequence[int],
    read: str,
    max_mismatches: int,
    strand: str,
    max_seed_hits: int,
) -> List[AlignmentResult]:
    #Proširenje seeda: provjeri svakog kandidata iz `_candidate_start_positions` naspram cijelog reada
    read_len = len(read)
    seq_len = len(sequence)
    results: List[AlignmentResult] = []
    candidates = _candidate_start_positions(sequence, suffix_array, read, max_mismatches, max_seed_hits)
    for start in candidates:
        if start < 0 or start + read_len > seq_len:
            continue  # ispada s jednog ili drugog kraja sekvence
        window = sequence[start:start + read_len]
        mismatches = _count_mismatches(read, window)
        if mismatches <= max_mismatches:
            results.append(AlignmentResult(start, mismatches, strand))
    return results


# API

def align_read(
    read: str,
    sequence: str,
    suffix_array: Sequence[int],
    max_mismatches: int = 2,
    max_seed_hits: int = 1000,
) -> List[AlignmentResult]:
    """Poravnaj jedan read naspram `sequence` koristeći njen suffix array

    Pokušava obje orijentacije: `read` kakav je zadan (strand="+") i njegov reverse complement (strand="-", pozicija se i dalje prijavljuje na
    forward lancu)
    Vraća svako poravnanje s najviše `max_mismatches` razlika, sortirano po poziciji.

    `max_seed_hits` ograničava koliko egzaktnih pogodaka seed smije imati prije nego što se preskoči kao neinformativan (repeati / N-runovi)
    (`_candidate_start_positions`)"""
    read = _validate_read(read, sequence, max_mismatches)

    forward_hits = _align_one_strand(
        sequence, suffix_array, read, max_mismatches, "+", max_seed_hits
    )
    reverse_hits = _align_one_strand(
        sequence, suffix_array, reverse_complement(read), max_mismatches, "-", max_seed_hits
    )

    results = forward_hits + reverse_hits
    results.sort(key=lambda r: (r.position, r.strand))
    return results


def align_reads(
    reads: List[str],
    sequence: str,
    suffix_array: Sequence[int],
    max_mismatches: int = 2,
    max_seed_hits: int = 1000,
) -> Dict[str, List[AlignmentResult]]:
    #Poravnaj skup readova, ključano po stringu reada. Duplicirani readovi se poravnaju samo jednom (ubrzanje za skupove s puno PCR duplikata)
    results: Dict[str, List[AlignmentResult]] = {}
    for read in reads:
        if read in results:
            continue
        results[read] = align_read(read, sequence, suffix_array, max_mismatches, max_seed_hits)
    return results


def load_reads(path: str) -> List[Tuple[str, str]]:
    """Učitaj readove iz FASTA ili FASTQ datoteke (format se automatski prepoznaje iz prve neprazne linije)
    Vraća parove (read_id, sequence)"""
    with open(path, "r") as handle:
        lines = [line.rstrip("\n\r") for line in handle]

    non_blank = next((line for line in lines if line.strip()), "")
    if non_blank.startswith(">"):
        return _parse_fasta_reads(lines)
    if non_blank.startswith("@"):
        return _parse_fastq_reads(lines)
    raise ValueError(
        f"Unrecognized reads file format in {path!r}: expected a FASTA "
        f"file (starting with '>') or a FASTQ file (starting with '@')"
    )


def _parse_fasta_reads(lines: List[str]) -> List[Tuple[str, str]]:
    reads: List[Tuple[str, str]] = []
    read_id = None
    chunks: List[str] = []
    for line in lines:
        if not line.strip():
            continue
        if line.startswith(">"):
            if read_id is not None:
                reads.append((read_id, "".join(chunks).upper()))
            read_id = line[1:].split(None, 1)[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if read_id is not None:
        reads.append((read_id, "".join(chunks).upper()))
    return reads


def _parse_fastq_reads(lines: List[str]) -> List[Tuple[str, str]]:
    reads: List[Tuple[str, str]] = []
    records = [line for line in lines if line.strip()]
    for i in range(0, len(records), 4):
        header, sequence = records[i], records[i + 1]
        if not header.startswith("@"):
            raise ValueError(
                f"Malformed FASTQ record at line {i}: expected a '@' header, "
                f"got {header!r}"
            )
        read_id = header[1:].split(None, 1)[0]
        reads.append((read_id, sequence.upper()))
    return reads


def align_reads_from_file(
    path: str,
    sequence: str,
    suffix_array: Sequence[int],
    max_mismatches: int = 2,
    max_seed_hits: int = 1000,
) -> Dict[str, List[AlignmentResult]]:
    #Učitaj readove iz FASTA/FASTQ datoteke i poravnaj ih, ključano po read ID-u (za razliku od `align_reads`, čuva duplikate s različitim ID-ovima)
    return {
        read_id: align_read(seq, sequence, suffix_array, max_mismatches, max_seed_hits)
        for read_id, seq in load_reads(path)
    }


# Primjer

def _sample_read_window(sequence: str, length: int, rng: "random.Random") -> Tuple[int, str]:
    """Odaberi nasumičan prozor duljine `length` bez "N" u sebi -- pravi readovi
    ne dolaze iz nerazriješenih assembly-gap regija, a hg38 ih ima
    (npr. kratki krak chr21) koje se protežu kroz milijune baza."""
    while True:
        start = rng.randint(0, len(sequence) - length)
        window = sequence[start:start + length]
        if "N" not in window:
            return start, window


if __name__ == "__main__":
    import random

    fasta_path = "hg38/hg38.u.fa"
    rng = random.Random(0)

    # demo na chrM: izreži readove, dodaj par točkastih mutacija, provjeri
    # poravnavaju li se natrag na mjesto odakle su uzeti.
    print("Building suffix array for chrM ...")
    chrm_sequence = parse_fasta_chromosome(fasta_path, "chrM")
    chrm_sa = build_suffix_array(fasta_path, "chrM", method="sais")

    read_length = 120
    demo_reads: Dict[str, Tuple[int, int]] = {}  # read -> (stvarna_pozicija, broj_uvedenih_mutacija)
    for _ in range(5):
        true_pos, window = _sample_read_window(chrm_sequence, read_length, rng)
        bases = list(window)
        num_mutations = rng.randint(0, 2)
        for idx in rng.sample(range(read_length), num_mutations):
            bases[idx] = rng.choice([b for b in "ACGT" if b != bases[idx]])
        demo_reads["".join(bases)] = (true_pos, num_mutations)

    print(f"\nAligning {len(demo_reads)} synthetic reads (chrM, max_mismatches=2):")
    for read, (true_pos, num_mutations) in demo_reads.items():
        alignments = align_read(read, chrm_sequence, chrm_sa, max_mismatches=2)
        best = alignments[0] if alignments else None
        status = "OK" if best and best.position == true_pos else "MISSED"
        print(
            f"  true_pos={true_pos:>6}  introduced_mutations={num_mutations}  "
            f"-> found={best}  [{status}]"
        )

    # Mjerenje performansi na cijelom kromosomu
    print("\nBuilding suffix array for chr21 (divsufsort backend) ...")
    chr21_sequence = parse_fasta_chromosome(fasta_path, "chr21")
    chr21_sa = build_suffix_array(fasta_path, "chr21", method="divsufsort")

    num_reads = 1000
    batch_reads = []
    for _ in range(num_reads):
        length = rng.randint(100, 150)
        _, window = _sample_read_window(chr21_sequence, length, rng)
        bases = list(window)
        for idx in rng.sample(range(length), rng.randint(0, 2)):
            bases[idx] = rng.choice([b for b in "ACGT" if b != bases[idx]])
        batch_reads.append("".join(bases))

    print(f"Aligning {num_reads} reads (100-150bp, max_mismatches=2) against chr21 "
          f"({len(chr21_sequence):,} bases) ...")
    t0 = time.time()
    batch_results = align_reads(batch_reads, chr21_sequence, chr21_sa, max_mismatches=2)
    elapsed = time.time() - t0

    aligned = sum(1 for hits in batch_results.values() if hits)
    print(f"Done in {elapsed:.2f} s total, {elapsed / num_reads * 1000:.2f} ms/read on average.")
    print(f"{aligned}/{len(batch_results)} distinct reads aligned at least once.")
