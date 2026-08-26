"""
Izgradnja suffix arraya
(Compeau & Pevzner pogl. 9)
SUFFIXARRAY(Text): sortirane početne pozicije svih sufiksa od Text ("$" je najmanji), npr. SUFFIXARRAY("panamabananas$") = [13,5,3,1,7,9,11,6,4,2,8,10,0,12]
`build_suffix_array_naive` sortira sufikse izravno O(n^2 log n) - koristi se samo kao oracle za provjeru ispravnosti na kratkim testovima

SA-IS (Nong, Zhang & Chen, DCC 2009), O(n) - odabran umjesto DC3/skew jer treba samo jedan rekurzivni korak (sortiraj LMS sufikse indukcijom, 
imenuj ih, rekurziraj na sudarima) umjesto DC3-ovog "split-then-merge" pristupa, a to je i ono što interno implementira sama libdivsufsort

Dva backenda: čisti Python "sais" je O(n), ali spor u praksi (~0.6-1 Mbp/s, satima za cijeli kromosom); "divsufsort" delegira na
pydivsufsort/libdivsufsort (~25-40 Mbp/s, sekunde po kromosomu)
"sais" se drži kao referenca koja se objašnjava u radu; "divsufsort" se koristi na skali stvarnog kromosoma
"""
from __future__ import annotations

import sys
import time
from typing import List, Sequence, Tuple

import numpy as np

from genome_index.fasta_io import parse_fasta_chromosome

SENTINEL = "$"


# Naivna izgradnja - referentna

def build_suffix_array_naive(text: str) -> List[int]:
    """Izravno sortiraj sve sufikse od `text` i vrati njihove početne pozicije - O(n^2 log n) u najgorem slučaju
    Oracle za provjeru ispravnosti samo na kratkim test stringovima, nikad na podacima stvarne kromosomske skale."""
    n = len(text)
    return sorted(range(n), key=lambda i: text[i:])


# SA-IS
# izgradnja u linearnom vremenu putem induciranog sortiranja

def _classify_suffix_types(s: Sequence[int]) -> List[bool]:
    """
    Klasificiraj svaki sufiks s[i:] kao S-tip (True) ili L-tip (False):

        sufiks i je S-tip  <==>  s[i] < s[i+1], ili
                                  (s[i] == s[i+1] i sufiks i+1 je S-tip)
        sufiks i je L-tip  <==>  s[i] > s[i+1], ili
                                  (s[i] == s[i+1] i sufiks i+1 je L-tip)

    sam sentinel je po definiciji S-tip
    Jedan prolaz zdesna nalijevo, jer tip(i) ovisi samo o tipu(i+1)
    """
    n = len(s)
    t = [False] * n
    t[n - 1] = True
    for i in range(n - 2, -1, -1):
        if s[i] < s[i + 1]:
            t[i] = True
        elif s[i] > s[i + 1]:
            t[i] = False
        else:
            t[i] = t[i + 1]
    return t


def _sa_is(s: List[int], alphabet_size: int) -> List[int]:
    """
    Izgradi suffix array cijelobrojno-kodiranog stringa `s` u O(n) vremenu koristeći SA-IS (Nong, Zhang & Chen, 2009)

    Preduvjet: s[-1] == 0 i 0 se ne pojavljuje nigdje drugdje u `s` (jedinstven, strogo najmanji sentinel na kraju)

    Struktura na visokoj razini:
      1. Klasificiraj svaku poziciju kao S-tip ili L-tip (`_classify_suffix_types`)
      2. LMS pozicije: i > 0 je "leftmost S-type" ako je S-tip, a i-1 je L-tip
         Uzastopne LMS pozicije omeđuju "LMS substringove".
      3. Nabaci LMS sufikse na kraj njihovih kanti u proizvoljnom redoslijedu, zatim jedan L-indukcijski prolaz + jedan S-indukcijski prolaz
         Ključni rezultat SA-IS-a:
         ovaj "round-trip" već ispravno sortira LMS *substringove* jedne prema drugima, iako puni sufiksi još nisu na konačnim mjestima
      4. Imenuj svaki LMS substring (rangiraj po tom poretku, uspoređujući susjedne radi jednakosti)
      5. Ako su sva imena različita, poredak LMS sufiksa je poznat
         Inače rekurziraj: izgradi reducirani string s1 (imena u izvornom redoslijedu slijeva nadesno, najviše pola duljine s) i rekurziraj na njemu
      6. Nabaci sad ispravno poredane LMS sufikse na kraj kanti i pokreni oba indukcijska prolaza još jednom
         ovo postavlja sve ostalo i daje konačni suffix array
    """
    n = len(s)
    if n == 1:
        return [0]
    if n == 2:
        return [1, 0]

    t = _classify_suffix_types(s)

    def is_lms(i: int) -> bool:
        return i > 0 and t[i] and not t[i - 1]

    bucket_sizes = [0] * alphabet_size
    for c in s:
        bucket_sizes[c] += 1

    def bucket_heads() -> List[int]:
        heads = [0] * alphabet_size
        total = 0
        for c in range(alphabet_size):
            heads[c] = total
            total += bucket_sizes[c]
        return heads

    def bucket_tails() -> List[int]:
        tails = [0] * alphabet_size
        total = 0
        for c in range(alphabet_size):
            total += bucket_sizes[c]
            tails[c] = total
        return tails

    def induce_sort_l(sa: List[int]) -> None:
        """Skenirajući sa slijeva nadesno
        postavi svakog L-tip prethodnika j-1 već postavljenog sufiksa j na početak kante od j-1"""
        heads = bucket_heads()
        for i in range(n):
            j = sa[i]
            if j <= 0:
                continue
            if not t[j - 1]:
                c = s[j - 1]
                sa[heads[c]] = j - 1
                heads[c] += 1

    def induce_sort_s(sa: List[int]) -> None:
        """Skenirajući sa zdesna nalijevo
        postavi svakog S-tip prethodnika j-1 već postavljenog sufiksa j na kraj kante od j-1"""
        tails = bucket_tails()
        for i in range(n - 1, -1, -1):
            j = sa[i]
            if j <= 0:
                continue
            if t[j - 1]:
                c = s[j - 1]
                tails[c] -= 1
                sa[tails[c]] = j - 1


    sa: List[int] = [-1] * n
    tails = bucket_tails()
    for i in range(n - 1, -1, -1):
        if is_lms(i):
            c = s[i]
            tails[c] -= 1
            sa[tails[c]] = i

    induce_sort_l(sa)
    induce_sort_s(sa)

    # sortirane LMS pozicije na početak sa, imenuj svaku
    n1 = 0
    for i in range(n):
        if sa[i] > 0 and is_lms(sa[i]):
            sa[n1] = sa[i]
            n1 += 1

    for i in range(n1, n):
        sa[i] = -1

    name = -1
    prev_pos = -1
    for i in range(n1):
        pos = sa[i]
        different = prev_pos == -1
        if not different:
            d = 0
            while True:
                if (prev_pos + d >= n or pos + d >= n
                        or s[prev_pos + d] != s[pos + d]
                        or t[prev_pos + d] != t[pos + d]):
                    different = True
                    break
                if d > 0 and (is_lms(prev_pos + d) or is_lms(pos + d)):
                    break  # oba substringa su završila zajedno: jednaki
                d += 1
        if different:
            name += 1
        prev_pos = pos
        # LMS pozicije su udaljene barem 2, pa je pos // 2 jedinstveno, kompaktno mjesto za spremiti ime ovog LMS sufiksa
        sa[n1 + pos // 2] = name

    s1 = [x for x in sa[n1:] if x != -1]


    if name + 1 == n1:
        sa1 = [0] * n1
        for i, rank in enumerate(s1):
            sa1[rank] = i
    else:
        sa1 = _sa_is(s1, name + 1)

    lms_positions = [i for i in range(1, n) if is_lms(i)]


    sa = [-1] * n
    tails = bucket_tails()
    for i in range(n1 - 1, -1, -1):
        pos = lms_positions[sa1[i]]
        c = s[pos]
        tails[c] -= 1
        sa[tails[c]] = pos

    induce_sort_l(sa)
    induce_sort_s(sa)

    return sa


# API

def _encode_with_sentinel(text: str) -> Tuple[List[int], int]:
    """Mapiraj svaki različit znak u `text` na mali cijeli broj (1, 2, ... po ASCII redu -- za SA-IS radi bilo koji potpuni uređaj), rezervirajući 0 za
    dodani sentinel
    Vraća kodiranu sekvencu i veličinu abecede"""
    alphabet = sorted(set(text))
    if SENTINEL in alphabet:
        raise ValueError(
            f"Input text must not already contain the sentinel {SENTINEL!r}."
        )
    code_of = {ch: i + 1 for i, ch in enumerate(alphabet)}
    encoded = [code_of[ch] for ch in text]
    encoded.append(0)  # sentinel, strogo manji od svakog pravog simbola
    return encoded, len(alphabet) + 1

# Bira int32 ili int64 ovisno o duljini sekvence, radi uštede memorije
def _suffix_array_dtype(n: int) -> type:
    return np.int32 if n < np.iinfo(np.int32).max else np.int64


def build_suffix_array_from_sequence(
    sequence: str, method: str = "sais"
) -> np.ndarray:
    """
    Izgradi suffix array za string `sequence` koji je već u memoriji
    Izdvojeno iz `build_suffix_array` kako bi ga testovi mogli pozvati nad proizvoljnim stringovima, ne samo onima iz FASTA datoteka

    `sequence` ne smije sadržavati sentinel "$"
    `method`: "sais" (zadano, čisti Python, u redu do nekoliko milijuna baza) ili "divsufsort" (pydivsufsort/libdivsufsort, za stvarnu kromosomsku skalu)

    Vraća numpy array (int32, ili int64 iznad 2**31) umjesto obične Python liste -- na skali hg38 obična lista intova troši nekoliko puta više memorije
    Pozovi `.tolist()` za običnu listu
    """
    n = len(sequence)
    dtype = _suffix_array_dtype(n + 1)

    if method == "sais":
        encoded, alphabet_size = _encode_with_sentinel(sequence)
        sa = _sa_is(encoded, alphabet_size)
        # Pozicija n (dodani sentinel) uvijek se sortira prva i nije dio sufiksa izvorne sekvence; izbaci je
        sa = [p for p in sa if p != n]
        return np.array(sa, dtype=dtype)

    if method == "divsufsort":
        if SENTINEL in sequence:
            raise ValueError(
                f"Input text must not already contain the sentinel {SENTINEL!r}."
            )
        try:
            import pydivsufsort
        except ImportError as exc:
            raise ImportError(
                "method='divsufsort' requires the optional 'pydivsufsort' "
                "package. Install it with: pip install pydivsufsort"
            ) from exc
        sa = pydivsufsort.divsufsort(sequence.encode("ascii"))
        return sa.astype(dtype, copy=False)

    raise ValueError(
        f"Unknown method {method!r}; expected 'sais' or 'divsufsort'."
    )


def build_suffix_array(
    fasta_path: str, chromosome_name: str, method: str = "sais"
) -> np.ndarray:
    """Izvuci `chromosome_name` iz FASTA datoteke na `fasta_path` i izgradi njegov suffix array (detalji u `parse_fasta_chromosome` i
    `build_suffix_array_from_sequence`; FileNotFoundError/ValueError se propagiraju iz prve)"""
    sequence = parse_fasta_chromosome(fasta_path, chromosome_name)
    return build_suffix_array_from_sequence(sequence, method=method)


# Primjer

if __name__ == "__main__":
    fasta_path = sys.argv[1] if len(sys.argv) > 1 else "hg38/hg38.u.fa"
    chromosome = sys.argv[2] if len(sys.argv) > 2 else "chr21"
    method = sys.argv[3] if len(sys.argv) > 3 else "divsufsort"

    print(f"Building suffix array for {chromosome!r} from {fasta_path!r} "
          f"using method={method!r} ...")

    t0 = time.time()
    suffix_array = build_suffix_array(fasta_path, chromosome, method=method)
    elapsed = time.time() - t0

    print(f"Done in {elapsed:.2f} s.")
    print(f"Sequence length: {suffix_array.shape[0]:,} bases")
    print(f"First 20 suffix array entries: {suffix_array[:20].tolist()}")
