"""Spajanje svih kontiga referentnog genoma u jedan string radi izgradnje jednog suffix arraya nad cijelim genomom
Kontizi se razdvajaju N-spacerima kako se read ne bi mogao lažno poravnati preko granice dva kontiga
Izgradnja indeksa i poravnanje moraju koristiti isti raspored, inače vraćene pozicije ne odgovaraju
"""
from __future__ import annotations

import bisect
from typing import List, Optional, Tuple

SPACER = "N" * 200


def load_genome_concat(fasta_path: str) -> Tuple[List[str], List[int], List[int], str]:
    """Učitaj multi-FASTA i spoji sve kontige
    Vraća (imena, offseti, duljine, genom), gdje je offseti[i] 0-bazirani početak kontiga imena[i] unutar `genom`"""
    names: List[str] = []
    seqs: List[str] = []
    chunks: List[str] = []
    current_name = None

    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_name is not None:
                    seqs.append("".join(chunks).upper())
                current_name = line[1:].split(None, 1)[0]
                names.append(current_name)
                chunks = []
            else:
                chunks.append(line.rstrip())
    if current_name is not None:
        seqs.append("".join(chunks).upper())

    offsets: List[int] = []
    lengths: List[int] = []
    position = 0
    for seq in seqs:
        offsets.append(position)
        lengths.append(len(seq))
        position += len(seq) + len(SPACER)

    return names, offsets, lengths, SPACER.join(seqs)


def global_pos_to_contig(
    global_pos: int, names: List[str], offsets: List[int], lengths: List[int]
) -> Optional[Tuple[str, int]]:
    """Preslikaj 0-baziranu poziciju u spojenom genomu natrag u (kontig, 1-bazirana lokalna pozicija)
    Vraća None ako pozicija pada unutar N-spacera između kontiga"""
    i = bisect.bisect_right(offsets, global_pos) - 1
    if i < 0:
        return None
    local = global_pos - offsets[i]
    if local >= lengths[i]:
        return None
    return names[i], local + 1
