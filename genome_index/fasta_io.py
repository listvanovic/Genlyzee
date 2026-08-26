"""
Pomoćne funkcije za parsiranje FASTA datoteka
namjerno odvojene od suffix_array.py kako bi ih faze poravnanja i pozivanja mutacija mogle koristiti bez izmjena
"""
from __future__ import annotations

import os
from typing import List


def parse_fasta_chromosome(fasta_path: str, chromosome_name: str) -> str:
    """
    Izvuci jedan kromosom iz (potencijalno višegigabajtne) multi-FASTA datoteke (poput hg38), uspoređujući `chromosome_name` bez obzira na 
    velika/mala slova s prvim tokenom svakog zaglavlja
    Čita liniju po liniju i staje čim se pojavi sljedeće zaglavlje pa kromosom blizu početka datoteke ne zahtijeva čitanje ostatka

    Vraća sekvencu velikim slovima, bez newlineova
    Baca FileNotFoundError ako `fasta_path` ne postoji, ValueError ako nijedno zaglavlje ne odgovara
    """
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path!r}")

    target = chromosome_name.strip().lower()
    sequence_chunks: List[str] = []
    found = False
    reading_target = False

    with open(fasta_path, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                if reading_target:
                    break  # ciljani kromosom je potpuno prikupljen - ne čitaj ostatak
                header_name = line[1:].split(None, 1)[0].strip().lower()
                reading_target = header_name == target
                if reading_target:
                    found = True
                continue

            if reading_target:
                sequence_chunks.append(line.rstrip())  # rstrip uklanja i '\r'

    if not found:
        raise ValueError(
            f"Chromosome {chromosome_name!r} not found in FASTA file "
            f"{fasta_path!r}."
        )

    return "".join(sequence_chunks).upper()
