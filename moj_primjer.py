# Primjer upotrebe rada: od FASTA datoteke do pozvanih mutacija

import random

from genome_index.fasta_io import parse_fasta_chromosome
from genome_index.suffix_array import build_suffix_array_from_sequence
from genome_index.alignment import align_read, align_reads
from genome_index.mutation_calling import call_mutations

FASTA_PATH = "hg38/hg38.u.fa"

# Faza 1: učitaj kromosom i izgradi njegov indeks
print("Učitavam chrM iz FASTA datoteke ...")
sequence = parse_fasta_chromosome(FASTA_PATH, "chrM")
print(f"Duljina chrM: {len(sequence):,} slova")

print("Gradim suffix array ...")
sa = build_suffix_array_from_sequence(sequence, method="sais")
print("Gotovo.\n")

# Faza 2: uzmi pravi komad genoma, ubaci jednu "mutaciju", poravnaj ga natrag
true_position = 5000
original_piece = sequence[true_position:true_position + 80]

mutated_piece = list(original_piece)
mutated_piece[30] = "A" if mutated_piece[30] != "A" else "T"
my_read = "".join(mutated_piece)

print(f"Moj 'read' (80 slova, s 1 ubačenom razlikom na mjestu 30):\n  {my_read}\n")

results = align_read(my_read, sequence, sa, max_mismatches=2)

print("Rezultati poravnanja:")
for r in results:
    print(f"  pozicija={r.position}  razlike={r.mismatches}  lanac={r.strand}")

best = results[0] if results else None
success = best is not None and best.position == true_position and best.mismatches == 1
print(f"\nStvarna pozicija odakle je read uzet: {true_position}")
print(f"Provjera Faze 2: {'USPJEH' if success else 'NEUSPJEH'} (najbolji pogodak: {best})")
assert success, "Poravnanje nije pronašlo očekivanu poziciju/broj razlika."

# Faza 3: poravnaj skup readova (dio s pravim mutacijama, dio sa slučajnim greškama sekvenciranja) i pozovi mutacije iz pileupa
print("\nPripremam skup readova s dvije zasađene mutacije i par slučajnih grešaka ...")
rng = random.Random(0)

true_mutation_positions = [3800, 5500]
true_mutations = {}
for pos in true_mutation_positions:
    ref_base = sequence[pos]
    alt_base = next(b for b in "ACGT" if b != ref_base)
    true_mutations[pos] = (ref_base, alt_base)

read_length = 100
num_reads = 200
sample_start, sample_end = 3000, 7000  # readovi koncentrirani oko mutacija radi dovoljnog coveragea

reads = []
for _ in range(num_reads):
    start = rng.randint(sample_start, sample_end - read_length)
    bases = list(sequence[start:start + read_length])
    for pos, (ref_base, alt_base) in true_mutations.items():
        if start <= pos < start + read_length:
            bases[pos - start] = alt_base
    reads.append("".join(bases))

for i in rng.sample(range(num_reads), 5):  # izolirane greške sekvenciranja
    bases = list(reads[i])
    offset = rng.randrange(read_length)
    original = bases[offset]
    bases[offset] = rng.choice([b for b in "ACGT" if b != original])
    reads[i] = "".join(bases)

print(f"Poravnavam {num_reads} readova (duljina {read_length}, max_mismatches=2) ...")
alignment_results = align_reads(reads, sequence, sa, max_mismatches=2)

print("Pozivam mutacije (min_coverage=3, min_variant_fraction=0.5) ...\n")
calls = call_mutations(sequence, reads, alignment_results)
called_positions = {c.position for c in calls}

print(f"{len(calls)} mutacija(e) pozvano:")
for c in calls:
    print(
        f"  pozicija={c.position:>6}  {c.ref_base}->{c.alt_base}  "
        f"podrška={c.supporting_reads:>3}/{c.coverage:<3}  VAF={c.variant_allele_fraction:.2f}"
    )

all_planted_found = all(pos in called_positions for pos in true_mutations)
spurious = called_positions - set(true_mutations)

print("\nProvjera Faze 3:")
for pos, (ref_base, alt_base) in true_mutations.items():
    found = pos in called_positions
    print(f"  zasađena mutacija na {pos} ({ref_base}->{alt_base}): [{'OK' if found else 'PROPUŠTENO'}]")
print(
    "  slučajne greške ispravno filtrirane: "
    f"[{'OK' if not spurious else f'NEOČEKIVANO POZVANO: {spurious}'}]"
)

assert all_planted_found, "Nisu pronađene sve zasađene mutacije."
assert not spurious, "Slučajna greška je pogrešno pozvana kao mutacija."
