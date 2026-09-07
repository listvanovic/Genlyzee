#!/usr/bin/env bash
# Gradi STAR genome index za referentni genom (ako već ne postoji) i poravnava
# RNA-seq readove nad njim, s .sam izlazom. Vrijeme i peak memorija oba koraka
# mjere se s /usr/bin/time -l (BSD time na macOS-u).
#
# Koristi lokalno izgrađen binary ../bin/STAR (vidi build_star_macos.sh), a ne
# Homebrew rna-star: taj se linka na Appleov libc++, gdje je
# std::stringbuf::pubsetbuf() no-op, pa se STAR-ovi bufferi za čitanje i pisanje
# nikad ne popune i svaki run tiho obradi nula reada. Popravak je u ../patches/.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FLY_DIR="$ROOT_DIR/Fly"
STAR_DIR="$ROOT_DIR/star_comparison"
INDEX_DIR="$STAR_DIR/index"
RESULTS_DIR="$STAR_DIR/results/star"
STAR_BIN="$STAR_DIR/bin/STAR"

GENOME_FASTA="$FLY_DIR/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
READS_FASTQ="$FLY_DIR/IFM48h_1.fastq"

THREADS="${THREADS:-6}"

mkdir -p "$INDEX_DIR" "$RESULTS_DIR"

# Genom je malen (~144 Mb, ~1870 kontiga uklj. mnogo kratkih scaffolda), pa su
# genomeSAindexNbases i genomeChrBinNbits smanjeni ispod zadanih vrijednosti,
# kako STAR manual preporučuje za male/fragmentirane genome.
if [ ! -f "$INDEX_DIR/SAindex" ]; then
  echo "== Gradim STAR genome index =="
  /usr/bin/time -l "$STAR_BIN" \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --genomeSAindexNbases 12 \
    --genomeChrBinNbits 16 \
    2> "$RESULTS_DIR/genome_generate_time_mem.log"
else
  echo "== STAR genome index već postoji, preskačem genomeGenerate =="
fi

echo "== Pokrećem STAR alignment =="
/usr/bin/time -l "$STAR_BIN" \
  --runMode alignReads \
  --runThreadN "$THREADS" \
  --genomeDir "$INDEX_DIR" \
  --readFilesIn "$READS_FASTQ" \
  --outSAMtype SAM \
  --outFileNamePrefix "$RESULTS_DIR/IFM48h_1." \
  2> "$RESULTS_DIR/star_alignment_time_mem.log"

echo "== Gotovo. SAM izlaz: $RESULTS_DIR/IFM48h_1.Aligned.out.sam =="
