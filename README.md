# genome_index

Implementacija indeksiranja genoma suffix arrayem i pozivanja mutacija (SNP-ova) iz poravnatih readova, rađena kao dio završnog rada.

Algoritamska podloga je Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*, pogl. 3 i 9 — izgradnja suffix arraya (SA-IS), egzaktno/približno poravnanje preko binarne pretrage i seed-and-extend pristupa, te ideja "bubble"-a iz assembly grafova primijenjena na filtriranje mutacija od grešaka sekvenciranja.

## Struktura projekta

```
genome_index/
  fasta_io.py          parsiranje kromosoma iz (multi-)FASTA datoteke
  suffix_array.py       izgradnja suffix arraya (SA-IS i divsufsort backend)
  alignment.py           poravnanje reada/skupa readova nad suffix arrayem
  mutation_calling.py    pileup i pozivanje SNP-ova iz poravnatih readova
tests/                  jedinični i fuzz testovi (pytest)
moj_primjer.py           primjer upotrebe cijelog pipelinea, od FASTA do mutacija
```

## Instalacija

Potreban je Python 3.10+.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` uključuje `numpy`, `pydivsufsort` (brzi backend za suffix array na skali cijelog kromosoma) i `pytest`.

## Referentni genom

Repozitorij ne sadrži referentni genom (hg38 je nekoliko GB pa je isključen preko `.gitignore`). Za pokretanje `moj_primjer.py` ili primjera u `alignment.py`/`suffix_array.py` treba lokalno imati:

```
hg38/hg38.u.fa
```

FASTA datoteku s hg38 sklopom (npr. s UCSC-a, `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/`), raspakiranu i eventualno "unwrapped" (jedna linija po zaglavlju + sekvenca) — otud `.u.` u imenu.

## Pokretanje testova

```bash
pytest tests/ -v
```

Uključuje usporedbu SA-IS izgradnje s naivnim (izravno sortiranim) oracleom, te fuzz test koji poravnanje uspoređuje s brute-force pretragom kroz 560 nasumičnih slučajeva (egzaktni, mutirani, reverse-complement i potpuno nasumični readovi).

## Primjer

```bash
python moj_primjer.py
```

Učitava mitohondrijski kromosom (chrM) iz hg38, gradi mu suffix array, zatim demonstrira:

1. poravnanje jednog reada s ubačenom mutacijom natrag na njegovu pravu poziciju,
2. poravnanje skupa od 200 readova sa zasađenim mutacijama i slučajnim greškama sekvenciranja, te pozivanje mutacija iz pileupa — provjerava se da su zasađene mutacije pronađene i da su slučajne greške ispravno odbačene.
