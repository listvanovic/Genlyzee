# Usporedba: vlastiti aligner (suffix array) vs STAR

Podaci: Drosophila melanogaster RNA-seq, `Fly/IFM48h_1.fastq` (single-end, 76bp,
14.988.127 reada), referentni genom `Fly/Drosophila_melanogaster.BDGP6.dna.toplevel.fa`
(144 Mb, 1870 kontiga).

## Metodologija

1. STAR (2.7.11b) poravnava sve readove nad genomom, izlaz kao `.sam`.
2. Iz SAM-a se izdvajaju "čisti" readovi: primarni, ne-supplementary
   alignmenti (`FLAG & 0x100 == 0`, `FLAG & 0x800 == 0`, `FLAG & 0x4 == 0`)
   čiji je CIGAR točno `<duljina>M` (bez insercija, delecija, clippinga,
   splicanja) — 9.299.988 reada, 63.60% primarnih alignmenata. Ovo je i
   test-skup i STAR-ova referentna pozicija (RNAME + POS) za usporedbu. STAR
   se ovdje koristi kao referentni aligner, ne kao dokazana ground-truth
   vrijednost — pravi biološki izvor reada nije poznat jer podaci nisu
   simulirani, pa se rezultati u nastavku čitaju kao podudarnost s STAR-om,
   ne kao apsolutna točnost.
3. Isti readovi (kao FASTA) poravnavaju se vlastitim alignerom (suffix array,
   seed-and-extend, `max_mismatches=3`) nad jednim suffix arrayem izgrađenim
   nad cijelim genomom — svi kontigi spojeni s N-spacerima od 200bp da se
   spriječe lažna poklapanja preko granica kontiga.
4. Za svaki read se uspoređuje (RNAME, POS) iz STAR-a s (RNAME, POS) vlastitog
   alignera; podudaranje zahtijeva istoimeni RNAME i identičan POS.

STAR je poravnavao svih 14.988.127 reada, s punom podrškom za splicing i
indele; vlastiti aligner je testiran samo na podskupu od 9.299.988 "čistih"
reada, jer ne podržava splice/indel poravnanja. To je namjerno — mjeri se
točnost pozicioniranja na readovima koje bi vlastiti aligner uopće trebao
moći obraditi, ne cjelokupna RNA-seq funkcionalnost.

## Struktura SAM filea i CIGAR string

SAM je tekstualni format: zaglavlje (linije koje počinju s `@`; ovdje `@HD` s
verzijom formata i 1870 `@SQ` linija, po jedna za svaki kontig referentnog
genoma), a zatim jedan red po alignmentu s 11 obaveznih tab-odvojenih polja.
Primjer iz našeg outputa:

```
ILLUMINA-075005_0053_FC:3:1:17967:18033#0  0  2L  3515966  255  76M  *  0  0  GAGAACTTCGCC...  hhhgehghhggg...  NH:i:1  HI:i:1  AS:i:74  nM:i:0
```

| # | Polje | Vrijednost u primjeru | Značenje |
|---|---|---|---|
| 1 | QNAME | `ILLUMINA-...:18033#0` | ID reada |
| 2 | FLAG | `0` | bitmaska (0x4 = unmapped, 0x10 = reverse strand, 0x100 = sekundarni alignment, 0x800 = supplementary) |
| 3 | RNAME | `2L` | referentni kontig na koji je read poravnat |
| 4 | POS | `3515966` | 1-bazirana pozicija početka poravnanja |
| 5 | MAPQ | `255` | kvaliteta mapiranja (STAR koristi 255 za unique) |
| 6 | CIGAR | `76M` | kako read leži na referenci |
| 7-9 | RNEXT/PNEXT/TLEN | `*  0  0` | podaci o paru (nerelevantno, single-end) |
| 10 | SEQ | `GAGAACTTCGCC...` | sekvenca reada, orijentirana prema forward lancu reference |
| 11 | QUAL | `hhhgehghhggg...` | Phred kvalitete |
| 12+ | opcionalni tagovi | `NH:i:1 nM:i:0` | broj lokusa (`NH`), broj nepodudaranja (`nM`) |

CIGAR opisuje poravnanje kao niz parova `<broj><operacija>`, čitano slijeva
nadesno duž reada:

| Op | Značenje | Troši baze reada | Troši baze reference |
|---|---|---|---|
| `M` | poravnata baza (podudaranje ili nepodudaranje) | da | da |
| `I` | insercija u readu (baze kojih nema u referenci) | da | ne |
| `D` | delecija (baze reference kojih nema u readu) | ne | da |
| `N` | preskočena regija reference — kod RNA-seq to je intron (splice) | ne | da |
| `S` | soft clipping — kraj reada koji nije poravnat, ali je prisutan u SEQ | da | ne |
| `H` | hard clipping — odrezani kraj, nije ni prisutan u SEQ | ne | ne |

`M` ne znači da se baze podudaraju, nego samo da su poravnate bez pomaka —
nepodudaranja se broje odvojeno, u `nM` tagu. Read s CIGAR-om `76M` može
imati poneku point-mutaciju ili grešku sekvenciranja, što je upravo scenarij
koji vlastiti aligner podržava (mismatch-tolerant, bez indela).

Kriterij "poravnan bez prekida" je stoga: cijeli CIGAR mora biti točno jedna
`M` operacija čiji broj odgovara duljini reada — za naše readove od 76bp to
je `76M` (za read duljine 80 bio bi `80M`). Filter u
`scripts/filter_clean_reads.py` to provjerava regexom `^(\d+)M$` uz dodatnu
provjeru da se broj podudara s `len(SEQ)`, pa svaki `S`, `H`, `I`, `D` ili `N`
automatski ispada. Primjeri koji se odbacuju:

```
...:17532:18022#0  16  X   2001069   255  75M1S       ...   -> 1 baza soft-clipana
...:19344:18030#0  16  2R  13969542  255  61M69N15M   ...   -> read prelazi preko introna (69bp splice)
```

Raspodjela CIGAR operacija po svih 14.621.937 primarnih alignmenata
(izračunato nezavisnim `awk` prolazom, poklapa se s Python filterom):

| Kategorija | Broj reada | Udio |
|---|---|---|
| čisti `<duljina>M` (uzeti u analizu) | 9.299.988 | 63.60% |
| sadrži `S` (soft clipping) | 4.236.843 | 28.98% |
| sadrži `N` (splice/intron) | 1.338.755 | 9.16% |
| sadrži `D` (delecija) | 51.549 | 0.35% |
| sadrži `I` (insercija) | 46.401 | 0.32% |
| sadrži `H` (hard clipping) | 0 | 0.00% |

(kategorije se preklapaju — jedan read može imati i splice i clipping)

## Rezultati

### Točnost

| Mjera | Vrijednost | Što mjeri |
|---|---|---|
| Stopa jedinstvenog pozicioniranja (coverage) | 8.966.941 / 9.299.988 = 96.42% | koliko čistih reada vlastiti aligner uopće jedinstveno smjesti, bez obzira je li točno |
| — od toga multi-mapped | 179.897 / 9.299.988 = 1.93% | aligner našao >1 jednako dobru poziciju unutar `max_mismatches=3` |
| — od toga unmapped | 153.150 / 9.299.988 = 1.65% | aligner nije našao nijednu poziciju unutar `max_mismatches=3` |
| Podudarnost sa STAR-om među unique pozivima | 8.966.932 / 8.966.941 = 99.9999% | kad se aligner odluči za jednu poziciju, poklapa li se sa STAR-ovom |
| Stroga stopa identičnog RNAME+POS podudaranja | 8.966.932 / 9.299.988 = 96.42% | podudaranja / svi čisti readovi — multi-mapped i unmapped ovdje broje se kao promašaj |

Ove mjere se namjerno ne svode na jedan broj, jer odgovaraju na različita
pitanja. Gubitak u strogoj stopi podudaranja (96.42% umjesto blizu 100%)
dolazi gotovo isključivo iz nedostatka pokrivenosti — među 8.966.941 reada
gdje se aligner odlučio za jednu poziciju, s STAR-om se ne slaže samo 9.

### Podjela po STAR-ovom NH

STAR sam prijavljuje `NH:i:N` kad je read podjednako dobro poravnat na N
mjesta i jedno bira kao primarno. Od 9.299.988 čistih reada, 109.411 (1.18%)
ima STAR `NH≥2` — dakle i STAR je za njih neodlučan, pa stroga usporedba s
jednom STAR pozicijom tu nije baš najbolja mjera (`scripts/extract_nh.py`):

| Skup | Broj reada | Unique | Multi | Unmapped | Podudarnost sa STAR-om |
|---|---|---|---|---|---|
| STAR `NH:i:1` | 9.190.577 | 8.962.784 (97.52%) | 78.951 (0.86%) | 148.842 (1.62%) | 97.52% strogo, 99.9999% od unique poziva |
| STAR `NH≥2` | 109.411 | 4.157 (3.80%) | 100.946 (92.26%) | 4.308 (3.94%) | 3.80% strogo, 99.9759% od unique poziva |

Na podskupu gdje je i STAR siguran (`NH:i:1`) vlastiti aligner postiže
97.52% strogu podudarnost, više od ukupnog prosjeka, jer se isključuju
readovi koji su dvosmisleni po definiciji. Na podskupu gdje je STAR sam
prijavio više mogućih lokusa, vlastiti aligner uglavnom (92.26%) također
vraća multi-mapping umjesto jedne pozicije — dva neovisna alata se slažu da
su ti readovi dvosmisleni, umjesto da aligner nasumično bira jedno mjesto.

Jedini "unique" poziv unutar `NH≥2` skupa koji se ne slaže sa STAR-om (read
`...3:14:17033:20323#0`) provjeren je i protiv STAR-ovog sekundarnog lokusa,
ne samo primarnog: STAR je prijavio 3L:28.103.368 (primarni) i X:5.752.974
(sekundarni), a vlastiti aligner treću poziciju, 3L:18.780.015 — ne poklapa
se ni s jednom od te dvije.

### Tih 9 neslaganja

Za svih 9 spornih reada izbrojano je koliko nepodudaranja read ima naspram
genoma na STAR-ovoj i na poziciji vlastitog alignera
(`scripts/inspect_mismatches.py`, rezultat u
`results/comparison/mismatches_verified.csv`):

| Read | STAR pozicija (nepodudaranja) | Pozicija vlastitog alignera (nepodudaranja) |
|---|---|---|
| `...3:8:10065:3458#0` | 3L:3.901.348 (5) | 3L:3.901.576 (3) |
| `...3:14:17033:20323#0` | 3L:28.103.368 (9) | 3L:18.780.015 (3) |
| preostalih 7 reada | 2L:14.743.463 (4) | 2L:14.743.493 (3) |

Na svih 9 pozicija vlastiti aligner ima manju Hammingovu udaljenost od
STAR-a — brojanje se poklapa i sa STAR-ovim vlastitim `nM` tagom (`nM:i:5`,
`nM:i:9`, `nM:i:4`), pa je isključena greška u samom brojanju. Osam od 9
reada ima STAR `NH:i:1` (STAR ih smatra jedinstveno mapiranima, bez ikakve
dvosmislenosti na svojoj strani).

Ovo pokazuje da vlastiti aligner na tih 9 pozicija daje bolje rangirano
poravnanje po broju nepodudaranja — ne dokazuje se time da je ta pozicija i
stvarni biološki izvor reada, jer prava koordinata podrijetla nije poznata
(readovi nisu simulirani). Niža Hammingova udaljenost može doći od prave
lokacije, ali i od slučajnog poklapanja s ponovljenom sekvencom negdje
drugdje u genomu. Uzrok razlike vjerojatno leži u različitim strategijama
generiranja kandidata i bodovanja — STAR maksimizira svoj alignment score
(`AS`) i ima drugačiju seed-pretragu od pigeonhole seed-and-extend pristupa
korištenog ovdje. Budući da svih 9 pozicija ima CIGAR `76M` na obje strane
usporedbe, razlika se ne može pripisati STAR-ovom splice modelu — ovim
testom uzrok STAR-ovog izbora nije izoliran.

Brojanje nepodudaranja provjereno je i na 2000 nasumičnih reada iz skupa na
kojem se alati slažu: 73% ima 0 nepodudaranja na prijavljenoj poziciji, 94%
najviše 1 — što je sanity check metode, ne dokaz o cijelom skupu.

### Off-by-one provjera na cijelom skupu

`compare_coordinates.py` ima parametar `--tolerance`. Pokrenuto na svih
9.299.988 reada s tolerancijom 0 i ±1 baza: oba daju identičnih 8.966.932
podudaranja. Proširenje tolerancije ništa ne mijenja, čime je isključena
sustavna off-by-one greška u pretvorbi koordinata (0-bazirane pozicije
suffix arraya u 1-bazirane SAM pozicije).

### Vrijeme izvršavanja

Vrijeme i memorija mjereni su alatom `/usr/bin/time -l` (macOS BSD `time`).
"real" je stvarno proteklo vrijeme, "user" zbroj CPU vremena preko svih
threadova, "sys" vrijeme u kernelu.

| Korak | Alat | Threadova | Broj reada | real | user | sys | reads/s |
|---|---|---|---|---|---|---|---|
| Genome index | STAR | 6 | — | 53.49 s | 189.83 s | 15.47 s | — |
| Genome index | vlastiti aligner | 1 | — | 9.06 s | 8.35 s | 0.48 s | — |
| Alignment | STAR | 6 | 14.988.127 | 97.08 s | 414.57 s | 52.20 s | ~154.400 |
| Alignment | vlastiti aligner | 1 | 9.299.988 | 1009.52 s (16m50s) | 968.57 s | 24.54 s | ~9.210 |

STAR je pokrenut s `--runThreadN 6`; vlastiti aligner je single-threaded
(nema threadinga ni multiprocessinga u `run_my_aligner.py`), što se vidi i
u omjeru user/real vremena — STAR-ov ≈4.27, vlastiti ≈0.96. STAR je ~16.8x
brži po readu (6.5 µs naspram 108.6 µs), ali dio te razlike dolazi od
nejednakog broja threadova, ne samo od algoritma — brojka je throughput u
ovoj konkretnoj konfiguraciji, ne normalizirana usporedba po jezgri. STAR je
uz to industrijski C++ alat sa suffix-array pretragom optimiziranom
godinama; vlastiti aligner je čist Python (`divsufsort` C backend koristi se
samo za izgradnju suffix arraya, sama pretraga po readu je Python petlja).

### Memory usage (peak RSS)

Vrijednost je `maximum resident set size` iz `/usr/bin/time -l`, što
odgovara `ru_maxrss` iz `getrusage(2)`. Na macOS-u to polje je u bajtovima
(na Linuxu bi bilo u kilobajtima):

| Korak | Alat | Peak RSS (bajtova) | Peak RSS (GB) |
|---|---|---|---|
| Genome index | STAR | 1.535.000.576 | 1.54 GB |
| Genome index | vlastiti aligner | 891.928.576 | 0.89 GB |
| Alignment | STAR | 1.978.810.368 | 1.98 GB |
| Alignment | vlastiti aligner | 1.562.148.864 | 1.56 GB |

Memorija je usporediva unatoč razlici u implementaciji, a vlastiti aligner
je u alignment koraku čak nešto štedljiviji. Oba alata dominantno troše
memoriju na genomski indeks (STAR-ov suffix-array indeks vs. naš suffix
array + string genoma), koji je za ovaj genom (144 Mb) reda veličine 1-2GB
kod oba pristupa.

## Komentar

Podudarnost pozicioniranja vlastitog alignera sa STAR-om je vrlo visoka na
readovima koje je i dizajniran obraditi: kad pronađe jedinstvenu poziciju,
poklapa se sa STAR-om u 99.9999% slučajeva, a na preostalih 9 neslaganja
daje poravnanje s manjom Hammingovom udaljenošću od STAR-ovog. To je dobar
znak da temeljni algoritam (suffix array binarna pretraga + pigeonhole
seed-and-extend) radi ispravno, iako se time ne dokazuje koja je pozicija
stvarni biološki izvor reada.

Glavni nedostatak nije netočnost nego pokrivenost: 3.58% čistih reada
vlastiti aligner ne uspije jedinstveno smjestiti. To je očekivano s obzirom
da je genom pun ponavljajućih sekvenci (1870 kontiga, mnogo kratkih
scaffolda), a `max_mismatches=3` seed-and-extend pristup namjerno preskače
seedove s prekomjernim brojem pogodaka radi performansi — što se poklapa s
mehanizmom multi-mapping/unmapped ishoda, i s time da STAR na istim readovima
(`NH≥2`) također pokazuje neodlučnost.

Brzina po readu je očekivano nepovoljnija za vlastiti aligner (~16.8x
sporije u izmjerenom throughputu), dijelom zbog nejednakog broja threadova
(STAR 6, vlastiti aligner 1). Apsolutno gledano vrijeme je i dalje praktično
za skup ove veličine (16m50s za 9.3M reada), a memory footprint usporediv sa
STAR-om zahvaljujući `divsufsort` C backendu za izgradnju indeksa.

Usporedba nije potpuno "fer" u smislu opsega funkcionalnosti — STAR radi
cijeli RNA-seq pipeline sa splicingom, vlastiti aligner rješava samo ungapped
mismatch-tolerant slučaj — ali je informativna za pitanje koje ovaj rad
postavlja: koliko je pouzdana pozicija koju vrati vlastiti suffix-array
aligner, mjereno naspram etabliranog alata na istim readovima.

## Poznati problem i rješenje: STAR na macOS-u (Apple Silicon)

Konfiguracija na kojoj je problem opažen: STAR 2.7.11b, Homebrew build
(`rna-star`), macOS 26.5 (build 25F71), Apple Silicon (arm64), Apple clang
15.0.0 (`arm64-apple-darwin25.5.0`) s Appleovom implementacijom standardne
C++ biblioteke (libc++). Na toj konfiguraciji `alignReads` tiho vraća 0
poravnatih reada — exit code 0, izgleda kao uspjeh. Nije provjereno javlja
li se isti problem na drugim verzijama macOS-a ili STAR-a.

Uzrok: `std::stringbuf::pubsetbuf()` je standardna C++ funkcija, ali njeno
konkretno ponašanje za `stringbuf` standard ostavlja implementation-defined
— ne postoji jamstvo da će predani `char*` postati stvarni interni buffer
streama. STAR 2.7.11b se oslanja upravo na to da hoće. U testiranom buildu
to je funkcioniralo s GNU-ovim libstdc++ (ponašanje koje ta implementacija
bira, ne dio standarda), dok Appleov libc++ isti poziv tretira kao no-op i
ne povezuje vanjski buffer. Rezultat: STAR-ov buffer za čitanje readova i
buffer za pisanje SAM zapisa ostaju prazni bez obzira na stvarni sadržaj
streamova, a program ipak završava s exit kodom 0.

Riješeno kompajliranjem STAR-a iz izvornog koda uz patch koji zamjenjuje
`pubsetbuf` pristup eksplicitnim `std::stringbuf::str()` pozivom — kopira
sadržaj u stringbuf umjesto da se oslanja na aliasing vanjskog buffera (vidi
`star_comparison/patches/star_2.7.11b_macos_libcxx_fix.patch` i
`star_comparison/scripts/build_star_macos.sh`). `stringbuf::str()` je
standardno, prenosivo ponašanje, pa se na Linuxu/libstdc++-u očekuje isti
rezultat, iako to nije testirano u sklopu ovog rada.

## Reproducibilnost

Ulazni podaci (`Fly/`) i međurezultati veći od nekoliko stotina KB nisu u
gitu (vidi `.gitignore`) — commitane su `.sample.*` verzije (header + prvih
~2000 redaka) i skripte kojima se puni fileovi regeneriraju lokalno, redom:

```bash
star_comparison/scripts/build_star_macos.sh                     # build patchanog STAR-a -> bin/STAR
star_comparison/scripts/run_star.sh                              # genome index + SAM     -> results/star/
star_comparison/scripts/filter_clean_reads.py  <sam> <out>       # -> results/filtered/clean_reads.tsv
star_comparison/scripts/make_fasta_input.py    <tsv> <out>       # -> results/my_aligner/clean_reads_input.fasta
star_comparison/scripts/build_my_index.py                        # -> my_index/
star_comparison/scripts/run_my_aligner.py      <fasta> <out>     # -> results/my_aligner/my_alignments.tsv
star_comparison/scripts/compare_coordinates.py <star> <my> <out> --mismatches-out <out2>
star_comparison/scripts/inspect_mismatches.py  <out2> <clean> <out3>
star_comparison/scripts/extract_nh.py          <sam> <out>       # -> results/filtered/clean_reads_nh.tsv
```

Pipeline je deterministički, pa ponovno pokretanje na istim ulazima daje
identične brojke.
