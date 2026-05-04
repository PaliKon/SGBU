# Analyzátor mRNA sekvencií

Semestrálna práca z predmetu **Biomolekulárne základy v genetike a biomedicíne**.
Projekt je interaktívna konzolová aplikácia v jazyku Python na prácu s mRNA sekvenciami vo formáte FASTA.

## Implementované funkcionality

Program pokrýva povinné úlohy 1–16 a voliteľnú úlohu 17:

1. kontrola identity sekvencií a validácia abecedy,
2. preklad kodónu na aminokyselinu podľa štandardnej genetickej tabuľky,
3. preklad mRNA na proteín od štart kodónu `AUG` po stop kodón,
4. načítanie a uloženie FASTA súborov,
5. načítanie viacerých mRNA sekvencií do pamäte,
6. výpis zoznamu načítaných sekvencií,
7. výpis konkrétnej sekvencie podľa poradového čísla,
8. frekvenčná analýza mRNA,
9. generovanie náhodnej mRNA sekvencie,
10. hľadanie ORF kandidátov vrátane dĺžky a GC obsahu,
11. označenie ORF ako CDS a preklad na proteín,
12. bodová mutácia,
13. Hammingova vzdialenosť,
14. dot plot a export do PNG,
15. Needleman-Wunsch globálne zarovnanie,
16. Smith-Waterman lokálne zarovnanie,
17. voliteľne index k-tíc pre `k = 3..9`.

Úlohy BLAST a MSA nie sú v tejto verzii implementované.

## Štruktúra projektu

```text
bioinformatika_semestralka/
├── main.py
├── sequence_utils.py
├── fasta_utils.py
├── sequence_database.py
├── mrna_analysis.py
├── mrna_analysis_menu.py
├── alignment.py
├── alignment_menu.py
├── dotplot.py
├── dotplot_menu.py
├── kmer_index.py
├── kmer_index_menu.py
├── experiments.py
├── plot_experiments.py
├── requirements.txt
├── data/
│   └── real/
│       ├── actb_human.fasta
│       ├── brca1_human.fasta
│       ├── gapdh_human.fasta
│       ├── hbb_human.fasta
│       └── tp53_human.fasta
└── output/
    └── experiments/
```

## Požiadavky

Odporúčaná verzia Pythonu:

```bash
python 3.10+
```

Použité externé knižnice:

```bash
pillow
matplotlib
```

## Inštalácia

V koreňovom priečinku projektu je možné vytvoriť virtuálne prostredie:

### Windows

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

### Linux / macOS

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Spustenie hlavnej aplikácie

```bash
python main.py
```

Po spustení sa zobrazí interaktívne menu. Hodnota `0` ukončí program.

Odporúčaný postup pri práci s dátami:

1. vybrať možnosť `2` a načítať sekvencie z priečinka `data/real`,
2. možnosťou `5` skontrolovať načítané sekvencie,
3. následne spúšťať požadované analýzy a algoritmy.

## Spustenie experimentov

Experimenty sa spúšťajú samostatne:

```bash
python experiments.py
```

Skript načíta sekvencie z priečinka `data/real`, porovná algoritmy Needleman-Wunsch a Smith-Waterman, zmeria čas výpočtu a vytvorí súhrn indexu k-tíc.

Výstupy sa ukladajú do priečinka:

```text
output/experiments/
```

Vytvorené CSV súbory:

```text
dataset_summary.csv
nw_sw_scores.csv
nw_sw_timing.csv
kmer_index_summary.csv
experiment_info.txt
```

## Vytvorenie grafov

Po spustení experimentov je možné vytvoriť grafy:

```bash
python plot_experiments.py
```

Grafy sa uložia do:

```text
output/experiments/plots/
```

## FASTA dáta

Použité sú reálne mRNA sekvencie z databázy NCBI GenBank/RefSeq:

| Súbor | Accession ID | Gén |
|---|---|---|
| `actb_human.fasta` | `NM_001101.5` | ACTB |
| `brca1_human.fasta` | `NM_007294.4` | BRCA1 |
| `gapdh_human.fasta` | `NM_002046.7` | GAPDH |
| `hbb_human.fasta` | `NM_000518.5` | HBB |
| `tp53_human.fasta` | `NM_000546.6` | TP53 |

Pri načítaní FASTA súborov program konvertuje znak `T` na `U`, aby sa ďalej pracovalo s mRNA abecedou `A, U, C, G`.

## Výstupy programu

- FASTA export: `output/exported_all_from_data.fasta` alebo `output/individual_sequences/`
- dot plot PNG: `output/dotplot_manual.png`, `output/dotplot_seqX_seqY_wK.png`
- experimentálne CSV a grafy: `output/experiments/`
- ostatné výsledky sa vypisujú do konzoly

## Poznámky

- Pri Hammingovej vzdialenosti musia mať sekvencie rovnakú dĺžku.
- Pri reálnych dlhších sekvenciách je vhodné pri NW/SW používať prefixy, keďže dynamické programovanie má časovú aj pamäťovú zložitosť `O(n*m)`.
- Pri experimentoch je počet rekonštruovaných optimálnych zarovnaní obmedzený na 1, aby meranie času nebolo ovplyvnené veľkým počtom traceback ciest. Interaktívna aplikácia však podporuje výpis všetkých najlepších zarovnaní do nastaveného limitu.
