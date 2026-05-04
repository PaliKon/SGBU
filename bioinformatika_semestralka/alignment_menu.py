# alignment_menu.py

from sequence_database import SequenceDatabase

from alignment import (
    hamming_comparison_details,
    needleman_wunsch,
    smith_waterman,
    alignment_match_line,
)


def hamming_distance_manual_menu():
    """
    Úloha 13 - manuálne zadanie dvoch sekvencií.
    """
    print("\n=== 13. Hammingova vzdialenost - manualne zadanie ===")

    seq1 = input("Zadaj prvu sekvenciu: ").strip()
    seq2 = input("Zadaj druhu sekvenciu: ").strip()

    try:
        result = hamming_comparison_details(seq1, seq2)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    print_hamming_result(result)


def hamming_distance_database_menu(database: SequenceDatabase):
    """
    Úloha 13 - výber dvoch sekvencií z databázy.
    """
    print("\n=== 13. Hammingova vzdialenost - sekvencie z databazy ===")

    if database.count() < 2:
        print("Na porovnanie potrebujes aspon 2 sekvencie v databaze.")
        print("Najprv nacitaj viac FASTA suborov.")
        return

    print("\nDostupne sekvencie:")
    database.list_sequences()

    try:
        first_number = int(input("\nZadaj cislo prvej sekvencie: "))
        second_number = int(input("Zadaj cislo druhej sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    if first_number < 1 or first_number > database.count():
        print("Neplatne cislo prvej sekvencie.")
        return

    if second_number < 1 or second_number > database.count():
        print("Neplatne cislo druhej sekvencie.")
        return

    seq1 = database.get_sequence(first_number)
    seq2 = database.get_sequence(second_number)

    print("\nPoznamka:")
    print("Hammingova vzdialenost vyzaduje rovnako dlhe sekvencie.")
    print("Ak su realne mRNA sekvencie roznej dlzky, porovnanie skonci chybou.")

    try:
        result = hamming_comparison_details(seq1, seq2)
    except ValueError as error:
        print(f"Chyba: {error}")
        print("Tip: pre realne sekvencie roznej dlzky bude vhodnejsie globalne alebo lokalne zarovnanie.")
        return

    print_hamming_result(result)


def print_hamming_result(result: dict):
    """
    Pekne vypíše výsledok Hammingovej vzdialenosti.
    """
    print("\nVysledok porovnania:")
    print(f"Dlzka sekvencii: {result['length']}")
    print(f"Hammingova vzdialenost: {result['distance']}")
    print(f"Pocet rovnakych pozicii: {result['matches']}")
    print(f"Pocet rozdielnych pozicii: {result['distance']}")

    if result["length"] > 0:
        percent_identity = (result["matches"] / result["length"]) * 100
        percent_difference = (result["distance"] / result["length"]) * 100

        print(f"Percentualna zhoda: {percent_identity:.2f}%")
        print(f"Percentualny rozdiel: {percent_difference:.2f}%")

    if result["differences"]:
        print("\nRozdielne pozicie:")

        max_to_print = 30

        for diff in result["differences"][:max_to_print]:
            print(
                f"pozicia {diff['position']}: "
                f"{diff['seq1_base']} -> {diff['seq2_base']}"
            )

        if len(result["differences"]) > max_to_print:
            remaining = len(result["differences"]) - max_to_print
            print(f"... a dalsich {remaining} rozdielov")
    else:
        print("\nSekvencie su identicke.")

def needleman_wunsch_manual_menu():
    """
    Úloha 15 - Needleman-Wunsch pre ručne zadané sekvencie.
    """
    print("\n=== 15. Needleman-Wunsch - manualne zadanie ===")

    seq1 = input("Zadaj prvu sekvenciu: ").strip()
    seq2 = input("Zadaj druhu sekvenciu: ").strip()

    scoring = read_scoring_parameters()

    if scoring is None:
        return

    match_score, mismatch_score, gap_score = scoring

    result = needleman_wunsch(
        seq1,
        seq2,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score
    )

    print_needleman_wunsch_result(result)


def needleman_wunsch_database_menu(database: SequenceDatabase):
    """
    Úloha 15 - Needleman-Wunsch pre dve sekvencie z databázy.
    """
    print("\n=== 15. Needleman-Wunsch - sekvencie z databazy ===")

    if database.count() < 2:
        print("Na zarovnanie potrebujes aspon 2 sekvencie v databaze.")
        print("Najprv nacitaj FASTA subory.")
        return

    print("\nDostupne sekvencie:")
    database.list_sequences()

    try:
        first_number = int(input("\nZadaj cislo prvej sekvencie: "))
        second_number = int(input("Zadaj cislo druhej sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    if first_number < 1 or first_number > database.count():
        print("Neplatne cislo prvej sekvencie.")
        return

    if second_number < 1 or second_number > database.count():
        print("Neplatne cislo druhej sekvencie.")
        return

    seq1 = database.get_sequence(first_number)
    seq2 = database.get_sequence(second_number)

    print("\nPozor:")
    print("Needleman-Wunsch ma casovu a pamatovu narocnost O(n*m).")
    print("Pri velmi dlhych realnych mRNA sekvenciach moze byt vypocet pomalsi.")
    print("Na rychly test mozes pouzit iba prvych N baz.")

    use_prefix = input("\nChces pouzit iba prvych N baz zo sekvencii? (a/n): ").strip().lower()

    if use_prefix == "a":
        try:
            n = int(input("Zadaj N: "))
        except ValueError:
            print("Neplatny vstup. Zadaj cele cislo.")
            return

        seq1 = seq1[:n]
        seq2 = seq2[:n]

    scoring = read_scoring_parameters()

    if scoring is None:
        return

    match_score, mismatch_score, gap_score = scoring

    result = needleman_wunsch(
        seq1,
        seq2,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score
    )

    print_needleman_wunsch_result(result)


def read_scoring_parameters():
    """
    Načíta skórovacie parametre od používateľa.
    """
    print("\nZadaj skore.")
    print("Ak nechas prazdne, pouziju sa predvolene hodnoty:")
    print("match = 1, mismatch = -1, gap = -2")

    match_input = input("Match score [1]: ").strip()
    mismatch_input = input("Mismatch score [-1]: ").strip()
    gap_input = input("Gap score [-2]: ").strip()

    try:
        match_score = int(match_input) if match_input else 1
        mismatch_score = int(mismatch_input) if mismatch_input else -1
        gap_score = int(gap_input) if gap_input else -2
    except ValueError:
        print("Neplatny vstup. Skore musi byt cele cislo.")
        return None

    return match_score, mismatch_score, gap_score


def print_needleman_wunsch_result(result: dict):
    """
    Vypíše výsledok Needleman-Wunsch algoritmu.
    """
    print("\nVysledok Needleman-Wunsch:")
    print(f"Najlepsie globalne skore: {result['score']}")
    print(f"Pocet vypisanych najlepsich zarovnani: {result['alignment_count']}")

    if result["max_alignments_reached"]:
        print("Upozornenie: Dosiahnuty limit poctu zarovnani, nevypisuju sa dalsie.")

    print("\nNajlepsie globalne zarovnania:")

    for index, alignment in enumerate(result["alignments"], start=1):
        aligned_seq1 = alignment["seq1"]
        aligned_seq2 = alignment["seq2"]
        match_line = alignment_match_line(aligned_seq1, aligned_seq2)

        print(f"\nZarovnanie {index}:")
        print(aligned_seq1)
        print(match_line)
        print(aligned_seq2)

def smith_waterman_manual_menu():
    """
    Úloha 16 - Smith-Waterman pre ručne zadané sekvencie.
    """
    print("\n=== 16. Smith-Waterman - manualne zadanie ===")

    seq1 = input("Zadaj prvu sekvenciu: ").strip()
    seq2 = input("Zadaj druhu sekvenciu: ").strip()

    scoring = read_sw_scoring_parameters()

    if scoring is None:
        return

    match_score, mismatch_score, gap_score = scoring

    result = smith_waterman(
        seq1,
        seq2,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score
    )

    print_smith_waterman_result(result)


def smith_waterman_database_menu(database: SequenceDatabase):
    """
    Úloha 16 - Smith-Waterman pre dve sekvencie z databázy.
    """
    print("\n=== 16. Smith-Waterman - sekvencie z databazy ===")

    if database.count() < 2:
        print("Na zarovnanie potrebujes aspon 2 sekvencie v databaze.")
        print("Najprv nacitaj FASTA subory.")
        return

    print("\nDostupne sekvencie:")
    database.list_sequences()

    try:
        first_number = int(input("\nZadaj cislo prvej sekvencie: "))
        second_number = int(input("Zadaj cislo druhej sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    if first_number < 1 or first_number > database.count():
        print("Neplatne cislo prvej sekvencie.")
        return

    if second_number < 1 or second_number > database.count():
        print("Neplatne cislo druhej sekvencie.")
        return

    seq1 = database.get_sequence(first_number)
    seq2 = database.get_sequence(second_number)

    print("\nSmith-Waterman hlada najlepsie lokalne zarovnanie.")
    print("Pri dlhych realnych sekvenciach moze byt vypocet pomalsi.")
    print("Na rychly test mozes pouzit iba prvych N baz.")

    use_prefix = input("\nChces pouzit iba prvych N baz zo sekvencii? (a/n): ").strip().lower()

    if use_prefix == "a":
        try:
            n = int(input("Zadaj N: "))
        except ValueError:
            print("Neplatny vstup. Zadaj cele cislo.")
            return

        seq1 = seq1[:n]
        seq2 = seq2[:n]

    scoring = read_sw_scoring_parameters()

    if scoring is None:
        return

    match_score, mismatch_score, gap_score = scoring

    result = smith_waterman(
        seq1,
        seq2,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score
    )

    print_smith_waterman_result(result)


def read_sw_scoring_parameters():
    """
    Načíta skórovacie parametre pre Smith-Waterman.
    """
    print("\nZadaj skore.")
    print("Ak nechas prazdne, pouziju sa predvolene hodnoty:")
    print("match = 2, mismatch = -1, gap = -2")

    match_input = input("Match score [2]: ").strip()
    mismatch_input = input("Mismatch score [-1]: ").strip()
    gap_input = input("Gap score [-2]: ").strip()

    try:
        match_score = int(match_input) if match_input else 2
        mismatch_score = int(mismatch_input) if mismatch_input else -1
        gap_score = int(gap_input) if gap_input else -2
    except ValueError:
        print("Neplatny vstup. Skore musi byt cele cislo.")
        return None

    return match_score, mismatch_score, gap_score


def print_smith_waterman_result(result: dict):
    """
    Vypíše výsledok Smith-Waterman algoritmu.
    """
    print("\nVysledok Smith-Waterman:")
    print(f"Najlepsie lokalne skore: {result['score']}")
    print(f"Pocet vypisanych najlepsich lokalnych zarovnani: {result['alignment_count']}")

    if result["score"] == 0:
        print("Nenaslo sa ziadne pozitivne lokalne zarovnanie.")
        return

    if result["max_alignments_reached"]:
        print("Upozornenie: Dosiahnuty limit poctu zarovnani, nevypisuju sa dalsie.")

    print("\nNajlepsie lokalne zarovnania:")

    for index, alignment in enumerate(result["alignments"], start=1):
        aligned_seq1 = alignment["seq1"]
        aligned_seq2 = alignment["seq2"]
        match_line = alignment_match_line(aligned_seq1, aligned_seq2)

        print(f"\nZarovnanie {index}:")
        print(
            f"Seq1 pozicie: {alignment['seq1_start']} - {alignment['seq1_end']}"
        )
        print(
            f"Seq2 pozicie: {alignment['seq2_start']} - {alignment['seq2_end']}"
        )
        print(aligned_seq1)
        print(match_line)
        print(aligned_seq2)