# alignment_menu.py

from sequence_database import SequenceDatabase
from alignment import hamming_comparison_details


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