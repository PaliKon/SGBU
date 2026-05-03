# dotplot_menu.py

from sequence_database import SequenceDatabase
from dotplot import save_dotplot_png, dotplot_statistics


def dotplot_manual_menu():
    """
    Úloha 14 - dot plot pre ručne zadané sekvencie.
    """
    print("\n=== 14. Dot plot - manualne zadanie sekvencii ===")

    seq1 = input("Zadaj prvu sekvenciu: ").strip()
    seq2 = input("Zadaj druhu sekvenciu: ").strip()

    try:
        word_size = int(input("Zadaj word size: "))
    except ValueError:
        print("Neplatny vstup. Word size musi byt cele cislo.")
        return

    output_file = input("Zadaj vystupny PNG subor [output/dotplot_manual.png]: ").strip()

    if output_file == "":
        output_file = "output/dotplot_manual.png"

    try:
        save_dotplot_png(seq1, seq2, output_file, word_size=word_size)
        stats = dotplot_statistics(seq1, seq2, word_size=word_size)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    print("\nDot plot bol ulozeny.")
    print(f"Subor: {output_file}")
    print(f"Rozmery matice: {stats['matrix_width']} x {stats['matrix_height']}")
    print(f"Pocet porovnanych buniek: {stats['total_cells']}")
    print(f"Pocet zhod: {stats['matches']}")
    print(f"Percento zhod: {stats['match_percent']:.2f}%")


def dotplot_database_menu(database: SequenceDatabase):
    """
    Úloha 14 - dot plot pre dve sekvencie z databázy.
    """
    print("\n=== 14. Dot plot - sekvencie z databazy ===")

    if database.count() < 2:
        print("Na dot plot potrebujes aspon 2 sekvencie v databaze.")
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

    try:
        word_size = int(input("Zadaj word size: "))
    except ValueError:
        print("Neplatny vstup. Word size musi byt cele cislo.")
        return

    default_output = f"output/dotplot_seq{first_number}_seq{second_number}_w{word_size}.png"

    output_file = input(f"Zadaj vystupny PNG subor [{default_output}]: ").strip()

    if output_file == "":
        output_file = default_output

    try:
        save_dotplot_png(seq1, seq2, output_file, word_size=word_size)
        stats = dotplot_statistics(seq1, seq2, word_size=word_size)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    print("\nDot plot bol ulozeny.")
    print(f"Subor: {output_file}")
    print(f"Prva sekvencia: {first_number}")
    print(f"Druha sekvencia: {second_number}")
    print(f"Word size: {word_size}")
    print(f"Rozmery matice: {stats['matrix_width']} x {stats['matrix_height']}")
    print(f"Pocet porovnanych buniek: {stats['total_cells']}")
    print(f"Pocet zhod: {stats['matches']}")
    print(f"Percento zhod: {stats['match_percent']:.2f}%")