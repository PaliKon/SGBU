# main.py

from sequence_utils import (
    MRNA_ALPHABET,
    are_identical_and_valid,
    codon_to_amino_acid,
    translate_mrna_to_protein,
)


from fasta_utils import write_fasta, write_fasta_individual
from sequence_database import SequenceDatabase

from mrna_analysis_menu import (
    nucleotide_frequency_menu,
    generate_random_mrna_menu,
    find_orfs_menu,
    cds_translation_menu,
    point_mutation_menu,
)

from alignment_menu import (
    hamming_distance_manual_menu,
    hamming_distance_database_menu,
)

from dotplot_menu import (
    dotplot_manual_menu,
    dotplot_database_menu,
)


def test_basic_operations():
    seq1 = "AUGGCUUAA"
    seq2 = "AUGGCUUAA"

    print("\n=== Test 1-3: Zakladne operacie ===")

    print("Test identity a abecedy:")
    print(are_identical_and_valid(seq1, seq2, MRNA_ALPHABET))

    print("\nTest prevodu kodonu:")
    print("AUG ->", codon_to_amino_acid("AUG"))
    print("GCU ->", codon_to_amino_acid("GCU"))
    print("UAA ->", codon_to_amino_acid("UAA"))

    print("\nTest prekladu mRNA na protein:")
    print(seq1, "->", translate_mrna_to_protein(seq1))


def load_sequences_menu(database: SequenceDatabase):
    directory = input("Zadaj priecinok so FASTA subormi [data]: ").strip()

    if directory == "":
        directory = "data"

    database.load_from_directory(directory)

    print(f"\nPocet sekvencii v databaze: {database.count()}")


def list_sequences_menu(database: SequenceDatabase):
    print("\n=== Zoznam sekvencii ===")
    database.list_sequences()


def print_sequence_menu(database: SequenceDatabase):
    if database.count() == 0:
        print("Databaza je prazdna. Najprv nacitaj sekvencie.")
        return

    database.list_sequences()

    try:
        number = int(input("\nZadaj cislo sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    if number < 1 or number > database.count():
        print("Neplatne cislo sekvencie.")
        return

    print("\nAko chces vypisat sekvenciu?")
    print("1. Celu sekvenciu")
    print("2. Len prvych N baz")

    choice = input("Vyber moznost: ").strip()

    if choice == "1":
        database.print_sequence(number)

    elif choice == "2":
        try:
            length = int(input("Kolko prvych baz chces vypisat? "))
        except ValueError:
            print("Neplatny vstup. Zadaj cele cislo.")
            return

        sequence = database.get_sequence(number)
        record = database.records[number - 1]

        print(f">{record.get('header') or f'sequence_{number}'}")
        print(sequence[:length])

    else:
        print("Neplatna moznost.")


def export_all_to_one_file_menu(database: SequenceDatabase):
    if database.count() == 0:
        print("Databaza je prazdna. Najprv nacitaj sekvencie.")
        return

    output_file = input(
        "Zadaj vystupny subor [output/exported_all_from_data.fasta]: "
    ).strip()

    if output_file == "":
        output_file = "output/exported_all_from_data.fasta"

    write_fasta(database.records, output_file)

    print(f"Sekvencie boli ulozene do jedneho suboru: {output_file}")


def export_individual_files_menu(database: SequenceDatabase):
    if database.count() == 0:
        print("Databaza je prazdna. Najprv nacitaj sekvencie.")
        return

    output_directory = input(
        "Zadaj vystupny priecinok [output/individual_sequences]: "
    ).strip()

    if output_directory == "":
        output_directory = "output/individual_sequences"

    write_fasta_individual(database.records, output_directory)

    print(f"Kazda sekvencia bola ulozena samostatne do priecinka: {output_directory}")


def show_database_summary(database: SequenceDatabase):
    print("\n=== Kontrola nacitanych sekvencii ===")

    if database.count() == 0:
        print("Databaza je prazdna.")
        return

    for i in range(1, database.count() + 1):
        record = database.records[i - 1]
        sequence = record["sequence"]

        print(f"\nSekvencia {i}:")
        print("Header:", record.get("header") or f"sequence_{i}")
        print("Dlzka:", len(sequence))
        print("Prvych 80 baz:", sequence[:80])




def show_menu():
    print("\n=== Bioinformatika semestralna praca ===")
    print("1. Test zakladnych operacii 1-3")
    print("2. Nacitat vsetky FASTA subory z priecinka")
    print("3. Vypisat zoznam nacitanych sekvencii")
    print("4. Vypisat konkretnu sekvenciu")
    print("5. Kontrola nacitanych sekvencii")
    print("6. Exportovat vsetky sekvencie do jedneho FASTA suboru")
    print("7. Exportovat kazdu sekvenciu zvlast")
    print("8. Frekvencna analyza mRNA sekvencie")
    print("9. Generovat nahodnu mRNA sekvenciu")
    print("10. Hladat ORF kandidatov")
    print("11. Oznacit ORF ako CDS a prelozit na protein")
    print("12. Vytvorit bodovu mutaciu")
    print("13. Hammingova vzdialenost - manualne zadanie")
    print("14. Hammingova vzdialenost - sekvencie z databazy")
    print("15. Dot plot - manualne zadanie sekvencii")
    print("16. Dot plot - sekvencie z databazy")
    print("0. Koniec")


def main():
    database = SequenceDatabase()

    while True:
        show_menu()
        choice = input("Vyber moznost: ").strip()

        if choice == "1":
            test_basic_operations()

        elif choice == "2":
            load_sequences_menu(database)

        elif choice == "3":
            list_sequences_menu(database)

        elif choice == "4":
            print_sequence_menu(database)

        elif choice == "5":
            show_database_summary(database)

        elif choice == "6":
            export_all_to_one_file_menu(database)

        elif choice == "7":
            export_individual_files_menu(database)

        elif choice == "8":
            nucleotide_frequency_menu(database)

        elif choice == "9":
            generate_random_mrna_menu()

        elif choice == "10":
            find_orfs_menu(database)

        elif choice == "11":
            cds_translation_menu(database)

        elif choice == "12":
            point_mutation_menu(database)

        elif choice == "13":
            hamming_distance_manual_menu()

        elif choice == "14":
            hamming_distance_database_menu(database)

        elif choice == "15":
            dotplot_manual_menu()

        elif choice == "16":
            dotplot_database_menu(database)

        elif choice == "0":
            print("Koniec programu.")
            break

        else:
            print("Neplatna moznost. Skus znova.")


if __name__ == "__main__":
    main()