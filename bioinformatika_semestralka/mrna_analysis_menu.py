# mrna_analysis_menu.py

from sequence_database import SequenceDatabase

from mrna_analysis import (
    nucleotide_frequency,
    generate_random_mrna,
    find_orfs,
    cds_to_protein,
    point_mutation,
)


def select_sequence(database: SequenceDatabase) -> tuple[int, dict] | None:
    """
    Pomocná funkcia na výber sekvencie z databázy.
    """
    if database.count() == 0:
        print("Databaza je prazdna. Najprv nacitaj sekvencie.")
        return None

    database.list_sequences()

    try:
        number = int(input("\nZadaj cislo sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return None

    if number < 1 or number > database.count():
        print("Neplatne cislo sekvencie.")
        return None

    return number, database.records[number - 1]


def nucleotide_frequency_menu(database: SequenceDatabase):
    """
    Menu pre úlohu 8.
    """
    print("\n=== 8. Frekvencna analyza mRNA sekvencie ===")

    selected = select_sequence(database)

    if selected is None:
        return

    number, record = selected
    sequence = record["sequence"]

    frequencies = nucleotide_frequency(sequence)

    print(f"\nSekvencia {number}: {record.get('header') or f'sequence_{number}'}")
    print(f"Dlzka sekvencie: {len(sequence)}")

    for base, data in frequencies.items():
        print(
            f"{base}: pocet = {data['count']}, "
            f"relativne = {data['relative']:.4f}, "
            f"percent = {data['percent']:.2f}%"
        )


def generate_random_mrna_menu():
    """
    Menu pre úlohu 9.
    """
    print("\n=== 9. Generovanie nahodnej mRNA sekvencie ===")

    try:
        length = int(input("Zadaj dlzku nahodnej mRNA sekvencie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    try:
        sequence = generate_random_mrna(length)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    print("\nVygenerovana mRNA sekvencia:")
    print(sequence)


def find_orfs_menu(database: SequenceDatabase):
    """
    Menu pre úlohu 10.
    """
    print("\n=== 10. Hladanie ORF kandidatov ===")

    selected = select_sequence(database)

    if selected is None:
        return

    number, record = selected
    sequence = record["sequence"]

    orfs = find_orfs(sequence)

    print(f"\nSekvencia {number}: {record.get('header') or f'sequence_{number}'}")
    print(f"Dlzka sekvencie: {len(sequence)}")
    print(f"Pocet najdenych ORF kandidatov: {len(orfs)}")

    if not orfs:
        print("Nenasiel sa ziadny ORF.")
        return

    for index, orf in enumerate(orfs, start=1):
        print(
            f"\nORF {index}: "
            f"start = {orf['start']}, "
            f"end = {orf['end']}, "
            f"frame = {orf['frame']}, "
            f"dlzka = {orf['length']}, "
            f"GC = {orf['gc_content']:.2f}%"
        )

        preview = orf["sequence"][:90]

        if len(orf["sequence"]) > 90:
            preview += "..."

        print(f"Zaciatok ORF: {preview}")


def cds_translation_menu(database: SequenceDatabase):
    """
    Menu pre úlohu 11.
    """
    print("\n=== 11. Oznacenie ORF ako CDS a preklad na protein ===")

    selected = select_sequence(database)

    if selected is None:
        return

    number, record = selected
    sequence = record["sequence"]

    orfs = find_orfs(sequence)

    if not orfs:
        print("Pre zvolenu sekvenciu sa nenasiel ziadny ORF.")
        return

    print(f"\nNajdene ORF kandidaty pre sekvenciu {number}:")

    for index, orf in enumerate(orfs, start=1):
        print(
            f"{index}. start = {orf['start']}, "
            f"end = {orf['end']}, "
            f"frame = {orf['frame']}, "
            f"dlzka = {orf['length']}, "
            f"GC = {orf['gc_content']:.2f}%"
        )

    try:
        selected_orf_index = int(input("\nVyber cislo ORF, ktory oznacis ako CDS: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    if selected_orf_index < 1 or selected_orf_index > len(orfs):
        print("Neplatne cislo ORF.")
        return

    selected_orf = orfs[selected_orf_index - 1]
    protein = cds_to_protein(selected_orf["sequence"])

    print("\nVybrany CDS:")
    print(f"Start: {selected_orf['start']}")
    print(f"End: {selected_orf['end']}")
    print(f"Frame: {selected_orf['frame']}")
    print(f"Dlzka CDS: {selected_orf['length']}")
    print(f"GC obsah: {selected_orf['gc_content']:.2f}%")

    print("\nProtein:")
    print(protein)
    print(f"Dlzka proteinu: {len(protein)} aminokyselin")


def point_mutation_menu(database: SequenceDatabase):
    """
    Menu pre úlohu 12.
    """
    print("\n=== 12. Bodova mutacia mRNA sekvencie ===")

    selected = select_sequence(database)

    if selected is None:
        return

    number, record = selected
    sequence = record["sequence"]

    print(f"\nZvolena sekvencia: {record.get('header') or f'sequence_{number}'}")
    print(f"Dlzka sekvencie: {len(sequence)}")

    try:
        position = int(input("Zadaj poziciu mutacie: "))
    except ValueError:
        print("Neplatny vstup. Zadaj cele cislo.")
        return

    new_base = input("Zadaj novu bazu A/U/C/G: ").strip().upper()

    try:
        mutated_sequence = point_mutation(sequence, position, new_base)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    old_base = sequence[position - 1]

    print("\nMutacia bola vytvorena.")
    print(f"Pozicia: {position}")
    print(f"Povodna baza: {old_base}")
    print(f"Nova baza: {new_base}")

    start = max(0, position - 11)
    end = min(len(sequence), position + 10)

    print("\nOkolie mutacie v povodnej sekvencii:")
    print(sequence[start:end])

    print("\nOkolie mutacie v mutovanej sekvencii:")
    print(mutated_sequence[start:end])

    save_choice = input("\nChces pridat mutovanu sekvenciu do databazy? (a/n): ").strip().lower()

    if save_choice == "a":
        mutated_header = (
            f"{record.get('header') or f'sequence_{number}'} "
            f"| mutation_pos_{position}_{old_base}_to_{new_base}"
        )

        database.add_sequence(mutated_sequence, mutated_header)

        print("Mutovana sekvencia bola pridana do databazy.")
        print(f"Novy pocet sekvencii v databaze: {database.count()}")