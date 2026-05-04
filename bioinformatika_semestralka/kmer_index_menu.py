# kmer_index_menu.py

from sequence_database import SequenceDatabase


def build_kmer_index_menu(database: SequenceDatabase):
    """
    Úloha 17 - vybudovanie indexu k-tíc.
    """
    print("\n=== 17. Vybudovanie indexu k-tic ===")

    if database.count() == 0:
        print("Databaza je prazdna. Najprv nacitaj FASTA sekvencie.")
        return

    try:
        k = int(input("Zadaj k v rozsahu 3..9: "))
    except ValueError:
        print("Neplatny vstup. k musi byt cele cislo.")
        return

    try:
        database.build_kmer_index(k)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    index = database.kmer_index

    print("\nIndex bol uspesne vytvoreny.")
    print(f"k = {index.k}")
    print(f"Pocet sekvencii v databaze: {database.count()}")
    print(f"Pocet unikatnych k-tic: {index.number_of_unique_kmers()}")
    print(f"Celkovy pocet vyskytov k-tic: {index.number_of_occurrences()}")

    print("\nNajcastejsie k-tice:")
    for kmer, count in index.most_frequent_kmers(limit=10):
        print(f"{kmer}: {count}x")


def search_kmer_menu(database: SequenceDatabase):
    """
    Vyhľadanie konkrétnej k-tice v indexe.
    """
    print("\n=== Vyhladanie k-tice v indexe ===")

    if not database.has_kmer_index():
        print("Index k-tic este nie je vytvoreny.")
        print("Najprv pouzi moznost na vybudovanie indexu.")
        return

    index = database.kmer_index

    kmer = input(f"Zadaj k-ticu dlzky {index.k}: ").strip().upper()

    try:
        results = index.search(kmer)
    except ValueError as error:
        print(f"Chyba: {error}")
        return

    print(f"\nK-tica: {kmer}")
    print(f"Pocet vyskytov: {len(results)}")

    if not results:
        print("K-tica sa v databaze nenasla.")
        return

    max_to_print = 30

    print("\nVyskyty:")
    for result in results[:max_to_print]:
        print(
            f"sekvencia {result['sequence_number']}, "
            f"pozicia {result['position']}, "
            f"header: {result['header']}"
        )

    if len(results) > max_to_print:
        remaining = len(results) - max_to_print
        print(f"... a dalsich {remaining} vyskytov")


def kmer_index_summary_menu(database: SequenceDatabase):
    """
    Vypíše súhrn aktuálneho indexu.
    """
    print("\n=== Suhrn indexu k-tic ===")

    if not database.has_kmer_index():
        print("Index k-tic este nie je vytvoreny.")
        return

    index = database.kmer_index

    print(f"k = {index.k}")
    print(f"Pocet sekvencii v databaze: {database.count()}")
    print(f"Pocet unikatnych k-tic: {index.number_of_unique_kmers()}")
    print(f"Celkovy pocet vyskytov k-tic: {index.number_of_occurrences()}")

    print("\nNajcastejsie k-tice:")
    for kmer, count in index.most_frequent_kmers(limit=10):
        print(f"{kmer}: {count}x")