# main.py

from sequence_utils import (
    MRNA_ALPHABET,
    are_identical_and_valid,
    codon_to_amino_acid,
    translate_mrna_to_protein,
)

from fasta_utils import read_fasta, write_fasta
from sequence_database import SequenceDatabase


def test_basic_operations():
    seq1 = "AUGGCUUAA"
    seq2 = "AUGGCUUAA"

    print("=== Test 1-3: Zakladne operacie ===")

    print("Test identity a abecedy:")
    print(are_identical_and_valid(seq1, seq2, MRNA_ALPHABET))

    print("\nTest prevodu kodonu:")
    print("AUG ->", codon_to_amino_acid("AUG"))
    print("GCU ->", codon_to_amino_acid("GCU"))
    print("UAA ->", codon_to_amino_acid("UAA"))

    print("\nTest prekladu mRNA na protein:")
    print(seq1, "->", translate_mrna_to_protein(seq1))


def test_fasta_and_database():
    print("\n=== Test 4-7: FASTA a databaza sekvencii ===")

    records = read_fasta("data/actb_human.fasta")

    print("\nNacitane FASTA zaznamy:")
    for record in records:
        print("Header:", record["header"])
        print("Dlzka sekvencie:", len(record["sequence"]))
        print("Prvych 80 baz:", record["sequence"][:80])

    write_fasta(records, "output/exported_sample.fasta")
    print("\nSekvencie boli ulozene do output/exported_actb_human.fasta")

    database = SequenceDatabase()
    database.load_from_fasta("data/actb_human.fasta")

    print("\nZoznam sekvencii v databaze:")
    database.list_sequences()

    print("\nPrvych 200 baz sekvencie cislo 1:")
    seq = database.get_sequence(1)
    print(seq[:200])


def main():
    test_basic_operations()
    test_fasta_and_database()


if __name__ == "__main__":
    main()