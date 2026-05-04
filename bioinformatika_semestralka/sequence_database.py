# sequence_database.py

from fasta_utils import read_fasta
from pathlib import Path
from kmer_index import KmerIndex


class SequenceDatabase:
    """
    Jednoduchá databáza mRNA sekvencií v pamäti.
    """

    def __init__(self):
        self.records = []
        self.kmer_index = None

    def load_from_fasta(self, file_path: str) -> None:
        """
        Načíta viacero mRNA sekvencií zo súboru FASTA do pamäte.
        Ak už existuje k-mer index, automaticky ho aktualizuje.
        """
        loaded_records = read_fasta(file_path)

        for record in loaded_records:
            self.add_sequence(record["sequence"], record.get("header"))

    def load_from_directory(self, directory_path: str) -> None:
        """
        Načíta všetky FASTA súbory z priečinka do databázy.
        Podporované prípony: .fasta, .fa, .fna
        """
        directory = Path(directory_path)

        if not directory.exists():
            raise FileNotFoundError(f"Priečinok neexistuje: {directory_path}")

        if not directory.is_dir():
            raise NotADirectoryError(f"Zadaná cesta nie je priečinok: {directory_path}")

        fasta_extensions = {".fasta", ".fa", ".fna"}

        fasta_files = sorted(
            file_path for file_path in directory.iterdir()
            if file_path.is_file() and file_path.suffix.lower() in fasta_extensions
        )

        if not fasta_files:
            print(f"V priečinku {directory_path} sa nenašli žiadne FASTA súbory.")
            return

        for file_path in fasta_files:
            print(f"Načítavam: {file_path}")

            try:
                self.load_from_fasta(str(file_path))
            except ValueError as error:
                print(f"Chyba pri načítaní súboru {file_path}: {error}")

    def add_sequence(self, sequence: str, header: str | None = None) -> None:
        """
        Pridá jednu sekvenciu manuálne do databázy.
        Ak už existuje k-mer index, automaticky ho aktualizuje.
        """
        record = {
            "header": header,
            "sequence": sequence
        }

        self.records.append(record)

        if self.kmer_index is not None:
            sequence_number = len(self.records)
            self.kmer_index.add_sequence(record, sequence_number)

    def list_sequences(self) -> None:
        """
        Vypíše zoznam načítaných mRNA sekvencií.
        Ak hlavička nie je k dispozícii, vypíše iba poradové číslo.
        """
        if not self.records:
            print("Databáza je prázdna.")
            return

        for index, record in enumerate(self.records, start=1):
            header = record.get("header")

            if header:
                print(f"{index}. {header}")
            else:
                print(f"{index}.")

    def print_sequence(self, sequence_number: int) -> None:
        """
        Vypíše konkrétnu mRNA sekvenciu podľa poradového čísla.
        Číslovanie začína od 1.
        """
        if sequence_number < 1 or sequence_number > len(self.records):
            print("Neplatné číslo sekvencie.")
            return

        record = self.records[sequence_number - 1]
        header = record.get("header") or f"sequence_{sequence_number}"
        sequence = record["sequence"]

        print(f">{header}")
        print(sequence)

    def count(self) -> int:
        """
        Vráti počet sekvencií v databáze.
        """
        return len(self.records)

    def get_sequence(self, sequence_number: int) -> str:
        """
        Vráti sekvenciu podľa poradového čísla.
        Číslovanie začína od 1.
        """
        if sequence_number < 1 or sequence_number > len(self.records):
            raise IndexError("Neplatné číslo sekvencie.")

        return self.records[sequence_number - 1]["sequence"]
    
    def build_kmer_index(self, k: int) -> None:
        """
        Vybuduje index k-tíc nad aktuálnou databázou.
        """
        self.kmer_index = KmerIndex(k)
        self.kmer_index.build(self.records)

    def has_kmer_index(self) -> bool:
        """
        Zistí, či je index k-tíc vytvorený.
        """
        return self.kmer_index is not None