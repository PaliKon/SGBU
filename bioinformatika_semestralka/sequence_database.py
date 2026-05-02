# sequence_database.py

from fasta_utils import read_fasta


class SequenceDatabase:
    """
    Jednoduchá databáza mRNA sekvencií v pamäti.
    """

    def __init__(self):
        self.records = []

    def load_from_fasta(self, file_path: str) -> None:
        """
        Načíta viacero mRNA sekvencií zo súboru FASTA do pamäte.
        """
        loaded_records = read_fasta(file_path)
        self.records.extend(loaded_records)

    def add_sequence(self, sequence: str, header: str | None = None) -> None:
        """
        Pridá jednu sekvenciu manuálne do databázy.
        """
        record = {
            "header": header,
            "sequence": sequence
        }
        self.records.append(record)

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

    def get_sequence(self, sequence_number: int) -> str:
        """
        Vráti sekvenciu podľa poradového čísla.
        Číslovanie začína od 1.
        """
        if sequence_number < 1 or sequence_number > len(self.records):
            raise IndexError("Neplatné číslo sekvencie.")

        return self.records[sequence_number - 1]["sequence"]