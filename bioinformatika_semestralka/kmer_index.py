# kmer_index.py

from sequence_utils import normalize_sequence


class KmerIndex:
    """
    Index k-tíc nad databázou mRNA sekvencií.

    Index má tvar:
    {
        "AUG": [
            {"sequence_number": 1, "position": 5, "header": "..."},
            {"sequence_number": 2, "position": 18, "header": "..."}
        ],
        ...
    }

    Pozície sú 1-based.
    """

    def __init__(self, k: int):
        if k < 3 or k > 9:
            raise ValueError("k musi byt v rozsahu 3..9")

        self.k = k
        self.index = {}

    def build(self, records: list[dict]) -> None:
        """
        Vybuduje index nad celou databázou sekvencií.
        """
        self.index = {}

        for sequence_number, record in enumerate(records, start=1):
            self.add_sequence(record, sequence_number)

    def add_sequence(self, record: dict, sequence_number: int) -> None:
        """
        Pridá jednu sekvenciu do existujúceho indexu.
        """
        sequence = normalize_sequence(record["sequence"])
        header = record.get("header") or f"sequence_{sequence_number}"

        if len(sequence) < self.k:
            return

        for i in range(0, len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            position = i + 1

            if kmer not in self.index:
                self.index[kmer] = []

            self.index[kmer].append({
                "sequence_number": sequence_number,
                "position": position,
                "header": header,
            })

    def search(self, kmer: str) -> list[dict]:
        """
        Vyhľadá k-ticu v indexe.
        """
        kmer = normalize_sequence(kmer)

        if len(kmer) != self.k:
            raise ValueError(f"Hladane slovo musi mat dlzku {self.k}.")

        return self.index.get(kmer, [])

    def contains(self, kmer: str) -> bool:
        """
        Zistí, či sa k-tica nachádza v indexe.
        """
        kmer = normalize_sequence(kmer)
        return kmer in self.index

    def number_of_unique_kmers(self) -> int:
        """
        Vráti počet unikátnych k-tíc v indexe.
        """
        return len(self.index)

    def number_of_occurrences(self) -> int:
        """
        Vráti celkový počet výskytov všetkých k-tíc.
        """
        return sum(len(positions) for positions in self.index.values())

    def most_frequent_kmers(self, limit: int = 10) -> list[tuple[str, int]]:
        """
        Vráti najčastejšie k-tice.
        """
        frequencies = [
            (kmer, len(positions))
            for kmer, positions in self.index.items()
        ]

        frequencies.sort(key=lambda item: item[1], reverse=True)

        return frequencies[:limit]