# fasta_utils.py

from sequence_utils import normalize_sequence, is_defined_over_alphabet, MRNA_ALPHABET


def read_fasta(file_path: str) -> list[dict]:
    """
    Načíta jeden alebo viac FASTA záznamov zo súboru.

    Výstup:
    [
        {"header": "sample_mRNA_1", "sequence": "AUGGCUUAA"},
        {"header": "sample_mRNA_2", "sequence": "AUGUUUGGCUAG"}
    ]
    """
    records = []
    current_header = None
    current_sequence_parts = []

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequence = normalize_sequence("".join(current_sequence_parts))
                    sequence = sequence.replace("T", "U")

                    if not is_defined_over_alphabet(sequence, MRNA_ALPHABET):
                        raise ValueError(f"Sekvencia '{current_header}' obsahuje neplatné znaky.")

                    records.append({
                        "header": current_header,
                        "sequence": sequence
                    })

                current_header = line[1:].strip()
                current_sequence_parts = []
            else:
                current_sequence_parts.append(line)

    if current_header is not None:
        sequence = normalize_sequence("".join(current_sequence_parts))
        sequence = sequence.replace("T", "U")

        if not is_defined_over_alphabet(sequence, MRNA_ALPHABET):
            raise ValueError(f"Sekvencia '{current_header}' obsahuje neplatné znaky.")

        records.append({
            "header": current_header,
            "sequence": sequence
        })

    return records


def write_fasta(records: list[dict], file_path: str, line_width: int = 70) -> None:
    """
    Uloží jeden alebo viac FASTA záznamov do súboru.
    """
    with open(file_path, "w", encoding="utf-8") as file:
        for index, record in enumerate(records, start=1):
            header = record.get("header") or f"sequence_{index}"
            sequence = normalize_sequence(record["sequence"])

            file.write(f">{header}\n")

            for i in range(0, len(sequence), line_width):
                file.write(sequence[i:i + line_width] + "\n")