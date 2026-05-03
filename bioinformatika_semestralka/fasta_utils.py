# fasta_utils.py

from sequence_utils import normalize_sequence, is_defined_over_alphabet, MRNA_ALPHABET
from pathlib import Path


def read_fasta(file_path: str) -> list[dict]:
    """
    Načíta jeden alebo viac FASTA záznamov zo súboru.

    Pri načítaní sa T konvertuje na U, aby sa so sekvenciami ďalej
    pracovalo ako s mRNA sekvenciami.
    """
    records = []
    current_header = None
    current_sequence_parts = []

    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Súbor je prázdny: {file_path}")

    for line in lines:
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

    if not records:
        raise ValueError(f"Súbor neobsahuje žiadny FASTA záznam: {file_path}")

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

def safe_filename(text: str) -> str:
    """
    Vytvorí bezpečný názov súboru z FASTA hlavičky.
    """
    allowed_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"

    result = []
    for char in text:
        if char in allowed_chars:
            result.append(char)
        elif char in {" ", ".", "|", ":", ",", "(", ")"}:
            result.append("_")

    filename = "".join(result)

    while "__" in filename:
        filename = filename.replace("__", "_")

    filename = filename.strip("_")

    if not filename:
        filename = "sequence"

    return filename[:80]


def write_fasta_individual(records: list[dict], output_directory: str, line_width: int = 70) -> None:
    """
    Uloží každú sekvenciu do samostatného FASTA súboru.
    """
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    for index, record in enumerate(records, start=1):
        header = record.get("header") or f"sequence_{index}"
        sequence = normalize_sequence(record["sequence"])

        filename = f"{index}_{safe_filename(header)}.fasta"
        file_path = output_path / filename

        write_fasta(
            [{"header": header, "sequence": sequence}],
            str(file_path),
            line_width=line_width
        )