# mrna_analysis.py

import random

from sequence_utils import (
    MRNA_ALPHABET,
    normalize_sequence,
    is_defined_over_alphabet,
    translate_mrna_to_protein,
)


START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}


def nucleotide_frequency(sequence: str) -> dict:
    """
    Vypočíta absolútnu a relatívnu frekvenciu báz A, U, C, G.
    """
    sequence = normalize_sequence(sequence)

    if not is_defined_over_alphabet(sequence, MRNA_ALPHABET):
        raise ValueError("Sekvencia obsahuje neplatné znaky. Povolené sú A, U, C, G.")

    if len(sequence) == 0:
        raise ValueError("Sekvencia je prázdna.")

    result = {}

    for base in ["A", "U", "C", "G"]:
        count = sequence.count(base)
        relative = count / len(sequence)

        result[base] = {
            "count": count,
            "relative": relative,
            "percent": relative * 100,
        }

    return result


def generate_random_mrna(length: int) -> str:
    """
    Vygeneruje náhodnú mRNA sekvenciu zadanej dĺžky.
    """
    if length <= 0:
        raise ValueError("Dĺžka sekvencie musí byť kladné číslo.")

    bases = ["A", "U", "C", "G"]
    return "".join(random.choice(bases) for _ in range(length))


def gc_content(sequence: str) -> float:
    """
    Vypočíta GC obsah sekvencie v percentách.
    """
    sequence = normalize_sequence(sequence)

    if len(sequence) == 0:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def find_orfs(sequence: str) -> list[dict]:
    """
    Nájde všetkých kandidátov na ORF v mRNA sekvencii.

    ORF začína štart kodónom AUG a končí stop-kodónom
    UAA, UAG alebo UGA v rovnakom čítacom rámci.

    Pozície start/end sú 1-based, teda prvá báza má pozíciu 1.
    """
    sequence = normalize_sequence(sequence)

    if not is_defined_over_alphabet(sequence, MRNA_ALPHABET):
        raise ValueError("Sekvencia obsahuje neplatné znaky. Povolené sú A, U, C, G.")

    orfs = []

    for frame in range(3):
        i = frame

        while i <= len(sequence) - 3:
            codon = sequence[i:i + 3]

            if codon == START_CODON:
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j + 3]

                    if stop_codon in STOP_CODONS:
                        orf_sequence = sequence[i:j + 3]

                        orfs.append({
                            "start": i + 1,
                            "end": j + 3,
                            "frame": frame + 1,
                            "length": len(orf_sequence),
                            "sequence": orf_sequence,
                            "gc_content": gc_content(orf_sequence),
                        })

                        break

            i += 3

    return orfs


def cds_to_protein(orf_sequence: str) -> str:
    """
    Označený ORF berieme ako CDS a preložíme ho na proteín.
    """
    return translate_mrna_to_protein(orf_sequence)


def point_mutation(sequence: str, position: int, new_base: str) -> str:
    """
    Vytvorí bodovú mutáciu v mRNA.

    position je 1-based, teda prvá báza má pozíciu 1.
    """
    sequence = normalize_sequence(sequence)
    new_base = normalize_sequence(new_base)

    if not is_defined_over_alphabet(sequence, MRNA_ALPHABET):
        raise ValueError("Sekvencia obsahuje neplatné znaky. Povolené sú A, U, C, G.")

    if new_base not in MRNA_ALPHABET:
        raise ValueError("Nová báza musí byť jedna z A, U, C, G.")

    if position < 1 or position > len(sequence):
        raise ValueError("Pozícia mutácie je mimo rozsahu sekvencie.")

    index = position - 1
    return sequence[:index] + new_base + sequence[index + 1:]