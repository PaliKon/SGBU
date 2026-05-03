# alignment.py

from sequence_utils import normalize_sequence


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Vypočíta Hammingovu vzdialenosť dvoch rovnako dlhých sekvencií.

    Hammingova vzdialenosť je počet pozícií, na ktorých sa dve sekvencie líšia.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    if len(seq1) != len(seq2):
        raise ValueError("Hammingova vzdialenost je definovana iba pre sekvencie rovnakej dlzky.")

    distance = 0

    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            distance += 1

    return distance


def hamming_comparison_details(seq1: str, seq2: str) -> dict:
    """
    Vráti detailné porovnanie dvoch sekvencií:
    - Hammingovu vzdialenosť
    - počet rovnakých pozícií
    - počet rozdielnych pozícií
    - zoznam rozdielnych pozícií

    Pozície sú 1-based.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    if len(seq1) != len(seq2):
        raise ValueError("Hammingova vzdialenost je definovana iba pre sekvencie rovnakej dlzky.")

    differences = []

    for index, (base1, base2) in enumerate(zip(seq1, seq2), start=1):
        if base1 != base2:
            differences.append({
                "position": index,
                "seq1_base": base1,
                "seq2_base": base2,
            })

    return {
        "length": len(seq1),
        "distance": len(differences),
        "matches": len(seq1) - len(differences),
        "differences": differences,
    }