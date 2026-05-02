# sequence_utils.py

MRNA_ALPHABET = {"A", "U", "C", "G"}

GENETIC_CODE = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C",
    "UGA": "*", "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def normalize_sequence(sequence: str) -> str:
    """
    Odstráni medzery, nové riadky a prevedie sekvenciu na veľké písmená.
    """
    return sequence.replace(" ", "").replace("\n", "").replace("\r", "").upper()


def is_defined_over_alphabet(sequence: str, alphabet: set[str]) -> bool:
    """
    Skontroluje, či je sekvencia definovaná iba nad povolenou abecedou.
    """
    sequence = normalize_sequence(sequence)
    return all(char in alphabet for char in sequence)


def are_identical_and_valid(seq1: str, seq2: str, alphabet: set[str]) -> bool:
    """
    Skontroluje, či sú dve sekvencie identické a zároveň platné nad danou abecedou.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    return (
        seq1 == seq2
        and is_defined_over_alphabet(seq1, alphabet)
        and is_defined_over_alphabet(seq2, alphabet)
    )


def codon_to_amino_acid(codon: str) -> str:
    """
    Prevedie kodón z 3 mRNA báz na aminokyselinu.
    Stop kodón je označený znakom '*'.
    """
    codon = normalize_sequence(codon)

    if len(codon) != 3:
        raise ValueError("Kodón musí mať presne 3 nukleotidy.")

    if not is_defined_over_alphabet(codon, MRNA_ALPHABET):
        raise ValueError("Kodón obsahuje neplatné znaky. Povolené sú A, U, C, G.")

    return GENETIC_CODE[codon]


def translate_mrna_to_protein(mrna: str) -> str:
    """
    Preloží mRNA na proteín od prvého štart kodónu AUG po stop-kodón.
    Stop kodón sa do výsledného proteínu nezapisuje.
    """
    mrna = normalize_sequence(mrna)

    if not is_defined_over_alphabet(mrna, MRNA_ALPHABET):
        raise ValueError("mRNA obsahuje neplatné znaky. Povolené sú A, U, C, G.")

    start_index = mrna.find("AUG")

    if start_index == -1:
        return ""

    protein = ""

    for i in range(start_index, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]
        amino_acid = codon_to_amino_acid(codon)

        if amino_acid == "*":
            break

        protein += amino_acid

    return protein