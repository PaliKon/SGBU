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

def needleman_wunsch(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_score: int = -2,
    max_alignments: int = 100
) -> dict:
    """
    Implementácia algoritmu Needleman-Wunsch pre globálne zarovnanie.

    Vracia:
    - najlepšie skóre,
    - skórovaciu maticu,
    - všetky najlepšie zarovnania, obmedzené parametrom max_alignments.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    rows = len(seq2) + 1
    cols = len(seq1) + 1

    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    traceback_matrix = [[[] for _ in range(cols)] for _ in range(rows)]

    # Inicializácia prvého riadku
    for col in range(1, cols):
        score_matrix[0][col] = score_matrix[0][col - 1] + gap_score
        traceback_matrix[0][col].append("left")

    # Inicializácia prvého stĺpca
    for row in range(1, rows):
        score_matrix[row][0] = score_matrix[row - 1][0] + gap_score
        traceback_matrix[row][0].append("up")

    # Vyplnenie matice
    for row in range(1, rows):
        for col in range(1, cols):
            char1 = seq1[col - 1]
            char2 = seq2[row - 1]

            if char1 == char2:
                diagonal_score = score_matrix[row - 1][col - 1] + match_score
            else:
                diagonal_score = score_matrix[row - 1][col - 1] + mismatch_score

            up_score = score_matrix[row - 1][col] + gap_score
            left_score = score_matrix[row][col - 1] + gap_score

            best_score = max(diagonal_score, up_score, left_score)
            score_matrix[row][col] = best_score

            if diagonal_score == best_score:
                traceback_matrix[row][col].append("diagonal")

            if up_score == best_score:
                traceback_matrix[row][col].append("up")

            if left_score == best_score:
                traceback_matrix[row][col].append("left")

    alignments = []

    def traceback(row: int, col: int, aligned_seq1: str, aligned_seq2: str):
        """
        Rekurzívne spätne prechádza maticu a skladá všetky optimálne zarovnania.
        """
        if len(alignments) >= max_alignments:
            return

        if row == 0 and col == 0:
            alignments.append({
                "seq1": aligned_seq1,
                "seq2": aligned_seq2,
            })
            return

        for direction in traceback_matrix[row][col]:
            if direction == "diagonal":
                traceback(
                    row - 1,
                    col - 1,
                    seq1[col - 1] + aligned_seq1,
                    seq2[row - 1] + aligned_seq2
                )

            elif direction == "up":
                traceback(
                    row - 1,
                    col,
                    "-" + aligned_seq1,
                    seq2[row - 1] + aligned_seq2
                )

            elif direction == "left":
                traceback(
                    row,
                    col - 1,
                    seq1[col - 1] + aligned_seq1,
                    "-" + aligned_seq2
                )

    traceback(rows - 1, cols - 1, "", "")

    return {
        "score": score_matrix[rows - 1][cols - 1],
        "score_matrix": score_matrix,
        "alignments": alignments,
        "alignment_count": len(alignments),
        "max_alignments_reached": len(alignments) >= max_alignments,
    }


def alignment_match_line(aligned_seq1: str, aligned_seq2: str) -> str:
    """
    Vytvorí pomocný riadok medzi zarovnanými sekvenciami.

    | znamená zhodu.
    medzera znamená nezhodu alebo gap.
    """
    line = []

    for char1, char2 in zip(aligned_seq1, aligned_seq2):
        if char1 == char2 and char1 != "-":
            line.append("|")
        else:
            line.append(" ")

    return "".join(line)

def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
    max_alignments: int = 100
) -> dict:
    """
    Implementácia algoritmu Smith-Waterman pre lokálne zarovnanie.

    Vracia:
    - najlepšie lokálne skóre,
    - skórovaciu maticu,
    - všetky najlepšie lokálne zarovnania, obmedzené parametrom max_alignments.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    rows = len(seq2) + 1
    cols = len(seq1) + 1

    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    traceback_matrix = [[[] for _ in range(cols)] for _ in range(rows)]

    best_score = 0
    best_positions = []

    for row in range(1, rows):
        for col in range(1, cols):
            char1 = seq1[col - 1]
            char2 = seq2[row - 1]

            if char1 == char2:
                diagonal_score = score_matrix[row - 1][col - 1] + match_score
            else:
                diagonal_score = score_matrix[row - 1][col - 1] + mismatch_score

            up_score = score_matrix[row - 1][col] + gap_score
            left_score = score_matrix[row][col - 1] + gap_score

            cell_score = max(0, diagonal_score, up_score, left_score)
            score_matrix[row][col] = cell_score

            if cell_score == 0:
                continue

            if diagonal_score == cell_score:
                traceback_matrix[row][col].append("diagonal")

            if up_score == cell_score:
                traceback_matrix[row][col].append("up")

            if left_score == cell_score:
                traceback_matrix[row][col].append("left")

            if cell_score > best_score:
                best_score = cell_score
                best_positions = [(row, col)]
            elif cell_score == best_score:
                best_positions.append((row, col))

    alignments = []

    def traceback(row: int, col: int, aligned_seq1: str, aligned_seq2: str, end_row: int, end_col: int):
        """
        Rekurzívny traceback pre lokálne zarovnanie.
        Končí, keď narazí na bunku so skóre 0.
        """
        if len(alignments) >= max_alignments:
            return

        if score_matrix[row][col] == 0:
            alignments.append({
                "seq1": aligned_seq1,
                "seq2": aligned_seq2,
                "seq1_start": col + 1,
                "seq1_end": end_col,
                "seq2_start": row + 1,
                "seq2_end": end_row,
            })
            return

        for direction in traceback_matrix[row][col]:
            if direction == "diagonal":
                traceback(
                    row - 1,
                    col - 1,
                    seq1[col - 1] + aligned_seq1,
                    seq2[row - 1] + aligned_seq2,
                    end_row,
                    end_col
                )

            elif direction == "up":
                traceback(
                    row - 1,
                    col,
                    "-" + aligned_seq1,
                    seq2[row - 1] + aligned_seq2,
                    end_row,
                    end_col
                )

            elif direction == "left":
                traceback(
                    row,
                    col - 1,
                    seq1[col - 1] + aligned_seq1,
                    "-" + aligned_seq2,
                    end_row,
                    end_col
                )

    for row, col in best_positions:
        if len(alignments) >= max_alignments:
            break

        traceback(
            row,
            col,
            "",
            "",
            end_row=row,
            end_col=col
        )

    unique_alignments = []
    seen = set()

    for alignment in alignments:
        key = (
            alignment["seq1"],
            alignment["seq2"],
            alignment["seq1_start"],
            alignment["seq1_end"],
            alignment["seq2_start"],
            alignment["seq2_end"],
        )

        if key not in seen:
            seen.add(key)
            unique_alignments.append(alignment)

    return {
        "score": best_score,
        "score_matrix": score_matrix,
        "alignments": unique_alignments,
        "alignment_count": len(unique_alignments),
        "max_alignments_reached": len(alignments) >= max_alignments,
    }