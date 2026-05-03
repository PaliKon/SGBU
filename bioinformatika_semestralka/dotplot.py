# dotplot.py

from pathlib import Path
from PIL import Image

from sequence_utils import normalize_sequence


def create_dotplot_matrix(seq1: str, seq2: str, word_size: int = 1) -> list[list[int]]:
    """
    Vytvorí bodovú maticu pre dve sekvencie.

    Hodnota 1 znamená zhodu slov danej dĺžky.
    Hodnota 0 znamená nezhodu.
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)

    if word_size <= 0:
        raise ValueError("Word size musi byt kladne cislo.")

    if word_size > len(seq1) or word_size > len(seq2):
        raise ValueError("Word size nemoze byt vacsi ako dlzka porovnavanych sekvencii.")

    width = len(seq1) - word_size + 1
    height = len(seq2) - word_size + 1

    matrix = []

    for y in range(height):
        row = []

        word2 = seq2[y:y + word_size]

        for x in range(width):
            word1 = seq1[x:x + word_size]

            if word1 == word2:
                row.append(1)
            else:
                row.append(0)

        matrix.append(row)

    return matrix


def save_dotplot_png(
    seq1: str,
    seq2: str,
    output_file: str,
    word_size: int = 1,
    pixel_size: int = 20
) -> None:
    """
    Vytvorí dot plot pre dve sekvencie a uloží ho do PNG súboru.

    Čierny bod = zhoda.
    Biely bod = nezhoda.
    """
    matrix = create_dotplot_matrix(seq1, seq2, word_size)

    matrix_height = len(matrix)
    matrix_width = len(matrix[0])

    image_width = matrix_width * pixel_size
    image_height = matrix_height * pixel_size

    image = Image.new("RGB", (image_width, image_height), "white")
    pixels = image.load()

    for y in range(matrix_height):
        for x in range(matrix_width):
            if matrix[y][x] == 1:
                for dy in range(pixel_size):
                    for dx in range(pixel_size):
                        pixels[
                            x * pixel_size + dx,
                            y * pixel_size + dy
                        ] = (0, 0, 0)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    image.save(output_path)


def dotplot_statistics(seq1: str, seq2: str, word_size: int = 1) -> dict:
    """
    Vráti základné štatistiky dot plotu.
    """
    matrix = create_dotplot_matrix(seq1, seq2, word_size)

    total_cells = len(matrix) * len(matrix[0])
    matches = sum(sum(row) for row in matrix)

    return {
        "matrix_width": len(matrix[0]),
        "matrix_height": len(matrix),
        "total_cells": total_cells,
        "matches": matches,
        "match_percent": (matches / total_cells) * 100 if total_cells > 0 else 0,
    }