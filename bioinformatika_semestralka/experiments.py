# experiments.py

import csv
import time
from itertools import combinations
from pathlib import Path

from sequence_database import SequenceDatabase
from alignment import needleman_wunsch, smith_waterman
from kmer_index import KmerIndex


DATA_DIR = "data"
OUTPUT_DIR = "output/experiments"

PREFIX_LENGTHS = [50, 100, 200, 500]
K_VALUES = [3, 5, 7, 9]

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

MAX_ALIGNMENTS_FOR_EXPERIMENT = 1


def ensure_output_dir() -> Path:
    output_path = Path(OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def load_experiment_database() -> SequenceDatabase:
    database = SequenceDatabase()
    database.load_from_directory(DATA_DIR)

    if database.count() < 5:
        print("Upozornenie: V databaze je menej ako 5 sekvencii.")
        print("Skontroluj, ci mas v data/real presne 5 platnych FASTA suborov.")

    return database


def write_dataset_summary(database: SequenceDatabase, output_path: Path) -> None:
    file_path = output_path / "dataset_summary.csv"

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow(["sequence_number", "header", "length", "first_80_bases"])

        for index, record in enumerate(database.records, start=1):
            sequence = record["sequence"]
            writer.writerow([
                index,
                record.get("header") or f"sequence_{index}",
                len(sequence),
                sequence[:80],
            ])

    print(f"Dataset summary ulozeny do: {file_path}")


def measure_runtime_ms(function, *args, **kwargs) -> tuple[float, object]:
    start = time.perf_counter()
    result = function(*args, **kwargs)
    end = time.perf_counter()

    elapsed_ms = (end - start) * 1000
    return elapsed_ms, result


def run_nw_sw_score_experiment(database: SequenceDatabase, output_path: Path, prefix_length: int = 500) -> None:
    """
    Porovna NW a SW skore na dvojiciach sekvencii pri jednom prefixe.
    """
    file_path = output_path / "nw_sw_scores.csv"

    records = database.records[:5]

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow([
            "pair",
            "seq1_header",
            "seq2_header",
            "prefix_length",
            "match",
            "mismatch",
            "gap",
            "nw_score",
            "sw_score",
            "nw_alignment_length",
            "sw_alignment_length",
        ])

        for (i, record1), (j, record2) in combinations(enumerate(records, start=1), 2):
            seq1 = record1["sequence"][:prefix_length]
            seq2 = record2["sequence"][:prefix_length]

            nw_result = needleman_wunsch(
                seq1,
                seq2,
                match_score=MATCH_SCORE,
                mismatch_score=MISMATCH_SCORE,
                gap_score=GAP_SCORE,
                max_alignments=MAX_ALIGNMENTS_FOR_EXPERIMENT,
            )

            sw_result = smith_waterman(
                seq1,
                seq2,
                match_score=MATCH_SCORE,
                mismatch_score=MISMATCH_SCORE,
                gap_score=GAP_SCORE,
                max_alignments=MAX_ALIGNMENTS_FOR_EXPERIMENT,
            )

            nw_alignment_length = 0
            sw_alignment_length = 0

            if nw_result["alignments"]:
                nw_alignment_length = len(nw_result["alignments"][0]["seq1"])

            if sw_result["alignments"]:
                sw_alignment_length = len(sw_result["alignments"][0]["seq1"])

            writer.writerow([
                f"{i}-{j}",
                record1.get("header") or f"sequence_{i}",
                record2.get("header") or f"sequence_{j}",
                prefix_length,
                MATCH_SCORE,
                MISMATCH_SCORE,
                GAP_SCORE,
                nw_result["score"],
                sw_result["score"],
                nw_alignment_length,
                sw_alignment_length,
            ])

    print(f"NW/SW skore ulozene do: {file_path}")


def run_timing_experiment(database: SequenceDatabase, output_path: Path) -> None:
    """
    Zmeria cas NW a SW pre rozne prefixy.
    """
    file_path = output_path / "nw_sw_timing.csv"

    records = database.records[:5]

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow([
            "pair",
            "prefix_length",
            "algorithm",
            "time_ms",
            "score",
            "matrix_cells",
        ])

        for prefix_length in PREFIX_LENGTHS:
            for (i, record1), (j, record2) in combinations(enumerate(records, start=1), 2):
                seq1 = record1["sequence"][:prefix_length]
                seq2 = record2["sequence"][:prefix_length]

                matrix_cells = (len(seq1) + 1) * (len(seq2) + 1)

                nw_time_ms, nw_result = measure_runtime_ms(
                    needleman_wunsch,
                    seq1,
                    seq2,
                    match_score=MATCH_SCORE,
                    mismatch_score=MISMATCH_SCORE,
                    gap_score=GAP_SCORE,
                    max_alignments=MAX_ALIGNMENTS_FOR_EXPERIMENT,
                )

                writer.writerow([
                    f"{i}-{j}",
                    prefix_length,
                    "Needleman-Wunsch",
                    round(nw_time_ms, 4),
                    nw_result["score"],
                    matrix_cells,
                ])

                sw_time_ms, sw_result = measure_runtime_ms(
                    smith_waterman,
                    seq1,
                    seq2,
                    match_score=MATCH_SCORE,
                    mismatch_score=MISMATCH_SCORE,
                    gap_score=GAP_SCORE,
                    max_alignments=MAX_ALIGNMENTS_FOR_EXPERIMENT,
                )

                writer.writerow([
                    f"{i}-{j}",
                    prefix_length,
                    "Smith-Waterman",
                    round(sw_time_ms, 4),
                    sw_result["score"],
                    matrix_cells,
                ])

                print(
                    f"Prefix {prefix_length}, par {i}-{j}: "
                    f"NW {nw_time_ms:.2f} ms, SW {sw_time_ms:.2f} ms"
                )

    print(f"Casove merania ulozene do: {file_path}")


def run_kmer_index_experiment(database: SequenceDatabase, output_path: Path) -> None:
    """
    Otestuje index k-tic pre rozne hodnoty k.
    """
    file_path = output_path / "kmer_index_summary.csv"

    records = database.records[:5]

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow([
            "k",
            "build_time_ms",
            "unique_kmers",
            "total_occurrences",
            "most_frequent_kmer",
            "most_frequent_count",
        ])

        for k in K_VALUES:
            index = KmerIndex(k)

            build_time_ms, _ = measure_runtime_ms(index.build, records)

            most_frequent = index.most_frequent_kmers(limit=1)

            if most_frequent:
                most_frequent_kmer, most_frequent_count = most_frequent[0]
            else:
                most_frequent_kmer, most_frequent_count = "", 0

            writer.writerow([
                k,
                round(build_time_ms, 4),
                index.number_of_unique_kmers(),
                index.number_of_occurrences(),
                most_frequent_kmer,
                most_frequent_count,
            ])

            print(
                f"k={k}: build {build_time_ms:.2f} ms, "
                f"unikatne={index.number_of_unique_kmers()}, "
                f"vyskyty={index.number_of_occurrences()}"
            )

    print(f"Vysledky indexu k-tic ulozene do: {file_path}")


def create_simple_summary_txt(output_path: Path) -> None:
    file_path = output_path / "experiment_info.txt"

    with open(file_path, "w", encoding="utf-8") as file:
        file.write("Experimenty pre semestralnu pracu BZGI\n")
        file.write("======================================\n\n")
        file.write(f"Data directory: {DATA_DIR}\n")
        file.write(f"Output directory: {OUTPUT_DIR}\n")
        file.write(f"Prefix lengths: {PREFIX_LENGTHS}\n")
        file.write(f"K values: {K_VALUES}\n")
        file.write(f"Scoring: match={MATCH_SCORE}, mismatch={MISMATCH_SCORE}, gap={GAP_SCORE}\n\n")
        file.write("Teoreticka zlozitost:\n")
        file.write("Needleman-Wunsch: O(n*m), pamat O(n*m)\n")
        file.write("Smith-Waterman: O(n*m), pamat O(n*m)\n")
        file.write("Index k-tic: O(N), kde N je celkova dlzka sekvencii\n")

    print(f"Info subor ulozeny do: {file_path}")


def main():
    output_path = ensure_output_dir()

    print("Nacitavam experimentálne data...")
    database = load_experiment_database()

    print(f"Pocet nacitanych sekvencii: {database.count()}")

    write_dataset_summary(database, output_path)
    run_nw_sw_score_experiment(database, output_path, prefix_length=500)
    run_timing_experiment(database, output_path)
    run_kmer_index_experiment(database, output_path)
    create_simple_summary_txt(output_path)

    print("\nExperimenty dokoncene.")
    print(f"Vysledky najdes v priecinku: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()