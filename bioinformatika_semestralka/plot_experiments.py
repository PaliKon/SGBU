# plot_experiments.py

import csv
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt


INPUT_DIR = Path("output/experiments")
OUTPUT_DIR = INPUT_DIR / "plots"


def ensure_output_dir():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def read_semicolon_csv(file_path: Path):
    with open(file_path, "r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter=";")
        return list(reader)


def plot_timing():
    rows = read_semicolon_csv(INPUT_DIR / "nw_sw_timing.csv")

    grouped = defaultdict(list)

    for row in rows:
        prefix = int(row["prefix_length"])
        algorithm = row["algorithm"]
        time_ms = float(row["time_ms"])
        grouped[(prefix, algorithm)].append(time_ms)

    prefixes = sorted(set(int(row["prefix_length"]) for row in rows))

    nw_avg = []
    sw_avg = []

    for prefix in prefixes:
        nw_values = grouped[(prefix, "Needleman-Wunsch")]
        sw_values = grouped[(prefix, "Smith-Waterman")]

        nw_avg.append(sum(nw_values) / len(nw_values))
        sw_avg.append(sum(sw_values) / len(sw_values))

    plt.figure(figsize=(8, 5))
    plt.plot(prefixes, nw_avg, marker="o", label="Needleman-Wunsch")
    plt.plot(prefixes, sw_avg, marker="o", label="Smith-Waterman")
    plt.xlabel("Dĺžka prefixu")
    plt.ylabel("Priemerný čas (ms)")
    plt.title("Porovnanie času NW a SW")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "timing_nw_sw.png", dpi=200)
    plt.close()


def plot_scores():
    rows = read_semicolon_csv(INPUT_DIR / "nw_sw_scores.csv")

    pairs = [row["pair"] for row in rows]
    nw_scores = [float(row["nw_score"]) for row in rows]
    sw_scores = [float(row["sw_score"]) for row in rows]

    x = list(range(len(pairs)))
    width = 0.4

    plt.figure(figsize=(10, 5))
    plt.bar([i - width / 2 for i in x], nw_scores, width=width, label="Needleman-Wunsch")
    plt.bar([i + width / 2 for i in x], sw_scores, width=width, label="Smith-Waterman")
    plt.xticks(x, pairs, rotation=45)
    plt.xlabel("Dvojica sekvencií")
    plt.ylabel("Skóre")
    plt.title("Porovnanie skóre NW a SW")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scores_nw_sw.png", dpi=200)
    plt.close()


def plot_kmer_unique():
    rows = read_semicolon_csv(INPUT_DIR / "kmer_index_summary.csv")

    k_values = [int(row["k"]) for row in rows]
    unique_kmers = [int(row["unique_kmers"]) for row in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(k_values, unique_kmers, marker="o")
    plt.xlabel("k")
    plt.ylabel("Počet unikátnych k-tíc")
    plt.title("Počet unikátnych k-tíc podľa hodnoty k")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "kmer_unique.png", dpi=200)
    plt.close()


def plot_kmer_build_time():
    rows = read_semicolon_csv(INPUT_DIR / "kmer_index_summary.csv")

    k_values = [int(row["k"]) for row in rows]
    build_times = [float(row["build_time_ms"]) for row in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(k_values, build_times, marker="o")
    plt.xlabel("k")
    plt.ylabel("Čas vybudovania indexu (ms)")
    plt.title("Čas vybudovania indexu k-tíc")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "kmer_build_time.png", dpi=200)
    plt.close()


def main():
    ensure_output_dir()

    plot_timing()
    plot_scores()
    plot_kmer_unique()
    plot_kmer_build_time()

    print("Grafy boli vytvorené v priečinku:")
    print(OUTPUT_DIR)


if __name__ == "__main__":
    main()