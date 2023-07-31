#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional
import argparse
import json
from pathlib import Path
import sys


class ValidationError(Exception):
    pass


def colored_print(text, color, newline=True):
    colors = {
        "red": "\033[91m",
        "green": "\033[92m",
        "yellow": "\033[93m",
        "blue": "\033[94m",
        "purple": "\033[95m",
        "cyan": "\033[96m",
        "white": "\033[97m",
        "black": "\033[98m",
        "bold": "\033[1m",
        "underline": "\033[4m",
        "end": "\033[0m"
    }
    if color not in colors:
        raise ValueError(f"Color {color} not found")
    print(colors[color] + text + colors["end"], end="\n" if newline else "")


def filter_cpg_rows(filepath, lines: list[str], amp: int) -> list[list[int]]:
    amp_cpg_sites = [16, 15, 11, 5, 8]
    char_to_num = {".": 0, "G": 1}
    filtered = []
    for i, line in enumerate(lines):
        if any([c != "" for c in line.strip("\n") if c not in char_to_num.keys()]):
            continue
        cpg_sites = [char_to_num[c] for c in line.strip("\n")]
        if len(cpg_sites) != amp_cpg_sites[amp - 1]:
            colored_print(f"WARNING: amplicon {amp} expects {amp_cpg_sites[amp - 1]} "
                          f"but got {len(cpg_sites)} CPG sites"
                          f"line {i+1} of {filepath}", "yellow")
            continue
        filtered.append(cpg_sites)
    return filtered


def validate_aligned_file(filepath: str, amp: int) -> Optional[list[list[int]]]:
    validation_strs = [
        "4181151171181191201201211221241241241241251261",
        "3141415181101131141141151161191191201231",
        "618181919191101131151221251",
        "317181191201",
        "314171101101131161191"
    ]
    expected_amplicon = amp
    with open(filepath, "r") as f:
        lines = f.readlines()
    if len(lines) == 0:
        raise ValidationError(f"{filepath} is empty")
    validated_line = lines[0].strip()
    if validation_strs[expected_amplicon - 1] != validated_line:
        raise ValidationError(f"error with amplicon {amp}: {validated_line} != "
                              f"{validation_strs[expected_amplicon - 1]} ")
    return filter_cpg_rows(filepath, lines, amp)


def create_heatmap(arr: np.ndarray, filepath: Path):
    width = 2
    height = 5
    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(arr, cmap="hot_r", aspect="auto", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    border_width = 1.5
    ax.spines["left"].set_linewidth(border_width)
    ax.spines["right"].set_linewidth(border_width)
    ax.spines["bottom"].set_linewidth(border_width)
    ax.spines["top"].set_linewidth(border_width)

    image_format = "svg"

    fig.savefig(filepath, format=image_format, dpi=1200)


def save_stats(arr: np.ndarray, filepath: Path, amp: int):
    methylated_count = np.sum(arr)
    data = {
        "percent_methylation": (methylated_count / (arr.shape[0] * arr.shape[1])) * 100,
        "count_methylation": int(methylated_count),
        "total_sites": int(arr.shape[0] * arr.shape[1]),

    }
    # fully methylated, partially methylated and un-methylated are only classifications
    # that exist on amp3. These terms don't mean strictly what they say.
    # fully methylated 11 / 11 CPG sites methylated
    # partially methylated 3-10 / 11 CPG sites methylated
    # un-methylated 0-2 / 11 CPG sites methylated
    if amp == 3:
        fully_methylated_count = int(arr[np.sum(arr, axis=1) == 11].shape[0])
        partially_methylated_count = int(arr[(np.sum(arr, axis=1) >= 3) &
                                         (np.sum(arr, axis=1) <= 10)].shape[0])
        unmethylated_count = int(arr[np.sum(arr, axis=1) <= 2].shape[0])
        data["percent_fully_methylated"] = (fully_methylated_count / arr.shape[0]) * 100
        data["percent_partially_methylated"] = (partially_methylated_count / arr.shape[0]) * 100
        data["percent_unmethylated"] = (unmethylated_count / arr.shape[0]) * 100
        data["total_rows"] = int(arr.shape[0])
        data["count_fully_methylated"] = int(fully_methylated_count)
        data["count_partially_methylated"] = int(partially_methylated_count)
        data["count_unmethylated"] = int(unmethylated_count)
    with open(filepath, "w") as f:
        json.dump(data, f)


def as_ndarray(arr: list[list[int]], count: int) -> Optional[np.ndarray]:
    arr = np.array(arr)
    arr = arr[:count]
    row_sums = np.sum(arr, axis=1)
    sorted_indices = np.argsort(row_sums)
    sorted_arr = arr[sorted_indices]
    return sorted_arr


def stats(args):
    try:
        arr = validate_aligned_file(args.aligned_file, args.amplicon)
    except ValidationError as e:
        colored_print("FAILED", "red")
        colored_print(str(e), "red")
        sys.exit(1)
    arr = as_ndarray(arr, args.selection_count)
    create_heatmap(arr, Path(args.output_path, f"figure{args.amplicon}.svg"))
    save_stats(arr, Path(args.output_path, f"stats{args.amplicon}.json"), args.amplicon)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("aligned_file", help="Alignment file to be processed")
    parser.add_argument("-n", "--selection-count", default=300, type=int, help="Selection count")
    parser.add_argument("-o", "--output-path", help="Directory for the stats")
    parser.add_argument("-a", "--amplicon", type=int,  help="Amplicon being processed")
    stats(parser.parse_args())
