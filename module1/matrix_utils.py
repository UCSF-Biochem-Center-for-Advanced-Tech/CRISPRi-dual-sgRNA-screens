import os
from glob import glob
import pandas as pd


def load_counts_file(path: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=None, names=["sgID_AB", "count"])
    return df.set_index("sgID_AB")["count"]


def merge_counts(count_files: list[str], sgid_order: list[str]) -> pd.DataFrame:
    samples = [os.path.basename(f).split(".")[0] for f in count_files]
    first_series = load_counts_file(count_files[0]).reindex(sgid_order).fillna(0)
    matrix = pd.DataFrame(index=sgid_order)
    matrix[samples[0]] = first_series.values
    for sample, path in zip(samples[1:], count_files[1:]):
        series = load_counts_file(path).reindex(sgid_order).fillna(0)
        matrix[sample] = series.values
    matrix.insert(0, "target", matrix.index.str.split("_").str[0])
    matrix = matrix.fillna(0)
    num_cols = matrix.columns[1:]
    matrix[num_cols] = matrix[num_cols].astype(int)
    return matrix


def write_matrices(output_dir: str, sgid_order: list[str]) -> None:
    patterns = {
        "AB.match.counts.txt": ".AB.match.counts.txt",
        "all.aligned.counts.txt": ".all.aligned.counts.txt",
    }
    for label, suffix in patterns.items():
        paths = sorted(glob(os.path.join(output_dir, f"*{suffix}")))
        if not paths:
            continue
        mat = merge_counts(paths, sgid_order)
        out_path = os.path.join(output_dir, f"counts_{label}.tsv")
        mat.to_csv(out_path, sep="\t")
