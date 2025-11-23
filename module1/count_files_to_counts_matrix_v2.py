import argparse
import os
from glob import glob
import pandas as pd


def load_counts_file(path: str) -> pd.Series:
    """
    Load a single counts file (either *.AB.match.counts.txt or *.all.aligned.counts.txt)
    into a Series keyed by sgID_AB. Assumes two tab-delimited columns: sgID_AB and count.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=["sgID_AB", "count"])
    return df.set_index("sgID_AB")["count"]


def make_count_matrix(count_files: list[str]) -> pd.DataFrame:
    """
    Merge multiple count files into a single matrix with rows = sgID_AB,
    columns = samples. Assumes every file contains the full guide list in the
    same order (as produced by dualguide_fastqgz_to_counts_v2.py).
    """
    if not count_files:
        raise ValueError("No count files found to merge.")

    samples = [os.path.basename(f).split(".")[0] for f in count_files]

    # Use the first file to establish row order/index
    first_series = load_counts_file(count_files[0])
    matrix = pd.DataFrame(index=first_series.index)
    matrix[samples[0]] = first_series.values

    for sample, path in zip(samples[1:], count_files[1:]):
        series = load_counts_file(path)
        # Align by index to be safe, but preserve the established order
        matrix[sample] = series.reindex(matrix.index).fillna(0).astype(int)

    matrix.insert(0, "target", matrix.index.str.split("_").str[0])
    matrix = matrix.fillna(0).astype(int)
    return matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Merge per-sample count files into a counts matrix. "
            "Works with either *.AB.match.counts.txt (UMI-deduped) or "
            "*.all.aligned.counts.txt (read-level) files."
        )
    )
    parser.add_argument("Count_Files_Path", help="Directory containing count files.")
    parser.add_argument("Out_File_Path", help="Directory to write the merged matrix.")
    parser.add_argument(
        "--suffix",
        default=".AB.match.counts.txt",
        help="Filename suffix to select count files (default: %(default)s).",
    )
    args = parser.parse_args()

    counts_dir = args.Count_Files_Path.rstrip("/")
    out_dir = args.Out_File_Path.rstrip("/")

    os.makedirs(out_dir, exist_ok=True)

    count_files = sorted(glob(os.path.join(counts_dir, f"*{args.suffix}")))
    if not count_files:
        raise FileNotFoundError(f"No files found matching *{args.suffix} in {counts_dir}")

    matrix = make_count_matrix(count_files)
    out_path = os.path.join(out_dir, "counts.txt")
    matrix.to_csv(out_path, sep="\t")
