import os
from glob import glob
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message=".*ImplicitModificationWarning.*")


def load_counts_file(path: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=None, names=["sgID_AB", "count"], dtype={"sgID_AB": str})
    return df.set_index("sgID_AB")["count"]


def merge_counts(count_files: list[str], sgid_order: list[str] | None, include_umi: bool = False) -> pd.DataFrame:
    samples = [os.path.basename(f).split(".")[0] for f in count_files]
    if sgid_order is not None:
        idx = pd.Index(sgid_order, dtype=str)
    else:
        # union of all indices
        all_idx = set()
        for path in count_files:
            all_idx.update(load_counts_file(path).index.tolist())
        idx = pd.Index(sorted(all_idx))

    if include_umi:
        # rows: sgID_AB, cols: sample-UMI
        matrix = pd.DataFrame(index=idx)
        for sample, path in zip(samples, count_files):
            series = load_counts_file(path)
            # split guide and UMI suffix if present
            guide_part = series.index.to_series().astype(str).str.split("++", n=1, expand=True, regex=False)
            if guide_part.shape[1] == 1:
                guides = guide_part[0].values
                umis = pd.Series([""] * len(guide_part))
            else:
                guides = guide_part[0].values
                umis = guide_part[1]
            col_names = [f"{sample}-{u}" if u else sample for u in umis]
            df_sample = pd.DataFrame({"guide": guides, "col": col_names, "count": series.values})
            for col, subdf in df_sample.groupby("col"):
                vals = pd.Series(subdf["count"].values, index=subdf["guide"]).reindex(idx).fillna(0).astype(int)
                matrix[col] = vals.values
        matrix = matrix.fillna(0)
        matrix.insert(0, "target", matrix.index.str.split("_").str[0])
        return matrix
    else:
        matrix = pd.DataFrame(index=idx)
        for sample, path in zip(samples, count_files):
            series = load_counts_file(path).reindex(idx).fillna(0).astype(int)
            matrix[sample] = series.values

        if sgid_order is not None:
            matrix.insert(0, "target", matrix.index.str.split("_").str[0])
        matrix = matrix.fillna(0)
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
        use_order = None if "AB.match" in label else sgid_order
        mat = merge_counts(paths, use_order, include_umi="AB.match" in label)
        out_path = os.path.join(output_dir, f"counts_{label}.tsv")
        mat.to_csv(out_path, sep="\t")
