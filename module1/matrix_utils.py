import os
from glob import glob
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message=".*ImplicitModificationWarning.*")
warnings.filterwarnings("ignore", message=".*Transforming to str index.*")
warnings.filterwarnings("ignore", message=".*PerformanceWarning: DataFrame is highly fragmented.*")


def load_counts_file(path: str) -> pd.Series:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["sgID_AB", "count"],
        dtype={"sgID_AB": str},
        low_memory=False,
    )
    return df.set_index("sgID_AB")["count"]


def merge_counts(count_files: list[str], sgid_order: list[str] | None, include_umi: bool = False) -> pd.DataFrame:
    samples = [os.path.basename(f).split(".")[0] for f in count_files]
    if sgid_order is not None:
        idx = pd.Index(sgid_order, dtype=str)
    else:
        # union of all indices
        all_idx = set()
        if include_umi:
            for path in count_files:
                col0 = pd.read_csv(
                    path,
                    sep="\t",
                    usecols=[0],
                    header=0,
                    dtype=str,
                    low_memory=False,
                ).iloc[:, 0]
                all_idx.update(col0.dropna().astype(str).tolist())
        else:
            for path in count_files:
                all_idx.update(load_counts_file(path).index.tolist())
        idx = pd.Index(sorted(all_idx, key=str))

    if include_umi:
        # rows: sgID_AB, cols: sample-UMI (wide)
        umi_cols = set()
        df_list = []
        for sample, path in zip(samples, count_files):
            df = pd.read_csv(path, sep="\t", low_memory=False)
            # Drop any accidental index columns
            df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed:")]]
            df = df.rename(columns={df.columns[0]: "sgID_AB"})
            df = df.set_index("sgID_AB")
            df.columns = [f"{sample}-{c}" for c in df.columns]
            umi_cols.update(df.columns.tolist())
            df_list.append(df)
        umi_cols = sorted(umi_cols)
        merged = pd.DataFrame(index=idx)
        # concat all columns at once to avoid fragmentation
        merged = pd.concat([df.reindex(idx) for df in df_list], axis=1)
        merged = merged.reindex(columns=umi_cols).fillna(0).astype(int)
        merged.insert(0, "target", merged.index.str.split("_").str[0])
        return merged
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
