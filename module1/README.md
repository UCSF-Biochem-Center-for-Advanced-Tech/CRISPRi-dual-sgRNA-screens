# Module 1: Dual sgRNA FASTQ to Counts (v2)

This module processes paired-end dual-sgRNA FASTQ files into count tables with optional UMI-aware deduplication. It supersedes the original scripts with improved validation, pair-aware disambiguation, UMI handling, and richer stats reporting.

## Key Features (v2)
- Pair-aware matching: rescues ambiguous A/B hits when the pair uniquely matches the library.
- UMI handling:
  - Uses a UMI whitelist when provided (`--UMI_Table`).
  - If UMIs are present in read headers but no UMI table is provided, accepts raw UMIs (no whitelist).
  - Produces read-level counts (`*.all.aligned.counts.txt`) and UMI-deduped molecule counts (`*.AB.match.counts.txt`).
- Robust input validation (guide table columns, file pairing, output directory).
- Detailed stats to stdout; progress to stderr; optional text/JSON stats files.
- Parameterized test subset, quality threshold, UMI mismatch option, and buffered FASTQ reading.

## Requirements
- Python 3.8+ (tested on 3.13)
- pandas
- No other external dependencies required; `numba` is optional but not used by default.

Install dependencies:
```bash
pip install -r module1/requirements.txt
```

## Scripts
- `dualguide_fastqgz_to_counts_v2.py` — main FASTQ→counts pipeline (v2).
- `count_files_to_counts_matrix_v2.py` — merges per-sample count files into a matrix (works with either read-level or UMI-deduped counts).
- Legacy scripts (`dualguide_fastqgz_to_counts.py`, `count_files_to_counts_matrix.py`, UMI variants) are retained but not recommended.

## Input FASTQ assumptions
- FASTQs are produced by `bcl-convert` (or equivalent) with R1/R2 only; I1 is not required as a separate file.
- The UMI (from Index 1) is included in the read IDs, e.g. `@...:1080:TCAGTCGA 1:N:0:ACTCGNTA`. The pipeline parses the UMI from the header; no separate I1 file is needed.

## Usage: FASTQ to Counts (v2)
```bash
python module1/dualguide_fastqgz_to_counts_v2.py \\
  [--UMI_Table UMI_sequences.csv] \\
  [--test] [--test-lines N] \\
  [--low-quality-threshold Q] \\
  [--umi-allow-one-mismatch] \\
  [--write-stats-file] [--write-stats-json] \\
  [--write-offlibrary] [--write-anndata] [--write-count-matrix] \\
  Guide_Table.csv Out_File_Path Seq_File_Names...
```

### Required inputs
- `Guide_Table.csv`: must contain columns `sgID_AB`, `sgID_A`, `protospacer_A`, `sgID_B`, `protospacer_B`. Protospacers are trimmed internally to match 19 bp reads (R1 uses positions 1:20; R2 is reverse-complemented/trimmed).
- `Seq_File_Names`: FASTQ/FA patterns; expects Illumina-style paired R1/R2 names. Multiple patterns allowed (wildcards ok).
- `Out_File_Path`: directory for outputs.

### Optional UMI whitelist
- `--UMI_Table`: CSV with column `UMI`. If omitted and UMIs are present in headers, raw UMIs are accepted without whitelist.

### Useful flags
- `--test`: process only the first `--test-lines` reads (default 100000).
- `--low-quality-threshold`: Phred cutoff for allowing mismatches at low-Q positions (default 20).
- `--umi-allow-one-mismatch`: allow 1-mismatch UMI matching when a UMI table is provided (default exact).
- `--write-stats-file`: write `{outprefix}.stats.txt` alongside stdout stats.
- `--write-stats-json`: write `{outprefix}.stats.json` with the same stats in JSON.
- `--write-offlibrary`: write off-library A/B combos to `{outprefix}.offlibrary.counts.txt`.
- `--write-anndata`: write AnnData `.h5ad` for mapped counts (and off-library if requested; requires `anndata`).
- `--write-count-matrix`: after all samples are processed, merge per-sample counts into `counts_AB.match.counts.txt` (UMI-deduped) and `counts_all.aligned.counts.txt` (read-level) in the output directory.

### Outputs (per sample)
- `{outprefix}.all.aligned.counts.txt` — read-level counts keyed by sgID_AB; zero-filled and ordered by the guide table.
- `{outprefix}.AB.match.counts.txt` — UMI-deduped molecule counts keyed by sgID_AB (or read-level if no UMIs); zero-filled and ordered by the guide table.
- `{outprefix}.offlibrary.counts.txt` (opt): off-library A/B combinations.
- `{outprefix}.h5ad` / `{outprefix}.offlibrary.h5ad` (opt): AnnData exports if `--write-anndata` and `anndata` is installed.
- `{outprefix}.stats.txt` / `{outprefix}.stats.json` (optional) — per-sample stats.
- stdout: human-readable stats; stderr: progress, collision summaries, counts written.

## Usage: Merge Counts into a Matrix (v2)
Merge per-sample count files into a matrix (`counts.txt` with a `target` column):
```bash
python module1/count_files_to_counts_matrix_v2.py \\
  counts_dir out_dir \\
  --suffix .AB.match.counts.txt   # or .all.aligned.counts.txt for read-level
```
Assumes all count files share the same sgID_AB order (as emitted by v2).

## Notes
- Off-library pairs: Reads where A and B each map but the combination isn’t in the guide table are counted and reported; they are not assigned.
- Pair rescue: If A or B is ambiguous alone, but the A×B intersection yields a single library pair, the read is rescued.
- UMI dedup: A UMI is counted once per sgID_AB; duplicate UMIs increment the “reads collapsed as UMI duplicates” stat.
