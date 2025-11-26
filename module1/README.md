# Module 1: Dual sgRNA FASTQ to Counts (v2)

This module processes paired-end dual-sgRNA FASTQ files into count tables with integration barcode (UIB; formerly “UMI”) tracking. It supersedes the original scripts with improved validation, pair-aware disambiguation, UIB handling, and richer stats reporting.

## Key Features (v2)
- Pair-aware matching: rescues ambiguous A/B hits when the pair uniquely matches the library.
- UIB handling:
  - Requires a UIB whitelist when UIBs are present in read headers (`--UMI_Table` with column `UMI`).
  - Counts are not deduplicated; each UIB is treated as a unique integration barcode and retained as a column.
  - Produces read-level counts (`*.all.aligned.counts.txt`) and UIB-tagged counts (`*.AB.match.counts.txt`, wide matrix with one column per observed UIB).
- Robust input validation (guide table columns, file pairing, output directory).
- Detailed stats to stdout; progress to stderr; optional text/JSON stats files.
- Parameterized test subset, quality threshold, UIB mismatch option, and buffered FASTQ reading.

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
- `count_files_to_counts_matrix_v2.py` — merges per-sample count files into a matrix (read-level or UIB-tagged).
- Legacy scripts (`dualguide_fastqgz_to_counts.py`, `count_files_to_counts_matrix.py`, UMI/UMI-dedup variants) are retained but not recommended.

## Input FASTQ assumptions
- FASTQs are produced by `bcl-convert` (or equivalent) with R1/R2 only; I1 is not required as a separate file.
- The UIB (often captured in Index 1) is included in the read IDs, e.g. `@...:1080:TCAGTCGA 1:N:0:ACTCGNTA`. The pipeline parses the UIB from the header; no separate I1 file is needed.

## Usage: FASTQ to Counts (v2)
```bash
python module1/dualguide_fastqgz_to_counts_v2.py \\
  [--UMI_Table UIB_sequences.csv] \\
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

### UIB whitelist
- `--UMI_Table`: CSV with column `UMI` (UIB whitelist). Required if UIBs are present in headers.

### Useful flags
- `--test`: process only the first `--test-lines` reads (default 100000).
- `--low-quality-threshold`: Phred cutoff for allowing mismatches at low-Q positions (default 20).
- `--umi-allow-one-mismatch`: allow 1-mismatch UIB matching when a UIB table is provided (default exact).
- `--write-stats-file`: write `{outprefix}.stats.txt` alongside stdout stats.
- `--write-stats-json`: write `{outprefix}.stats.json` with the same stats in JSON.
- `--write-offlibrary`: write off-library A/B combos to `{outprefix}.offlibrary.counts.txt`.
- `--write-anndata`: write AnnData `.h5ad` for mapped counts (and off-library if requested; requires `anndata`).
- `--write-count-matrix`: after all samples are processed, merge per-sample counts into `counts_AB.match.counts.txt` (UIB-tagged, wide) and `counts_all.aligned.counts.txt` (read-level) in the output directory.

### Outputs (per sample)
- `{outprefix}.all.aligned.counts.txt` — read-level counts keyed by sgID_AB; zero-filled and ordered by the guide table.
- `{outprefix}.AB.match.counts.txt` — wide, UIB-tagged counts: rows = sgID_AB, columns = UIBs observed in that sample, entries = read counts (no dedup).
- `{outprefix}.offlibrary.counts.txt` (opt): off-library A/B combinations (including UIB suffix).
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
For `.all.aligned` (read-level), rows follow the guide order. For `.AB.match` (UIB-tagged), rows = `sgID_AB` and columns = `sample-UIB` combinations, aggregating all observed UIBs across samples.

## Notes
- Off-library pairs: Reads where A and B each map but the combination isn’t in the guide table are counted and reported; they are not assigned.
- Pair rescue: If A or B is ambiguous alone, but the A×B intersection yields a single library pair, the read is rescued.
- UIBs are treated as integration barcodes: no deduplication is performed; per-UIB counts are preserved.
