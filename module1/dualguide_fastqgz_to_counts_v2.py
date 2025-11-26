import sys
import os
import re
import argparse
import gzip
import io
import multiprocessing
from glob import glob
from typing import List
import pandas as pd
# Silence pandas implicit modification warnings globally
import warnings
warnings.filterwarnings("ignore", message=".*ImplicitModificationWarning.*")
warnings.filterwarnings("ignore", message=".*Transforming to str index.*")
warnings.filterwarnings("ignore", message=".*PerformanceWarning: DataFrame is highly fragmented.*")

# Silence pandas implicit modification warnings from index coercion
try:
    from pandas.errors import ImplicitModificationWarning
    warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
except Exception:
    warnings.filterwarnings("ignore", message=".*ImplicitModificationWarning.*")
# Allow running as script by adding module path when not installed as a package
try:
    import matrix_utils
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    import matrix_utils
try:
    import anndata as ad
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False

DEFAULT_LOW_QUALITY_THRESHOLD = 20
DEFAULT_TEST_LINES = 100000

def getMismatchDict(table, column, trim_range=None, allowOneMismatch=True, id_column=None):
    mismatch_dict: dict[str, set[str]] = dict()
    clash_zero = 0
    clash_one = 0

    if trim_range:
        col = table[column].apply(lambda seq: seq[trim_range[0]:trim_range[1]])
    else:
        col = table[column]

    if id_column is not None:
        ids = table[id_column].astype(str)
    else:
        ids = col.index.astype(str)

    for sgRNA_id, seq in zip(ids, col):
        if seq in mismatch_dict:
            clash_zero += 1
            mismatch_dict[seq].add(sgRNA_id)
        else:
            mismatch_dict[seq] = {sgRNA_id}

        if allowOneMismatch:
            for position in range(len(seq)):
                mismatchSeq = seq[:position] + 'N' + seq[position + 1:]

                if mismatchSeq in mismatch_dict:
                    if sgRNA_id not in mismatch_dict[mismatchSeq]:
                        clash_one += 1
                        mismatch_dict[mismatchSeq].add(sgRNA_id)
                else:
                    mismatch_dict[mismatchSeq] = {sgRNA_id}

    if clash_zero or clash_one:
        summary = []
        if clash_zero:
            summary.append(f'{clash_zero} exact collisions')
        if clash_one:
            summary.append(f'{clash_one} one-mismatch collisions')
            print(
            f"{column}: {'; '.join(summary)}; storing all candidates, will require disambiguation at pair level.",
            file=sys.stderr,
        )

    return mismatch_dict


def matchBarcode(mismatch_dict, barcode, allowOneMismatch=True, quality=None, low_q_threshold=None):
    """
    Return a tuple: (call, used_mismatch, candidates)
      - call: 'none', 'multiple', or the single matching sgRNA id
      - used_mismatch: True if a 1-mismatch lookup was needed
      - candidates: set of all matching sgRNA ids (empty if none)
    """
    candidates: set[str] = set()
    used_mismatch = False

    direct = mismatch_dict.get(barcode)
    if direct:
        candidates.update(direct)

    # Only search 1-mismatch neighborhood if no exact match and mismatches allowed
    if not candidates and allowOneMismatch:
        mismatch_hits: set[str] = set()
        length = len(barcode)
        if quality is not None and low_q_threshold is not None:
            low_q_positions = [
                idx for idx, qchar in enumerate(quality)
                if (ord(qchar) - 33) < low_q_threshold
            ]
            for position in low_q_positions:
                mismatchSeq = barcode[:position] + 'N' + barcode[position + 1:]
                hits = mismatch_dict.get(mismatchSeq)
                if hits:
                    mismatch_hits.update(hits)
        else:
            for position in range(length):
                mismatchSeq = barcode[:position] + 'N' + barcode[position + 1:]
                hits = mismatch_dict.get(mismatchSeq)
                if hits:
                    mismatch_hits.update(hits)

        if mismatch_hits:
            used_mismatch = True
            candidates.update(mismatch_hits)

    if not candidates:
        return 'none', used_mismatch, set()
    if len(candidates) == 1:
        return next(iter(candidates)), used_mismatch, candidates
    return 'multiple', used_mismatch, candidates


def parse_umi_from_header(header_line):
    """
    Parse the UIB (integration barcode) from a FASTQ header line if present.
    Header formats:
      without UIB: @...:1080 1:N:0:ACTCGNTA
      with UIB   : @...:1080:TCAGTCGA 1:N:0:ACTCGNTA
    """
    header = header_line.decode('utf-8').strip()
    first_token = header.split(' ')[0]
    parts = first_token.split(':')
    if len(parts) >= 8:
        return parts[-1]
    return None


def peek_has_umi(fastq_path):
    with gzip.open(fastq_path) as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError(f'FASTQ file appears empty: {fastq_path}')
        umi = parse_umi_from_header(first_line)
        return umi is not None


def writeToCounts(fileTup):
    (
        r1file,
        r2file,
        outprefix,
        mismatchDicts,
        umiMismatchDict,
        has_umi,
        umi_allow_one_mismatch,
        use_umi_whitelist,
        testRun,
        valid_pairs,
        test_lines,
        low_q_threshold,
        sgid_order,
        write_stats_file,
        write_stats_json,
        write_offlibrary,
        write_anndata,
    ) = fileTup

    print(f'Processing files: R1={r1file}, R2={r2file}', file=sys.stderr)
    sys.stderr.flush()

    statsCounts = {
        'A sgRNA not mapped': 0,
        'B sgRNA not mapped': 0,
        'A sgRNA multiple mappings': 0,
        'B sgRNA multiple mappings': 0,
        'All sgRNAs uniquely map': 0,  # A and B (and UIB if present) unique, before pair validation
        'Valid library A/B pair': 0,
        'A sgRNA and B sgRNA do not match': 0,
        'A sgRNA and B sgRNA match': 0,
        'A sgRNA exact match': 0,
        'A sgRNA 1-mismatch match': 0,
        'B sgRNA exact match': 0,
        'B sgRNA 1-mismatch match': 0,
        'Reads rescued by pair intersection': 0,
        'Reads unresolved due to pair ambiguity': 0,
        'Reads with A/B not in library': 0,
    }

    if has_umi:
        statsCounts.update({
            'UIB not mapped': 0,
            'UIB multiple mappings': 0,
            'UIB present in reads': 0,
        })

    pairCounts_sgRNAs = dict()
    pairCounts_double = dict()
    # No dedup; UIB is part of the key
    offlibrary_counts: dict[str, int] = {}

    current_umi = None
    read_count = 0

    with gzip.open(r1file) as infile_r1_raw, gzip.open(r2file) as infile_r2_raw:
        infile_r1 = io.BufferedReader(infile_r1_raw, buffer_size=1024 * 1024)
        infile_r2 = io.BufferedReader(infile_r2_raw, buffer_size=1024 * 1024)
        r1_seq = None
        r2_seq = None
        r1_qual = None
        r2_qual = None
        for i, (r1, r2) in enumerate(zip(infile_r1, infile_r2)):
            read_idx = i // 4
            if testRun and read_idx >= test_lines:
                break
            if i % 4 == 0:
                current_umi = parse_umi_from_header(r1) if has_umi else None
                if has_umi and current_umi is not None:
                    statsCounts['UIB present in reads'] += 1
            elif i % 4 == 1:
                r1_seq = r1.strip().decode('utf-8')
                r2_seq = r2.strip().decode('utf-8')
            elif i % 4 == 3:
                r1_qual = r1.strip().decode('utf-8')
                r2_qual = r2.strip().decode('utf-8')

                # Reads are 19bp long after trimming; use the full read to match the
                # guide table (which has already been trimmed to positions 1:20).
                seq_a = r1_seq[:19]
                qual_a = r1_qual[:19]
                seq_b = r2_seq[:19]
                qual_b = r2_qual[:19]

                protospacer_a_r1, a_used_mismatch, a_candidates = matchBarcode(
                    mismatchDicts['protospacer_a_r1'],
                    seq_a,
                    allowOneMismatch=True,
                    quality=qual_a,
                    low_q_threshold=low_q_threshold,
                )
                protospacer_b_r2, b_used_mismatch, b_candidates = matchBarcode(
                    mismatchDicts['protospacer_b_r2'],
                    seq_b,
                    allowOneMismatch=True,
                    quality=qual_b,
                    low_q_threshold=low_q_threshold,
                )

                umi_candidates = set()
                if has_umi:
                    umi_call, umi_used_mismatch, umi_candidates = matchBarcode(
                        umiMismatchDict,
                        current_umi if current_umi else '',
                        allowOneMismatch=umi_allow_one_mismatch,
                    )

                def status(cands: set[str]) -> str:
                    if not cands:
                        return 'none'
                    if len(cands) > 1:
                        return 'multiple'
                    return 'unique'

                a_status = status(a_candidates)
                b_status = status(b_candidates)
                umi_status = status(umi_candidates) if has_umi else 'unique'

                # Per-guide stats
                if a_status == 'unique':
                    if a_used_mismatch:
                        statsCounts['A sgRNA 1-mismatch match'] += 1
                    else:
                        statsCounts['A sgRNA exact match'] += 1
                elif a_status == 'none':
                    statsCounts['A sgRNA not mapped'] += 1
                else:
                    statsCounts['A sgRNA multiple mappings'] += 1

                if b_status == 'unique':
                    if b_used_mismatch:
                        statsCounts['B sgRNA 1-mismatch match'] += 1
                    else:
                        statsCounts['B sgRNA exact match'] += 1
                elif b_status == 'none':
                    statsCounts['B sgRNA not mapped'] += 1
                else:
                    statsCounts['B sgRNA multiple mappings'] += 1

                if has_umi:
                    if umi_status == 'none':
                        statsCounts['UIB not mapped'] += 1
                    elif umi_status == 'multiple':
                        statsCounts['UIB multiple mappings'] += 1

                # Pair-level resolution (including rescue from ambiguous A/B)
                final_pair_id = None  # sgID_AB
                rescued = False

                can_attempt_pair = (
                    a_status != 'none'
                    and b_status != 'none'
                    and (not has_umi or umi_status == 'unique')
                )

                if can_attempt_pair:
                    statsCounts['All sgRNAs uniquely map'] += 1
                    candidate_pairs = []
                    for a_id in a_candidates:
                        for b_id in b_candidates:
                            sg_ab = valid_pairs.get((a_id, b_id))
                            if sg_ab:
                                candidate_pairs.append((a_id, b_id, sg_ab))

                    if len(candidate_pairs) == 1:
                        final_pair_id = candidate_pairs[0][2]
                        if len(a_candidates) > 1 or len(b_candidates) > 1:
                            rescued = True
                            statsCounts['Reads rescued by pair intersection'] += 1
                    elif len(candidate_pairs) > 1:
                        statsCounts['Reads unresolved due to pair ambiguity'] += 1
                    else:
                        # Both mapped but combination not in library
                        statsCounts['Reads with A/B not in library'] += 1
                        umi_suffix = f"++{current_umi}" if has_umi and current_umi else ""
                        offlibrary_key = f"{'|'.join(sorted(a_candidates))}++{'|'.join(sorted(b_candidates))}{umi_suffix}"
                        offlibrary_counts[offlibrary_key] = offlibrary_counts.get(offlibrary_key, 0) + 1

                # Only count uniquely resolved A+B (whether direct or rescued)
                if final_pair_id is not None:
                    statsCounts['Valid library A/B pair'] += 1
                    combinedSgId = final_pair_id
                    if combinedSgId not in pairCounts_sgRNAs:
                        pairCounts_sgRNAs[combinedSgId] = 0
                    pairCounts_sgRNAs[combinedSgId] += 1

                    statsCounts['A sgRNA and B sgRNA match'] += 1

                    if has_umi:
                        umi_str = str(next(iter(umi_candidates)) if umi_candidates else "")
                        combined_key = combinedSgId + '++' + umi_str
                        pairCounts_double[combined_key] = pairCounts_double.get(combined_key, 0) + 1
                    else:
                        pairCounts_double[combinedSgId] = pairCounts_double.get(combinedSgId, 0) + 1
                elif can_attempt_pair:
                    statsCounts['A sgRNA and B sgRNA do not match'] += 1

                read_count += 1

    with open(outprefix + '.all.aligned.counts.txt', 'w') as outfile:
        for pair in sgid_order:
            outfile.write(pair + '\t' + str(pairCounts_sgRNAs.get(pair, 0)) + '\n')

    with open(outprefix + '.AB.match.counts.txt', 'w') as outfile:
        if has_umi:
            # Write as wide matrix: rows = sgID_AB, cols = UIBs
            umi_values = sorted({key.split('++', 1)[1] if '++' in key else '' for key in pairCounts_double.keys()})
            header = ['sgID_AB'] + umi_values
            outfile.write('\t'.join(header) + '\n')
            for guide in sgid_order:
                row = [guide]
                for umi in umi_values:
                    key = guide + ('++' + umi if umi else '')
                    row.append(str(pairCounts_double.get(key, 0)))
                outfile.write('\t'.join(row) + '\n')
        else:
            for pair in sgid_order:
                outfile.write(pair + '\t' + str(pairCounts_double.get(pair, 0)) + '\n')
    # stderr summary for file outputs
    written_all = sum(pairCounts_sgRNAs.values())
    written_ab = sum(pairCounts_double.values())
    print(
        f"{os.path.basename(outprefix)}: wrote {written_all} counts to .all.aligned.counts.txt "
        f"and {written_ab} counts to .AB.match.counts.txt",
        file=sys.stderr,
    )

    numReads = read_count if read_count else 1  # avoid zero-division in degenerate cases
    stats_lines = []
    stats_dict = {}
    stats_lines.append(f"{os.path.basename(outprefix)} {numReads} reads")
    stats_dict["sample"] = os.path.basename(outprefix)
    stats_dict["total_reads_processed"] = numReads

    mapped_a = numReads - statsCounts['A sgRNA not mapped']
    mapped_b = numReads - statsCounts['B sgRNA not mapped']
    multi_a = statsCounts['A sgRNA multiple mappings']
    multi_b = statsCounts['B sgRNA multiple mappings']
    stats_lines.append(f"A sgRNAs mapping {mapped_a} ({mapped_a * 100.0 / numReads:.2f}%)")
    stats_lines.append(f"  A sgRNAs exact matches {statsCounts['A sgRNA exact match']}")
    stats_lines.append(f"  A sgRNAs via 1 mismatch {statsCounts['A sgRNA 1-mismatch match']}")
    stats_lines.append(f"  A sgRNAs multi-match {multi_a}")
    stats_lines.append(f"B sgRNAs mapping {mapped_b} ({mapped_b * 100.0 / numReads:.2f}%)")
    stats_lines.append(f"  B sgRNAs exact matches {statsCounts['B sgRNA exact match']}")
    stats_lines.append(f"  B sgRNAs via 1 mismatch {statsCounts['B sgRNA 1-mismatch match']}")
    stats_lines.append(f"  B sgRNAs multi-match {multi_b}")
    stats_dict.update({
        "A_mapped": mapped_a,
        "A_exact": statsCounts['A sgRNA exact match'],
        "A_1mm": statsCounts['A sgRNA 1-mismatch match'],
        "A_multi": multi_a,
        "B_mapped": mapped_b,
        "B_exact": statsCounts['B sgRNA exact match'],
        "B_1mm": statsCounts['B sgRNA 1-mismatch match'],
        "B_multi": multi_b,
    })
    if has_umi:
        mapped_umi = numReads - statsCounts['UIB not mapped']
        multi_umi = statsCounts['UIB multiple mappings']
        stats_lines.append(f"Total UIBs parsed from headers {statsCounts['UIB present in reads']}")
        stats_lines.append(f"UIBs mapping {mapped_umi} ({mapped_umi * 100.0 / numReads:.2f}%)")
        stats_lines.append(f"  UIBs multi-match {multi_umi}")
        stats_dict.update({
            "UIBs_parsed": statsCounts['UIB present in reads'],
            "UIBs_mapped": mapped_umi,
            "UIBs_multi": multi_umi,
        })
    both_unique = statsCounts['All sgRNAs uniquely map']
    valid_pairs_assigned = statsCounts['Valid library A/B pair']
    rescued = statsCounts['Reads rescued by pair intersection']
    unresolved = statsCounts['Reads unresolved due to pair ambiguity']
    not_in_library = statsCounts['Reads with A/B not in library']

    stats_lines.append(f"Reads with A, B, and UIB mapped (pre-pair validation) {both_unique}")
    stats_lines.append(f"  Assigned to valid library A/B pair {valid_pairs_assigned} (includes {rescued} rescued by pair intersection)")
    stats_lines.append(f"  Unresolved due to pair ambiguity {unresolved}")
    stats_lines.append(f"  A/B combo not present in library {not_in_library}")
    stats_dict.update({
        "pre_pair_unique": both_unique,
        "assigned_valid_pair": valid_pairs_assigned,
        "rescued_by_pair": rescued,
        "unresolved_ambiguity": unresolved,
        "ab_not_in_library": not_in_library,
    })
    if has_umi:
        uib_events = sum(pairCounts_double.values())
        guides_with_uib_events = len(pairCounts_double)
        stats_lines.append(f"UIB-tagged events (A/B/UIB) {uib_events}")
        stats_lines.append(f"Guides with ≥1 UIB-tagged event {guides_with_uib_events}")
        stats_dict.update({
            "uib_tagged_events": uib_events,
            "guides_with_uib_events": guides_with_uib_events,
        })

    if not (write_stats_file or write_stats_json):
        for line in stats_lines:
            print(line)
    if write_stats_file:
        stats_path = outprefix + '.stats.txt'
        with open(stats_path, 'w') as sf:
            sf.write('\n'.join(stats_lines))
    if write_stats_json:
        stats_path_json = outprefix + '.stats.json'
        with open(stats_path_json, 'w') as jf:
            import json
            json.dump(stats_dict, jf, indent=2)

    print()

    if write_offlibrary:
        offlib_path = outprefix + '.offlibrary.counts.txt'
        with open(offlib_path, 'w') as off:
            for key, val in sorted(offlibrary_counts.items()):
                off.write(f"{key}\t{val}\n")
        print(
            f"{os.path.basename(outprefix)}: wrote {len(offlibrary_counts)} off-library A/B combos to .offlibrary.counts.txt",
            file=sys.stderr,
        )

    if write_anndata:
        if not ANNDATA_AVAILABLE:
            print("anndata not installed; skipping AnnData export", file=sys.stderr)
        else:
            if has_umi:
                # pivot guide x UIB
                rows = {}
                for key, val in pairCounts_double.items():
                    if '++' in key:
                        guide, umi = key.split('++', 1)
                    else:
                        guide, umi = key, ''
                    rows.setdefault(guide, {})[umi] = val
                adata_df = pd.DataFrame.from_dict(rows, orient='index').reindex(pd.Index(sgid_order, dtype=str)).fillna(0)
                adata = ad.AnnData(X=adata_df)
            else:
                adata = ad.AnnData(
                    X=pd.DataFrame(pairCounts_double, index=[0]).T.reindex(pd.Index(sgid_order, dtype=str)).fillna(0),
                )
            adata.write_h5ad(outprefix + '.h5ad')
            if write_offlibrary:
                off_df = pd.DataFrame(offlibrary_counts, index=[0]).T if offlibrary_counts else pd.DataFrame()
                ad_off = ad.AnnData(X=off_df) if not off_df.empty else ad.AnnData(X=pd.DataFrame())
                ad_off.write_h5ad(outprefix + '.offlibrary.h5ad')

    sys.stdout.flush()


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Process raw R1/R2 sequencing data from screens to counts files, "
            "with optional UIB (integration barcode) detection from read IDs."
        )
    )
    parser.add_argument("Guide_Table", help="Table of sgRNA pairs in the library (CSV).")
    parser.add_argument("Out_File_Path", help="Directory where output files should be written.")
    parser.add_argument(
        "Seq_File_Names",
        nargs="+",
        help=(
            "Name(s) of sequencing file(s). Unix wildcards can be used to select "
            "multiple files at once. The script will search for all *.fastq.gz, "
            "*.fastq, *.fq.gz, *.fq, *.fa, *.fasta, and *.fna files matching the pattern(s)."
        ),
    )
    parser.add_argument(
        "--UMI_Table",
        dest="UMI_Table",
        default=None,
        help="Table listing valid UIBs (column name: UMI). Required if UIBs are present in read IDs.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help=(
            "Run the entire script on only the first %d lines of each file. "
            "Be sure to delete or move all test files before re-running script "
            "as they will not be overwritten."
        )
        % DEFAULT_TEST_LINES,
    )
    parser.add_argument(
        "--test-lines",
        type=int,
        default=DEFAULT_TEST_LINES,
        help="Number of reads to process when --test is supplied (default: %(default)s).",
    )
    parser.add_argument(
        "--low-quality-threshold",
        type=int,
        default=DEFAULT_LOW_QUALITY_THRESHOLD,
        help=(
            "Phred quality threshold used to decide where mismatches are allowed "
            "when matching protospacers (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--umi-allow-one-mismatch",
        action="store_true",
        help="Allow 1-mismatch matching for UIBs (default: exact only).",
    )
    parser.add_argument(
        "--write-stats-file",
        action="store_true",
        help="Also write per-sample stats to {outprefix}.stats.txt in addition to printing to stdout.",
    )
    parser.add_argument(
        "--write-stats-json",
        action="store_true",
        help="Also write per-sample stats to {outprefix}.stats.json in addition to printing to stdout.",
    )
    parser.add_argument(
        "--write-offlibrary",
        action="store_true",
        help="Write counts for off-library A/B combinations to {outprefix}.offlibrary.counts.txt.",
    )
    parser.add_argument(
        "--write-anndata",
        action="store_true",
        help="Write AnnData objects for mapped counts (and off-library if requested).",
    )
    parser.add_argument(
        "--write-count-matrix",
        action="store_true",
        help="After processing, merge per-sample count files into counts matrices (.AB.match and .all.aligned) in the output directory.",
    )
    return parser.parse_args()


def ensure_dir_for_output(path: str) -> str:
    """
    Validate/create output directory early to avoid confusing downstream errors.
    Ensures the directory exists and is writable.
    """
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise NotADirectoryError(f"Output path exists but is not a directory: {path}")
    else:
        try:
            os.makedirs(path)
        except OSError as e:
            raise OSError(f"Failed to create output directory '{path}': {e}") from e

    # Check write permission by trying to create a temp file
    try:
        test_path = os.path.join(path, ".write_test.tmp")
        with open(test_path, "w") as fh:
            fh.write("")
        os.remove(test_path)
    except OSError as e:
        raise PermissionError(f"Output directory is not writable: {path} ({e})") from e

    return path


def validate_file_exists(path: str, description: str) -> None:
    """Common helper to make errors more readable."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{description} not found: {path}")
    if not os.path.isfile(path):
        raise IsADirectoryError(f"{description} is not a regular file: {path}")


def read_and_validate_guide_table(path: str, required_cols: list[str]) -> pd.DataFrame:
    """Read guide table CSV and ensure required columns are present."""
    validate_file_exists(path, "Guide table")
    try:
        df = pd.read_csv(path, sep=",", header=0)
    except Exception as e:
        raise ValueError(f"Failed to read guide table '{path}': {e}") from e

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            "Guide table is missing required column(s): "
            f"{', '.join(missing)}\n"
            f"Columns found: {', '.join(df.columns)}"
        )
    if df.empty:
        raise ValueError(f"Guide table '{path}' is empty after parsing.")

    # Keep only required columns and drop any rows with missing values
    df = df[required_cols].dropna()
    # Keep explicit columns; do not set index so we can map sgID_A/B separately
    # convert protospacers to mimic reads according to sequencing strategy, assuming R1=19, R2=19
    df['protospacer_A'] = df['protospacer_A'].str[1:20].str.upper()
    trans = str.maketrans('ATGC', 'TACG')
    df['protospacer_B'] = df['protospacer_B'].str[1:20].str.upper().str.translate(trans).str[::-1]

    return df


def resolve_seq_files(patterns: List[str]) -> List[str]:
    """
    Expand Unix wildcards, filter by allowed extensions, and ensure we
    actually find paired R1/R2 files for each sample.

    Assumes Illumina-style filenames like:
        sample_L001_R1_001.fastq.gz
        sample_L001_R2_001.fastq.gz
        sample_L001_I1_001.fastq.gz  (ignored)
        sample_L001_I2_001.fastq.gz  (ignored)
    """
    allowed_exts = (
        ".fastq", ".fastq.gz", ".fq", ".fq.gz",
        ".fa", ".fasta", ".fna",
    )

    resolved: List[str] = []
    unmatched_patterns: List[str] = []

    # 1) Expand patterns and filter by extension
    for pattern in patterns:
        matches = glob(pattern)
        matches = [m for m in matches if m.lower().endswith(allowed_exts)]
        if not matches:
            unmatched_patterns.append(pattern)
        else:
            resolved.extend(matches)

    resolved = sorted(set(resolved))  # dedupe + deterministic order

    if unmatched_patterns:
        raise FileNotFoundError(
            "No sequencing files found matching the following pattern(s): "
            + ", ".join(unmatched_patterns)
        )

    if not resolved:
        raise FileNotFoundError(
            "No sequencing files found. "
            "Check your Seq_File_Names patterns and file extensions."
        )

    # 2) Group by sample and read type using Illumina-style pattern
    #
    # Captures:
    #   group(1) -> sample+lane prefix
    #   group(2) -> read label (R1, R2, I1, I2, etc.)
    #   group(3) -> optional chunk/index (_001, _002, etc.)
    read_regex = re.compile(
        r"(.+)_R([12I][12]?)(_?\d+)?\.(fastq|fq|fa|fasta|fna)(\.gz)?$",
        re.IGNORECASE,
    )

    groups: dict[str, dict[str, str]] = {}
    unparsed: list[str] = []

    for path in resolved:
        fname = os.path.basename(path)
        m = read_regex.match(fname)
        if not m:
            unparsed.append(path)
            continue

        sample_key = m.group(1)
        # include lane/chunk if present, so different lanes are separate samples
        if m.group(3):
            sample_key = f"{sample_key}{m.group(3)}"

        read_label = "R" + m.group(2).upper()  # ensures R1/R2/I1/I2-style

        if read_label not in ("R1", "R2", "I1", "I2"):
            unparsed.append(path)
            continue

        slot = groups.setdefault(sample_key, {})
        if read_label in slot and slot[read_label] != path:
            raise ValueError(
                f"Duplicate {read_label} files for sample '{sample_key}':\n"
                f"  {slot[read_label]}\n"
                f"  {path}"
            )
        slot[read_label] = path

    # 3) Enforce R1/R2 pairing (ignore I1/I2 presence)
    bad_samples = [
        s for s, g in groups.items()
        if "R1" not in g or "R2" not in g
    ]

    if bad_samples:
        raise ValueError(
            "Each sample must have both R1 and R2 files, but the following "
            "sample(s) are missing one of them:\n  "
            + "\n  ".join(bad_samples)
        )

    # 4) Build final list: return only R1/R2, ignore I1/I2
    paired_files: List[str] = []
    for sample_key, g in sorted(groups.items()):
        # deterministic: R1 then R2 for each sample
        paired_files.append(g["R1"])
        paired_files.append(g["R2"])

    # 5) Optional: check readability
    for path in paired_files:
        if not os.access(path, os.R_OK):
            raise PermissionError(f"Sequencing file is not readable: {path}")

    return paired_files


def read_and_validate_umi_table(path: str | None, allow_one_mismatch: bool) -> pd.Series | None:
    """
    If a UIB table is provided, validate that it exists and has a 'UMI' column.
    Returns a Series of UIBs or None.
    """
    if path is None:
        return None

    validate_file_exists(path, "UMI table")
    try:
        df = pd.read_csv(path, sep=",", header=0)
    except Exception as e:
        raise ValueError(f"Failed to read UMI table '{path}': {e}") from e

    if "UMI" not in df.columns:
        raise ValueError(
            f"UMI table '{path}' must contain a column named 'UMI'. "
            f"Columns found: {', '.join(df.columns)}"
        )
    if df["UMI"].isna().all():
        raise ValueError(f"UMI table '{path}' has no non-empty UMIs in column 'UMI'.")

    print('UIBs in table', df.shape[0])

    umiMismatchDict = getMismatchDict(df, 'UMI', allowOneMismatch=allow_one_mismatch, id_column='UMI')

    return umiMismatchDict

if __name__ == '__main__':
    args = parse_args()

    outputDirectory = ensure_dir_for_output(args.Out_File_Path)

    inputFileList = resolve_seq_files(args.Seq_File_Names)

    required_cols = ["sgID_AB", "sgID_A", "protospacer_A", "sgID_B", "protospacer_B"]
    guideTable = read_and_validate_guide_table(args.Guide_Table, required_cols)

    print('sgRNAs in library', len(guideTable), file=sys.stderr)
    sgid_order = list(guideTable['sgID_AB'])

    # Map of valid A/B combos to sgID_AB for downstream pair disambiguation
    valid_pairs = {
        (str(row['sgID_A']), str(row['sgID_B'])): str(row['sgID_AB'])
        for _, row in guideTable.iterrows()
    }

    combinedMismatchDicts = {
        'protospacer_a_r1': getMismatchDict(guideTable, 'protospacer_A', allowOneMismatch=True, id_column='sgID_A'),
        'protospacer_b_r2': getMismatchDict(guideTable, 'protospacer_B', allowOneMismatch=True, id_column='sgID_B'),
    }

    umiMismatchDict = read_and_validate_umi_table(args.UMI_Table, allow_one_mismatch=args.umi_allow_one_mismatch)

    fileTups = []
    for i, fastqfile in enumerate(inputFileList):
        if i % 2 == 0:
            r1file = fastqfile
        elif i % 2 == 1:
            r2file = fastqfile
            sample_has_umi = peek_has_umi(r1file)
            if sample_has_umi and umiMismatchDict is None:
                raise ValueError("UIBs detected in read IDs but no UIB table was provided. Supply --UMI_Table.")
            use_umi_whitelist = umiMismatchDict is not None
            sample_name = os.path.split(fastqfile)[-1].split('_R')[0]
            # Drop trailing Illumina chunk like _S##_L00# if present for consistency
            sample_name = re.sub(r"_S\d+_L\d+$", "", sample_name)
            outputfile = os.path.join(outputDirectory, sample_name)
            fileTups.append(
                (
                    r1file,
                    r2file,
                    outputfile,
                    combinedMismatchDicts,
                    umiMismatchDict,
                    sample_has_umi,
                    args.umi_allow_one_mismatch,
                    use_umi_whitelist,
                    args.test,
                    valid_pairs,
                    args.test_lines,
                    args.low_quality_threshold,
                    sgid_order,
                    args.write_stats_file,
                    args.write_stats_json,
                    args.write_offlibrary,
                    args.write_anndata,
                )
            )

    pool = multiprocessing.Pool(len(fileTups))
    pool.map(writeToCounts, fileTups)
    pool.close()
    pool.join()

    if args.write_count_matrix:
        # Build matrices for both AB.match and all.aligned
        matrix_utils.write_matrices(outputDirectory, sgid_order)
