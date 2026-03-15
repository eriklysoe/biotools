"""Shared utilities: file parsing, validation, rarefaction."""

import io
import numpy as np
import pandas as pd


def parse_abundance_table(file_bytes):
    """
    Parse abundance table from bytes.
    Returns (DataFrame with samples as columns, feature IDs as index, warnings).
    """
    warnings = []

    # Try UTF-8 then latin-1
    for enc in ("utf-8", "latin-1"):
        try:
            text = file_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError("Could not decode file as UTF-8 or latin-1")

    # Auto-detect delimiter
    first_line = text.split("\n")[0]
    sep = "\t" if "\t" in first_line else ","

    df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)

    # Remove completely empty rows/cols
    df = df.dropna(how="all").dropna(axis=1, how="all")

    if df.shape[1] < 2:
        raise ValueError(f"Need at least 2 samples, found {df.shape[1]}")

    # Ensure numeric
    try:
        df = df.apply(pd.to_numeric)
    except (ValueError, TypeError):
        raise ValueError("Abundance table must contain only numeric values")

    # Check for negative values
    if (df.values < 0).any():
        raise ValueError("Abundance table contains negative values — counts must be non-negative")

    # Check for all-zero samples
    zero_samples = df.columns[df.sum(axis=0) == 0].tolist()
    if zero_samples:
        warnings.append(f"All-zero samples removed: {', '.join(zero_samples)}")
        df = df.drop(columns=zero_samples)

    if df.shape[1] < 2:
        raise ValueError("Fewer than 2 samples remain after removing all-zero samples")

    return df, warnings


def parse_metadata(file_bytes):
    """Parse metadata file. Returns DataFrame with sample IDs as index."""
    for enc in ("utf-8", "latin-1"):
        try:
            text = file_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError("Could not decode metadata file")

    first_line = text.split("\n")[0]
    sep = "\t" if "\t" in first_line else ","

    df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0, dtype=str)
    df.index = df.index.astype(str)
    return df


def match_metadata(abundance_df, metadata_df):
    """
    Match metadata to abundance table samples.
    Returns (matched_metadata_df, warnings).
    """
    warnings = []
    samples = set(abundance_df.columns.astype(str))
    meta_samples = set(metadata_df.index.astype(str))

    common = samples & meta_samples
    only_abundance = samples - meta_samples
    only_meta = meta_samples - samples

    if not common:
        raise ValueError("No sample IDs match between abundance table and metadata")

    if only_abundance:
        warnings.append(
            f"{len(only_abundance)} sample(s) in abundance table not in metadata: "
            f"{', '.join(sorted(list(only_abundance)[:5]))}"
        )
    if only_meta:
        warnings.append(
            f"{len(only_meta)} sample(s) in metadata not in abundance table"
        )

    # Reindex metadata to match abundance table column order
    matched = metadata_df.loc[metadata_df.index.isin(common)]
    return matched, warnings


def get_sample_depths(abundance_df):
    """Get read depth statistics per sample."""
    depths = abundance_df.sum(axis=0)
    return {
        "min": int(depths.min()),
        "max": int(depths.max()),
        "median": int(depths.median()),
        "per_sample": {str(s): int(d) for s, d in depths.items()},
    }


def rarefy_table(abundance_df, depth, seed=42):
    """
    Rarefy abundance table to given depth.
    Drops samples below depth. Returns (rarefied_df, dropped_samples, warnings).
    """
    warnings = []
    rng = np.random.RandomState(seed)

    sample_depths = abundance_df.sum(axis=0)
    keep = sample_depths[sample_depths >= depth].index.tolist()
    dropped = sample_depths[sample_depths < depth].index.tolist()

    if dropped:
        warnings.append(
            f"{len(dropped)} sample(s) dropped (below rarefaction depth {depth}): "
            f"{', '.join(str(s) for s in dropped[:5])}"
        )

    if len(keep) < 2:
        raise ValueError(
            f"Fewer than 2 samples have depth >= {depth}. "
            f"Min depth is {int(sample_depths.min())}."
        )

    rarefied = pd.DataFrame(index=abundance_df.index, columns=keep, dtype=int)

    for sample in keep:
        counts = abundance_df[sample].values.astype(int)
        total = counts.sum()

        # Create pool of feature indices
        pool = np.repeat(np.arange(len(counts)), counts)
        chosen = rng.choice(pool, size=depth, replace=False)
        rarefied[sample] = np.bincount(chosen, minlength=len(counts))

    return rarefied, dropped, warnings


def parse_newick_tree(file_bytes):
    """Parse Newick tree from bytes. Returns skbio TreeNode."""
    from skbio import TreeNode

    for enc in ("utf-8", "latin-1"):
        try:
            text = file_bytes.decode(enc).strip()
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError("Could not decode tree file")

    tree = TreeNode.read(io.StringIO(text))
    return tree


def validate_tree_tips(tree, feature_ids):
    """Check overlap between tree tip labels and feature IDs."""
    tip_names = {tip.name for tip in tree.tips() if tip.name}
    feature_set = set(str(f) for f in feature_ids)

    common = tip_names & feature_set
    if not common:
        return False, "No tree tip labels match feature IDs in abundance table"

    if len(common) < len(feature_set):
        missing = len(feature_set) - len(common)
        return True, f"{missing} feature(s) not found in tree — they will be excluded from phylogenetic metrics"

    return True, None
