"""PCoA pipeline: parsing, distance matrix, ordination."""

import csv
import io

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

DISTANCE_METRICS = {
    "braycurtis": "Bray-Curtis dissimilarity",
    "jaccard": "Jaccard distance",
    "euclidean": "Euclidean distance",
    "correlation": "Correlation distance",
}

NORMALISATION_METHODS = {
    "none": "None",
    "tss": "Total-sum scaling",
    "hellinger": "Hellinger transformation",
}


# ── Parsing ──────────────────────────────────────────────────────────

def _read_text(file_bytes):
    """Decode bytes to text, trying UTF-8 then latin-1."""
    try:
        return file_bytes.decode("utf-8")
    except UnicodeDecodeError:
        return file_bytes.decode("latin-1")


def _detect_delimiter(text):
    """Auto-detect CSV vs TSV."""
    first_line = text.split("\n", 1)[0]
    if "\t" in first_line:
        return "\t"
    return ","


def parse_abundance_table(file_bytes):
    """
    Parse an abundance table (CSV/TSV).

    Returns (DataFrame, warnings) where:
    - DataFrame: features (rows) x samples (columns), all numeric
    - warnings: list of warning strings
    """
    warnings = []
    text = _read_text(file_bytes)
    delimiter = _detect_delimiter(text)

    df = pd.read_csv(io.StringIO(text), sep=delimiter, index_col=0)

    # Strip whitespace from column/index names
    df.columns = df.columns.str.strip()
    df.index = df.index.astype(str).str.strip()

    # Coerce to numeric, dropping non-numeric columns
    numeric_df = df.apply(pd.to_numeric, errors="coerce")
    non_numeric_cols = numeric_df.columns[numeric_df.isna().all()].tolist()
    if non_numeric_cols:
        warnings.append(
            f"Dropped non-numeric columns: {', '.join(non_numeric_cols)}"
        )
        numeric_df = numeric_df.drop(columns=non_numeric_cols)

    if numeric_df.shape[1] < 3:
        raise ValueError(
            f"Need at least 3 samples, found {numeric_df.shape[1]}"
        )

    # Fill remaining NaN with 0
    numeric_df = numeric_df.fillna(0)

    # Check for all-zero samples
    zero_samples = numeric_df.columns[numeric_df.sum(axis=0) == 0].tolist()
    if zero_samples:
        warnings.append(
            f"Removed all-zero samples: {', '.join(zero_samples)}"
        )
        numeric_df = numeric_df.drop(columns=zero_samples)
        if numeric_df.shape[1] < 3:
            raise ValueError(
                "Fewer than 3 non-zero samples remain after removing "
                "all-zero samples"
            )

    return numeric_df, warnings


def parse_metadata(file_bytes):
    """
    Parse a metadata file (CSV/TSV).

    Returns DataFrame with sample IDs as index.
    """
    text = _read_text(file_bytes)
    delimiter = _detect_delimiter(text)

    df = pd.read_csv(io.StringIO(text), sep=delimiter, index_col=0)
    df.columns = df.columns.str.strip()
    df.index = df.index.astype(str).str.strip()
    # Convert all values to strings for JSON safety
    df = df.astype(str)
    return df


def match_metadata(abundance_df, metadata_df):
    """
    Match metadata samples to abundance table samples.

    Returns (matched_metadata_df, warnings).
    """
    warnings = []
    abundance_samples = set(abundance_df.columns)
    metadata_samples = set(metadata_df.index)

    only_abundance = abundance_samples - metadata_samples
    only_metadata = metadata_samples - abundance_samples

    if only_abundance:
        warnings.append(
            f"Samples in abundance table but not metadata: "
            f"{', '.join(sorted(only_abundance))}"
        )
    if only_metadata:
        warnings.append(
            f"Samples in metadata but not abundance table: "
            f"{', '.join(sorted(only_metadata))}"
        )

    common = abundance_samples & metadata_samples
    if not common:
        raise ValueError(
            "No matching sample IDs between abundance table and metadata"
        )

    matched = metadata_df.loc[metadata_df.index.isin(common)]
    return matched, warnings


# ── Preprocessing ────────────────────────────────────────────────────

def normalise(abundance_df, method):
    """Apply normalisation to the abundance table."""
    if method == "none":
        return abundance_df

    if method == "tss":
        # Total-sum scaling: divide each sample by its total
        totals = abundance_df.sum(axis=0)
        totals = totals.replace(0, 1)  # avoid div by zero
        return abundance_df.div(totals, axis=1)

    if method == "hellinger":
        # Hellinger: sqrt of relative abundance
        totals = abundance_df.sum(axis=0)
        totals = totals.replace(0, 1)
        relative = abundance_df.div(totals, axis=1)
        return np.sqrt(relative)

    raise ValueError(f"Unknown normalisation method: {method}")


# ── Distance matrix ──────────────────────────────────────────────────

def compute_distance_matrix(abundance_df, metric):
    """
    Compute distance/dissimilarity matrix.

    abundance_df: features (rows) x samples (columns)
    Returns skbio DistanceMatrix.
    """
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric: {metric}")

    # Transpose: samples as rows for distance computation
    sample_matrix = abundance_df.T.values
    sample_ids = list(abundance_df.columns)

    if metric in ("braycurtis", "jaccard"):
        dm = beta_diversity(metric, sample_matrix, ids=sample_ids)
    else:
        condensed = pdist(sample_matrix, metric=metric)
        dm = DistanceMatrix(squareform(condensed), ids=sample_ids)

    # Check for NaN
    if np.any(np.isnan(dm.data)):
        raise ValueError(
            "Distance matrix contains NaN values. This often happens when "
            "samples have all-zero counts. Remove empty samples and retry."
        )

    return dm


# ── PCoA ─────────────────────────────────────────────────────────────

def run_pcoa(dm):
    """
    Run PCoA ordination.

    Returns dict with samples, coordinates, variance explained.
    """
    result = pcoa(dm)

    # Extract first 3 axes (or fewer if not enough)
    n_axes = min(3, result.samples.shape[1])
    coords = result.samples.iloc[:, :n_axes]

    variance = result.proportion_explained.iloc[:n_axes].values * 100

    return {
        "samples": list(coords.index),
        "pc1": coords.iloc[:, 0].tolist() if n_axes >= 1 else [],
        "pc2": coords.iloc[:, 1].tolist() if n_axes >= 2 else [],
        "pc3": coords.iloc[:, 2].tolist() if n_axes >= 3 else [],
        "variance_explained": [round(v, 2) for v in variance],
    }


# ── Main pipeline ────────────────────────────────────────────────────

def run_pipeline(abundance_bytes, metadata_bytes, metric, norm_method):
    """
    Run the full PCoA pipeline.

    Returns (result_dict, warnings).
    """
    all_warnings = []

    # Parse abundance table
    abundance_df, parse_warnings = parse_abundance_table(abundance_bytes)
    all_warnings.extend(parse_warnings)

    # Parse metadata if provided
    metadata = None
    metadata_columns = []
    if metadata_bytes:
        metadata_df = parse_metadata(metadata_bytes)
        metadata_df, meta_warnings = match_metadata(abundance_df, metadata_df)
        all_warnings.extend(meta_warnings)

        # Reorder metadata to match abundance columns
        common_samples = [
            s for s in abundance_df.columns if s in metadata_df.index
        ]
        metadata_df = metadata_df.loc[common_samples]
        # Also filter abundance to common samples if metadata was provided
        abundance_df = abundance_df[common_samples]

        metadata_columns = list(metadata_df.columns)
        metadata = {
            col: metadata_df[col].tolist() for col in metadata_columns
        }

    # Normalise
    abundance_df = normalise(abundance_df, norm_method)

    # Distance matrix
    dm = compute_distance_matrix(abundance_df, metric)

    # PCoA
    pcoa_result = run_pcoa(dm)

    # Distance matrix as nested list for download
    dm_data = {
        "ids": list(dm.ids),
        "matrix": dm.data.tolist(),
    }

    return {
        **pcoa_result,
        "distance_metric": metric,
        "normalisation": norm_method,
        "num_samples": len(pcoa_result["samples"]),
        "num_features": abundance_df.shape[0],
        "metadata": metadata,
        "metadata_columns": metadata_columns,
        "distance_matrix": dm_data,
        "warnings": all_warnings,
    }
