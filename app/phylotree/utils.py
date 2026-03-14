"""PhyloTree pipeline: alignment and tree building."""

import io
import shutil
import subprocess
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceMatrix,
    DistanceTreeConstructor,
)

SUBPROCESS_TIMEOUT = 300  # 5 minutes


# ── Binary availability ──────────────────────────────────────────────

def check_binary(name):
    """Return True if *name* is on PATH."""
    return shutil.which(name) is not None


def get_available_methods():
    """Return which alignment / tree methods are available."""
    return {
        "alignment": {
            "muscle": check_binary("muscle"),
            "clustalo": check_binary("clustalo"),
            "builtin": True,
        },
        "tree": {
            "nj": True,
            "upgma": True,
            "ml": check_binary("iqtree2") or check_binary("iqtree"),
        },
    }


def _iqtree_bin():
    for name in ("iqtree2", "iqtree"):
        if check_binary(name):
            return name
    return None


# ── Alignment ────────────────────────────────────────────────────────

def align_muscle(fasta_path, output_path):
    """Align with MUSCLE (v5 uses -align/-output)."""
    r = subprocess.run(
        ["muscle", "-align", str(fasta_path), "-output", str(output_path)],
        capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT,
    )
    if r.returncode != 0:
        raise RuntimeError(f"MUSCLE failed:\n{r.stderr.strip()}")


def align_clustalo(fasta_path, output_path):
    """Align with Clustal Omega."""
    r = subprocess.run(
        ["clustalo", "-i", str(fasta_path), "-o", str(output_path), "--force"],
        capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT,
    )
    if r.returncode != 0:
        raise RuntimeError(f"ClustalOmega failed:\n{r.stderr.strip()}")


def _pairwise_distance_matrix(records):
    """Compute a distance matrix from unaligned sequences (percent identity)."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    names = [r.id for r in records]
    n = len(records)
    matrix = [[0.0]]
    for i in range(1, n):
        row = []
        for j in range(i):
            score = aligner.score(records[i].seq, records[j].seq)
            max_len = max(len(records[i].seq), len(records[j].seq))
            identity = score / max_len if max_len > 0 else 1.0
            row.append(max(0.001, 1.0 - identity))
        row.append(0.0)
        matrix.append(row)
    return DistanceMatrix(names, matrix)


# ── Tree building ────────────────────────────────────────────────────

EXPORT_FORMATS = {
    "newick":   {"ext": ".nwk",  "label": "Newick (.nwk)",        "mime": "text/plain"},
    "nexus":    {"ext": ".nex",  "label": "Nexus (.nex)",         "mime": "text/plain"},
    "phyloxml": {"ext": ".xml",  "label": "PhyloXML (.xml)",      "mime": "application/xml"},
    "nexml":    {"ext": ".xml",  "label": "NeXML (.xml)",         "mime": "application/xml"},
}


def newick_to_format(newick_str, fmt):
    """Convert a Newick string to another tree format via Biopython."""
    if fmt not in EXPORT_FORMATS:
        raise ValueError(f"Unsupported export format: {fmt}")

    tree = Phylo.read(io.StringIO(newick_str), "newick")
    buf = io.StringIO()
    Phylo.write(tree, buf, fmt)
    return buf.getvalue()


def _tree_to_newick(tree):
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()


def _build_distance_tree(dm, method):
    constructor = DistanceTreeConstructor()
    if method == "nj":
        tree = constructor.nj(dm)
    else:
        tree = constructor.upgma(dm)
    # Strip Biopython's default "Inner*" names from internal nodes
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None
    return tree


def _run_iqtree(alignment_path, work_dir):
    """Run IQ-TREE 2 and return (newick, best_model)."""
    binary = _iqtree_bin()
    if not binary:
        raise RuntimeError("IQ-TREE binary not found in container")

    prefix = str(work_dir / "iqtree")
    r = subprocess.run(
        [binary, "-s", str(alignment_path), "-m", "MFP", "-B", "1000",
         "--prefix", prefix],
        capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT,
    )
    if r.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed:\n{r.stderr.strip()}")

    treefile = Path(f"{prefix}.treefile")
    if not treefile.exists():
        raise RuntimeError("IQ-TREE did not produce a tree file")
    newick = treefile.read_text().strip()

    model = None
    logfile = Path(f"{prefix}.iqtree")
    if logfile.exists():
        for line in logfile.read_text().splitlines():
            if "Best-fit model:" in line:
                model = line.split("Best-fit model:")[1].strip().split()[0]
                break

    return newick, model


# ── Main pipeline ────────────────────────────────────────────────────

def run_pipeline(fasta_path, align_method, tree_method, work_dir):
    """
    Run full phylogenetic pipeline.

    Returns dict with: newick, method, model, bootstrap,
                       num_sequences, alignment_length.
    """
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    num_sequences = len(records)

    if num_sequences < 3:
        raise ValueError(
            "At least 3 sequences are required for phylogenetic analysis"
        )

    if tree_method == "ml" and align_method == "builtin":
        raise ValueError(
            "Maximum Likelihood requires MUSCLE or ClustalOmega for alignment"
        )

    # ── Alignment ──
    alignment_path = work_dir / "aligned.fasta"
    alignment_length = None

    if align_method == "builtin":
        # Skip MSA; compute distances directly from unaligned sequences
        dm = _pairwise_distance_matrix(records)
        alignment_length = None
    else:
        if align_method == "muscle":
            align_muscle(fasta_path, alignment_path)
        elif align_method == "clustalo":
            align_clustalo(fasta_path, alignment_path)
        else:
            raise ValueError(f"Unknown alignment method: {align_method}")

        alignment = AlignIO.read(str(alignment_path), "fasta")
        alignment_length = alignment.get_alignment_length()

    # ── Tree ──
    if tree_method in ("nj", "upgma"):
        if align_method != "builtin":
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)

        tree = _build_distance_tree(dm, tree_method)
        newick = _tree_to_newick(tree)

        return {
            "newick": newick,
            "method": tree_method,
            "model": None,
            "bootstrap": False,
            "num_sequences": num_sequences,
            "alignment_length": alignment_length,
        }

    elif tree_method == "ml":
        newick, model = _run_iqtree(alignment_path, work_dir)
        return {
            "newick": newick,
            "method": "ml",
            "model": model,
            "bootstrap": True,
            "num_sequences": num_sequences,
            "alignment_length": alignment_length,
        }

    else:
        raise ValueError(f"Unknown tree method: {tree_method}")
