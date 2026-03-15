"""Shared utilities: IUPAC validation, sequence parsing, complement."""

import re

# Valid IUPAC nucleotide characters
IUPAC_DNA = set("ACGTURYKMSWBDHVN")
IUPAC_RNA = set("ACGURYKMSWBDHVN")

# IUPAC ambiguity codes → regex character classes
IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "[AG]", "Y": "[CT]", "K": "[GT]", "M": "[AC]",
    "S": "[GC]", "W": "[AT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

DNA_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKWSVHDBN")
RNA_COMPLEMENT = str.maketrans("ACGURYKMSWBDHVN", "UGCAYRMKWSVHDBN")


def validate_sequence(seq, oligo_type="DNA"):
    """Validate and clean an oligo sequence. Returns (cleaned, errors)."""
    if not seq or not seq.strip():
        return "", ["Empty sequence"]

    cleaned = re.sub(r"\s+", "", seq.upper())
    allowed = IUPAC_DNA if oligo_type == "DNA" else IUPAC_RNA

    invalid = set(cleaned) - allowed
    if invalid:
        return cleaned, [f"Invalid characters: {', '.join(sorted(invalid))}"]

    return cleaned, []


def convert_oligo_type(seq, from_type, to_type):
    """Convert between DNA (T) and RNA (U)."""
    if from_type == to_type:
        return seq
    if to_type == "RNA":
        return seq.replace("T", "U")
    return seq.replace("U", "T")


def reverse_complement(seq, oligo_type="DNA"):
    """Return reverse complement of a sequence."""
    table = DNA_COMPLEMENT if oligo_type == "DNA" else RNA_COMPLEMENT
    return seq.translate(table)[::-1]


def iupac_to_regex(seq):
    """Convert an IUPAC sequence to a regex pattern."""
    return "".join(IUPAC_EXPAND.get(c, c) for c in seq.upper())


def gc_content(seq):
    """Calculate GC content as a percentage."""
    seq = seq.upper().replace("U", "T")
    gc = sum(1 for c in seq if c in "GC")
    return round(gc / len(seq) * 100, 1) if seq else 0


def find_runs(seq, min_length=4):
    """Find runs of identical bases >= min_length."""
    runs = []
    for m in re.finditer(r"(.)\1{" + str(min_length - 1) + r",}", seq):
        runs.append({"base": m.group(1), "start": m.start(), "length": len(m.group())})
    return runs


def gc_clamp_check(seq, window=5, min_gc=2, max_gc=3):
    """Check 3' GC clamp. Returns (ok, gc_count_in_window)."""
    tail = seq[-window:] if len(seq) >= window else seq
    gc = sum(1 for c in tail if c in "GC")
    ok = min_gc <= gc <= max_gc
    return ok, gc


def parse_fasta_text(text):
    """Parse FASTA-formatted text, return list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []

    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences.append((current_name, "".join(current_seq)))
            current_name = line[1:].strip() or "unnamed"
            current_seq = []
        else:
            current_seq.append(re.sub(r"\s+", "", line.upper()))

    if current_name is not None:
        sequences.append((current_name, "".join(current_seq)))
    elif current_seq or text.strip():
        # Plain sequence without FASTA header
        seq = re.sub(r"\s+", "", text.upper())
        sequences.append(("sequence", seq))

    return sequences
