"""Tm mismatch heatmap: calculate Tm for every single-nucleotide mismatch."""

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False

BASES = ["A", "T", "G", "C"]


def calc_tm_mismatch_heatmap(seq, na_mM=50, mg_mM=0, dntp_mM=0,
                              oligo_nM=250):
    """
    For each position in the oligo, substitute each possible base and
    calculate the resulting Tm.

    Returns:
        dict with 'positions', 'bases', 'heatmap' (2D list), 'reference_tm',
        'table' (list of dicts for tabular display)
    """
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq = seq.upper().replace("U", "T")

    # Reference Tm (no mismatch)
    ref_tm = primer3.calc_tm(
        seq, mv_conc=na_mM, dv_conc=mg_mM,
        dntp_conc=dntp_mM, dna_conc=oligo_nM,
    )
    ref_tm = round(ref_tm, 1)

    positions = list(range(1, len(seq) + 1))
    heatmap = []  # bases × positions
    table = []

    for base in BASES:
        row = []
        for i, orig_base in enumerate(seq):
            if base == orig_base:
                tm = ref_tm
                delta = 0.0
            else:
                mutant = seq[:i] + base + seq[i+1:]
                tm = primer3.calc_tm(
                    mutant, mv_conc=na_mM, dv_conc=mg_mM,
                    dntp_conc=dntp_mM, dna_conc=oligo_nM,
                )
                tm = round(tm, 1)
                delta = round(tm - ref_tm, 1)

            row.append(tm)
            table.append({
                "position": i + 1,
                "original_base": orig_base,
                "substituted_base": base,
                "tm": tm,
                "delta_tm": delta,
                "is_match": base == orig_base,
            })
        heatmap.append(row)

    return {
        "sequence": seq,
        "positions": positions,
        "bases": BASES,
        "heatmap": heatmap,
        "reference_tm": ref_tm,
        "table": table,
    }
