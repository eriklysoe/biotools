"""Structure utilities: hairpin, self-dimer, hetero-dimer via primer3-py."""

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


def _thermo_to_dict(result):
    """Convert primer3 ThermoResult to dict."""
    return {
        "dg": round(result.dg / 1000, 2),      # cal → kcal
        "dh": round(result.dh / 1000, 2),
        "ds": round(result.ds, 2),
        "tm": round(result.tm, 1),
        "structure_found": result.structure_found,
    }


def calc_hairpin(seq, na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """Calculate hairpin thermodynamics."""
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq = seq.upper().replace("U", "T")
    result = primer3.calc_hairpin(
        seq,
        mv_conc=na_mM,
        dv_conc=mg_mM,
        dntp_conc=dntp_mM,
        dna_conc=oligo_nM,
    )
    data = _thermo_to_dict(result)

    # Warning if dG < -2.0 kcal/mol
    data["warning"] = data["dg"] < -2.0
    data["warning_message"] = (
        f"Significant hairpin: \u0394G = {data['dg']} kcal/mol"
        if data["warning"] else None
    )

    return data


def calc_homodimer(seq, na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """Calculate self-dimer (homodimer) thermodynamics."""
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq = seq.upper().replace("U", "T")
    result = primer3.calc_homodimer(
        seq,
        mv_conc=na_mM,
        dv_conc=mg_mM,
        dntp_conc=dntp_mM,
        dna_conc=oligo_nM,
    )
    data = _thermo_to_dict(result)

    # Warning if dG < -6.0 kcal/mol
    data["warning"] = data["dg"] < -6.0
    data["warning_message"] = (
        f"Significant self-dimer: \u0394G = {data['dg']} kcal/mol"
        if data["warning"] else None
    )

    # Check if 3' end is involved
    data["three_prime_involved"] = _check_3prime_dimer(seq, seq)

    return data


def calc_heterodimer(seq1, seq2, na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """Calculate hetero-dimer thermodynamics between two oligos."""
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq1 = seq1.upper().replace("U", "T")
    seq2 = seq2.upper().replace("U", "T")

    result = primer3.calc_heterodimer(
        seq1, seq2,
        mv_conc=na_mM,
        dv_conc=mg_mM,
        dntp_conc=dntp_mM,
        dna_conc=oligo_nM,
    )
    data = _thermo_to_dict(result)

    # Warning if dG < -6.0 kcal/mol
    data["warning"] = data["dg"] < -6.0
    data["warning_message"] = (
        f"Significant hetero-dimer: \u0394G = {data['dg']} kcal/mol"
        if data["warning"] else None
    )

    # Check if 3' end of either oligo is involved
    data["three_prime_involved_seq1"] = _check_3prime_dimer(seq1, seq2)
    data["three_prime_involved_seq2"] = _check_3prime_dimer(seq2, seq1)

    return data


def _check_3prime_dimer(seq1, seq2, tail_len=5):
    """
    Heuristic check if the 3' end of seq1 can base-pair with seq2.
    Returns True if >=3 of the last tail_len bases can pair.
    """
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    tail = seq1[-tail_len:]
    tail_rc = "".join(complement.get(c, "N") for c in reversed(tail))

    # Check if tail_rc appears anywhere in seq2
    for i in range(len(seq2) - len(tail_rc) + 1):
        matches = sum(1 for a, b in zip(tail_rc, seq2[i:i+len(tail_rc)]) if a == b)
        if matches >= 3:
            return True

    return False


def generate_dimer_ascii(seq1, seq2, is_homodimer=False):
    """
    Generate an ASCII representation of dimer alignment.
    Finds the best alignment and shows base pairing.
    """
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    seq2_rc = seq2[::-1]

    best_score = 0
    best_offset = 0

    # Try all offsets
    for offset in range(-(len(seq2_rc) - 1), len(seq1)):
        score = 0
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(seq2_rc):
                if seq1[i] == complement.get(seq2_rc[j], ""):
                    score += 1
        if score > best_score:
            best_score = score
            best_offset = offset

    if best_score == 0:
        return {"top": seq1, "middle": " " * len(seq1), "bottom": seq2_rc}

    # Build alignment at best offset
    top_line = ""
    mid_line = ""
    bot_line = ""

    start = min(0, best_offset)
    end = max(len(seq1), best_offset + len(seq2_rc))

    for pos in range(start, end):
        i = pos  # position in seq1
        j = pos - best_offset  # position in seq2_rc

        c1 = seq1[i] if 0 <= i < len(seq1) else " "
        c2 = seq2_rc[j] if 0 <= j < len(seq2_rc) else " "

        top_line += c1
        bot_line += c2

        if c1 != " " and c2 != " " and c1 == complement.get(c2, ""):
            mid_line += "|"
        elif c1 != " " and c2 != " ":
            mid_line += "\u00b7"
        else:
            mid_line += " "

    label_top = "5'" if not is_homodimer else "5'"
    label_bot = "3'" if not is_homodimer else "3'"

    return {
        "top": f"{label_top} {top_line} 3'",
        "middle": f"   {mid_line}",
        "bottom": f"3' {bot_line} 5'",
    }
