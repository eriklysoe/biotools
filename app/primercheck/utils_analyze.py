"""Analyze utilities: Tm calculation, molecular weight, extinction coefficient."""

import math

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False

from .utils_shared import gc_content, find_runs, gc_clamp_check


# ── Molecular weights (g/mol) ──────────────────────────────────────
DNA_MW = {"A": 313.21, "T": 304.19, "G": 329.21, "C": 289.18}
RNA_MW = {"A": 329.21, "U": 306.17, "G": 345.21, "C": 305.18}


# ── Nearest-neighbour extinction coefficients at 260 nm ───────────
# Published NN parameters for calculating molar extinction coefficients
DNA_NN_EXTINCTION = {
    "AA": 27400, "AT": 22800, "AG": 25000, "AC": 21200,
    "TA": 23400, "TT": 27400, "TG": 19000, "TC": 21200,
    "GA": 25200, "GT": 22800, "GG": 21600, "GC": 17600,
    "CA": 21200, "CT": 25200, "CG": 18000, "CC": 21600,
}
DNA_MONO_EXTINCTION = {"A": 15400, "T": 8700, "G": 11500, "C": 7400}

RNA_NN_EXTINCTION = {
    "AA": 27400, "AU": 24000, "AG": 25000, "AC": 21200,
    "UA": 24600, "UU": 19600, "UG": 20000, "UC": 21200,
    "GA": 25200, "GU": 21200, "GG": 21600, "GC": 17600,
    "CA": 21200, "CU": 25200, "CG": 18000, "CC": 21600,
}
RNA_MONO_EXTINCTION = {"A": 15400, "U": 9900, "G": 11500, "C": 7200}


# ── SantaLucia 1998 nearest-neighbour parameters ──────────────────
# delta-H (cal/mol) and delta-S (cal/mol·K) for DNA/DNA
NN_DH = {
    "AA": -7900, "AT": -7200, "AG": -7800, "AC": -8400,
    "TA": -7200, "TT": -7900, "TG": -8500, "TC": -8200,
    "GA": -8200, "GT": -8400, "GG": -8000, "GC": -9800,
    "CA": -8500, "CT": -7800, "CG": -10600, "CC": -8000,
}
NN_DS = {
    "AA": -22.2, "AT": -20.4, "AG": -21.0, "AC": -22.4,
    "TA": -21.3, "TT": -22.2, "TG": -22.7, "TC": -22.2,
    "GA": -22.2, "GT": -22.4, "GG": -19.9, "GC": -24.4,
    "CA": -22.7, "CT": -21.0, "CG": -27.2, "CC": -19.9,
}

# Initiation parameters
INIT_DH = {"G": 100, "A": 2300}  # 5' base
INIT_DS = {"G": -2.8, "A": 4.1}

R = 1.987  # gas constant cal/(mol·K)


def calc_tm_wallace(seq):
    """Wallace rule: Tm = 2(A+T) + 4(G+C). Only valid for <20 bp."""
    seq = seq.upper().replace("U", "T")
    at = sum(1 for c in seq if c in "AT")
    gc = sum(1 for c in seq if c in "GC")
    tm = 2 * at + 4 * gc
    warn = len(seq) > 20
    return round(tm, 1), warn


def calc_tm_basic_nn(seq, oligo_conc_nM=250):
    """Basic nearest-neighbour Tm (SantaLucia 1998), no salt correction."""
    seq = seq.upper().replace("U", "T")
    if len(seq) < 2:
        return 0.0

    # Sum NN parameters
    dH = 0.0
    dS = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in NN_DH:
            dH += NN_DH[dinuc]
            dS += NN_DS[dinuc]

    # Initiation
    first = seq[0]
    last = seq[-1]
    dH += INIT_DH.get(first, INIT_DH["A"]) + INIT_DH.get(last, INIT_DH["A"])
    dS += INIT_DS.get(first, INIT_DS["A"]) + INIT_DS.get(last, INIT_DS["A"])

    # Self-complementary correction
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    if seq == rc:
        dS += -1.4

    oligo_conc_M = oligo_conc_nM * 1e-9
    # Tm = dH / (dS + R * ln(Ct/4))   for non-self-complementary
    ct = oligo_conc_M
    if seq != rc:
        tm = dH / (dS + R * math.log(ct / 4)) - 273.15
    else:
        tm = dH / (dS + R * math.log(ct)) - 273.15

    return round(tm, 1)


def calc_tm_salt_corrected(seq, na_mM=50, mg_mM=0, dntp_mM=0,
                           oligo_nM=250):
    """Salt-corrected Tm using primer3-py."""
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq = seq.upper().replace("U", "T")
    tm = primer3.calc_tm(
        seq,
        mv_conc=na_mM,
        dv_conc=mg_mM,
        dntp_conc=dntp_mM,
        dna_conc=oligo_nM,
    )
    return round(tm, 1)


def calc_tm_owczarzy_mg(seq, na_mM=50, mg_mM=3, dntp_mM=0, oligo_nM=250):
    """
    Owczarzy 2008 Mg2+ correction.
    Uses primer3-py which internally applies divalent salt correction.
    """
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    seq = seq.upper().replace("U", "T")
    tm = primer3.calc_tm(
        seq,
        mv_conc=na_mM,
        dv_conc=mg_mM,
        dntp_conc=dntp_mM,
        dna_conc=oligo_nM,
    )
    return round(tm, 1)


def calc_molecular_weight(seq, oligo_type="DNA"):
    """Calculate molecular weight in g/mol."""
    seq = seq.upper()
    weights = DNA_MW if oligo_type == "DNA" else RNA_MW

    mw = sum(weights.get(c, 0) for c in seq)

    if oligo_type == "DNA":
        # Subtract water for phosphodiester bonds and add 5'-OH
        mw -= 61.96
    else:
        mw -= 61.96

    return round(mw, 1)


def calc_extinction_coefficient(seq, oligo_type="DNA"):
    """
    Calculate molar extinction coefficient at 260 nm using NN method.
    Returns epsilon in L/(mol·cm).
    """
    seq = seq.upper()
    if oligo_type == "DNA":
        seq = seq.replace("U", "T")
        nn_ext = DNA_NN_EXTINCTION
        mono_ext = DNA_MONO_EXTINCTION
    else:
        seq = seq.replace("T", "U")
        nn_ext = RNA_NN_EXTINCTION
        mono_ext = RNA_MONO_EXTINCTION

    if len(seq) < 2:
        return mono_ext.get(seq, 0) if seq else 0

    # NN method: sum of NN - sum of individual (except first and last)
    nn_sum = 0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        nn_sum += nn_ext.get(dinuc, 0)

    # Subtract internal individual contributions
    mono_sum = 0
    for i in range(1, len(seq) - 1):
        mono_sum += mono_ext.get(seq[i], 0)

    epsilon = nn_sum - mono_sum
    return epsilon


def analyze_oligo(seq, oligo_type="DNA", target_type="DNA",
                  na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """
    Full property analysis of an oligo.
    Returns dict with all calculated properties and warnings.
    """
    seq_for_calc = seq.upper()
    if oligo_type == "RNA":
        seq_for_calc = seq_for_calc.replace("T", "U")
    else:
        seq_for_calc = seq_for_calc.replace("U", "T")

    # For Tm calculations, always use DNA (T) form
    seq_dna = seq_for_calc.replace("U", "T")

    length = len(seq_for_calc)
    gc = gc_content(seq_for_calc)
    mw = calc_molecular_weight(seq_for_calc, oligo_type)
    ext_coeff = calc_extinction_coefficient(seq_for_calc, oligo_type)

    # Tm calculations
    tm_wallace, wallace_warn = calc_tm_wallace(seq_dna)

    tm_basic_nn = calc_tm_basic_nn(seq_dna, oligo_nM)

    tm_salt = None
    tm_mg = None
    if HAS_PRIMER3:
        tm_salt = calc_tm_salt_corrected(seq_dna, na_mM, mg_mM, dntp_mM, oligo_nM)
        if mg_mM > 0:
            tm_mg = calc_tm_owczarzy_mg(seq_dna, na_mM, mg_mM, dntp_mM, oligo_nM)

    recommended_tm = tm_salt if tm_salt is not None else tm_basic_nn

    # Derived properties
    ug_per_od = round(mw / ext_coeff, 1) if ext_coeff > 0 else None
    nmol_per_od = round(1e6 / ext_coeff, 1) if ext_coeff > 0 else None

    # Warnings
    warnings = []
    if gc < 40:
        warnings.append({"type": "gc_low", "message": f"GC content low ({gc}% < 40%)"})
    elif gc > 60:
        warnings.append({"type": "gc_high", "message": f"GC content high ({gc}% > 60%)"})

    if length < 18:
        warnings.append({"type": "length_short", "message": f"Length short ({length} < 18 nt)"})
    elif length > 30:
        warnings.append({"type": "length_long", "message": f"Length long ({length} > 30 nt)"})

    runs = find_runs(seq_for_calc, 4)
    if runs:
        for r in runs:
            warnings.append({
                "type": "base_run",
                "message": f"Run of {r['length']}x {r['base']} at position {r['start']+1}",
            })

    clamp_ok, clamp_gc = gc_clamp_check(seq_for_calc)
    if not clamp_ok:
        if clamp_gc < 2:
            warnings.append({"type": "gc_clamp_weak", "message": f"3' GC clamp weak ({clamp_gc} G/C in last 5 bases)"})
        else:
            warnings.append({"type": "gc_clamp_strong", "message": f"3' GC clamp too strong ({clamp_gc} G/C in last 5 bases)"})

    return {
        "sequence": seq_for_calc,
        "oligo_type": oligo_type,
        "length": length,
        "gc_content": gc,
        "molecular_weight": mw,
        "extinction_coefficient": ext_coeff,
        "ug_per_od260": ug_per_od,
        "nmol_per_od260": nmol_per_od,
        "tm_wallace": tm_wallace,
        "tm_wallace_warning": wallace_warn,
        "tm_basic_nn": tm_basic_nn,
        "tm_salt_corrected": tm_salt,
        "tm_mg_corrected": tm_mg,
        "recommended_tm": recommended_tm,
        "warnings": warnings,
        "gc_clamp_ok": clamp_ok,
        "gc_clamp_count": clamp_gc,
    }
