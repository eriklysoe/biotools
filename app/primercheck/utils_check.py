"""In-silico PCR: check primers against a template sequence."""

import re

from .utils_shared import reverse_complement, iupac_to_regex, validate_sequence
from .utils_analyze import analyze_oligo
from .utils_structure import calc_hairpin, calc_homodimer, calc_heterodimer


def find_binding_sites(primer_seq, template_seq, max_mismatches=2, strand="+"):
    """
    Search for primer binding sites on a template, allowing mismatches.
    Returns list of dicts with position, mismatches, strand info.
    """
    primer_seq = primer_seq.upper().replace("U", "T")
    template_seq = template_seq.upper().replace("U", "T")
    plen = len(primer_seq)
    sites = []

    for i in range(len(template_seq) - plen + 1):
        segment = template_seq[i:i+plen]
        mismatches = sum(1 for a, b in zip(primer_seq, segment) if a != b)
        if mismatches <= max_mismatches:
            sites.append({
                "position": i + 1,  # 1-based
                "strand": strand,
                "mismatches": mismatches,
                "matched_sequence": segment,
            })

    return sites


def find_amplicons(fwd_sites, rev_sites, template_seq, fwd_len, rev_len,
                   max_amplicon=10000):
    """
    Given binding sites for forward and reverse primers,
    find valid amplicon combinations.
    """
    amplicons = []

    for fs in fwd_sites:
        for rs in rev_sites:
            # Forward on + strand, reverse on - strand
            # Forward binds at fs['position'], reverse binds at rs['position']
            # Amplicon goes from fwd start to rev end (inclusive)

            if fs["strand"] == "+" and rs["strand"] == "-":
                amp_start = fs["position"] - 1  # 0-based
                amp_end = rs["position"] - 1 + rev_len  # 0-based exclusive
                amp_size = amp_end - amp_start

                if 0 < amp_size <= max_amplicon:
                    amp_seq = template_seq[amp_start:amp_end]
                    amplicons.append({
                        "fwd_position": fs["position"],
                        "rev_position": rs["position"],
                        "fwd_mismatches": fs["mismatches"],
                        "rev_mismatches": rs["mismatches"],
                        "fwd_strand": "+",
                        "rev_strand": "-",
                        "amplicon_size": amp_size,
                        "amplicon_sequence": amp_seq,
                        "fwd_primer_end": amp_start + fwd_len,
                        "rev_primer_start": rs["position"] - 1,
                    })

    # Sort by amplicon size
    amplicons.sort(key=lambda a: a["amplicon_size"])
    return amplicons


def run_check(fwd_seq, rev_seq, template_text, max_mismatches=2,
              oligo_type="DNA", na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """
    Full primer check pipeline:
    1. Analyze both primers
    2. Check hairpin, self-dimer for each
    3. Check hetero-dimer between them
    4. Find binding sites on template
    5. Calculate amplicons
    """
    fwd_seq = fwd_seq.upper().replace("U", "T")
    rev_seq = rev_seq.upper().replace("U", "T")
    template_seq = re.sub(r"[^A-Za-z]", "", template_text).upper().replace("U", "T")

    if len(template_seq) < 3:
        raise ValueError("Template sequence too short (< 3 bases)")

    salt_args = dict(na_mM=na_mM, mg_mM=mg_mM, dntp_mM=dntp_mM, oligo_nM=oligo_nM)

    # Analyze primers
    fwd_analysis = analyze_oligo(fwd_seq, oligo_type=oligo_type, **salt_args)
    rev_analysis = analyze_oligo(rev_seq, oligo_type=oligo_type, **salt_args)

    # Structure analysis
    fwd_hairpin = calc_hairpin(fwd_seq, **salt_args)
    rev_hairpin = calc_hairpin(rev_seq, **salt_args)
    fwd_homodimer = calc_homodimer(fwd_seq, **salt_args)
    rev_homodimer = calc_homodimer(rev_seq, **salt_args)
    heterodimer = calc_heterodimer(fwd_seq, rev_seq, **salt_args)

    # Tm difference
    fwd_tm = fwd_analysis["recommended_tm"] or 0
    rev_tm = rev_analysis["recommended_tm"] or 0
    tm_diff = round(abs(fwd_tm - rev_tm), 1)

    # Find binding sites on both strands
    template_rc = reverse_complement(template_seq, "DNA")

    fwd_sites_plus = find_binding_sites(fwd_seq, template_seq, max_mismatches, "+")
    rev_rc = reverse_complement(rev_seq, "DNA")
    rev_sites_minus = find_binding_sites(rev_rc, template_seq, max_mismatches, "-")

    # Also try reverse complement orientations for completeness
    # Forward binding on minus strand
    fwd_rc = reverse_complement(fwd_seq, "DNA")
    fwd_sites_minus = find_binding_sites(fwd_rc, template_seq, max_mismatches, "-")
    # Reverse binding on plus strand
    rev_sites_plus = find_binding_sites(rev_seq, template_seq, max_mismatches, "+")

    # Calculate amplicons (fwd on +, rev on -)
    amplicons = find_amplicons(
        fwd_sites_plus, rev_sites_minus, template_seq,
        len(fwd_seq), len(rev_seq),
    )

    # Warnings
    warnings = []
    if tm_diff > 3:
        warnings.append(f"Tm difference between primers is {tm_diff}\u00b0C (> 3\u00b0C)")
    if heterodimer["warning"]:
        warnings.append(heterodimer["warning_message"])
    if not amplicons:
        warnings.append(
            "No amplicons found. Try increasing the mismatch threshold."
        )

    return {
        "forward": {
            "sequence": fwd_seq,
            "analysis": fwd_analysis,
            "hairpin": fwd_hairpin,
            "homodimer": fwd_homodimer,
        },
        "reverse": {
            "sequence": rev_seq,
            "analysis": rev_analysis,
            "hairpin": rev_hairpin,
            "homodimer": rev_homodimer,
        },
        "heterodimer": heterodimer,
        "tm_difference": tm_diff,
        "template_length": len(template_seq),
        "fwd_binding_sites": fwd_sites_plus,
        "rev_binding_sites": rev_sites_minus,
        "amplicons": amplicons,
        "warnings": warnings,
    }
