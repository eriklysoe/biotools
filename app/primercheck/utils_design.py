"""Primer3 design: design optimal primer pairs from a target sequence."""

import re

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


def design_primers(template, target_start=None, target_length=None,
                   excluded_regions=None,
                   min_size=100, max_size=500,
                   min_length=18, opt_length=20, max_length=25,
                   min_tm=57.0, opt_tm=60.0, max_tm=63.0,
                   min_gc=40.0, max_gc=60.0,
                   max_tm_diff=3.0, num_pairs=5,
                   na_mM=50, mg_mM=0, dntp_mM=0, oligo_nM=250):
    """
    Design primer pairs using Primer3.
    Returns dict with primer pairs and explain text.
    """
    if not HAS_PRIMER3:
        raise RuntimeError("primer3-py is not installed")

    template = re.sub(r"[^A-Za-z]", "", template).upper().replace("U", "T")

    if len(template) < 20:
        raise ValueError("Template too short for primer design (< 20 bp)")

    num_pairs = min(max(1, num_pairs), 10)

    seq_args = {
        "SEQUENCE_ID": "target",
        "SEQUENCE_TEMPLATE": template,
    }

    if target_start is not None and target_length is not None:
        seq_args["SEQUENCE_TARGET"] = [int(target_start), int(target_length)]

    if excluded_regions:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = excluded_regions

    global_args = {
        "PRIMER_OPT_SIZE": opt_length,
        "PRIMER_MIN_SIZE": min_length,
        "PRIMER_MAX_SIZE": max_length,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": min_tm,
        "PRIMER_MAX_TM": max_tm,
        "PRIMER_MIN_GC": min_gc,
        "PRIMER_MAX_GC": max_gc,
        "PRIMER_MAX_DIFF_TM": max_tm_diff,
        "PRIMER_NUM_RETURN": num_pairs,
        "PRIMER_PRODUCT_SIZE_RANGE": [[min_size, max_size]],
        "PRIMER_MAX_HAIRPIN_TH": 24.0,
        "PRIMER_MAX_SELF_ANY_TH": 45.0,
        "PRIMER_MAX_SELF_END_TH": 35.0,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_SALT_MONOVALENT": float(na_mM),
        "PRIMER_SALT_DIVALENT": float(mg_mM),
        "PRIMER_DNTP_CONC": float(dntp_mM),
        "PRIMER_DNA_CONC": float(oligo_nM),
        "PRIMER_EXPLAIN_FLAG": 1,
    }

    result = primer3.design_primers(seq_args, global_args)

    # Parse results
    num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    pairs = []

    for i in range(num_returned):
        fwd_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        rev_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        fwd_pos = result.get(f"PRIMER_LEFT_{i}", (0, 0))
        rev_pos = result.get(f"PRIMER_RIGHT_{i}", (0, 0))
        product_size = result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)

        pair = {
            "rank": i + 1,
            "penalty": round(result.get(f"PRIMER_PAIR_{i}_PENALTY", 0), 3),
            "forward": {
                "sequence": fwd_seq,
                "start": fwd_pos[0] + 1,  # 1-based
                "length": fwd_pos[1],
                "tm": round(result.get(f"PRIMER_LEFT_{i}_TM", 0), 1),
                "gc_percent": round(result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0), 1),
                "self_any_th": round(result.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH", 0), 2),
                "self_end_th": round(result.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0), 2),
                "hairpin_th": round(result.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH", 0), 2),
            },
            "reverse": {
                "sequence": rev_seq,
                "start": rev_pos[0] + 1,  # 1-based
                "length": rev_pos[1],
                "tm": round(result.get(f"PRIMER_RIGHT_{i}_TM", 0), 1),
                "gc_percent": round(result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0), 1),
                "self_any_th": round(result.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH", 0), 2),
                "self_end_th": round(result.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0), 2),
                "hairpin_th": round(result.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH", 0), 2),
            },
            "product_size": product_size,
            "compl_any_th": round(result.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0), 2),
            "compl_end_th": round(result.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0), 2),
        }
        pairs.append(pair)

    # Explain text
    explain = {}
    for key in ["PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN",
                "PRIMER_PAIR_EXPLAIN", "PRIMER_INTERNAL_EXPLAIN"]:
        if key in result:
            explain[key] = result[key]

    # Suggestions if no pairs found
    suggestions = []
    if num_returned == 0:
        suggestions.append("No primer pairs found. Try:")
        suggestions.append("- Widening the Tm range (e.g., 55-65\u00b0C)")
        suggestions.append("- Expanding the amplicon size range")
        suggestions.append("- Relaxing GC content limits")
        suggestions.append("- Increasing primer length range")
        suggestions.append("- Removing target region constraints")

    return {
        "pairs": pairs,
        "num_returned": num_returned,
        "explain": explain,
        "suggestions": suggestions,
        "template_length": len(template),
    }
