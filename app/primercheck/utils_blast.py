"""NCBI BLAST: submit oligo to BLAST and retrieve results."""

import time
import re
import xml.etree.ElementTree as ET

import requests


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
MAX_POLL_TIME = 120  # seconds
POLL_INTERVAL = 5  # seconds


def submit_blast(seq, database="nt", program="blastn"):
    """Submit a BLAST query and return the RID (Request ID)."""
    seq = seq.upper().replace("U", "T")

    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": seq,
        "WORD_SIZE": "7",
        "EXPECT": "1000",
        "REWARD": "1",
        "PENALTY": "-3",
        "FILTER": "L",
        "FORMAT_TYPE": "XML",
    }

    resp = requests.post(BLAST_URL, data=params, timeout=30)
    resp.raise_for_status()

    # Extract RID from response
    rid_match = re.search(r"RID = (\S+)", resp.text)
    if not rid_match:
        raise RuntimeError("Failed to get BLAST RID from NCBI")

    rid = rid_match.group(1)
    return rid


def poll_blast(rid, max_time=MAX_POLL_TIME, interval=POLL_INTERVAL):
    """Poll BLAST for results. Returns XML string when ready."""
    start = time.time()

    while time.time() - start < max_time:
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
        }
        resp = requests.get(BLAST_URL, params=params, timeout=30)
        resp.raise_for_status()

        if "Status=WAITING" in resp.text:
            time.sleep(interval)
            continue
        elif "Status=FAILED" in resp.text:
            raise RuntimeError("BLAST search failed at NCBI")
        elif "Status=UNKNOWN" in resp.text:
            raise RuntimeError("BLAST RID not found or expired")
        elif "Status=READY" in resp.text:
            # Fetch actual results
            result_params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
            }
            result_resp = requests.get(BLAST_URL, params=result_params, timeout=60)
            result_resp.raise_for_status()
            return result_resp.text

        time.sleep(interval)

    raise TimeoutError(
        f"BLAST timed out after {max_time}s. RID: {rid} — "
        f"retrieve manually at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}"
    )


def parse_blast_xml(xml_text):
    """Parse BLAST XML output and return structured hits."""
    root = ET.fromstring(xml_text)

    hits = []
    iterations = root.findall(".//Iteration")

    for iteration in iterations:
        for hit in iteration.findall(".//Hit"):
            hit_def = hit.findtext("Hit_def", "")
            hit_accession = hit.findtext("Hit_accession", "")

            for hsp in hit.findall(".//Hsp"):
                identity = int(hsp.findtext("Hsp_identity", "0"))
                align_len = int(hsp.findtext("Hsp_align-len", "1"))
                mismatches = align_len - identity
                gaps = int(hsp.findtext("Hsp_gaps", "0"))
                evalue = float(hsp.findtext("Hsp_evalue", "999"))
                bit_score = float(hsp.findtext("Hsp_bit-score", "0"))
                pct_identity = round(identity / align_len * 100, 1) if align_len > 0 else 0

                # Extract organism from Hit_def (usually in brackets)
                org_match = re.search(r"\[(.+?)\]", hit_def)
                organism = org_match.group(1) if org_match else ""

                hits.append({
                    "description": hit_def,
                    "accession": hit_accession,
                    "organism": organism,
                    "pct_identity": pct_identity,
                    "alignment_length": align_len,
                    "mismatches": mismatches,
                    "gaps": gaps,
                    "evalue": evalue,
                    "bit_score": round(bit_score, 1),
                    "ncbi_url": f"https://www.ncbi.nlm.nih.gov/nucleotide/{hit_accession}",
                })

    # Sort by bit score descending
    hits.sort(key=lambda h: h["bit_score"], reverse=True)
    return hits


def run_blast(seq):
    """Full BLAST pipeline: submit, poll, parse."""
    rid = submit_blast(seq)
    xml_text = poll_blast(rid)
    hits = parse_blast_xml(xml_text)
    return {
        "rid": rid,
        "hits": hits,
        "num_hits": len(hits),
    }
