"""
Sequence format definitions for Bio.SeqIO conversion.
Based on Biopython's supported formats.
"""

# Formats that can be read (input)
INPUT_FORMATS = [
    "abi", "abi-trim", "ace",
    "cif-atom", "cif-seqres",
    "clustal", "embl",
    "fasta", "fasta-2line",
    "fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina",
    "gck", "genbank", "gb", "ig", "imgt",
    "nexus",
    "pdb-seqres", "pdb-atom",
    "phd", "phylip", "pir",
    "seqxml", "sff", "sff-trim", "snapgene", "stockholm",
    "swiss", "tab", "qual",
    "uniprot-xml", "xdna",
]

# Formats that can be written (output)
OUTPUT_FORMATS = [
    "clustal", "embl",
    "fasta", "fasta-2line",
    "fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina",
    "genbank", "gb", "imgt",
    "nexus", "phd", "phylip", "pir",
    "seqxml", "sff", "stockholm",
    "tab", "qual", "xdna",
]

# Molecule type annotation values (replaces deprecated Alphabet)
MOLECULE_TYPES = {
    "none": None,
    "DNA": "DNA",
    "RNA": "RNA",
    "protein": "protein",
}

# Human-readable descriptions
FORMAT_INFO = {
    "abi": "ABI Sanger capillary trace files",
    "abi-trim": "ABI with quality trimming (Mott's algorithm)",
    "ace": "ACE assembly contigs",
    "cif-atom": "mmCIF structure-based sequence",
    "cif-seqres": "mmCIF header sequence",
    "clustal": "Clustal X/W alignment",
    "embl": "EMBL flat file",
    "fasta": "FASTA format",
    "fasta-2line": "FASTA (no line wrapping, 2 lines per record)",
    "fastq": "FASTQ (Sanger, ASCII offset 33)",
    "fastq-sanger": "FASTQ Sanger (ASCII offset 33)",
    "fastq-solexa": "FASTQ Solexa",
    "fastq-illumina": "FASTQ Illumina",
    "gck": "Gene Construction Kit",
    "genbank": "GenBank / GenPept flat file",
    "gb": "GenBank (alias)",
    "ig": "IntelliGenetics / MASE",
    "imgt": "IMGT variant of EMBL",
    "nexus": "NEXUS / PAUP alignment",
    "pdb-seqres": "PDB header sequence",
    "pdb-atom": "PDB structure-based sequence",
    "phd": "PHD (PHRED output)",
    "phylip": "PHYLIP alignment (names truncated to 10 chars)",
    "pir": "PIR / NBRF format",
    "seqxml": "SeqXML",
    "sff": "Standard Flowgram Format (454)",
    "sff-trim": "SFF with trimming",
    "snapgene": "SnapGene native format",
    "stockholm": "Stockholm / PFAM alignment",
    "swiss": "Swiss-Prot / UniProt text",
    "tab": "Tab-separated (ID + sequence)",
    "qual": "QUAL (space-separated PHRED scores)",
    "uniprot-xml": "UniProt XML",
    "xdna": "DNA Strider / Serial Cloner",
}
