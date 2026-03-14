# BioTools

A self-hosted web application that combines multiple bioinformatics tools into a single service. Built with Flask and served by Gunicorn in Docker.

## Tools

### SeqConvert
Biological sequence format converter powered by [Biopython](https://biopython.org/). Converts between ~35 bioinformatics file formats including FASTA, GenBank, FASTQ, EMBL, Clustal, PHYLIP, and more. Supports drag-and-drop file upload and provides both a web UI and a REST API.

### Venn Diagram
Interactive Venn diagram tool for comparing 2–4 ID lists. Paste or upload your lists, and get an interactive SVG diagram with clickable regions. Supports exporting to SVG, PNG, PDF, and CSV.

### PhyloTree
Phylogenetic tree builder with interactive visualization. Upload or paste FASTA sequences, choose an alignment method (MUSCLE, ClustalOmega, or built-in pairwise) and tree-building algorithm (Neighbor-Joining, UPGMA, or Maximum Likelihood via IQ-TREE). Results are displayed as an interactive D3.js phylogram with branch length support. Export trees in Newick, Nexus, PhyloXML, or NeXML format.

### PCoA
Principal Coordinates Analysis for microbial ecology and community analysis. Upload an abundance table (CSV/TSV, features × samples) with optional metadata. Choose from four distance metrics (Bray-Curtis, Jaccard, Euclidean, Correlation) and three normalisation methods (None, TSS, Hellinger). Results are displayed as interactive 2D/3D Plotly.js scatter plots with metadata-based colouring and shaping. Download coordinates and distance matrices as CSV.

## Quick Start

### Using Docker Hub (recommended)

```bash
docker run -d -p 5590:5590 eriklysoe/biotools:v1.0.1
```

### Building from source

```bash
docker compose up -d --build
```

The app runs at **http://localhost:5590**.

## Running Without Docker

```bash
pip install -r requirements.txt
# PhyloTree also needs: apt-get install muscle clustalo iqtree
mkdir -p tmp
TMP_DIR=./tmp gunicorn -b 0.0.0.0:5590 -w 2 --timeout 120 app.main:app
```

## Configuration

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
| `PORT` | `5590` | Host port (docker-compose only) |
| `MAX_UPLOAD_MB` | `50` | Maximum upload file size in MB |
| `TMP_DIR` | `/app/tmp` | Temporary file directory |
| `BASE_URL` | `http://localhost:5590` | Base URL shown in API usage hints |
| `ADMIN_USER` | _(empty)_ | HTTP Basic Auth username (optional) |
| `ADMIN_PASS` | _(empty)_ | HTTP Basic Auth password (optional) |

## API

### SeqConvert

```bash
# Convert a FASTA file to GenBank
curl -F "file=@input.fasta" \
     -F "input_format=fasta" \
     -F "output_format=genbank" \
     http://localhost:5590/api/convert \
     -o output.gb

# List available formats
curl http://localhost:5590/api/formats
```

### Venn Diagram

```bash
# Compare two sets
curl -X POST http://localhost:5590/api/compare \
     -H "Content-Type: application/json" \
     -d '{"sets": {"Set A": ["gene1","gene2","gene3"], "Set B": ["gene2","gene3","gene4"]}}'
```

### PhyloTree

```bash
# Build a tree from a FASTA file
curl -F "sequence_file=@sequences.fasta" \
     -F "alignment_method=muscle" \
     -F "tree_method=nj" \
     http://localhost:5590/api/phylotree

# Check available methods
curl http://localhost:5590/api/phylotree/methods
```

### PCoA

```bash
# Run PCoA on an abundance table
curl -F "abundance_file=@abundance.csv" \
     -F "distance_metric=braycurtis" \
     -F "normalisation=none" \
     http://localhost:5590/api/pcoa

# With optional metadata
curl -F "abundance_file=@abundance.csv" \
     -F "metadata_file=@metadata.csv" \
     -F "distance_metric=braycurtis" \
     -F "normalisation=tss" \
     http://localhost:5590/api/pcoa
```

### Health Check

```
GET /health
```

## License

MIT
