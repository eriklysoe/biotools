# BioTools

A self-hosted web application that combines multiple bioinformatics tools into a single service. Built with Flask and served by Gunicorn in Docker.

## Tools

### SeqConvert
Biological sequence format converter powered by [Biopython](https://biopython.org/). Converts between ~35 bioinformatics file formats including FASTA, GenBank, FASTQ, EMBL, Clustal, PHYLIP, and more. Supports drag-and-drop file upload and provides both a web UI and a REST API.

### Venn Diagram
Interactive Venn diagram tool for comparing 2–4 ID lists. Paste or upload your lists, and get an interactive SVG diagram with clickable regions. Supports exporting to SVG, PNG, PDF, and CSV.

## Quick Start

```bash
docker compose up -d --build
```

The app runs at **http://localhost:5590**.

## Running Without Docker

```bash
pip install -r requirements.txt
mkdir -p tmp
TMP_DIR=./tmp gunicorn -b 0.0.0.0:5590 -w 2 --timeout 120 app.main:app
```

## Configuration

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
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

### Health Check

```
GET /health
```
