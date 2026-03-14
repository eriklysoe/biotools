# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioTools is a self-hosted web app combining multiple bioinformatics tools behind a single Flask service served by Gunicorn in Docker. Currently includes three tools:

- **SeqConvert** — Biological sequence format converter wrapping Biopython's `Bio.SeqIO` (~35 formats: FASTA, GenBank, FASTQ, EMBL, etc.)
- **VennApp** — Interactive Venn diagram tool for comparing 2–4 ID lists with SVG rendering and export
- **PhyloTree** — Phylogenetic tree builder: FASTA input → MSA (MUSCLE/ClustalOmega/built-in) → tree (NJ/UPGMA/ML via IQ-TREE) → interactive D3 visualization

## Build & Run

```bash
docker compose up -d --build
# App runs at http://localhost:5590
# Health check: GET /health
```

### Running Locally Without Docker

```bash
pip install -r requirements.txt
# PhyloTree also needs: apt-get install muscle clustalo iqtree
mkdir -p tmp
TMP_DIR=./tmp gunicorn -b 0.0.0.0:5590 -w 2 --timeout 120 app.main:app
```

There are no tests, linter, or CI configured.

## Architecture

Single Flask app in `app/main.py` with SeqConvert and VennApp routes inline. PhyloTree is a Flask Blueprint registered from `app/phylotree/`. All UI is in self-contained Jinja2 templates with inline CSS/JS (no build step). The Flask app configures a `static/` folder but it is currently unused.

- **`app/main.py`** — SeqConvert/VennApp routes, conversion logic, auth decorator. Registers the PhyloTree blueprint.
- **`app/formats.py`** — Static data: `INPUT_FORMATS`, `OUTPUT_FORMATS`, `MOLECULE_TYPES`, `FORMAT_INFO`.
- **`app/phylotree/`** — PhyloTree Blueprint package:
  - `__init__.py` — Blueprint definition + `before_request` auth check
  - `routes.py` — `/phylotree` page and `/api/phylotree` endpoint
  - `utils.py` — Alignment (MUSCLE/ClustalOmega/built-in pairwise) and tree building (NJ/UPGMA via Biopython, ML via IQ-TREE subprocess)
- **`templates/seqconvert.html`** — SeqConvert UI with drag-and-drop upload.
- **`templates/venn.html`** — VennApp UI with client-side SVG Venn rendering.
- **`templates/phylotree.html`** — PhyloTree UI with D3.js tree visualization, file upload + paste input, method selectors.

All templates include a shared navigation bar (`<nav class="topnav">`). The nav bar CSS is duplicated in each template (no shared base template).

### Routes

| Route | Method | Tool | Purpose |
|-------|--------|------|---------|
| `/`, `/seqconvert` | GET | SeqConvert | Web UI |
| `/api/convert` | POST | SeqConvert | Multipart file conversion |
| `/api/formats` | GET | SeqConvert | Format listing JSON |
| `/venn` | GET | VennApp | Web UI |
| `/api/compare` | POST | VennApp | JSON set comparison |
| `/phylotree` | GET | PhyloTree | Web UI |
| `/api/phylotree` | POST | PhyloTree | FASTA → Newick tree (multipart or text) |
| `/api/phylotree/methods` | GET | PhyloTree | Available alignment/tree methods |
| `/health` | GET | — | Health check |

### Conversion Flow (SeqConvert)

Upload saved to `TMP_DIR` → `Bio.SeqIO.convert()` (or `SeqIO.parse` + `SeqIO.write` when molecule type annotation needed) → result streamed back as download → temp files cleaned up via `@response.call_on_close`.

### Venn Comparison

Set operations computed both client-side (JS, for rendering) and server-side (`/api/compare`, bitmask-based intersection of 2–4 sets). The UI primarily uses client-side computation; the API exists for programmatic access.

### PhyloTree Pipeline

FASTA input → alignment (MUSCLE/ClustalOmega subprocess or built-in pairwise distance) → tree building (NJ/UPGMA via `Bio.Phylo.TreeConstruction`, or ML via `iqtree2` subprocess with `-m MFP -B 1000`) → Newick string + metadata returned as JSON. Each job runs in a unique temp directory under `TMP_DIR`, cleaned up after response. External binaries (muscle, clustalo, iqtree) are gracefully detected; unavailable methods are disabled in the UI.

## Key Configuration

| Env Variable | Default | Purpose |
|-------------|---------|---------|
| `PORT` | `5590` | Host port (docker-compose only) |
| `MAX_UPLOAD_MB` | `50` | Max upload size in MB |
| `TMP_DIR` | `/app/tmp` | Temp file directory (mounted to `data/tmp/` on host) |
| `BASE_URL` | `http://localhost:5590` | Used in API usage hints |
| `ADMIN_USER` | _(empty)_ | HTTP Basic Auth username (optional) |
| `ADMIN_PASS` | _(empty)_ | HTTP Basic Auth password (optional) |

## Dependencies

Python 3.12, Flask >= 3.0, Gunicorn >= 22.0, Biopython >= 1.84, NumPy >= 1.26. System binaries (Docker only): muscle, clustalo, iqtree.
