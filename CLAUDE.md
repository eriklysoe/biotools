# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioTools is a self-hosted web app combining multiple bioinformatics tools behind a single Flask service served by Gunicorn in Docker. Currently includes two tools:

- **SeqConvert** — Biological sequence format converter wrapping Biopython's `Bio.SeqIO` (~35 formats: FASTA, GenBank, FASTQ, EMBL, etc.)
- **VennApp** — Interactive Venn diagram tool for comparing 2–4 ID lists with SVG rendering and export

## Build & Run

```bash
docker compose up -d --build
# App runs at http://localhost:5590
# Health check: GET /health
```

### Running Locally Without Docker

```bash
pip install -r requirements.txt
mkdir -p tmp
TMP_DIR=./tmp gunicorn -b 0.0.0.0:5590 -w 2 --timeout 120 app.main:app
```

There are no tests, linter, or CI configured.

## Architecture

Single Flask app in `app/main.py` serving both tools. All UI is in self-contained Jinja2 templates with inline CSS/JS (no build step, no static assets directory).

- **`app/main.py`** — All routes, SeqConvert conversion logic, VennApp set comparison logic, auth decorator.
- **`app/formats.py`** — Static data: `INPUT_FORMATS`, `OUTPUT_FORMATS`, `MOLECULE_TYPES`, `FORMAT_INFO` dicts/lists consumed by routes and templates.
- **`templates/seqconvert.html`** — SeqConvert single-page UI with drag-and-drop upload, format selectors, and fetch-based conversion.
- **`templates/venn.html`** — VennApp single-page UI with set input, client-side SVG Venn diagram rendering (2–4 sets), and export (SVG/PNG/PDF/CSV).

Both templates include a shared navigation bar (`<nav class="topnav">`) for switching between tools. The nav bar CSS is duplicated in each template (no shared base template).

### Routes

| Route | Method | Tool | Purpose |
|-------|--------|------|---------|
| `/`, `/seqconvert` | GET | SeqConvert | Web UI |
| `/api/convert` | POST | SeqConvert | Multipart file conversion |
| `/api/formats` | GET | SeqConvert | Format listing JSON |
| `/venn` | GET | VennApp | Web UI |
| `/api/compare` | POST | VennApp | JSON set comparison |
| `/health` | GET | — | Health check |

### Conversion Flow (SeqConvert)

Upload saved to `TMP_DIR` → `Bio.SeqIO.convert()` (or `SeqIO.parse` + `SeqIO.write` when molecule type annotation needed) → result streamed back as download → temp files cleaned up via `@response.call_on_close`.

### Venn Comparison

Set operations computed both client-side (JS, for rendering) and server-side (`/api/compare`, bitmask-based intersection of 2–4 sets). The UI primarily uses client-side computation; the API exists for programmatic access.

## Key Configuration

| Env Variable | Default | Purpose |
|-------------|---------|---------|
| `MAX_UPLOAD_MB` | `50` | Max upload size in MB |
| `TMP_DIR` | `/app/tmp` | Temp file directory |
| `BASE_URL` | `http://localhost:5590` | Used in API usage hints |
| `ADMIN_USER` | _(empty)_ | HTTP Basic Auth username (optional) |
| `ADMIN_PASS` | _(empty)_ | HTTP Basic Auth password (optional) |

## Dependencies

Python 3.12, Flask >= 3.0, Gunicorn >= 22.0, Biopython >= 1.84, NumPy >= 1.26.
