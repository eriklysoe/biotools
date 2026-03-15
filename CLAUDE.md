# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioTools is a self-hosted web app combining multiple bioinformatics tools behind a single Flask service served by Gunicorn in Docker. Currently includes six tools:

- **SeqConvert** — Biological sequence format converter wrapping Biopython's `Bio.SeqIO` (~35 formats: FASTA, GenBank, FASTQ, EMBL, etc.)
- **VennApp** — Interactive Venn diagram tool for comparing 2–4 ID lists with SVG rendering and export
- **PhyloTree** — Phylogenetic tree builder: FASTA input → MSA (MUSCLE/ClustalOmega/built-in) → tree (NJ/UPGMA/ML via IQ-TREE) → interactive D3 visualization
- **PCoA** — Principal Coordinates Analysis: abundance table → distance matrix (Bray-Curtis/Jaccard/Euclidean/Correlation) → PCoA ordination → interactive Plotly.js scatter plots (2D/3D) with optional metadata colouring
- **PrimerCheck** — Oligonucleotide analysis toolkit: oligo properties (Tm/MW/GC%/hairpin/dimers), in-silico PCR, Primer3 design, NCBI BLAST
- **Diversity** — Alpha diversity (10 metrics, rarefaction, stats) and beta diversity (7 distance metrics, 4 ordination methods, PERMANOVA/PERMDISP/ANOSIM)

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

Single Flask app in `app/main.py` with SeqConvert and VennApp routes inline. PhyloTree, PCoA, PrimerCheck, and Diversity are Flask Blueprints registered from `app/phylotree/`, `app/pcoa/`, `app/primercheck/`, and `app/diversity/`. All UI is in self-contained Jinja2 templates with inline CSS/JS (no build step). The Flask app configures a `static/` folder but it is currently unused.

- **`app/main.py`** — SeqConvert/VennApp routes, conversion logic, auth decorator. Registers all blueprints.
- **`app/formats.py`** — Static data: `INPUT_FORMATS`, `OUTPUT_FORMATS`, `MOLECULE_TYPES`, `FORMAT_INFO`.
- **`app/phylotree/`** — PhyloTree Blueprint package:
  - `__init__.py` — Blueprint definition + `before_request` auth check
  - `routes.py` — `/phylotree` page and `/api/phylotree` endpoint
  - `utils.py` — Alignment (MUSCLE/ClustalOmega/built-in pairwise) and tree building (NJ/UPGMA via Biopython, ML via IQ-TREE subprocess)
- **`app/pcoa/`** — PCoA Blueprint package:
  - `__init__.py` — Blueprint definition + `before_request` auth check
  - `routes.py` — `/pcoa` page and `/api/pcoa` endpoint
  - `utils.py` — Abundance table parsing, normalisation (TSS/Hellinger), distance matrix (scikit-bio/scipy), PCoA ordination
- **`app/primercheck/`** — PrimerCheck Blueprint package:
  - `__init__.py` — Blueprint definition + `before_request` auth check
  - `routes.py` — `/primercheck` page and 8 API endpoints (analyze, hairpin, selfdimer, heterodimer, tmmismatch, check, design, blast)
  - `utils_shared.py` — IUPAC validation, sequence parsing, complement, GC content
  - `utils_analyze.py` — Tm calculations (Wallace, basic NN, salt-corrected, Mg2+), MW, extinction coefficient
  - `utils_structure.py` — Hairpin/dimer prediction via primer3-py
  - `utils_mismatch.py` — Tm mismatch heatmap generation
  - `utils_check.py` — In-silico PCR engine
  - `utils_design.py` — Primer3 wrapper
  - `utils_blast.py` — NCBI BLAST submit/poll/parse
- **`app/diversity/`** — Diversity Blueprint package:
  - `__init__.py` — Blueprint definition + `before_request` auth check
  - `routes.py` — `/diversity` page and `/api/diversity/alpha`, `/api/diversity/beta`, `/api/diversity/mantel` endpoints
  - `utils_shared.py` — Abundance table/metadata parsing, rarefaction, tree validation
  - `utils_alpha_metrics.py` — 10 alpha diversity metrics
  - `utils_alpha_stats.py` — Statistical testing (Mann-Whitney, Kruskal-Wallis, Dunn's post-hoc)
  - `utils_alpha_rarefaction.py` — Rarefaction curve calculation
  - `utils_beta_distances.py` — 7 distance metrics (Bray-Curtis, Jaccard, Canberra, UniFrac, Euclidean, Correlation, Aitchison)
  - `utils_beta_ordination.py` — PCoA, NMDS, t-SNE, PCA ordination
  - `utils_beta_stats.py` — PERMANOVA, PERMDISP, ANOSIM, pairwise PERMANOVA, Mantel test
- **`templates/seqconvert.html`** — SeqConvert UI with drag-and-drop upload.
- **`templates/venn.html`** — VennApp UI with client-side SVG Venn rendering.
- **`templates/phylotree.html`** — PhyloTree UI with D3.js tree visualization, file upload + paste input, method selectors.
- **`templates/pcoa.html`** — PCoA UI with Plotly.js interactive scatter plots (2D/3D), metadata colouring/shaping, distance matrix display.
- **`templates/primercheck.html`** — PrimerCheck UI with 4 main tabs (Analyze/Check/Design/BLAST) and 5 sub-tabs in Analyze.
- **`templates/diversity.html`** — Diversity UI with Alpha/Beta tabs, Plotly.js visualizations, statistical output.

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
| `/pcoa` | GET | PCoA | Web UI |
| `/api/pcoa` | POST | PCoA | Abundance table → PCoA results (multipart) |
| `/primercheck` | GET | PrimerCheck | Web UI |
| `/api/primercheck/analyze` | POST | PrimerCheck | Oligo property analysis |
| `/api/primercheck/hairpin` | POST | PrimerCheck | Hairpin prediction |
| `/api/primercheck/selfdimer` | POST | PrimerCheck | Self-dimer prediction |
| `/api/primercheck/heterodimer` | POST | PrimerCheck | Hetero-dimer prediction |
| `/api/primercheck/tmmismatch` | POST | PrimerCheck | Tm mismatch heatmap |
| `/api/primercheck/check` | POST | PrimerCheck | In-silico PCR |
| `/api/primercheck/design` | POST | PrimerCheck | Primer3 primer design |
| `/api/primercheck/blast` | POST | PrimerCheck | NCBI BLAST search |
| `/diversity` | GET | Diversity | Web UI |
| `/api/diversity/alpha` | POST | Diversity | Alpha diversity metrics |
| `/api/diversity/beta` | POST | Diversity | Beta diversity + ordination |
| `/api/diversity/mantel` | POST | Diversity | Mantel test between distance matrices |
| `/health` | GET | — | Health check |

### Conversion Flow (SeqConvert)

Upload saved to `TMP_DIR` → `Bio.SeqIO.convert()` (or `SeqIO.parse` + `SeqIO.write` when molecule type annotation needed) → result streamed back as download → temp files cleaned up via `@response.call_on_close`.

### Venn Comparison

Set operations computed both client-side (JS, for rendering) and server-side (`/api/compare`, bitmask-based intersection of 2–4 sets). The UI primarily uses client-side computation; the API exists for programmatic access.

### PhyloTree Pipeline

FASTA input → alignment (MUSCLE/ClustalOmega subprocess or built-in pairwise distance) → tree building (NJ/UPGMA via `Bio.Phylo.TreeConstruction`, or ML via `iqtree2` subprocess with `-m MFP -B 1000`) → Newick string + metadata returned as JSON. Each job runs in a unique temp directory under `TMP_DIR`, cleaned up after response. External binaries (muscle, clustalo, iqtree) are gracefully detected; unavailable methods are disabled in the UI.

### PCoA Pipeline

Abundance table (CSV/TSV, features × samples) → optional metadata file → normalisation (none/TSS/Hellinger) → distance matrix (Bray-Curtis/Jaccard via `skbio.diversity.beta_diversity`, Euclidean/Correlation via `scipy.spatial.distance`) → PCoA ordination via `skbio.stats.ordination.pcoa()` → JSON with PC1-3 coordinates, variance explained, distance matrix, and optional metadata. The UI renders interactive 2D/3D scatter plots with Plotly.js.

### PrimerCheck Pipeline

Oligo sequence → validation (IUPAC) → property calculation (Tm via multiple methods, MW, extinction coefficient, GC%) → structure prediction (hairpin/dimer via primer3-py) → results as JSON. In-silico PCR: template + primers → binding site search (with mismatch tolerance) → amplicon prediction. Primer design: template → Primer3 wrapper → candidate primer pairs. BLAST: sequence → NCBI BLAST API submit → poll → parse XML results.

### Diversity Pipeline

**Alpha:** Abundance table → optional rarefaction → 10 diversity metrics per sample → optional statistical testing with metadata grouping (Mann-Whitney U for 2 groups, Kruskal-Wallis + Dunn's post-hoc for 3+) → rarefaction curves (subsampling at increasing depths). **Beta:** Abundance table → optional rarefaction → distance matrix (7 metrics including UniFrac with optional tree) → ordination (PCoA via scikit-bio, NMDS/t-SNE/PCA via scikit-learn) → statistical tests (PERMANOVA, PERMDISP, ANOSIM, pairwise PERMANOVA with BH correction, Mantel test).

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

Python 3.12, Flask >= 3.0, Gunicorn >= 22.0, Biopython >= 1.84, NumPy >= 1.26, pandas >= 2.0, scipy >= 1.11, scikit-bio >= 0.6, primer3-py >= 2.0, requests >= 2.28, scikit-learn >= 1.3, scikit-posthocs >= 0.9. System binaries (Docker only): muscle, clustalo, iqtree.
