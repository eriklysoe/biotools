# BioTools

A self-hosted bioinformatics web app combining SeqConvert (sequence format converter) and Venn Diagram (ID list comparison) in a single service.

## Quick Start

```bash
docker compose up -d
```

Open `http://localhost:5590` in your browser.

## Docker Hub

```
docker pull eriklysoe/biotools:latest
```

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `MAX_UPLOAD_MB` | `50` | Max upload file size in MB |
| `BASE_URL` | `http://localhost:5590` | Public URL (used in API usage hints) |
| `TMP_DIR` | `/app/tmp` | Temporary file directory |
| `ADMIN_USER` | _(empty)_ | HTTP Basic Auth username (optional) |
| `ADMIN_PASS` | _(empty)_ | HTTP Basic Auth password (optional) |

When using a reverse proxy (Pangolin, Cloudflare Tunnel, nginx, Traefik, etc.), set `BASE_URL=https://biotools.yourdomain.com` so API hints show the correct public URL.
