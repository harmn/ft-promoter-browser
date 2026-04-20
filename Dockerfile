# ═══════════════════════════════════════════════════════════════════════════════
#  Four-stage build
#
#  Stage 1  rust-build   compile the data-pipeline binary
#  Stage 2  data-gen     run the binary to fetch + process all 27 genomes
#                        (~10 min, requires internet, result is Docker-cached)
#  Stage 3  node-build   npm build the React/Vite frontend
#  Stage 4  serve        nginx serves the static bundle
#
#  Build:   docker build -t ft-promoter-browser .
#  Run:     docker run -d -p 8080:80 ft-promoter-browser
#
#  If you have already generated public/data/ft-promoter.json locally you can
#  skip stages 1-2 with:
#    docker build --target serve-prebuilt -t ft-promoter-browser .
#  (see the alias target at the bottom of this file)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Stage 1: compile Rust data-pipeline ────────────────────────────────────
FROM rust:1.82-alpine AS rust-build

RUN apk add --no-cache musl-dev

# Mirror the on-disk layout so that CARGO_MANIFEST_DIR resolves correctly.
# env!("CARGO_MANIFEST_DIR") is baked in at compile time as /app/rust/ft-promoter-builder
# so project_root = /app  →  output = /app/public/data/ft-promoter.json
WORKDIR /app
COPY rust/ ./rust/

RUN cd rust/ft-promoter-builder && cargo build --release

# ── Stage 2: generate dataset (needs internet, ~10 min) ─────────────────────
FROM alpine:3.20 AS data-gen

RUN apk add --no-cache ca-certificates

COPY --from=rust-build \
  /app/rust/ft-promoter-builder/target/release/ft-promoter-builder \
  /usr/local/bin/ft-promoter-builder

# Create the output directory that the binary expects to write into
RUN mkdir -p /app/public/data

RUN ft-promoter-builder
# → writes /app/public/data/ft-promoter.json

# ── Stage 3: build React/Vite frontend ─────────────────────────────────────
FROM node:20-alpine AS node-build
WORKDIR /app

COPY package.json package-lock.json ./
RUN npm ci

COPY tsconfig.json tsconfig.app.json tsconfig.node.json ./
COPY vite.config.ts index.html ./
COPY src/ ./src/
COPY public/ ./public/

# Overwrite/add the freshly-built dataset from stage 2
COPY --from=data-gen /app/public/data/ft-promoter.json ./public/data/

RUN npm run build

# ── Stage 4: serve ──────────────────────────────────────────────────────────
FROM nginx:1.27-alpine AS serve

COPY --from=node-build /app/dist /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf

EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]

# ── Alias target: skip Rust if JSON was pre-generated locally ───────────────
# Usage: docker build --target serve-prebuilt -t ft-promoter-browser .
FROM serve AS serve-prebuilt
