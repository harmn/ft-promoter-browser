# FT Promoter Browser

A Vite + React browser for visually exploring the 27 full-length Arabidopsis thaliana genomes from 1001 Genomes Plus Phase 1, with a specific focus on TFBS variation in the FT promoter.

## What this includes

- A promoter-centric interface inspired by the inspection workflow of Acacia.
- Real accession metadata for the 27 phase-1 genomes.
- A motif score heatmap across all accessions.
- Per-accession promoter tracks for the FT upstream window.
- A focused comparison panel against Col-0.
- A real sequence-derived FT promoter dataset generated from the public phase-1 assemblies and annotations.

## Included data

The app now includes a generated `public/data/ft-promoter.json` built from the public 1001 Genomes Plus Phase 1 assemblies and gene annotations. The default generator is a Rust CLI for faster sequence processing. For each accession it:

- locates `AT1G65480` in the annotation GFF,
- extracts the FT promoter window from `-2000` to `+200` around the TSS,
- compares that window against the Col-0 window, and
- scans it with JASPAR PWM motifs for sequence-derived TFBS matches.

These motif scores are real sequence-derived match scores, not experimental ChIP occupancy or calibrated binding probabilities.

## Run locally

```bash
npm install
npm run build:data
npm run dev
```

## Rebuild the FT promoter dataset

Run `npm run build:data` to regenerate `public/data/ft-promoter.json` from the public release using the Rust builder.

The Rust builder is in `rust/ft-promoter-builder`.
On Windows, the launcher auto-detects Visual Studio Build Tools and loads the MSVC environment before invoking Cargo.
The JavaScript reference implementation remains available via `npm run build:data:js`.

## Notes on interpretation

1. The included motif library uses JASPAR PWMs for NFYC2, AGL15, PIF4, TCP20, CDF2, and CCA1.
2. The score shown in the heatmap is a normalized promoter-level PWM score derived from the strongest windows in the extracted promoter sequence.
3. The promoter comparison uses fixed TSS-relative windows, so variant counts are position-wise differences against Col-0 within that window.

## Source context

- Reference UI inspiration: https://github.com/wur-bioinformatics/acacia
- Arabidopsis phase-1 project: https://1001genomes.org/1001Gp/phase1/introduction
- Phase-1 accessions: https://1001genomes.org/1001Gp/phase1/accessions
