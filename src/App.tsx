import { useEffect, useState } from 'react';
import { SequenceAlignment } from './components/SequenceAlignment';
import { seedDataset } from './data/seedDataset';
import type { AccessionProfile, Dataset } from './types';

export default function App() {
  const [dataset, setDataset] = useState<Dataset>(seedDataset);
  const [datasetState, setDatasetState] = useState('Using bundled starter dataset.');
  const [highlightMotif, setHighlightMotif] = useState(seedDataset.motifs[0].id);

  useEffect(() => {
    let active = true;
    (async () => {
      try {
        const res = await fetch(`${import.meta.env.BASE_URL}data/ft-promoter.json`);
        if (!res.ok) throw new Error('not found');
        const ext = (await res.json()) as Dataset;
        if (!active) return;
        setDataset(ext);
        setDatasetState('Loaded /public/data/ft-promoter.json');
        if (ext.motifs[0]) setHighlightMotif(ext.motifs[0].id);
      } catch {
        if (active) setDatasetState('Using bundled starter dataset. Add /public/data/ft-promoter.json to replace it.');
      }
    })();
    return () => { active = false; };
  }, []);

  const profiles = Object.fromEntries(dataset.profiles.map((p) => [p.accessionName, p]));

  // All accessions — Col-0 first as reference row
  const allRows = [
    dataset.accessions.find((a) => a.accessionName === 'Col-0'),
    ...dataset.accessions.filter((a) => a.accessionName !== 'Col-0'),
  ]
    .filter((a): a is NonNullable<typeof a> => a != null)
    .flatMap((a) => {
      const p = profiles[a.accessionName];
      return p ? [{ accession: a, profile: p as AccessionProfile }] : [];
    });

  return (
    <div className="app-shell">

      {/* ── Alignment (full view) ─────────────────────────────────────── */}
      <SequenceAlignment
        rows={allRows}
        motifs={dataset.motifs}
        promoterWindow={dataset.promoterWindow}
        highlightMotif={highlightMotif}
        onHighlightMotif={setHighlightMotif}
      />

      {/* ── Hero ──────────────────────────────────────────────────────── */}
      <header className="hero">
        <div className="hero-copy">
          <p className="eyebrow">Arabidopsis promoter browser</p>
          <h1>FT promoter TFBS browser</h1>
          <p className="hero-text">
            FT (FLOWERING LOCUS T) is the central floral integrator in Arabidopsis. This browser
            shows the FT promoter ({dataset.promoterWindow.start} to +{dataset.promoterWindow.end} bp) aligned across {dataset.accessions.length} accessions
            from the 1001 Genomes Plus Phase 1 dataset, with binding site predictions for key
            flowering-time transcription factors: MADS-box repressors (FLC, SVP, AGL15), the
            CO–NF-Y photoperiod complex (NFYC2), G-box bHLH factors (PIF4, MYC3), circadian
            regulator CCA1, DOF factor CDF2, TCP20, AP2/ERF repressor SMZ, and RAV/B3 repressor
            TEM1. Scores are derived from JASPAR PWM scans; accessions are clustered by promoter
            sequence similarity (UPGMA).
          </p>
        </div>
        <div className="hero-panel">
          <span className="status-pill">{datasetState}</span>
          <dl>
            <div><dt>Gene target</dt><dd>{dataset.gene}</dd></div>
            <div><dt>Accessions</dt><dd>{dataset.accessions.length}</dd></div>
            <div><dt>Promoter window</dt><dd>{dataset.promoterWindow.start} to +{dataset.promoterWindow.end} bp</dd></div>
            <div><dt>Data mode</dt><dd>{dataset.isDemo ? 'Demo + real metadata' : 'External dataset'}</dd></div>
          </dl>
        </div>
      </header>

    </div>
  );
}
