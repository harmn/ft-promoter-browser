import { Fragment } from 'react';
import type { AccessionMeta, AccessionProfile, MotifDefinition } from '../types';

interface MotifHeatmapProps {
  accessions: AccessionMeta[];
  profiles: Record<string, AccessionProfile>;
  motifs: MotifDefinition[];
  selectedNames: string[];
  highlightMotif: string;
  onHighlightMotif: (motifId: string) => void;
}

function getCellColor(score: number): string {
  const hue = 28 + Math.round(score * 140);
  const saturation = 72;
  const lightness = 92 - Math.round(score * 42);
  return `hsl(${hue} ${saturation}% ${lightness}%)`;
}

export function MotifHeatmap({
  accessions,
  profiles,
  motifs,
  selectedNames,
  highlightMotif,
  onHighlightMotif,
}: MotifHeatmapProps) {
  return (
    <section className="panel heatmap-panel">
      <div className="panel-header">
        <div>
          <p className="eyebrow">TFBS matrix</p>
          <h2>Motif match scores across accessions</h2>
        </div>
        <p className="panel-copy">Click a motif to spotlight sequence-derived hits in the promoter track view.</p>
      </div>

      <div className="heatmap-grid" style={{ gridTemplateColumns: `minmax(150px, 1.2fr) repeat(${motifs.length}, minmax(84px, 1fr))` }}>
        <div className="heatmap-corner">Accession</div>
        {motifs.map((motif) => (
          <button
            className={motif.id === highlightMotif ? 'heatmap-heading is-active' : 'heatmap-heading'}
            key={motif.id}
            onClick={() => onHighlightMotif(motif.id)}
            type="button"
          >
            <span>{motif.label}</span>
            <small>{motif.element ?? motif.family}</small>
          </button>
        ))}

        {accessions.map((accession) => {
          const profile = profiles[accession.accessionName];
          const selected = selectedNames.includes(accession.accessionName);
          return (
            <Fragment key={accession.accessionName}>
              <div className={selected ? 'heatmap-row-label is-selected' : 'heatmap-row-label'} key={`${accession.accessionName}-label`}>
                <span>{accession.accessionName}</span>
                <small>{accession.countryCode}</small>
              </div>
              {motifs.map((motif) => {
                const score = profile?.motifScores[motif.id] ?? 0;
                const active = motif.id === highlightMotif;
                return (
                  <button
                    key={`${accession.accessionName}-${motif.id}`}
                    type="button"
                    className={active ? 'heatmap-cell is-highlighted' : 'heatmap-cell'}
                    style={{ background: getCellColor(score), borderColor: active ? motif.color : 'transparent' }}
                    title={`${accession.accessionName} · ${motif.label}: ${score.toFixed(2)}`}
                    onClick={() => onHighlightMotif(motif.id)}
                  >
                    {score.toFixed(2)}
                  </button>
                );
              })}
            </Fragment>
          );
        })}
      </div>
    </section>
  );
}
