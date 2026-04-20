import type { AccessionMeta, AccessionProfile, MotifDefinition, PromoterWindow, StructuralVariant } from '../types';

interface PromoterTracksProps {
  rows: Array<{ accession: AccessionMeta; profile: AccessionProfile }>;
  motifs: MotifDefinition[];
  promoterWindow: PromoterWindow;
  highlightMotif: string;
}

function scalePosition(position: number, promoterWindow: PromoterWindow): string {
  const width = promoterWindow.end - promoterWindow.start;
  return `${((position - promoterWindow.start) / width) * 100}%`;
}

function variantClass(variant: StructuralVariant): string {
  return `sv-chip sv-${variant.type}`;
}

export function PromoterTracks({ rows, motifs, promoterWindow, highlightMotif }: PromoterTracksProps) {
  const highlighted = motifs.find((motif) => motif.id === highlightMotif);

  return (
    <section className="panel track-panel">
      <div className="panel-header">
        <div>
          <p className="eyebrow">Promoter tracks</p>
          <h2>FT promoter architecture</h2>
        </div>
        <p className="panel-copy">
          Window {promoterWindow.start} to +{promoterWindow.end} bp around the TSS.
          {highlighted ? ` Highlighting ${highlighted.label}.` : ''}
        </p>
      </div>

      <div className="axis-row">
        <span>{promoterWindow.start}</span>
        <span>-1500</span>
        <span>-1000</span>
        <span>-500</span>
        <span>TSS</span>
        <span>+200</span>
      </div>

      <div className="track-stack">
        {rows.map(({ accession, profile }) => (
          <article className="track-row" key={accession.accessionName}>
            <div className="track-meta">
              <strong>{accession.accessionName}</strong>
              <span>{accession.country}</span>
              <small>{profile.sequenceIdentity}% identity · {profile.variantCount} local variants</small>
            </div>
            <div className="track-visual">
              <div className="track-line" />
              <div className="track-tss" style={{ left: scalePosition(promoterWindow.tss, promoterWindow) }} />
              {profile.structuralVariants.map((variant) => (
                <span
                  className={variantClass(variant)}
                  key={variant.id}
                  style={{
                    left: scalePosition(variant.start, promoterWindow),
                    width: `calc(${scalePosition(variant.end, promoterWindow)} - ${scalePosition(variant.start, promoterWindow)})`,
                  }}
                  title={`${variant.type} (${variant.impact}) ${variant.start}..${variant.end}`}
                />
              ))}
              {profile.motifInstances.map((instance, index) => {
                const motif = motifs.find((entry) => entry.id === instance.motifId);
                if (!motif) return null;
                const className = instance.motifId === highlightMotif ? 'motif-chip is-highlighted' : 'motif-chip';
                return (
                  <span
                    className={className}
                    key={`${accession.accessionName}-${instance.motifId}-${instance.start}-${instance.strand}-${index}`}
                    style={{
                      left: scalePosition(instance.start, promoterWindow),
                      width: `${(instance.length / (promoterWindow.end - promoterWindow.start)) * 100}%`,
                      background: motif.color,
                    }}
                    title={`${motif.label} ${instance.start}..${instance.start + instance.length} (${instance.score.toFixed(2)})`}
                  >
                    {motif.label}
                  </span>
                );
              })}
            </div>
          </article>
        ))}
      </div>
    </section>
  );
}
