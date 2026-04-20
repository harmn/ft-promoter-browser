import type { Dataset, AccessionMeta, AccessionProfile } from '../types';

interface SummaryStripProps {
  dataset: Dataset;
  displayedRows: Array<{ accession: AccessionMeta; profile: AccessionProfile }>;
  selectedCount: number;
  highlightMotifLabel: string;
}

export function SummaryStrip({ dataset, displayedRows, selectedCount, highlightMotifLabel }: SummaryStripProps) {
  const identities = displayedRows.map(({ profile }) => profile.sequenceIdentity);
  const variants = displayedRows.map(({ profile }) => profile.variantCount);
  const meanIdentity = identities.length
    ? (identities.reduce((total, value) => total + value, 0) / identities.length).toFixed(1)
    : '0.0';
  const meanVariants = variants.length
    ? Math.round(variants.reduce((total, value) => total + value, 0) / variants.length)
    : 0;
  const representedCountries = new Set(displayedRows.map(({ accession }) => accession.country)).size;

  const cards = [
    { label: 'Selected accessions', value: String(selectedCount), note: 'Pinned for promoter tracks' },
    { label: 'Mean promoter identity', value: `${meanIdentity}%`, note: `Across ${displayedRows.length || dataset.accessions.length} visible genomes` },
    { label: 'Mean local variant count', value: String(meanVariants), note: 'Predicted differences in the FT window' },
    { label: 'Highlighted TFBS', value: highlightMotifLabel, note: `${representedCountries} countries represented` },
  ];

  return (
    <section className="summary-strip">
      {cards.map((card) => (
        <article className="summary-card" key={card.label}>
          <p>{card.label}</p>
          <strong>{card.value}</strong>
          <span>{card.note}</span>
        </article>
      ))}
    </section>
  );
}
