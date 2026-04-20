import type { AccessionMeta, AccessionProfile } from '../types';

interface AccessionTableProps {
  accessions: AccessionMeta[];
  profiles: Record<string, AccessionProfile>;
}

export function AccessionTable({ accessions, profiles }: AccessionTableProps) {
  return (
    <section className="panel table-panel">
      <div className="panel-header">
        <div>
          <p className="eyebrow">Accessions</p>
          <h2>Phase-1 metadata</h2>
        </div>
        <p className="panel-copy">{accessions.length} accessions shown.</p>
      </div>

      <div className="table-scroll">
        <table>
          <thead>
            <tr>
              <th>Name</th>
              <th>Country</th>
              <th>Lat / Lon</th>
              <th>ID</th>
              <th>Identity vs Col-0</th>
              <th>SNPs</th>
            </tr>
          </thead>
          <tbody>
            {accessions.map((accession) => {
              const profile = profiles[accession.accessionName];
              return (
                <tr key={accession.accessionName}>
                  <td>
                    <strong>{accession.accessionName}</strong>
                    <small>{accession.stockId !== '—' ? `Stock ${accession.stockId}` : ''}</small>
                  </td>
                  <td>{accession.country}</td>
                  <td style={{ fontVariantNumeric: 'tabular-nums', fontSize: '0.82rem', color: 'var(--muted)' }}>
                    {accession.latitude.toFixed(2)}, {accession.longitude.toFixed(2)}
                  </td>
                  <td>{accession.accessionNumber}</td>
                  <td>{profile?.sequenceIdentity ?? '—'}%</td>
                  <td>{profile?.variantCount ?? '—'}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </section>
  );
}
