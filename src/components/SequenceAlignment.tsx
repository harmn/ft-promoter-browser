import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { AccessionMeta, AccessionProfile, MotifDefinition, MotifInstance, PromoterWindow } from '../types';

// ── Sequence-logo (PWM) renderer ─────────────────────────────────────────────
const LOGO_BASE_COLORS: Record<string, string> = {
  A: '#4ade80', C: '#60a5fa', G: '#fbbf24', T: '#f87171',
};

function MotifLogo({ motif }: { motif: MotifDefinition }) {
  const pfm = motif.pfm;
  if (!pfm) return null;

  const POS_W = 16;
  const MAX_H = 52;
  const len = pfm.A.length;

  const positions = Array.from({ length: len }, (_, i) => {
    const raw = [pfm.A[i], pfm.C[i], pfm.G[i], pfm.T[i]];
    const total = raw.reduce((s, v) => s + v, 0);
    const freqs = raw.map((v) => (total > 0 ? v / total : 0.25));
    const H = freqs.reduce((s, f) => (f > 0 ? s - f * Math.log2(f) : s), 0);
    const ic = Math.max(0, 2 - H); // 0–2 bits
    return [
      { b: 'A', f: freqs[0], c: LOGO_BASE_COLORS.A },
      { b: 'C', f: freqs[1], c: LOGO_BASE_COLORS.C },
      { b: 'G', f: freqs[2], c: LOGO_BASE_COLORS.G },
      { b: 'T', f: freqs[3], c: LOGO_BASE_COLORS.T },
    ]
      .sort((a, b) => a.f - b.f) // smallest first, stacks upward
      .map((entry) => ({ ...entry, h: entry.f * (ic / 2) * MAX_H }));
  });

  return (
    <div className="motif-logo-wrap">
      <div className="motif-logo-label">
        <strong style={{ color: motif.color }}>{motif.label}</strong>
        <span>{motif.family}{motif.element ? ` · ${motif.element}` : ''} · {motif.source}</span>
        <em>{motif.description}</em>
      </div>
      <svg width={len * POS_W} height={MAX_H + 18} style={{ display: 'block', overflow: 'visible' }}>
        {/* IC bars */}
        {positions.map((bases, i) => {
          let stackY = MAX_H;
          return (
            <g key={i}>
              {bases.map(({ b, c, h }) => {
                stackY -= h;
                if (h <= 0.4) return null;
                const midY = stackY + h / 2;
                const fontSize = Math.min(h * 0.72, POS_W - 3);
                return (
                  <g key={b}>
                    <rect x={i * POS_W + 1} y={stackY} width={POS_W - 2} height={h} fill={c} rx={1} />
                    {fontSize >= 5 && (
                      <text
                        x={i * POS_W + POS_W / 2} y={midY + fontSize * 0.35}
                        textAnchor="middle" fontSize={fontSize}
                        fill="white" fontWeight="bold"
                        fontFamily="IBM Plex Mono, monospace"
                        style={{ pointerEvents: 'none', userSelect: 'none' }}
                      >
                        {b}
                      </text>
                    )}
                  </g>
                );
              })}
            </g>
          );
        })}
        {/* Position ticks */}
        {Array.from({ length: len }, (_, i) =>
          (i === 0 || (i + 1) % 5 === 0 || i === len - 1) ? (
            <text key={i} x={(i + 0.5) * POS_W} y={MAX_H + 13} textAnchor="middle"
              fontSize={8} fill="#9ca3af" fontFamily="IBM Plex Mono, monospace">
              {i + 1}
            </text>
          ) : null
        )}
        {/* IC axis label */}
        <text x={-2} y={MAX_H / 2} fontSize={7} fill="#9ca3af" textAnchor="middle"
          transform={`rotate(-90, -2, ${MAX_H / 2})`} fontFamily="IBM Plex Mono, monospace">
          IC (bits)
        </text>
      </svg>
    </div>
  );
}

// ── UPGMA hierarchical clustering by pairwise SNP distance ──────────────────
type DendroNode =
  | { type: 'leaf'; idx: number }
  | { type: 'inner'; left: DendroNode; right: DendroNode; dist: number };

interface ClusterResult {
  order: number[];
  tree: DendroNode;
  maxDist: number;
  idxToPos: Map<number, number>;
}

function buildCluster(seqs: string[]): ClusterResult | null {
  const n = seqs.length;
  if (n <= 1) return null;
  const d: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let diff = 0;
      const len = Math.min(seqs[i].length, seqs[j].length);
      for (let k = 0; k < len; k++) {
        const a = seqs[i][k], b = seqs[j][k];
        if (a === '-' || b === '-' || a === 'N' || b === 'N') continue;
        if (a !== b) diff++;
      }
      d[i][j] = d[j][i] = diff;
    }
  }
  const nodes: DendroNode[] = seqs.map((_, i) => ({ type: 'leaf' as const, idx: i }));
  const dist = d.map((row) => [...row]);
  const sizes = new Array<number>(n).fill(1);
  let active = Array.from({ length: n }, (_, i) => i);
  let maxDist = 0;
  while (active.length > 1) {
    let minD = Infinity, mi = -1, mj = -1;
    for (let ai = 0; ai < active.length; ai++) {
      for (let aj = ai + 1; aj < active.length; aj++) {
        const i = active[ai], j = active[aj];
        if (dist[i][j] < minD) { minD = dist[i][j]; mi = i; mj = j; }
      }
    }
    maxDist = Math.max(maxDist, minD);
    nodes[mi] = { type: 'inner', left: nodes[mi], right: nodes[mj], dist: minD };
    const ni = sizes[mi], nj = sizes[mj];
    for (const k of active) {
      if (k === mi || k === mj) continue;
      dist[mi][k] = dist[k][mi] = (dist[mi][k] * ni + dist[mj][k] * nj) / (ni + nj);
    }
    sizes[mi] = ni + nj;
    active = active.filter((k) => k !== mj);
  }
  function leaves(node: DendroNode): number[] {
    return node.type === 'leaf' ? [node.idx] : [...leaves(node.left), ...leaves(node.right)];
  }
  const tree = nodes[active[0]];
  const order = leaves(tree);
  const idxToPos = new Map(order.map((idx, pos) => [idx, pos]));
  return { order, tree, maxDist, idxToPos };
}

// Flowering time under long-day conditions (days) — populate from AraPheno study #1
const FT_LD: Record<number, number> = {
  6909: 24,   // Col-0
  10015: 43,  // Shahdara (strong FLC)
  6966: 38,   // Sq-1 (Sweden)
};

// ── Dendrogram panel ─────────────────────────────────────────────────────────
const TREE_W = 110;

function DendroPanel({ result, totalRows, rowH }: {
  result: ClusterResult; totalRows: number; rowH: number;
}) {
  const { tree, maxDist, idxToPos } = result;
  const totalH = RULER_H + GENE_H + totalRows * rowH;
  const paths: string[] = [];

  function drawNode(node: DendroNode): { midY: number; x: number } {
    if (node.type === 'leaf') {
      const pos = idxToPos.get(node.idx) ?? 0;
      const y = RULER_H + GENE_H + (pos + 1.5) * rowH; // +1 row for Col-0, after ruler+gene track
      return { midY: y, x: TREE_W };
    }
    const L = drawNode(node.left);
    const R = drawNode(node.right);
    const midY = (L.midY + R.midY) / 2;
    const nx = maxDist > 0 ? TREE_W * (1 - node.dist / maxDist) : 0;
    paths.push(`M ${nx.toFixed(1)} ${L.midY.toFixed(1)} L ${nx.toFixed(1)} ${R.midY.toFixed(1)}`);
    paths.push(`M ${nx.toFixed(1)} ${L.midY.toFixed(1)} L ${L.x.toFixed(1)} ${L.midY.toFixed(1)}`);
    paths.push(`M ${nx.toFixed(1)} ${R.midY.toFixed(1)} L ${R.x.toFixed(1)} ${R.midY.toFixed(1)}`);
    return { midY, x: nx };
  }

  drawNode(tree);

  return (
    <svg width={TREE_W} height={totalH} style={{ flexShrink: 0, display: 'block' }}>
      {paths.map((d, i) => (
        <path key={i} d={d} stroke="var(--muted)" strokeWidth={1} fill="none" opacity={0.7} />
      ))}
    </svg>
  );
}

interface Props {
  rows: Array<{ accession: AccessionMeta; profile: AccessionProfile }>;
  motifs: MotifDefinition[];
  promoterWindow: PromoterWindow;
  highlightMotif: string;
  onHighlightMotif?: (id: string) => void;
}

// IGV-style base palette
const BASE_BG: Record<string, string> = {
  A: '#4ade80', T: '#f87171', C: '#60a5fa', G: '#fbbf24',
  N: '#d1d5db', '-': '#e5e7eb',
};
const BASE_TEXT: Record<string, string> = {
  A: '#14532d', T: '#7f1d1d', C: '#1e3a8a', G: '#78350f',
  N: '#6b7280', '-': '#9ca3af',
};
const INS_BG = 'rgba(32,113,120,0.45)';  // teal — base in query at a ref-gap column

const MIN_CW = 0.3;
const MAX_CW = 30;
const DEFAULT_CW = 9;
const ROW_H = 18;
const TFBS_H = 11;
const RULER_H = 36;   // taller to fit two-line tick (position + chr coord)
const GENE_H = 20;    // gene model track height

export function SequenceAlignment({ rows, motifs, promoterWindow, highlightMotif, onHighlightMotif }: Props) {
  const refRow = rows.find((r) => r.accession.accessionName === 'Col-0') ?? rows[0];
  const useAligned = rows.some((r) => r.profile.alignedSequence);
  if (!useAligned && !rows.some((r) => r.profile.sequence)) return null;

  const [hiddenMotifs, setHiddenMotifs] = useState<Set<string>>(new Set());
  const [minScore, setMinScore] = useState(0.8);

  const toggleHide = useCallback((id: string) => {
    setHiddenMotifs((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id); else next.add(id);
      return next;
    });
  }, []);

  const getSeq = (p: AccessionProfile) => (useAligned ? p.alignedSequence : p.sequence) ?? '';

  const refSeq = getSeq(refRow.profile);
  const alignedLen = refSeq.length;

  // Map 0-based unaligned ref position → aligned column index
  const refPosToCol = useMemo(() => {
    const map: number[] = [];
    for (let col = 0; col < alignedLen; col++) {
      if (refSeq[col] !== '-') map.push(col);
    }
    return map;
  }, [refSeq, alignedLen]);

  // Clustered row order (UPGMA by SNP distance; all accessions including Col-0)
  const [clustered, setClustered] = useState(true);
  const { orderedRows, clusterResult } = useMemo(() => {
    const refName = refRow.accession.accessionName;
    const others = rows.filter((r) => r.accession.accessionName !== refName);
    if (!clustered || rows.length <= 2) {
      return { orderedRows: [refRow, ...others], clusterResult: null as ClusterResult | null };
    }
    const seqs = rows.map((r) => getSeq(r.profile));
    const result = buildCluster(seqs);
    if (!result) return { orderedRows: [refRow, ...others], clusterResult: null as ClusterResult | null };
    return {
      orderedRows: result.order.map((i: number) => rows[i]),
      clusterResult: result,
    };
  }, [rows, clustered, useAligned]); // eslint-disable-line react-hooks/exhaustive-deps

  // Per-accession: 0-based unaligned position → aligned column
  const qPosMap = useMemo(() => {
    const m = new Map<string, number[]>();
    for (const { accession, profile } of rows) {
      const seq = getSeq(profile);
      const map: number[] = [];
      for (let col = 0; col < seq.length; col++) if (seq[col] !== '-') map.push(col);
      m.set(accession.accessionName, map);
    }
    return m;
  }, [rows, useAligned]);

  const [cw, setCw] = useState(DEFAULT_CW); // char width in px
  const [scrollX, setScrollX] = useState(0);
  const cwRef = useRef(cw);
  cwRef.current = cw;
  const scrollXRef = useRef(scrollX);
  scrollXRef.current = scrollX;

  const vpRef = useRef<HTMLDivElement>(null);
  const [vpW, setVpW] = useState(900);
  const vpWRef = useRef(vpW);
  vpWRef.current = vpW;

  useEffect(() => {
    const obs = new ResizeObserver((entries) => setVpW(entries[0].contentRect.width));
    if (vpRef.current) obs.observe(vpRef.current);
    return () => obs.disconnect();
  }, []);

  // Aligned column index of the TSS (ready only after sequence data loads)
  const tssCol = useMemo(() => {
    const tssBp = promoterWindow.tss - promoterWindow.start;
    return refPosToCol[tssBp] ?? 0;
  }, [refPosToCol, promoterWindow.tss, promoterWindow.start]);

  // One-time: zoom so ~1000 bp are visible and center on TSS
  const scrollInitialized = useRef(false);
  useEffect(() => {
    if (scrollInitialized.current || vpW <= 0 || tssCol <= 0 || alignedLen <= 0) return;
    scrollInitialized.current = true;
    const promoterLen = promoterWindow.end - promoterWindow.start + 1;
    const gapFactor = alignedLen / promoterLen;
    const newCw = Math.max(MIN_CW, Math.min(MAX_CW, vpW / (1000 * gapFactor)));
    setCw(newCw);
    setScrollX(Math.max(0, tssCol * newCw - vpW / 2));
  }, [tssCol, vpW, alignedLen, promoterWindow]);

  const totalW = alignedLen * cw;

  const clamp = (x: number) => Math.max(0, Math.min(x, totalW - vpW));

  // ── Drag to pan ──────────────────────────────────────────────────────────
  const dragRef = useRef<{ x: number; s0: number } | null>(null);
  const rafRef = useRef<number | null>(null);

  const onMouseDown = useCallback((e: React.MouseEvent) => {
    if (e.button !== 0) return;
    dragRef.current = { x: e.clientX, s0: scrollXRef.current };
    if (vpRef.current) vpRef.current.dataset.drag = '1';
    e.preventDefault();
  }, []);

  useEffect(() => {
    const move = (e: MouseEvent) => {
      if (!dragRef.current) return;
      const next = clamp(dragRef.current.s0 + (dragRef.current.x - e.clientX));
      if (rafRef.current) cancelAnimationFrame(rafRef.current);
      rafRef.current = requestAnimationFrame(() => setScrollX(next));
    };
    const up = () => {
      dragRef.current = null;
      if (vpRef.current) delete vpRef.current.dataset.drag;
    };
    window.addEventListener('mousemove', move);
    window.addEventListener('mouseup', up);
    return () => { window.removeEventListener('mousemove', move); window.removeEventListener('mouseup', up); };
  }, [totalW, vpW]);

  // ── Zoom + scroll via wheel (non-passive) ────────────────────────────────
  useEffect(() => {
    const el = vpRef.current;
    if (!el) return;
    const handler = (e: WheelEvent) => {
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        const rect = el.getBoundingClientRect();
        const pivot = (scrollXRef.current + e.clientX - rect.left) / cwRef.current;
        const factor = e.deltaY > 0 ? 0.85 : 1.18;
        const newCw = Math.max(MIN_CW, Math.min(MAX_CW, cwRef.current * factor));
        setCw(newCw);
        setScrollX(Math.max(0, pivot * newCw - (e.clientX - rect.left)));
      }
    };
    el.addEventListener('wheel', handler, { passive: false });
    return () => el.removeEventListener('wheel', handler);
  }, [totalW, vpW]);

  // ── Visible column window ─────────────────────────────────────────────────
  const startCol = Math.max(0, Math.floor(scrollX / cw) - 2);
  const endCol   = Math.min(alignedLen, Math.ceil((scrollX + vpW) / cw) + 3);
  const showText = cw >= 7;

  // ── Ruler ticks (interval scales with zoom) ───────────────────────────────
  const tickInterval = cw >= 3 ? 200 : cw >= 0.8 ? 500 : 1000;
  const ticks: number[] = [];
  for (let p = promoterWindow.start; p <= promoterWindow.end; p += tickInterval) ticks.push(p);

  return (
    <section className="panel aln-hero">
      <div className="aln-hero-header">
        <div>
          <p className="eyebrow">Sequence alignment · {rows.length} accessions</p>
          <h2>FT promoter — NW global alignment</h2>
          <p className="panel-copy" style={{ marginBottom: 0 }}>
            {useAligned ? 'Needleman–Wunsch star topology (Col-0 reference).' : 'Coordinate-aligned to TSS.'}{' '}
            Ctrl+scroll to zoom · drag to pan.
          </p>
        </div>
        <div className="zoom-btns">
          <button
            className={`cluster-btn${clustered ? ' is-active' : ''}`}
            onClick={() => setClustered((v) => !v)}
            title="Cluster sequences by SNP distance (UPGMA)"
          >
            cluster
          </button>
          <button onClick={() => { const n = Math.min(MAX_CW, cw * 1.5); setCw(n); }}>+</button>
          <span>{cw < 1 ? cw.toFixed(1) : Math.round(cw)}px</span>
          <button onClick={() => { const n = Math.max(MIN_CW, cw / 1.5); setCw(n); }}>−</button>
        </div>
      </div>

      {/* ── Main alignment body ────────────────────────────────────────── */}
      <div className="aln-body">

        {/* Dendrogram (when clustered) */}
        {clusterResult && (
          <DendroPanel result={clusterResult} totalRows={orderedRows.length} rowH={ROW_H + TFBS_H} />
        )}

        {/* Fixed label column */}
        <div className="aln-labels">
          <div style={{ height: RULER_H + GENE_H }} />
          {orderedRows.map(({ accession, profile }) => (
            <div key={accession.accessionName} className="aln-label-row" style={{ height: ROW_H + TFBS_H }}>
              <span className="aln-name">{accession.accessionName}</span>
              <span className="aln-meta">
                {profile.sequenceIdentity}%
                {FT_LD[accession.accessionNumber] != null && (
                  <> · {FT_LD[accession.accessionNumber]}d</>
                )}
              </span>
            </div>
          ))}
        </div>

        {/* Scrollable viewport */}
        <div ref={vpRef} className="aln-viewport" onMouseDown={onMouseDown}>
          <div
            className="aln-inner"
            style={{ width: totalW, transform: `translateX(${-scrollX}px)` }}
          >
            {/* ── Ruler ──────────────────────────────────────────────────── */}
            <div style={{ position: 'relative', height: RULER_H }}>
              {/* TSS marker */}
              {(() => {
                const col = refPosToCol[0 - promoterWindow.start] ?? 0;
                return (
                  <div style={{
                    position: 'absolute', left: col * cw, top: 0, bottom: 0,
                    width: 2, background: 'var(--warning)', opacity: 0.75,
                  }} />
                );
              })()}
              {ticks.map((pos) => {
                const col = refPosToCol[pos - promoterWindow.start] ?? (pos - promoterWindow.start);
                const x = col * cw;
                if (x < scrollX - 200 || x > scrollX + vpW + 200) return null;
                const chrCoord = promoterWindow.tssCoord != null
                  ? (promoterWindow.strand === '-'
                      ? promoterWindow.tssCoord - pos
                      : promoterWindow.tssCoord + pos)
                  : null;
                return (
                  <div key={pos} className="aln-tick" style={{ left: x }}>
                    <div className="aln-tick-line" />
                    <span>{pos === 0 ? 'TSS' : pos}</span>
                    {chrCoord != null && (
                      <span className="aln-tick-chr">
                        {(promoterWindow.chromosome ?? '').replace(/^\d+_/, '')}:{chrCoord.toLocaleString()}
                      </span>
                    )}
                  </div>
                );
              })}
            </div>

            {/* ── Gene model track (always shown) ────────────────────────── */}
            {(() => {
              const tssUnaligned = -promoterWindow.start;           // index of TSS in unaligned ref
              const tssCol  = refPosToCol[tssUnaligned] ?? tssUnaligned;
              const tssPx   = tssCol * cw;
              const midY    = GENE_H / 2;
              const isPlus  = (promoterWindow.strand ?? '+') === '+';

              // direction chevrons: one every ~36 screen-px along gene body
              const features = promoterWindow.geneFeatures ?? [];
              const lastFeat = features.length ? features[features.length - 1] : null;
              const geneEndPx = lastFeat
                ? (refPosToCol[lastFeat.end - promoterWindow.start] ?? 0) * cw
                : tssPx + 40;
              const chevronStep = Math.max(36, 36);
              const chevrons: number[] = [];
              if (geneEndPx > tssPx) {
                const cStart2 = Math.max(tssPx, scrollX);
                const cEnd2   = Math.min(geneEndPx, scrollX + vpW);
                for (let px = cStart2 + chevronStep / 2; px < cEnd2; px += chevronStep) {
                  chevrons.push(px);
                }
              }

              return (
                <div style={{ position: 'relative', height: GENE_H }}>
                  {/* backbone — full width */}
                  <div style={{
                    position: 'absolute', top: midY - 0.5, left: 0, right: 0,
                    height: 1, background: 'var(--muted)', opacity: 0.35,
                  }} />

                  {/* direction chevrons along gene body */}
                  {cw > 0.6 && chevrons.map((px, i) => (
                    <div key={i} style={{
                      position: 'absolute', left: px - 4, top: midY - 5,
                      fontSize: 10, lineHeight: '10px', color: 'var(--muted)',
                      opacity: 0.5, pointerEvents: 'none', userSelect: 'none',
                    }}>
                      {isPlus ? '›' : '‹'}
                    </div>
                  ))}

                  {/* TSS bent arrow */}
                  {tssPx >= scrollX - 40 && tssPx <= scrollX + vpW + 40 && (
                    <svg
                      style={{ position: 'absolute', left: tssPx - 1, top: 0, overflow: 'visible', pointerEvents: 'none' }}
                      width={1} height={GENE_H}
                    >
                      {/* vertical tick */}
                      <line x1={0} y1={2} x2={0} y2={midY} stroke="var(--accent)" strokeWidth={1.5} strokeLinecap="round" />
                      {/* horizontal arrow */}
                      <line x1={0} y1={midY} x2={isPlus ? 14 : -14} y2={midY} stroke="var(--accent)" strokeWidth={1.5} strokeLinecap="round" />
                      {/* arrowhead */}
                      {isPlus
                        ? <polygon points={`12,${midY - 3} 18,${midY} 12,${midY + 3}`} fill="var(--accent)" />
                        : <polygon points={`-12,${midY - 3} -18,${midY} -12,${midY + 3}`} fill="var(--accent)" />
                      }
                      {/* TSS label */}
                      <text x={isPlus ? 20 : -20} y={midY - 3} fontSize={7} fill="var(--accent)"
                        textAnchor={isPlus ? 'start' : 'end'}
                        fontFamily="IBM Plex Mono, monospace" fontWeight="600">
                        TSS
                      </text>
                    </svg>
                  )}

                  {/* exon / UTR blocks */}
                  {features.map((feat, fi) => {
                    const si = feat.start - promoterWindow.start;
                    const ei = feat.end   - promoterWindow.start;
                    const cS = refPosToCol[si] ?? si;
                    const cE = refPosToCol[ei] ?? ei;
                    const fx = cS * cw;
                    const fw = Math.max((cE - cS) * cw, 1);
                    if (fx + fw < scrollX - 50 || fx > scrollX + vpW + 50) return null;
                    const isExon = feat.type === 'exon';
                    const fh = isExon ? GENE_H - 4 : GENE_H / 2;
                    return (
                      <div key={fi} title={`${feat.type} ${feat.start}..${feat.end} bp`} style={{
                        position: 'absolute',
                        left: fx, width: fw,
                        top: (GENE_H - fh) / 2, height: fh,
                        background: isExon ? '#1d4ed8' : '#94a3b8',
                        borderRadius: isExon ? 2 : 1,
                        opacity: 0.85,
                      }} />
                    );
                  })}
                </div>
              );
            })()}

            {/* ── Sequence rows ───────────────────────────────────────────── */}
            {orderedRows.map(({ accession, profile }) => {
              const seq = getSeq(profile);
              if (!seq) return null;
              const isRef = accession.accessionName === refRow.accession.accessionName;
              const myPosToCol = qPosMap.get(accession.accessionName) ?? [];

              return (
                <div
                  key={accession.accessionName}
                  style={{ position: 'relative', height: ROW_H + TFBS_H }}
                >
                  {/* ── Row background strip at low zoom ───────────────── */}
                  {cw < 1.5 && (
                    <div style={{
                      position: 'absolute', left: startCol * cw, top: 0,
                      width: (endCol - startCol) * cw, height: ROW_H,
                      background: isRef ? 'rgba(250,204,21,0.08)' : 'rgba(209,213,219,0.18)',
                      borderBottom: '1px solid var(--border)',
                    }} />
                  )}

                  {/* ── Base cells (virtualized) ───────────────────────── */}
                  {cw >= 1.5 && Array.from({ length: endCol - startCol }, (_, i) => {
                    const col = startCol + i;
                    const base = seq[col] ?? 'N';
                    const refBase = refSeq[col] ?? '-';
                    const isGapQ = base === '-';
                    const isGapR = refBase === '-';
                    const isSnp = !isRef && !isGapQ && !isGapR
                      && base !== refBase && base !== 'N' && refBase !== 'N';

                    let bg: string;
                    let alpha: number;
                    if (isGapQ) {
                      bg = '#e5e7eb'; alpha = 1;
                    } else if (!isRef && isGapR) {
                      bg = INS_BG; alpha = 1;
                    } else {
                      bg = BASE_BG[base] ?? '#d1d5db';
                      alpha = isRef ? 1 : (isSnp ? 0.95 : 0.28);
                    }

                    return (
                      <div
                        key={col}
                        style={{
                          position: 'absolute',
                          left: col * cw,
                          top: 0,
                          width: cw,
                          height: ROW_H,
                          background: bg,
                          opacity: alpha,
                          boxShadow: isSnp ? 'inset 0 0 0 1.5px rgba(185,28,28,0.85)' : undefined,
                        }}
                      >
                        {showText && !isGapQ && (
                          <span style={{
                            display: 'block',
                            textAlign: 'center',
                            lineHeight: `${ROW_H}px`,
                            fontSize: Math.min(cw * 0.72, 13),
                            fontFamily: 'IBM Plex Mono, Courier New, monospace',
                            color: BASE_TEXT[base] ?? '#374151',
                            fontWeight: isSnp ? 700 : 400,
                            userSelect: 'none',
                            pointerEvents: 'none',
                          }}>
                            {base}
                          </span>
                        )}
                      </div>
                    );
                  })}

                  {/* ── TFBS chips ─────────────────────────────────────── */}
                  {profile.motifInstances.map((inst: MotifInstance, idx: number) => {
                    if (hiddenMotifs.has(inst.motifId)) return null;
                    if (inst.score < minScore) return null;
                    const motif = motifs.find((m) => m.id === inst.motifId);
                    if (!motif) return null;
                    const si = inst.start - promoterWindow.start;
                    const ei = si + inst.length - 1;
                    const cStart = myPosToCol[si];
                    const cEnd   = myPosToCol[ei];
                    if (cStart === undefined || cEnd === undefined) return null;
                    const chipX = cStart * cw;
                    const chipW = Math.max((cEnd - cStart + 1) * cw, 1);
                    if (chipX + chipW < scrollX - 20 || chipX > scrollX + vpW + 20) return null;
                    const isHL  = inst.motifId === highlightMotif;

                    const deltaGlow =
                      inst.delta === 'variant' ? 'inset 0 0 0 1.5px rgba(245,158,11,0.9)' :
                      inst.delta === 'gained'  ? 'inset 0 0 0 1.5px rgba(34,197,94,0.9)'  :
                      undefined;
                    const hlGlow = isHL ? '0 0 0 1.5px rgba(23,33,33,0.55)' : undefined;
                    const boxShadow = [hlGlow, deltaGlow].filter(Boolean).join(', ') || undefined;

                    return (
                      <div
                        key={`${inst.motifId}-${inst.start}-${idx}`}
                        title={`${motif.label} ${inst.start}..${inst.start + inst.length - 1} · score ${inst.score.toFixed(2)} · ${inst.delta}`}
                        style={{
                          position: 'absolute',
                          left: chipX,
                          top: ROW_H + 1,
                          width: chipW,
                          height: TFBS_H - 2,
                          background: motif.color,
                          opacity: isHL ? 1 : 0.6,
                          borderRadius: 2,
                          overflow: 'hidden',
                          boxShadow,
                        }}
                      >
                        {chipW > 22 && (
                          <span style={{
                            fontSize: 7, color: 'white', fontWeight: 700,
                            paddingLeft: 3, lineHeight: `${TFBS_H - 2}px`,
                            userSelect: 'none', pointerEvents: 'none',
                            fontFamily: 'IBM Plex Mono, monospace',
                          }}>{motif.label}</span>
                        )}
                      </div>
                    );
                  })}
                </div>
              );
            })}
          </div>
        </div>
      </div>

      {/* ── Motif legend ───────────────────────────────────────────────────── */}
      <div className="aln-legend">
        {motifs.map((motif) => {
          const isHL     = motif.id === highlightMotif;
          const isHidden = hiddenMotifs.has(motif.id);
          return (
            <div
              key={motif.id}
              className={`aln-legend-chip${isHL ? ' is-active' : ''}${isHidden ? ' is-hidden' : ''}`}
              style={{ borderColor: motif.color }}
            >
              <button
                type="button"
                className="aln-legend-chip-label"
                onClick={() => !isHidden && onHighlightMotif?.(motif.id)}
              >
                <span style={{ background: motif.color }} />
                {motif.label}
              </button>
              <button
                type="button"
                className="aln-legend-chip-eye"
                title={isHidden ? `Show ${motif.label}` : `Hide ${motif.label}`}
                onClick={() => toggleHide(motif.id)}
              >
                {isHidden ? '○' : '●'}
              </button>
            </div>
          );
        })}
      </div>

      {/* ── Score filter ───────────────────────────────────────────────────── */}
      <div className="aln-score-filter">
        <span>Min score</span>
        <input
          type="range" min={0.7} max={1.0} step={0.01} value={minScore}
          onChange={(e) => setMinScore(Number(e.target.value))}
        />
        <span>{minScore.toFixed(2)}</span>
      </div>

      {/* ── PWM sequence logo for highlighted motif ───────────────────────── */}
      {(() => {
        const active = motifs.find((m) => m.id === highlightMotif);
        return active ? <MotifLogo motif={active} /> : null;
      })()}
    </section>
  );
}
