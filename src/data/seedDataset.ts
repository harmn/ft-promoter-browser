import type {
  AccessionMeta,
  AccessionProfile,
  Dataset,
  MotifDefinition,
  MotifInstance,
  StructuralVariant,
} from '../types';

const accessions: AccessionMeta[] = [
  { id: '100265', accessionNumber: 1741, accessionName: 'KBS-Mac-74', stockId: 'CS78969', countryCode: 'USA', country: 'United States', latitude: 42.405, longitude: -85.398 },
  { id: '100266', accessionNumber: 6024, accessionName: 'Fly2-2', stockId: 'CS76864', countryCode: 'SWE', country: 'Sweden', latitude: 55.7501271, longitude: 13.3720798 },
  { id: '100267', accessionNumber: 6069, accessionName: 'Nyl-7', stockId: 'CS77137', countryCode: 'SWE', country: 'Sweden', latitude: 62.95681, longitude: 18.2763 },
  { id: '100268', accessionNumber: 6124, accessionName: 'T690', stockId: 'CS77309', countryCode: 'SWE', country: 'Sweden', latitude: 55.8362185, longitude: 13.2995709 },
  { id: '100269', accessionNumber: 6244, accessionName: 'TRA 01', stockId: 'CS77384', countryCode: 'SWE', country: 'Sweden', latitude: 62.9169, longitude: 18.4728 },
  { id: '100270', accessionNumber: 6909, accessionName: 'Col-0', stockId: 'CS76778', countryCode: 'USA', country: 'United States', latitude: 38.3, longitude: -92.3 },
  { id: '100271', accessionNumber: 6966, accessionName: 'Sq-1', stockId: 'CS77266', countryCode: 'UK', country: 'United Kingdom', latitude: 51.4083, longitude: -0.6383 },
  { id: '100272', accessionNumber: 8236, accessionName: 'HSm', stockId: 'CS76941', countryCode: 'CZE', country: 'Czech Republic', latitude: 49.33, longitude: 15.76 },
  { id: '100273', accessionNumber: 9075, accessionName: 'Lerik1-4', stockId: 'CS77023', countryCode: 'AZE', country: 'Azerbaijan', latitude: 38.7406, longitude: 48.6131 },
  { id: '100274', accessionNumber: 9537, accessionName: 'IP-Cum-1', stockId: 'CS76787', countryCode: 'ESP', country: 'Spain', latitude: 38.07, longitude: -6.66 },
  { id: '100275', accessionNumber: 9543, accessionName: 'IP-Gra-0', stockId: 'CS76886', countryCode: 'ESP', country: 'Spain', latitude: 36.77, longitude: -5.39 },
  { id: '100276', accessionNumber: 9638, accessionName: 'Noveg-3', stockId: 'CS77133', countryCode: 'RUS', country: 'Russia', latitude: 51.73, longitude: 80.86 },
  { id: '100277', accessionNumber: 9728, accessionName: 'Stiav-1', stockId: 'CS77279', countryCode: 'SVK', country: 'Slovakia', latitude: 48.46, longitude: 18.9 },
  { id: '100278', accessionNumber: 9764, accessionName: 'Qar-8a', stockId: 'CS76581', countryCode: 'LBN', country: 'Lebanon', latitude: 34.1, longitude: 35.84 },
  { id: '100279', accessionNumber: 9888, accessionName: 'IP-Pva-1', stockId: 'CS77197', countryCode: 'ESP', country: 'Spain', latitude: 40.93, longitude: -3.31 },
  { id: '100280', accessionNumber: 9905, accessionName: 'IP-Ven-0', stockId: 'CS78840', countryCode: 'ESP', country: 'Spain', latitude: 40.76, longitude: -4.01 },
  { id: '100281', accessionNumber: 9981, accessionName: 'Angit-1', stockId: 'CS76366', countryCode: 'ITA', country: 'Italy', latitude: 38.76, longitude: 16.24 },
  { id: '100282', accessionNumber: 10002, accessionName: 'TueWa1-2', stockId: 'CS76405', countryCode: 'GER', country: 'Germany', latitude: 48.53, longitude: 9.04 },
  { id: '100283', accessionNumber: 10015, accessionName: 'Shahdara', stockId: '—', countryCode: 'TJK', country: 'Tajikistan', latitude: 38.35, longitude: 68.48 },
  { id: '100284', accessionNumber: 10024, accessionName: 'Tnz-1', stockId: '—', countryCode: 'TZA', country: 'Tanzania', latitude: -2.87389, longitude: 36.2128 },
  { id: '100286', accessionNumber: 22001, accessionName: '85-3', stockId: '—', countryCode: 'CHN', country: 'China', latitude: 32.14, longitude: 115.06 },
  { id: '100287', accessionNumber: 22002, accessionName: '35-1', stockId: '—', countryCode: 'CHN', country: 'China', latitude: 27.94, longitude: 108.61 },
  { id: '100288', accessionNumber: 22003, accessionName: 'Taz-0', stockId: 'CS799913', countryCode: 'MAR', country: 'Morocco', latitude: 34.09166, longitude: -4.10258 },
  { id: '100289', accessionNumber: 22004, accessionName: 'Elh-2', stockId: 'CS799925', countryCode: 'MAR', country: 'Morocco', latitude: 31.47197, longitude: -7.40644 },
  { id: '100290', accessionNumber: 22005, accessionName: 'Rabacal-1', stockId: 'CS2107642', countryCode: 'POR', country: 'Portugal', latitude: 32.7536, longitude: -17.1297 },
  { id: '100291', accessionNumber: 22006, accessionName: 'Areeiro-1', stockId: '635AV', countryCode: 'POR', country: 'Portugal', latitude: 32.74, longitude: -16.93 },
  { id: '100292', accessionNumber: 22007, accessionName: 'ET-86.4', stockId: '—', countryCode: 'ETH', country: 'Ethiopia', latitude: 13.235084, longitude: 38.061701 },
];

const motifs: MotifDefinition[] = [
  { id: 'co-core', label: 'CONSTANS', family: 'CCT', element: 'CO-resp.', color: '#207178', anchorStart: -1725, length: 26, strand: '+', description: 'Photoperiod-responsive CO-like core site.' },
  { id: 'svpcaar', label: 'SVP', family: 'MADS-box', element: 'CArG-box', color: '#d95d39', anchorStart: -1490, length: 22, strand: '-', description: 'MADS-box occupancy candidate linked to floral repression.' },
  { id: 'flc-box', label: 'FLC', family: 'MADS-box', element: 'CArG-box', color: '#8a3ffc', anchorStart: -1315, length: 28, strand: '+', description: 'Putative FLC-associated cis-regulatory window.' },
  { id: 'nfy-triad', label: 'NF-Y', family: 'NF-Y', element: 'CCAAT-box', color: '#4e79a7', anchorStart: -980, length: 20, strand: '+', description: 'CCAAT-like element supporting CO complex binding.' },
  { id: 'tem1-site', label: 'TEM1', family: 'RAV / B3', element: 'RAV1-A', color: '#9c6f19', anchorStart: -760, length: 18, strand: '-', description: 'Repressor-associated TEM-like binding site.' },
  { id: 'spl9-site', label: 'SPL9', family: 'SBP', element: 'SBP-box', color: '#5f0f40', anchorStart: -510, length: 18, strand: '+', description: 'Age-dependent SPL candidate site near proximal promoter.' },
  { id: 'pif4-gbox', label: 'PIF4', family: 'bHLH', element: 'G-box', color: '#008148', anchorStart: -295, length: 16, strand: '+', description: 'G-box-like light and temperature responsive site.' },
  { id: 'tcp-box', label: 'TCP', family: 'TCP', element: 'TCP-box', color: '#fb8500', anchorStart: -115, length: 14, strand: '-', description: 'Proximal TCP-like motif near the FT TSS.' },
];

const clamp = (value: number, min: number, max: number) => Math.min(max, Math.max(min, value));

function hashName(value: string): number {
  return Array.from(value).reduce((acc, char, index) => acc + char.charCodeAt(0) * (index + 17), 0);
}

function buildStructuralVariants(seed: number): StructuralVariant[] {
  const variants: StructuralVariant[] = [];
  const variantCount = seed % 3;

  for (let index = 0; index < variantCount; index += 1) {
    const start = -1850 + ((seed + index * 97) % 1600);
    const size = 24 + ((seed + index * 41) % 110);
    const variantTypes: StructuralVariant['type'][] = ['insertion', 'deletion', 'duplication'];
    const impactLevels: StructuralVariant['impact'][] = ['low', 'moderate', 'high'];

    variants.push({
      id: `sv-${seed}-${index}`,
      start,
      end: start + size,
      type: variantTypes[(seed + index) % variantTypes.length],
      impact: impactLevels[(seed + index * 2) % impactLevels.length],
    });
  }

  return variants;
}

function buildInstances(accessionName: string): { motifScores: Record<string, number>; motifInstances: MotifInstance[] } {
  const seed = hashName(accessionName);
  const motifScores: Record<string, number> = {};
  const motifInstances: MotifInstance[] = [];

  motifs.forEach((motif, index) => {
    const shiftedSeed = seed + index * 131;
    const base = 0.36 + ((shiftedSeed % 59) / 100);
    const score = accessionName === 'Col-0' ? 0.82 - index * 0.04 : clamp(base, 0.18, 0.94);
    motifScores[motif.id] = Number(score.toFixed(2));

    if (score < 0.27) {
      return;
    }

    const offset = accessionName === 'Col-0' ? 0 : ((shiftedSeed % 31) - 15);
    const delta: MotifInstance['delta'] = score > 0.76 ? 'gained' : score < 0.38 ? 'lost' : 'stable';
    motifInstances.push({
      motifId: motif.id,
      start: motif.anchorStart + offset,
      length: motif.length,
      strand: motif.strand,
      score: Number(score.toFixed(2)),
      delta,
    });
  });

  return { motifScores, motifInstances };
}

function buildProfile(meta: AccessionMeta): AccessionProfile {
  if (meta.accessionName === 'Col-0') {
    const { motifScores, motifInstances } = buildInstances(meta.accessionName);
    return {
      accessionName: meta.accessionName,
      accessionNumber: meta.accessionNumber,
      promoterLength: 2200,
      sequenceIdentity: 100,
      variantCount: 0,
      motifScores,
      motifInstances,
      structuralVariants: [],
      notes: ['Reference-like FT promoter used as the baseline for comparisons.'],
    };
  }

  const seed = hashName(meta.accessionName);
  const { motifScores, motifInstances } = buildInstances(meta.accessionName);
  const variantCount = 8 + (seed % 34);
  const sequenceIdentity = Number((96.1 + ((seed % 33) / 10)).toFixed(1));
  const notes: string[] = [];

  if (motifScores['co-core'] > 0.75) {
    notes.push('High CO-like occupancy in the distal FT promoter.');
  }
  if (motifScores['flc-box'] < 0.4) {
    notes.push('Weak FLC-like signal relative to Col-0.');
  }
  if (variantCount > 28) {
    notes.push('Elevated local polymorphism load across the scanned promoter window.');
  }
  if (notes.length === 0) {
    notes.push('Promoter architecture remains broadly similar to the Col-0 reference.');
  }

  return {
    accessionName: meta.accessionName,
    accessionNumber: meta.accessionNumber,
    promoterLength: 2200,
    sequenceIdentity,
    variantCount,
    motifScores,
    motifInstances,
    structuralVariants: buildStructuralVariants(seed),
    notes,
  };
}

export const seedDataset: Dataset = {
  datasetLabel: '1001 Genomes Plus Phase 1 FT promoter browser',
  gene: 'FT / AT1G65480',
  promoterWindow: {
    start: -2000,
    end: 200,
    tss: 0,
  },
  accessions,
  motifs,
  profiles: accessions.map(buildProfile),
  citations: [
    'Igolkina et al. (2025), 1001 Genomes Plus Phase 1: 27 Arabidopsis thaliana genomes.',
    '1001genomes.org phase-1 accession table and linked JBrowse2 sessions.',
  ],
  isDemo: true,
};
