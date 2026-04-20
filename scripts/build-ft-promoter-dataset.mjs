import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import https from 'node:https';
import zlib from 'node:zlib';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const projectRoot = path.resolve(__dirname, '..');
const outputPath = path.join(projectRoot, 'public', 'data', 'ft-promoter.json');

const promoterWindow = { start: -2000, end: 200, tss: 0 };
const promoterLength = promoterWindow.end - promoterWindow.start + 1;
const ftGeneName = 'AT1G65480';
const baseUrl = 'https://1001genomes.org/data/1001Gp/27genomes/releases/current';

const accessions = [
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

const fileStemOverrides = {
  22001: '22001f',
};

const motifs = [
  {
    id: 'nfyc2',
    matrixId: 'UN0690.2',
    label: 'NFYC2',
    family: 'NF-Y / CO complex',
    element: 'CCAAT-box',
    color: '#207178',
    strongHitThreshold: 0.86,
    description: 'JASPAR UN0690.2 Arabidopsis NFYC2 PWM, relevant to the CO-NF-Y FT activation complex.',
    pfm: {
      A: [3455, 169, 116, 79, 120, 236, 343, 3557],
      C: [259, 4264, 31, 105, 53, 103, 4053, 184],
      G: [443, 87, 4428, 23, 4421, 3438, 174, 508],
      T: [523, 160, 105, 4473, 86, 903, 110, 431],
    },
  },
  {
    id: 'agl15',
    matrixId: 'MA0548.3',
    label: 'AGL15',
    family: 'MADS-box',
    element: 'CArG-box',
    color: '#d95d39',
    strongHitThreshold: 0.8,
    description: 'JASPAR MA0548.3 Arabidopsis AGL15 PWM as a MADS-box proxy for FT promoter repression-like sites.',
    pfm: {
      A: [48, 16, 84, 0, 0, 200, 88, 162, 7, 89, 138, 125, 0, 325, 543, 475],
      C: [33, 5, 8, 596, 435, 80, 85, 0, 8, 45, 42, 1, 0, 50, 13, 7],
      G: [5, 3, 21, 0, 0, 50, 53, 4, 0, 31, 168, 473, 585, 28, 5, 29],
      T: [513, 575, 486, 3, 164, 269, 373, 433, 584, 434, 251, 0, 14, 196, 38, 88],
    },
  },
  {
    id: 'pif4',
    matrixId: 'MA0561.1',
    label: 'PIF4',
    family: 'bHLH',
    element: 'G-box',
    color: '#008148',
    strongHitThreshold: 0.84,
    description: 'JASPAR MA0561.1 Arabidopsis PIF4 PWM for light and temperature responsive binding.',
    pfm: {
      A: [0, 335, 0, 49, 0, 0, 0, 64],
      C: [335, 0, 335, 0, 0, 0, 99, 206],
      G: [0, 0, 0, 286, 0, 335, 183, 65],
      T: [0, 0, 0, 0, 335, 0, 53, 0],
    },
  },
  {
    id: 'tcp20',
    matrixId: 'MA1065.3',
    label: 'TCP20',
    family: 'TCP',
    element: 'TCP-box',
    color: '#fb8500',
    strongHitThreshold: 0.84,
    description: 'JASPAR MA1065.3 Arabidopsis TCP20 PWM for class-I TCP promoter signals.',
    pfm: {
      A: [106, 19, 77, 4, 37, 379, 13, 4, 9, 887, 35],
      C: [52, 53, 7, 1, 3, 158, 915, 962, 876, 8, 769],
      G: [767, 12, 877, 962, 914, 159, 4, 0, 4, 61, 60],
      T: [45, 886, 9, 3, 16, 274, 38, 4, 81, 14, 106],
    },
  },
  {
    id: 'cdf2',
    matrixId: 'MA0973.2',
    label: 'CDF2',
    family: 'DOF',
    element: 'DOF-site',
    color: '#5f0f40',
    strongHitThreshold: 0.84,
    description: 'JASPAR MA0973.2 Arabidopsis CDF2 PWM, relevant to photoperiodic FT regulation.',
    pfm: {
      A: [608, 693, 980, 979, 868, 17, 53],
      C: [92, 16, 10, 10, 29, 7, 295],
      G: [150, 41, 4, 2, 81, 969, 154],
      T: [150, 250, 5, 9, 22, 6, 498],
    },
  },
  {
    id: 'cca1',
    matrixId: 'MA0972.1',
    label: 'CCA1',
    family: 'Circadian MYB-related',
    element: 'EE/LBS',
    color: '#8a3ffc',
    strongHitThreshold: 0.86,
    description: 'JASPAR MA0972.1 Arabidopsis CCA1 PWM for circadian-associated promoter signals.',
    pfm: {
      A: [934, 960, 985, 35, 992, 1, 5, 21],
      C: [1, 1, 4, 0, 5, 12, 991, 91],
      G: [10, 37, 10, 5, 1, 1, 0, 14],
      T: [54, 2, 1, 960, 2, 986, 4, 874],
    },
  },
  {
    id: 'myc3',
    matrixId: 'MA0568.2',
    label: 'MYC3',
    family: 'bHLH',
    element: 'G-box',
    color: '#3a86ff',
    strongHitThreshold: 0.82,
    description: 'JASPAR MA0568.2 Arabidopsis MYC3 bHLH — G-box (CACGTG) binding; integrates JA signalling with flowering.',
    pfm: {
      A: [10, 97, 0, 16, 1, 0],
      C: [89, 0, 84, 0, 2, 0],
      G: [0, 2, 0, 84, 0, 89],
      T: [0, 1, 16, 0, 97, 10],
    },
  },
  {
    id: 'smz',
    matrixId: 'UN0875.1',
    label: 'SMZ',
    family: 'AP2/ERF',
    element: 'GCC-box',
    color: '#ef233c',
    strongHitThreshold: 0.80,
    description: 'JASPAR UN0875.1 Arabidopsis SMZ (SCHLAFMÜTZE) AP2/ERF — binds GCC-box; represses FT in non-inductive conditions.',
    pfm: {
      A: [0, 0, 0, 0, 0, 0, 135, 0],
      C: [105, 119, 0, 147, 0, 0, 12, 147],
      G: [0, 28, 0, 0, 147, 25, 0, 0],
      T: [42, 0, 147, 0, 0, 122, 0, 0],
    },
  },
  {
    id: 'tem1',
    matrixId: 'MA1800.1',
    label: 'TEM1',
    family: 'RAV / B3',
    element: 'RAV1-A',
    color: '#bc6c25',
    strongHitThreshold: 0.82,
    description: 'JASPAR MA1800.1 Arabidopsis TEM1 (TEMPRANILLO 1) RAV/B3 — binds RAV1-A element; delays flowering by repressing FT.',
    pfm: {
      A: [219.542, 87.8331, 874.61, 999.532, 0.156152, 999.532, 0.172295, 599.65],
      C: [97.7096, 911.816, 15.7714, 0.156152, 999.532, 0.156152, 0.172295, 125.125],
      G: [682.505, 0.175316, 109.463, 0.156152, 0.156152, 0.156152, 241.385, 0.24975],
      T: [0.243665, 0.175316, 0.156152, 0.156152, 0.156152, 0.156152, 758.27, 274.975],
    },
  },
];

const complement = { A: 'T', C: 'G', G: 'C', T: 'A', W: 'W', S: 'S', M: 'K', K: 'M', R: 'Y', Y: 'R', B: 'V', D: 'H', H: 'D', V: 'B', N: 'N' };
const bases = ['A', 'C', 'G', 'T'];
const backgroundProbability = 0.25;

function buildPwm(motif) {
  const length = motif.pfm.A.length;
  const pwm = { A: [], C: [], G: [], T: [] };

  for (let index = 0; index < length; index += 1) {
    const total = bases.reduce((sum, base) => sum + motif.pfm[base][index], 0);
    for (const base of bases) {
      const probability = (motif.pfm[base][index] + 0.1) / (total + 0.4);
      pwm[base][index] = Math.log2(probability / backgroundProbability);
    }
  }

  const maxScore = Array.from({ length }, (_, index) => Math.max(...bases.map((base) => pwm[base][index]))).reduce((sum, value) => sum + value, 0);
  const minScore = Array.from({ length }, (_, index) => Math.min(...bases.map((base) => pwm[base][index]))).reduce((sum, value) => sum + value, 0);

  return { ...motif, pwm, length, maxScore, minScore };
}

function reverseComplementPwm(pwm) {
  return {
    A: [...pwm.T].reverse(),
    C: [...pwm.G].reverse(),
    G: [...pwm.C].reverse(),
    T: [...pwm.A].reverse(),
  };
}

const scoredMotifs = motifs.map(buildPwm).map((motif) => ({
  ...motif,
  reversePwm: reverseComplementPwm(motif.pwm),
}));

function fetchBuffer(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (response) => {
      if (response.statusCode !== 200) {
        reject(new Error(`Request failed for ${url}: ${response.statusCode}`));
        response.resume();
        return;
      }

      const chunks = [];
      response.on('data', (chunk) => chunks.push(chunk));
      response.on('end', () => resolve(Buffer.concat(chunks)));
      response.on('error', reject);
    }).on('error', reject);
  });
}

async function fetchGunzipText(url) {
  const gzBuffer = await fetchBuffer(url);
  return zlib.gunzipSync(gzBuffer).toString('utf8');
}

function parseAttributes(field) {
  return Object.fromEntries(field.split(';').filter(Boolean).map((entry) => {
    const [key, value = ''] = entry.split('=');
    return [key, value];
  }));
}

function findFtGene(gffText) {
  for (const line of gffText.split(/\r?\n/)) {
    if (!line || line.startsWith('#')) {
      continue;
    }

    const parts = line.split('\t');
    if (parts.length < 9 || parts[2] !== 'gene') {
      continue;
    }

    const attributes = parseAttributes(parts[8]);
    const name = attributes.Name?.split(',') ?? [];
    if (name.includes(ftGeneName)) {
      return {
        seqId: parts[0],
        start: Number(parts[3]),
        end: Number(parts[4]),
        strand: parts[6].startsWith('-') ? '-' : '+',
      };
    }
  }

  throw new Error(`Could not find ${ftGeneName} in annotation`);
}

function extractSequence(fastaText, seqId) {
  let active = false;
  const chunks = [];
  for (const line of fastaText.split(/\r?\n/)) {
    if (!line) {
      continue;
    }
    if (line.startsWith('>')) {
      const currentId = line.slice(1).trim().split(/\s+/)[0];
      active = currentId === seqId;
      continue;
    }
    if (active) {
      chunks.push(line.trim().toUpperCase());
    }
  }

  if (chunks.length === 0) {
    throw new Error(`Could not find sequence ${seqId} in FASTA`);
  }

  return chunks.join('');
}

function reverseComplement(sequence) {
  return Array.from(sequence).reverse().map((base) => complement[base] ?? 'N').join('');
}

function slicePromoter(sequence, gene) {
  if (gene.strand === '+') {
    const regionStart = gene.start + promoterWindow.start;
    const regionEnd = gene.start + promoterWindow.end;
    const leftPad = Math.max(0, 1 - regionStart);
    const rightPad = Math.max(0, regionEnd - sequence.length);
    const sliceStart = Math.max(1, regionStart);
    const sliceEnd = Math.min(sequence.length, regionEnd);
    const core = sequence.slice(sliceStart - 1, sliceEnd);
    return `${'N'.repeat(leftPad)}${core}${'N'.repeat(rightPad)}`.slice(0, promoterLength);
  }

  const regionStart = gene.end - promoterWindow.end;
  const regionEnd = gene.end - promoterWindow.start;
  const leftPad = Math.max(0, 1 - regionStart);
  const rightPad = Math.max(0, regionEnd - sequence.length);
  const sliceStart = Math.max(1, regionStart);
  const sliceEnd = Math.min(sequence.length, regionEnd);
  const core = sequence.slice(sliceStart - 1, sliceEnd);
  const padded = `${'N'.repeat(leftPad)}${core}${'N'.repeat(rightPad)}`.slice(0, promoterLength);
  return reverseComplement(padded);
}

function scoreWindow(windowSeq, pwm) {
  let score = 0;
  for (let index = 0; index < windowSeq.length; index += 1) {
    const base = windowSeq[index];
    const value = pwm[base]?.[index];
    if (value === undefined) {
      return Number.NEGATIVE_INFINITY;
    }
    score += value;
  }
  return score;
}

function normalizeScore(rawScore, motif) {
  if (motif.maxScore === motif.minScore) {
    return 0;
  }
  const normalized = (rawScore - motif.minScore) / (motif.maxScore - motif.minScore);
  return Math.max(0, Math.min(1, normalized));
}

function scanMotif(sequence, motif) {
  const hits = [];
  const windowScores = [];
  let bestScore = 0;

  for (let index = 0; index <= sequence.length - motif.length; index += 1) {
    const windowSeq = sequence.slice(index, index + motif.length);
    if (windowSeq.includes('N')) {
      continue;
    }
    const plusNormalized = normalizeScore(scoreWindow(windowSeq, motif.pwm), motif);
    const minusNormalized = normalizeScore(scoreWindow(windowSeq, motif.reversePwm), motif);
    const normalized = Math.max(plusNormalized, minusNormalized);
    const strand = minusNormalized > plusNormalized ? '-' : '+';
    bestScore = Math.max(bestScore, normalized);
    windowScores.push(normalized);

    if (normalized >= motif.strongHitThreshold) {
      hits.push({
        start: promoterWindow.start + index,
        length: motif.length,
        strand,
        score: Number(normalized.toFixed(3)),
      });
    }
  }

  hits.sort((left, right) => right.score - left.score || left.start - right.start);
  windowScores.sort((left, right) => right - left);
  const topWindowScores = windowScores.slice(0, 5);
  const promoterScore = topWindowScores.length
    ? Number((topWindowScores.reduce((sum, value) => sum + value, 0) / topWindowScores.length).toFixed(3))
    : 0;

  return {
    bestScore: Number(bestScore.toFixed(3)),
    promoterScore,
    hits: hits.slice(0, 8),
  };
}

function compareToReference(sequence, reference) {
  let informative = 0;
  let matches = 0;
  let differences = 0;

  for (let index = 0; index < Math.min(sequence.length, reference.length); index += 1) {
    const left = sequence[index];
    const right = reference[index];
    if (left === 'N' || right === 'N') {
      continue;
    }
    informative += 1;
    if (left === right) {
      matches += 1;
    } else {
      differences += 1;
    }
  }

  const sequenceIdentity = informative === 0 ? 0 : Number(((matches / informative) * 100).toFixed(2));
  return { sequenceIdentity, variantCount: differences };
}

function summariseNotes(profile, motifsById) {
  const notes = [];
  if ((profile.motifScores.nfyc2 ?? 0) > (profile.motifScores.agl15 ?? 0)) {
    notes.push('NF-Y-like promoter signal exceeds the AGL15-like MADS-box signal in the extracted FT promoter window.');
  }
  if ((profile.motifScores.pif4 ?? 0) >= 0.9) {
    notes.push('Contains a strong PIF4-like PWM hit in the scanned promoter window.');
  }
  if (profile.variantCount > 120) {
    notes.push('Promoter sequence is relatively diverged from Col-0 across the fixed -2000..+200 window.');
  }
  const topMotif = Object.entries(profile.motifScores).sort((left, right) => right[1] - left[1])[0];
  if (topMotif) {
    notes.push(`Best motif class by normalized match score: ${motifsById.get(topMotif[0])?.label ?? topMotif[0]}.`);
  }
  return notes;
}

async function buildProfile(accession, referenceSequence) {
  const accessionCode = String(accession.accessionNumber);
  const fileStem = fileStemOverrides[accession.accessionNumber] ?? accessionCode;
  const gffUrl = `${baseUrl}/annotations/genes_v05_${fileStem}.gff.gz`;
  const fastaUrl = `${baseUrl}/assemblies/${fileStem}.scaffolds_corrected.v2.1.fasta.gz`;

  const [gffText, fastaText] = await Promise.all([
    fetchGunzipText(gffUrl),
    fetchGunzipText(fastaUrl),
  ]);

  const ftGene = findFtGene(gffText);
  const chromosomeSequence = extractSequence(fastaText, ftGene.seqId);
  const promoterSequence = slicePromoter(chromosomeSequence, ftGene);
  const comparison = compareToReference(promoterSequence, referenceSequence ?? promoterSequence);

  const motifScores = {};
  const motifInstances = [];
  for (const motif of scoredMotifs) {
    const scan = scanMotif(promoterSequence, motif);
    motifScores[motif.id] = scan.promoterScore;
    const referenceHits = referenceSequence ? scanMotif(referenceSequence, motif).hits.length : scan.hits.length;
    for (const hit of scan.hits) {
      const delta = scan.hits.length > referenceHits ? 'gained' : scan.hits.length < referenceHits ? 'lost' : 'stable';
      motifInstances.push({ motifId: motif.id, ...hit, delta });
    }
  }

  const profile = {
    accessionName: accession.accessionName,
    accessionNumber: accession.accessionNumber,
    promoterLength: promoterSequence.length,
    sequenceIdentity: comparison.sequenceIdentity,
    variantCount: comparison.variantCount,
    motifScores,
    motifInstances,
    structuralVariants: [],
    notes: [],
  };

  profile.notes = summariseNotes(profile, new Map(scoredMotifs.map((motif) => [motif.id, motif])));
  return { profile, promoterSequence, ftGene };
}

async function main() {
  console.log('Building FT promoter dataset from 1001 Genomes Plus Phase 1...');
  const col0 = accessions.find((entry) => entry.accessionName === 'Col-0');
  if (!col0) {
    throw new Error('Col-0 metadata missing');
  }

  const reference = await buildProfile(col0, null);
  const profiles = [reference.profile];
  const sourceSummary = {
    [col0.accessionName]: {
      gene: reference.ftGene,
    },
  };

  for (const accession of accessions) {
    if (accession.accessionName === 'Col-0') {
      continue;
    }
    console.log(`Processing ${accession.accessionName} (${accession.accessionNumber})`);
    const result = await buildProfile(accession, reference.promoterSequence);
    profiles.push(result.profile);
    sourceSummary[accession.accessionName] = { gene: result.ftGene };
  }

  const dataset = {
    datasetLabel: '1001 Genomes Plus Phase 1 FT promoter sequence scan',
    gene: 'FT / AT1G65480',
    promoterWindow,
    accessions,
    motifs: scoredMotifs.map(({ pfm, pwm, reversePwm, maxScore, minScore, length, strongHitThreshold, matrixId, ...motif }) => ({
      ...motif,
      source: `JASPAR ${matrixId}`,
    })),
    profiles,
    citations: [
      '1001 Genomes Plus Phase 1 assemblies: https://1001genomes.org/data/1001Gp/27genomes/releases/current/assemblies/',
      '1001 Genomes Plus Phase 1 annotations: https://1001genomes.org/data/1001Gp/27genomes/releases/current/annotations/',
      'FT locus identified from gene annotations using Name=AT1G65480.',
      'Motif scores are promoter-level PWM scores derived from JASPAR matrices, not experimental binding measurements.',
    ],
    isDemo: false,
    provenance: {
      builtAt: new Date().toISOString(),
      sourceSummary,
    },
  };

  await mkdir(path.dirname(outputPath), { recursive: true });
  await writeFile(outputPath, JSON.stringify(dataset, null, 2));
  console.log(`Wrote ${outputPath}`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});