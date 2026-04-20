export type Strand = '+' | '-';

export interface AccessionMeta {
  id: string;
  accessionNumber: number;
  accessionName: string;
  stockId: string;
  countryCode: string;
  country: string;
  latitude: number;
  longitude: number;
}

export interface MotifDefinition {
  id: string;
  label: string;
  family: string;
  element?: string;
  color: string;
  description: string;
  source?: string;
  pfm?: { A: number[]; C: number[]; G: number[]; T: number[] };
  anchorStart?: number;
  length?: number;
  strand?: '+' | '-';
}

export interface MotifInstance {
  motifId: string;
  start: number;
  length: number;
  strand: Strand;
  score: number;
  delta: 'gained' | 'lost' | 'stable';
}

export interface StructuralVariant {
  id: string;
  start: number;
  end: number;
  type: 'insertion' | 'deletion' | 'duplication';
  impact: 'low' | 'moderate' | 'high';
}

export interface AccessionProfile {
  accessionName: string;
  accessionNumber: number;
  promoterLength: number;
  sequenceIdentity: number;
  variantCount: number;
  motifScores: Record<string, number>;
  motifInstances: MotifInstance[];
  structuralVariants: StructuralVariant[];
  notes: string[];
  sequence?: string;
  alignedSequence?: string;
}

export interface GeneFeature {
  type: 'exon' | '5utr' | '3utr';
  start: number;  // bp relative to TSS
  end: number;
}

export interface PromoterWindow {
  start: number;
  end: number;
  tss: number;
  chromosome?: string;
  tssCoord?: number;
  strand?: '+' | '-';
  geneFeatures?: GeneFeature[];
}

export interface Dataset {
  datasetLabel: string;
  gene: string;
  promoterWindow: PromoterWindow;
  accessions: AccessionMeta[];
  motifs: MotifDefinition[];
  profiles: AccessionProfile[];
  citations: string[];
  isDemo: boolean;
}
