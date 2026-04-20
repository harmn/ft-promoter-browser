use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::Read;
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow, bail};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use reqwest::blocking::Client;
use serde::Serialize;

const BASE_URL: &str = "https://1001genomes.org/data/1001Gp/27genomes/releases/current";
const FT_GENE_NAME: &str = "AT1G65480";
const PROMOTER_START: i32 = -5000;
const PROMOTER_END: i32 = 5000;
const PROMOTER_LENGTH: usize = (PROMOTER_END - PROMOTER_START + 1) as usize;
const BACKGROUND_PROBABILITY: f64 = 0.25;

#[derive(Clone, Serialize)]
struct AccessionMeta {
    id: &'static str,
    #[serde(rename = "accessionNumber")]
    accession_number: u32,
    #[serde(rename = "accessionName")]
    accession_name: &'static str,
    #[serde(rename = "stockId")]
    stock_id: &'static str,
    #[serde(rename = "countryCode")]
    country_code: &'static str,
    country: &'static str,
    latitude: f64,
    longitude: f64,
}

#[derive(Clone)]
struct MotifDefinition {
    id: &'static str,
    matrix_id: &'static str,
    label: &'static str,
    family: &'static str,
    element: &'static str,
    color: &'static str,
    strong_hit_threshold: f64,
    description: &'static str,
    pfm: [Vec<f64>; 4],
}

#[derive(Clone)]
struct ScoredMotif {
    definition: MotifDefinition,
    pwm: [Vec<f64>; 4],
    reverse_pwm: [Vec<f64>; 4],
    length: usize,
    max_score: f64,
    min_score: f64,
}

#[derive(Serialize)]
struct PfmOutput {
    #[serde(rename = "A")] a: Vec<f64>,
    #[serde(rename = "C")] c: Vec<f64>,
    #[serde(rename = "G")] g: Vec<f64>,
    #[serde(rename = "T")] t: Vec<f64>,
}

#[derive(Serialize)]
struct MotifOutput {
    id: String,
    label: String,
    family: String,
    element: String,
    color: String,
    description: String,
    source: String,
    pfm: PfmOutput,
}

#[derive(Clone, Serialize)]
struct MotifInstance {
    #[serde(rename = "motifId")]
    motif_id: String,
    start: i32,
    length: usize,
    strand: String,
    score: f64,
    delta: String,
}

#[derive(Serialize)]
struct StructuralVariant {
    id: String,
    start: i32,
    end: i32,
    #[serde(rename = "type")]
    variant_type: String,
    impact: String,
}

#[derive(Serialize)]
struct AccessionProfile {
    #[serde(rename = "accessionName")]
    accession_name: String,
    #[serde(rename = "accessionNumber")]
    accession_number: u32,
    #[serde(rename = "promoterLength")]
    promoter_length: usize,
    #[serde(rename = "sequenceIdentity")]
    sequence_identity: f64,
    #[serde(rename = "variantCount")]
    variant_count: usize,
    #[serde(rename = "motifScores")]
    motif_scores: BTreeMap<String, f64>,
    #[serde(rename = "motifInstances")]
    motif_instances: Vec<MotifInstance>,
    #[serde(rename = "structuralVariants")]
    structural_variants: Vec<StructuralVariant>,
    notes: Vec<String>,
    sequence: String,
    #[serde(rename = "alignedSequence")]
    aligned_sequence: String,
}

#[derive(Serialize, Clone)]
struct GeneFeature {
    #[serde(rename = "type")]
    feature_type: String,  // "exon" | "5utr" | "3utr"
    start: i32,  // relative to TSS
    end: i32,
}

#[derive(Serialize)]
struct PromoterWindow {
    start: i32,
    end: i32,
    tss: i32,
    chromosome: String,
    #[serde(rename = "tssCoord")]
    tss_coord: i64,
    strand: String,
    #[serde(rename = "geneFeatures")]
    gene_features: Vec<GeneFeature>,
}

#[derive(Serialize)]
struct GeneSummary {
    #[serde(rename = "seqId")]
    seq_id: String,
    start: usize,
    end: usize,
    strand: String,
}

#[derive(Serialize)]
struct Provenance {
    #[serde(rename = "builtAt")]
    built_at: String,
    #[serde(rename = "sourceSummary")]
    source_summary: BTreeMap<String, AccessionSourceSummary>,
}

#[derive(Serialize)]
struct AccessionSourceSummary {
    gene: GeneSummary,
}

#[derive(Serialize)]
struct Dataset {
    #[serde(rename = "datasetLabel")]
    dataset_label: &'static str,
    gene: &'static str,
    #[serde(rename = "promoterWindow")]
    promoter_window: PromoterWindow,
    accessions: Vec<AccessionMeta>,
    motifs: Vec<MotifOutput>,
    profiles: Vec<AccessionProfile>,
    citations: Vec<&'static str>,
    #[serde(rename = "isDemo")]
    is_demo: bool,
    provenance: Provenance,
}

#[derive(Clone)]
struct GeneRecord {
    seq_id: String,
    start: usize,
    end: usize,
    strand: char,
}

struct ScanResult {
    promoter_score: f64,
    hits: Vec<MotifInstance>,
}

struct BuildOutput {
    profile: AccessionProfile,
    promoter_sequence: String,
    gene: GeneSummary,
    gene_features: Vec<GeneFeature>,
}

// Needleman-Wunsch global alignment.
// Scoring: match +2, mismatch -1, gap -2 (linear).
// Returns (aligned_seq1, aligned_seq2) with '-' for gap characters.
fn needleman_wunsch(seq1: &[u8], seq2: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let m = seq1.len();
    let n = seq2.len();
    let cols = n + 1;
    let mut dp = vec![0i32; (m + 1) * cols];
    let mut tb = vec![0u8; (m + 1) * cols]; // 0=diag 1=up 2=left

    for i in 0..=m { dp[i * cols] = (i as i32) * -2; }
    for j in 0..=n { dp[j] = (j as i32) * -2; }

    for i in 1..=m {
        for j in 1..=n {
            let s = if seq1[i - 1] == seq2[j - 1] { 2i32 } else { -1i32 };
            let diag = dp[(i - 1) * cols + (j - 1)] + s;
            let up   = dp[(i - 1) * cols + j] - 2;
            let left = dp[i * cols + (j - 1)] - 2;
            let idx  = i * cols + j;
            if diag >= up && diag >= left { dp[idx] = diag; tb[idx] = 0; }
            else if up >= left            { dp[idx] = up;   tb[idx] = 1; }
            else                          { dp[idx] = left;  tb[idx] = 2; }
        }
    }

    let mut a1 = Vec::with_capacity(m + n);
    let mut a2 = Vec::with_capacity(m + n);
    let (mut i, mut j) = (m, n);
    while i > 0 || j > 0 {
        if i > 0 && j > 0 && tb[i * cols + j] == 0 {
            a1.push(seq1[i - 1]); a2.push(seq2[j - 1]); i -= 1; j -= 1;
        } else if i > 0 && (j == 0 || tb[i * cols + j] == 1) {
            a1.push(seq1[i - 1]); a2.push(b'-'); i -= 1;
        } else {
            a1.push(b'-'); a2.push(seq2[j - 1]); j -= 1;
        }
    }
    a1.reverse(); a2.reverse();
    (a1, a2)
}

// Star-topology MSA: align every query against ref (Col-0), merge all gap columns
// inserted into the reference into a single coordinate space, then project every
// sequence into that space.  Returns one aligned string per input sequence
// (index 0 = ref, index 1..N = queries in the same order as `query_seqs`).
fn build_star_msa(ref_seq: &str, query_seqs: &[String]) -> Vec<String> {
    let ref_bytes = ref_seq.as_bytes();
    let ref_len = ref_bytes.len();

    let pairs: Vec<(Vec<u8>, Vec<u8>)> = query_seqs
        .par_iter()
        .map(|q| needleman_wunsch(ref_bytes, q.as_bytes()))
        .collect();

    // gap_map[k] = maximum number of query bases inserted before ref position k,
    // across all pairwise alignments.
    let mut gap_map: HashMap<usize, usize> = HashMap::new();
    for (ref_aln, _) in &pairs {
        let mut rp = 0usize;
        let mut run = 0usize;
        for &b in ref_aln {
            if b == b'-' {
                run += 1;
            } else {
                if run > 0 {
                    let e = gap_map.entry(rp).or_insert(0);
                    *e = (*e).max(run);
                    run = 0;
                }
                rp += 1;
            }
        }
        if run > 0 {
            let e = gap_map.entry(rp).or_insert(0);
            *e = (*e).max(run);
        }
    }

    let build_ref = || -> String {
        let mut out = Vec::new();
        for k in 0..ref_len {
            for _ in 0..gap_map.get(&k).copied().unwrap_or(0) { out.push(b'-'); }
            out.push(ref_bytes[k]);
        }
        for _ in 0..gap_map.get(&ref_len).copied().unwrap_or(0) { out.push(b'-'); }
        String::from_utf8(out).expect("aligned ref is valid UTF-8")
    };

    let build_query = |pair_ref: &[u8], pair_query: &[u8]| -> String {
        let mut ins: HashMap<usize, Vec<u8>> = HashMap::new();
        let mut rp = 0usize;
        for (i, &b) in pair_ref.iter().enumerate() {
            if b == b'-' { ins.entry(rp).or_default().push(pair_query[i]); }
            else { rp += 1; }
        }

        let mut out = Vec::new();
        let mut pi = 0usize;
        for k in 0..ref_len {
            while pi < pair_ref.len() && pair_ref[pi] == b'-' { pi += 1; }
            let max_gaps = gap_map.get(&k).copied().unwrap_or(0);
            let my_ins = ins.get(&k).map(|v| v.as_slice()).unwrap_or(&[]);
            for g in 0..max_gaps { out.push(if g < my_ins.len() { my_ins[g] } else { b'-' }); }
            if pi < pair_query.len() { out.push(pair_query[pi]); pi += 1; }
        }
        let tail = gap_map.get(&ref_len).copied().unwrap_or(0);
        let my_tail = ins.get(&ref_len).map(|v| v.as_slice()).unwrap_or(&[]);
        for g in 0..tail { out.push(if g < my_tail.len() { my_tail[g] } else { b'-' }); }
        String::from_utf8(out).expect("aligned query is valid UTF-8")
    };

    let mut result = vec![build_ref()];
    for (pr, pq) in &pairs { result.push(build_query(pr, pq)); }
    result
}

fn main() -> Result<()> {
    println!("Building FT promoter dataset from 1001 Genomes Plus Phase 1 with Rust...");

    let project_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let output_path = project_root.join("public").join("data").join("ft-promoter.json");

    let client = Client::builder().build().context("failed to build HTTP client")?;
    let accessions = accessions();
    let motifs = scored_motifs();
    let col0 = accessions
        .iter()
        .find(|entry| entry.accession_name == "Col-0")
        .cloned()
        .ok_or_else(|| anyhow!("Col-0 metadata missing"))?;

    let reference = build_profile(&client, &col0, None, None, &motifs)?;

    let other_accessions: Vec<AccessionMeta> = accessions
        .iter()
        .filter(|entry| entry.accession_name != "Col-0")
        .cloned()
        .collect();

    let reference_sequence = reference.promoter_sequence.clone();
    let ref_instances = &reference.profile.motif_instances;
    let other_results: Vec<(AccessionMeta, BuildOutput)> = other_accessions
        .into_par_iter()
        .map(|accession| {
            println!("Processing {} ({})", accession.accession_name, accession.accession_number);
            let local_client = Client::builder().build().context("failed to build HTTP client")?;
            let result = build_profile(&local_client, &accession, Some(&reference_sequence), Some(ref_instances), &motifs)?;
            Ok::<_, anyhow::Error>((accession, result))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut source_summary = BTreeMap::new();

    let accession_order: HashMap<u32, usize> = accessions
        .iter()
        .enumerate()
        .map(|(index, accession)| (accession.accession_number, index))
        .collect();

    let mut sorted_results = other_results;
    sorted_results.sort_by_key(|(accession, _)| accession_order[&accession.accession_number]);

    let query_seqs: Vec<String> = sorted_results
        .iter()
        .map(|(_, output)| output.promoter_sequence.clone())
        .collect();

    println!("Running star MSA ({} sequences)...", query_seqs.len() + 1);
    let aligned = build_star_msa(&reference_sequence, &query_seqs);

    let ref_chromosome = reference.gene.seq_id.clone();
    let ref_tss_coord: i64 = if reference.gene.strand == "+" { reference.gene.start as i64 } else { reference.gene.end as i64 };
    let ref_strand = reference.gene.strand.clone();
    let ref_gene_features = reference.gene_features;
    source_summary.insert(col0.accession_name.to_string(), AccessionSourceSummary { gene: reference.gene });
    let mut reference_profile = reference.profile;
    reference_profile.aligned_sequence = aligned[0].clone();
    let mut profiles = vec![reference_profile];

    for (i, (accession, output)) in sorted_results.into_iter().enumerate() {
        source_summary.insert(accession.accession_name.to_string(), AccessionSourceSummary { gene: output.gene });
        let mut profile = output.profile;
        profile.aligned_sequence = aligned[i + 1].clone();
        profiles.push(profile);
    }

    let dataset = Dataset {
        dataset_label: "1001 Genomes Plus Phase 1 FT promoter sequence scan",
        gene: "FT / AT1G65480",
        promoter_window: PromoterWindow {
            start: PROMOTER_START,
            end: PROMOTER_END,
            tss: 0,
            chromosome: ref_chromosome,
            tss_coord: ref_tss_coord,
            strand: ref_strand,
            gene_features: ref_gene_features,
        },
        accessions,
        motifs: motifs
            .iter()
            .map(|motif| MotifOutput {
                id: motif.definition.id.to_string(),
                label: motif.definition.label.to_string(),
                family: motif.definition.family.to_string(),
                element: motif.definition.element.to_string(),
                color: motif.definition.color.to_string(),
                description: motif.definition.description.to_string(),
                source: format!("JASPAR {}", motif.definition.matrix_id),
                pfm: PfmOutput {
                    a: motif.definition.pfm[0].clone(),
                    c: motif.definition.pfm[1].clone(),
                    g: motif.definition.pfm[2].clone(),
                    t: motif.definition.pfm[3].clone(),
                },
            })
            .collect(),
        profiles,
        citations: vec![
            "1001 Genomes Plus Phase 1 assemblies: https://1001genomes.org/data/1001Gp/27genomes/releases/current/assemblies/",
            "1001 Genomes Plus Phase 1 annotations: https://1001genomes.org/data/1001Gp/27genomes/releases/current/annotations/",
            "FT locus identified from gene annotations using Name=AT1G65480.",
            "Motif scores are promoter-level PWM scores derived from JASPAR matrices, not experimental binding measurements.",
        ],
        is_demo: false,
        provenance: Provenance {
            built_at: iso_timestamp(),
            source_summary,
        },
    };

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).context("failed to create output directory")?;
    }
    fs::write(&output_path, serde_json::to_string_pretty(&dataset)?).context("failed to write output dataset")?;
    println!("Wrote {}", output_path.display());
    Ok(())
}

fn compute_instance_delta(
    hit: &MotifInstance,
    ref_instances: &[MotifInstance],
    sequence: &str,
    ref_sequence: &str,
) -> String {
    const TOLERANCE: i32 = 50;
    let best = ref_instances.iter()
        .filter(|rh| rh.motif_id == hit.motif_id)
        .filter(|rh| (rh.start - hit.start).abs() <= TOLERANCE)
        .min_by_key(|rh| (rh.start - hit.start).abs());
    match best {
        None => "gained".to_string(),
        Some(ref_hit) => {
            let si = (hit.start - PROMOTER_START) as usize;
            let ri = (ref_hit.start - PROMOTER_START) as usize;
            let len = hit.length.min(ref_hit.length);
            if si + len <= sequence.len() && ri + len <= ref_sequence.len() {
                if sequence[si..si + len] == ref_sequence[ri..ri + len] {
                    "conserved".to_string()
                } else {
                    "variant".to_string()
                }
            } else {
                "conserved".to_string()
            }
        }
    }
}

fn build_profile(client: &Client, accession: &AccessionMeta, reference_sequence: Option<&str>, ref_instances: Option<&[MotifInstance]>, motifs: &[ScoredMotif]) -> Result<BuildOutput> {
    let accession_code = accession.accession_number.to_string();
    let file_stem = file_stem(accession.accession_number);
    let gff_url = format!("{BASE_URL}/annotations/genes_v05_{file_stem}.gff.gz");
    let fasta_url = format!("{BASE_URL}/assemblies/{file_stem}.scaffolds_corrected.v2.1.fasta.gz");

    let gff_text = fetch_gunzip_text(client, &gff_url)?;
    let fasta_text = fetch_gunzip_text(client, &fasta_url)?;
    let gene = find_ft_gene(&gff_text).with_context(|| format!("failed to find FT gene in accession {accession_code}"))?;
    let chromosome_sequence = extract_sequence(&fasta_text, &gene.seq_id)
        .with_context(|| format!("failed to extract chromosome {} for accession {}", gene.seq_id, accession_code))?;
    let promoter_sequence = slice_promoter(&chromosome_sequence, &gene);
    let comparison = compare_to_reference(&promoter_sequence, reference_sequence.unwrap_or(&promoter_sequence));

    let mut motif_scores = BTreeMap::new();
    let mut motif_instances = Vec::new();
    for motif in motifs {
        let scan = scan_motif(&promoter_sequence, motif);
        motif_scores.insert(motif.definition.id.to_string(), scan.promoter_score);
        motif_instances.extend(scan.hits.into_iter().map(|mut hit| {
            hit.delta = match ref_instances {
                None => "conserved".to_string(),
                Some(refs) => compute_instance_delta(
                    &hit, refs, &promoter_sequence,
                    reference_sequence.unwrap_or(&promoter_sequence),
                ),
            };
            hit
        }));
    }

    let mut profile = AccessionProfile {
        accession_name: accession.accession_name.to_string(),
        accession_number: accession.accession_number,
        promoter_length: promoter_sequence.len(),
        sequence_identity: comparison.0,
        variant_count: comparison.1,
        motif_scores,
        motif_instances,
        structural_variants: Vec::new(),
        notes: Vec::new(),
        sequence: promoter_sequence.clone(),
        aligned_sequence: String::new(), // filled in after star MSA
    };
    profile.notes = summarise_notes(&profile, motifs);

    let gene_features = find_ft_gene_features(&gff_text, &gene);

    Ok(BuildOutput {
        profile,
        promoter_sequence,
        gene: GeneSummary {
            seq_id: gene.seq_id,
            start: gene.start,
            end: gene.end,
            strand: gene.strand.to_string(),
        },
        gene_features,
    })
}

fn fetch_gunzip_text(client: &Client, url: &str) -> Result<String> {
    let mut response = client.get(url).send().with_context(|| format!("request failed for {url}"))?;
    if !response.status().is_success() {
        bail!("request failed for {url}: {}", response.status());
    }
    let mut bytes = Vec::new();
    response.read_to_end(&mut bytes).context("failed to read response")?;
    let mut decoder = MultiGzDecoder::new(&bytes[..]);
    let mut text = String::new();
    decoder.read_to_string(&mut text).context("failed to gunzip response")?;
    Ok(text)
}

fn parse_attributes(field: &str) -> HashMap<&str, &str> {
    field
        .split(';')
        .filter(|entry| !entry.is_empty())
        .filter_map(|entry| entry.split_once('='))
        .collect()
}

fn find_ft_gene(gff_text: &str) -> Result<GeneRecord> {
    for line in gff_text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 || parts[2] != "gene" {
            continue;
        }
        let attributes = parse_attributes(parts[8]);
        let names = attributes.get("Name").copied().unwrap_or_default();
        if names.split(',').any(|name| name == FT_GENE_NAME) {
            return Ok(GeneRecord {
                seq_id: parts[0].to_string(),
                start: parts[3].parse()?,
                end: parts[4].parse()?,
                strand: if parts[6].starts_with('-') { '-' } else { '+' },
            });
        }
    }
    bail!("could not find {FT_GENE_NAME} in annotation")
}

fn find_ft_gene_features(gff_text: &str, gene: &GeneRecord) -> Vec<GeneFeature> {
    // Match features by chromosome + position overlap with the gene body, not by Parent name,
    // because these GFFs use accession-prefixed IDs (e.g. "6909_Chr1_gene_...").
    let tss = if gene.strand == '+' { gene.start as i32 } else { gene.end as i32 };
    let strand_char = gene.strand.to_string();
    let mut features = Vec::new();
    for line in gff_text.lines() {
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 { continue; }
        if parts[0] != gene.seq_id { continue; }
        if parts[6] != strand_char { continue; }
        let ftype = match parts[2] {
            "exon" => "exon",
            "five_prime_UTR" | "five_prime_utr" => "5utr",
            "three_prime_UTR" | "three_prime_utr" => "3utr",
            _ => continue,
        };
        let (Ok(fs), Ok(fe)): (Result<i32, _>, Result<i32, _>) = (parts[3].parse(), parts[4].parse()) else { continue };
        // Must lie within the gene boundaries (allow 1 bp tolerance for GFF conventions)
        if fs < gene.start as i32 - 1 || fe > gene.end as i32 + 1 { continue; }
        let (rel_start, rel_end) = if gene.strand == '+' {
            (fs - tss, fe - tss)
        } else {
            (tss - fe, tss - fs)
        };
        features.push(GeneFeature { feature_type: ftype.to_string(), start: rel_start, end: rel_end });
    }
    features.sort_by_key(|f| f.start);
    features
}

fn extract_sequence(fasta_text: &str, seq_id: &str) -> Result<String> {
    let mut active = false;
    let mut chunks = Vec::new();
    for line in fasta_text.lines() {
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            let current_id = rest.split_whitespace().next().unwrap_or_default();
            active = current_id == seq_id;
            continue;
        }
        if active {
            chunks.push(line.trim().to_uppercase());
        }
    }
    if chunks.is_empty() {
        bail!("could not find sequence {seq_id} in FASTA")
    }
    Ok(chunks.join(""))
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn slice_promoter(sequence: &str, gene: &GeneRecord) -> String {
    if gene.strand == '+' {
        let region_start = gene.start as i32 + PROMOTER_START;
        let region_end = gene.start as i32 + PROMOTER_END;
        slice_padded(sequence, region_start, region_end)
    } else {
        let region_start = gene.end as i32 - PROMOTER_END;
        let region_end = gene.end as i32 - PROMOTER_START;
        reverse_complement(&slice_padded(sequence, region_start, region_end))
    }
}

fn slice_padded(sequence: &str, region_start: i32, region_end: i32) -> String {
    let left_pad = (1 - region_start).max(0) as usize;
    let right_pad = (region_end - sequence.len() as i32).max(0) as usize;
    let slice_start = region_start.max(1) as usize;
    let slice_end = region_end.min(sequence.len() as i32) as usize;
    let core = if slice_start <= slice_end {
        sequence.get(slice_start - 1..slice_end).unwrap_or_default().to_string()
    } else {
        String::new()
    };
    let mut padded = String::with_capacity(PROMOTER_LENGTH);
    padded.push_str(&"N".repeat(left_pad));
    padded.push_str(&core);
    padded.push_str(&"N".repeat(right_pad));
    padded.truncate(PROMOTER_LENGTH);
    padded
}

fn score_window(window: &[u8], pwm: &[Vec<f64>; 4]) -> f64 {
    let mut score = 0.0;
    for (index, base) in window.iter().enumerate() {
        let base_index = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return f64::NEG_INFINITY,
        };
        score += pwm[base_index][index];
    }
    score
}

fn normalize_score(raw_score: f64, motif: &ScoredMotif) -> f64 {
    if (motif.max_score - motif.min_score).abs() < f64::EPSILON {
        return 0.0;
    }
    ((raw_score - motif.min_score) / (motif.max_score - motif.min_score)).clamp(0.0, 1.0)
}

fn scan_motif(sequence: &str, motif: &ScoredMotif) -> ScanResult {
    let bytes = sequence.as_bytes();
    let mut hits = Vec::new();
    let mut window_scores = Vec::new();

    if bytes.len() < motif.length {
        return ScanResult {
            promoter_score: 0.0,
            hits,
        };
    }

    for index in 0..=bytes.len() - motif.length {
        let window = &bytes[index..index + motif.length];
        if window.contains(&b'N') {
            continue;
        }
        let plus = normalize_score(score_window(window, &motif.pwm), motif);
        let minus = normalize_score(score_window(window, &motif.reverse_pwm), motif);
        let normalized = plus.max(minus);
        let strand = if minus > plus { "-" } else { "+" };
        window_scores.push(normalized);

        if normalized >= motif.definition.strong_hit_threshold {
            hits.push(MotifInstance {
                motif_id: motif.definition.id.to_string(),
                start: PROMOTER_START + index as i32,
                length: motif.length,
                strand: strand.to_string(),
                score: round(normalized, 3),
                delta: String::new(),
            });
        }
    }

    hits.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(Ordering::Equal)
            .then(left.start.cmp(&right.start))
    });
    window_scores.sort_by(|left, right| right.partial_cmp(left).unwrap_or(Ordering::Equal));
    let top_scores = &window_scores[..window_scores.len().min(5)];
    let promoter_score = if top_scores.is_empty() {
        0.0
    } else {
        round(top_scores.iter().sum::<f64>() / top_scores.len() as f64, 3)
    };

    ScanResult {
        promoter_score,
        hits,
    }
}

fn compare_to_reference(sequence: &str, reference: &str) -> (f64, usize) {
    let mut informative = 0usize;
    let mut matches = 0usize;
    let mut differences = 0usize;

    for (left, right) in sequence.bytes().zip(reference.bytes()) {
        if left == b'N' || right == b'N' {
            continue;
        }
        informative += 1;
        if left == right {
            matches += 1;
        } else {
            differences += 1;
        }
    }

    let identity = if informative == 0 {
        0.0
    } else {
        round(matches as f64 / informative as f64 * 100.0, 2)
    };
    (identity, differences)
}

fn summarise_notes(profile: &AccessionProfile, motifs: &[ScoredMotif]) -> Vec<String> {
    let mut notes = Vec::new();
    if profile.motif_scores.get("nfyc2").copied().unwrap_or(0.0) > profile.motif_scores.get("agl15").copied().unwrap_or(0.0) {
        notes.push("NF-Y-like promoter signal exceeds the AGL15-like MADS-box signal in the extracted FT promoter window.".to_string());
    }
    if profile.motif_scores.get("pif4").copied().unwrap_or(0.0) >= 0.9 {
        notes.push("Contains a strong PIF4-like PWM hit in the scanned promoter window.".to_string());
    }
    if profile.variant_count > 120 {
        notes.push("Promoter sequence is relatively diverged from Col-0 across the fixed -5000..+5000 window.".to_string());
    }
    if let Some((motif_id, _)) = profile
        .motif_scores
        .iter()
        .max_by(|left, right| left.1.partial_cmp(right.1).unwrap_or(Ordering::Equal))
    {
        let label = motifs
            .iter()
            .find(|motif| motif.definition.id == motif_id)
            .map(|motif| motif.definition.label)
            .unwrap_or(motif_id.as_str());
        notes.push(format!("Best motif class by normalized match score: {label}."));
    }
    notes
}

fn round(value: f64, decimals: i32) -> f64 {
    let factor = 10f64.powi(decimals);
    (value * factor).round() / factor
}

fn iso_timestamp() -> String {
    let output = std::process::Command::new("powershell")
        .args(["-NoProfile", "-Command", "[DateTime]::UtcNow.ToString('o').Replace('+00:00','Z')"])
        .output();
    if let Ok(output) = output {
        if output.status.success() {
            return String::from_utf8_lossy(&output.stdout).trim().to_string();
        }
    }
    "1970-01-01T00:00:00Z".to_string()
}

fn file_stem(accession_number: u32) -> String {
    match accession_number {
        22001 => "22001f".to_string(),
        value => value.to_string(),
    }
}

fn scored_motifs() -> Vec<ScoredMotif> {
    motifs()
        .into_iter()
        .map(|definition| {
            let length = definition.pfm[0].len();
            let mut pwm = [vec![0.0; length], vec![0.0; length], vec![0.0; length], vec![0.0; length]];
            for index in 0..length {
                let total: f64 = (0..4).map(|base| definition.pfm[base][index]).sum();
                for base in 0..4 {
                    let probability = (definition.pfm[base][index] + 0.1) / (total + 0.4);
                    pwm[base][index] = (probability / BACKGROUND_PROBABILITY).log2();
                }
            }
            let max_score: f64 = (0..length)
                .map(|index| (0..4).map(|base| pwm[base][index]).fold(f64::NEG_INFINITY, f64::max))
                .sum();
            let min_score: f64 = (0..length)
                .map(|index| (0..4).map(|base| pwm[base][index]).fold(f64::INFINITY, f64::min))
                .sum();
            let reverse_pwm = [
                pwm[3].iter().rev().copied().collect(),
                pwm[2].iter().rev().copied().collect(),
                pwm[1].iter().rev().copied().collect(),
                pwm[0].iter().rev().copied().collect(),
            ];
            ScoredMotif {
                definition,
                pwm,
                reverse_pwm,
                length,
                max_score,
                min_score,
            }
        })
        .collect()
}

fn motifs() -> Vec<MotifDefinition> {
    vec![
        MotifDefinition {
            id: "nfyc2",
            matrix_id: "UN0690.2",
            label: "NFYC2",
            family: "NF-Y / CO complex",
            element: "CCAAT-box",
            color: "#207178",
            strong_hit_threshold: 0.86,
            description: "JASPAR UN0690.2 Arabidopsis NFYC2 (CCAAT-box / CO-NF-Y element). The CCAAT boxes bound by NF-YB/NF-YC are the primary CO-responsive elements in the FT promoter.",
            pfm: [
                vec![3455.0, 169.0, 116.0, 79.0, 120.0, 236.0, 343.0, 3557.0],
                vec![259.0, 4264.0, 31.0, 105.0, 53.0, 103.0, 4053.0, 184.0],
                vec![443.0, 87.0, 4428.0, 23.0, 4421.0, 3438.0, 174.0, 508.0],
                vec![523.0, 160.0, 105.0, 4473.0, 86.0, 903.0, 110.0, 431.0],
            ],
        },
        MotifDefinition {
            id: "agl15",
            matrix_id: "MA0548.3",
            label: "AGL15",
            family: "MADS-box",
            element: "CArG-box",
            color: "#d95d39",
            strong_hit_threshold: 0.8,
            description: "JASPAR MA0548.3 Arabidopsis AGL15 PWM as a MADS-box proxy for FT promoter repression-like sites.",
            pfm: [
                vec![48.0, 16.0, 84.0, 0.0, 0.0, 200.0, 88.0, 162.0, 7.0, 89.0, 138.0, 125.0, 0.0, 325.0, 543.0, 475.0],
                vec![33.0, 5.0, 8.0, 596.0, 435.0, 80.0, 85.0, 0.0, 8.0, 45.0, 42.0, 1.0, 0.0, 50.0, 13.0, 7.0],
                vec![5.0, 3.0, 21.0, 0.0, 0.0, 50.0, 53.0, 4.0, 0.0, 31.0, 168.0, 473.0, 585.0, 28.0, 5.0, 29.0],
                vec![513.0, 575.0, 486.0, 3.0, 164.0, 269.0, 373.0, 433.0, 584.0, 434.0, 251.0, 0.0, 14.0, 196.0, 38.0, 88.0],
            ],
        },
        MotifDefinition {
            id: "pif4",
            matrix_id: "MA0561.1",
            label: "PIF4",
            family: "bHLH",
            element: "G-box",
            color: "#008148",
            strong_hit_threshold: 0.84,
            description: "JASPAR MA0561.1 Arabidopsis PIF4 (G-box / E-box CACGTG). bHLH; CO is recruited to the FT promoter G-box via FBH1-4 bHLH proteins with this motif.",
            pfm: [
                vec![0.0, 335.0, 0.0, 49.0, 0.0, 0.0, 0.0, 64.0],
                vec![335.0, 0.0, 335.0, 0.0, 0.0, 0.0, 99.0, 206.0],
                vec![0.0, 0.0, 0.0, 286.0, 0.0, 335.0, 183.0, 65.0],
                vec![0.0, 0.0, 0.0, 0.0, 335.0, 0.0, 53.0, 0.0],
            ],
        },
        MotifDefinition {
            id: "tcp20",
            matrix_id: "MA1065.3",
            label: "TCP20",
            family: "TCP",
            element: "TCP-box",
            color: "#fb8500",
            strong_hit_threshold: 0.84,
            description: "JASPAR MA1065.3 Arabidopsis TCP20 PWM for class-I TCP promoter signals.",
            pfm: [
                vec![106.0, 19.0, 77.0, 4.0, 37.0, 379.0, 13.0, 4.0, 9.0, 887.0, 35.0],
                vec![52.0, 53.0, 7.0, 1.0, 3.0, 158.0, 915.0, 962.0, 876.0, 8.0, 769.0],
                vec![767.0, 12.0, 877.0, 962.0, 914.0, 159.0, 4.0, 0.0, 4.0, 61.0, 60.0],
                vec![45.0, 886.0, 9.0, 3.0, 16.0, 274.0, 38.0, 4.0, 81.0, 14.0, 106.0],
            ],
        },
        MotifDefinition {
            id: "cdf2",
            matrix_id: "MA0973.2",
            label: "CDF2",
            family: "DOF",
            element: "DOF-site",
            color: "#5f0f40",
            strong_hit_threshold: 0.84,
            description: "JASPAR MA0973.2 Arabidopsis CDF2 PWM, relevant to photoperiodic FT regulation.",
            pfm: [
                vec![608.0, 693.0, 980.0, 979.0, 868.0, 17.0, 53.0],
                vec![92.0, 16.0, 10.0, 10.0, 29.0, 7.0, 295.0],
                vec![150.0, 41.0, 4.0, 2.0, 81.0, 969.0, 154.0],
                vec![150.0, 250.0, 5.0, 9.0, 22.0, 6.0, 498.0],
            ],
        },
        MotifDefinition {
            id: "cca1",
            matrix_id: "MA0972.1",
            label: "CCA1",
            family: "Circadian MYB-related",
            element: "EE/LBS",
            color: "#8a3ffc",
            strong_hit_threshold: 0.86,
            description: "JASPAR MA0972.1 Arabidopsis CCA1 PWM for circadian-associated promoter signals.",
            pfm: [
                vec![934.0, 960.0, 985.0, 35.0, 992.0, 1.0, 5.0, 21.0],
                vec![1.0, 1.0, 4.0, 0.0, 5.0, 12.0, 991.0, 91.0],
                vec![10.0, 37.0, 10.0, 5.0, 1.0, 1.0, 0.0, 14.0],
                vec![54.0, 2.0, 1.0, 960.0, 2.0, 986.0, 4.0, 874.0],
            ],
        },
        MotifDefinition {
            id: "flc",
            matrix_id: "MA0558.2",
            label: "FLC",
            family: "MADS-box",
            element: "CArG-box",
            color: "#6366f1",
            strong_hit_threshold: 0.82,
            description: "JASPAR MA0558.2 Arabidopsis FLC (FLOWERING LOCUS C) PWM — the primary vernalization-sensitive FT repressor.",
            pfm: [
                vec![54.0, 9.0, 135.0, 198.0, 261.0, 233.0, 193.0, 78.0, 141.0, 0.0, 217.0, 268.0, 257.0, 62.0],
                vec![193.0, 185.0, 70.0, 10.0, 0.0, 3.0, 3.0, 11.0, 0.0, 0.0, 21.0, 3.0, 1.0, 43.0],
                vec![12.0, 4.0, 37.0, 38.0, 3.0, 3.0, 18.0, 15.0, 134.0, 275.0, 6.0, 4.0, 15.0, 148.0],
                vec![16.0, 77.0, 33.0, 29.0, 11.0, 36.0, 61.0, 171.0, 0.0, 0.0, 31.0, 0.0, 2.0, 22.0],
            ],
        },
        MotifDefinition {
            id: "svp",
            matrix_id: "MA0555.2",
            label: "SVP",
            family: "MADS-box",
            element: "CArG-box",
            color: "#0e7490",
            strong_hit_threshold: 0.82,
            description: "JASPAR MA0555.2 Arabidopsis SVP (SHORT VEGETATIVE PHASE) PWM — MADS-box repressor acting with FLC to repress FT.",
            pfm: [
                vec![28.0, 28.0, 33.0, 7.0, 3.0, 58.0, 67.0, 86.0, 77.0, 78.0, 34.0, 16.0, 0.0, 71.0, 84.0, 77.0],
                vec![11.0, 1.0, 16.0, 85.0, 79.0, 16.0, 5.0, 0.0, 1.0, 8.0, 10.0, 0.0, 0.0, 10.0, 2.0, 7.0],
                vec![8.0, 10.0, 11.0, 0.0, 1.0, 7.0, 8.0, 0.0, 2.0, 3.0, 14.0, 74.0, 92.0, 0.0, 2.0, 2.0],
                vec![45.0, 53.0, 32.0, 0.0, 9.0, 11.0, 12.0, 6.0, 12.0, 3.0, 34.0, 2.0, 0.0, 11.0, 4.0, 6.0],
            ],
        },
        MotifDefinition {
            id: "myc3",
            matrix_id: "MA0568.2",
            label: "MYC3",
            family: "bHLH",
            element: "G-box",
            color: "#3a86ff",
            strong_hit_threshold: 0.82,
            description: "JASPAR MA0568.2 Arabidopsis MYC3 bHLH — G-box (CACGTG) binding; integrates JA signalling with flowering.",
            pfm: [
                vec![10.0, 97.0, 0.0, 16.0, 1.0, 0.0],
                vec![89.0, 0.0, 84.0, 0.0, 2.0, 0.0],
                vec![0.0, 2.0, 0.0, 84.0, 0.0, 89.0],
                vec![0.0, 1.0, 16.0, 0.0, 97.0, 10.0],
            ],
        },
        MotifDefinition {
            id: "smz",
            matrix_id: "UN0875.1",
            label: "SMZ",
            family: "AP2/ERF",
            element: "GCC-box",
            color: "#ef233c",
            strong_hit_threshold: 0.80,
            description: "JASPAR UN0875.1 Arabidopsis SMZ (SCHLAFMÜTZE) AP2/ERF — binds GCC-box; represses FT in non-inductive conditions.",
            pfm: [
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 135.0, 0.0],
                vec![105.0, 119.0, 0.0, 147.0, 0.0, 0.0, 12.0, 147.0],
                vec![0.0, 28.0, 0.0, 0.0, 147.0, 25.0, 0.0, 0.0],
                vec![42.0, 0.0, 147.0, 0.0, 0.0, 122.0, 0.0, 0.0],
            ],
        },
        MotifDefinition {
            id: "tem1",
            matrix_id: "MA1800.1",
            label: "TEM1",
            family: "RAV / B3",
            element: "RAV1-A",
            color: "#bc6c25",
            strong_hit_threshold: 0.82,
            description: "JASPAR MA1800.1 Arabidopsis TEM1 (TEMPRANILLO 1) RAV/B3 — binds RAV1-A element; delays flowering by repressing FT.",
            pfm: [
                vec![219.542, 87.8331, 874.61, 999.532, 0.156152, 999.532, 0.172295, 599.65],
                vec![97.7096, 911.816, 15.7714, 0.156152, 999.532, 0.156152, 0.172295, 125.125],
                vec![682.505, 0.175316, 109.463, 0.156152, 0.156152, 0.156152, 241.385, 0.24975],
                vec![0.243665, 0.175316, 0.156152, 0.156152, 0.156152, 0.156152, 758.27, 274.975],
            ],
        },
    ]
}

fn accessions() -> Vec<AccessionMeta> {
    vec![
        AccessionMeta { id: "100265", accession_number: 1741, accession_name: "KBS-Mac-74", stock_id: "CS78969", country_code: "USA", country: "United States", latitude: 42.405, longitude: -85.398 },
        AccessionMeta { id: "100266", accession_number: 6024, accession_name: "Fly2-2", stock_id: "CS76864", country_code: "SWE", country: "Sweden", latitude: 55.7501271, longitude: 13.3720798 },
        AccessionMeta { id: "100267", accession_number: 6069, accession_name: "Nyl-7", stock_id: "CS77137", country_code: "SWE", country: "Sweden", latitude: 62.95681, longitude: 18.2763 },
        AccessionMeta { id: "100268", accession_number: 6124, accession_name: "T690", stock_id: "CS77309", country_code: "SWE", country: "Sweden", latitude: 55.8362185, longitude: 13.2995709 },
        AccessionMeta { id: "100269", accession_number: 6244, accession_name: "TRA 01", stock_id: "CS77384", country_code: "SWE", country: "Sweden", latitude: 62.9169, longitude: 18.4728 },
        AccessionMeta { id: "100270", accession_number: 6909, accession_name: "Col-0", stock_id: "CS76778", country_code: "USA", country: "United States", latitude: 38.3, longitude: -92.3 },
        AccessionMeta { id: "100271", accession_number: 6966, accession_name: "Sq-1", stock_id: "CS77266", country_code: "UK", country: "United Kingdom", latitude: 51.4083, longitude: -0.6383 },
        AccessionMeta { id: "100272", accession_number: 8236, accession_name: "HSm", stock_id: "CS76941", country_code: "CZE", country: "Czech Republic", latitude: 49.33, longitude: 15.76 },
        AccessionMeta { id: "100273", accession_number: 9075, accession_name: "Lerik1-4", stock_id: "CS77023", country_code: "AZE", country: "Azerbaijan", latitude: 38.7406, longitude: 48.6131 },
        AccessionMeta { id: "100274", accession_number: 9537, accession_name: "IP-Cum-1", stock_id: "CS76787", country_code: "ESP", country: "Spain", latitude: 38.07, longitude: -6.66 },
        AccessionMeta { id: "100275", accession_number: 9543, accession_name: "IP-Gra-0", stock_id: "CS76886", country_code: "ESP", country: "Spain", latitude: 36.77, longitude: -5.39 },
        AccessionMeta { id: "100276", accession_number: 9638, accession_name: "Noveg-3", stock_id: "CS77133", country_code: "RUS", country: "Russia", latitude: 51.73, longitude: 80.86 },
        AccessionMeta { id: "100277", accession_number: 9728, accession_name: "Stiav-1", stock_id: "CS77279", country_code: "SVK", country: "Slovakia", latitude: 48.46, longitude: 18.9 },
        AccessionMeta { id: "100278", accession_number: 9764, accession_name: "Qar-8a", stock_id: "CS76581", country_code: "LBN", country: "Lebanon", latitude: 34.1, longitude: 35.84 },
        AccessionMeta { id: "100279", accession_number: 9888, accession_name: "IP-Pva-1", stock_id: "CS77197", country_code: "ESP", country: "Spain", latitude: 40.93, longitude: -3.31 },
        AccessionMeta { id: "100280", accession_number: 9905, accession_name: "IP-Ven-0", stock_id: "CS78840", country_code: "ESP", country: "Spain", latitude: 40.76, longitude: -4.01 },
        AccessionMeta { id: "100281", accession_number: 9981, accession_name: "Angit-1", stock_id: "CS76366", country_code: "ITA", country: "Italy", latitude: 38.76, longitude: 16.24 },
        AccessionMeta { id: "100282", accession_number: 10002, accession_name: "TueWa1-2", stock_id: "CS76405", country_code: "GER", country: "Germany", latitude: 48.53, longitude: 9.04 },
        AccessionMeta { id: "100283", accession_number: 10015, accession_name: "Shahdara", stock_id: "—", country_code: "TJK", country: "Tajikistan", latitude: 38.35, longitude: 68.48 },
        AccessionMeta { id: "100284", accession_number: 10024, accession_name: "Tnz-1", stock_id: "—", country_code: "TZA", country: "Tanzania", latitude: -2.87389, longitude: 36.2128 },
        AccessionMeta { id: "100286", accession_number: 22001, accession_name: "85-3", stock_id: "—", country_code: "CHN", country: "China", latitude: 32.14, longitude: 115.06 },
        AccessionMeta { id: "100287", accession_number: 22002, accession_name: "35-1", stock_id: "—", country_code: "CHN", country: "China", latitude: 27.94, longitude: 108.61 },
        AccessionMeta { id: "100288", accession_number: 22003, accession_name: "Taz-0", stock_id: "CS799913", country_code: "MAR", country: "Morocco", latitude: 34.09166, longitude: -4.10258 },
        AccessionMeta { id: "100289", accession_number: 22004, accession_name: "Elh-2", stock_id: "CS799925", country_code: "MAR", country: "Morocco", latitude: 31.47197, longitude: -7.40644 },
        AccessionMeta { id: "100290", accession_number: 22005, accession_name: "Rabacal-1", stock_id: "CS2107642", country_code: "POR", country: "Portugal", latitude: 32.7536, longitude: -17.1297 },
        AccessionMeta { id: "100291", accession_number: 22006, accession_name: "Areeiro-1", stock_id: "635AV", country_code: "POR", country: "Portugal", latitude: 32.74, longitude: -16.93 },
        AccessionMeta { id: "100292", accession_number: 22007, accession_name: "ET-86.4", stock_id: "—", country_code: "ETH", country: "Ethiopia", latitude: 13.235084, longitude: 38.061701 },
    ]
}
