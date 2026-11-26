use std::collections::HashSet;
use std::env;
use std::fs;
use std::process;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EmitTarget {
    Default,
    Ir,
    Python,
    Verilog,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum UncertaintyMode {
    Report,
    Error,
    Omit,
}

#[derive(Clone, Debug)]
struct Config {
    input_path: String,
    circular: bool,
    non_strict: bool,
    emit: EmitTarget,
    uncertainty: UncertaintyMode,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SwitchKind {
    Promoter(String),
    OperatorSiteCandidate(String),
    Rbs,
    LeaderlessTss,
    Coupled(usize),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SwitchEntry {
    kind: SwitchKind,
    count: usize,
    positions: Vec<usize>,
}

#[derive(Clone, Debug)]
struct CandidateGene {
    start: usize,
    end: usize,
    strand: char,
    prom_pos: usize,
    prom_end: usize,
    terminator_end: Option<usize>,
    prom_type: String,
    init_mode: String,
    promoter_matches: usize,
    rbs_matches: Option<usize>,
    terminator_dg: Option<f64>,
    unterminated: bool,
    uncertain_evidence: bool,
    uncertainty_reasons: Vec<&'static str>,
    coupled_to: Option<usize>,
}

#[derive(Clone, Debug)]
struct NoiseFloors {
    sigma70_prom_floor: Option<usize>,
    sigma32_prom_floor: Option<usize>,
    rbs_floor: Option<usize>,
    term_floor: Option<f64>,
}

const SHUFFLE_COUNT: usize = 10;
const SHUFFLE_SEED: u64 = 42;
// Empirical null quantiles, not tuned sequence-score weights.
const NULL_QUANTILE_PERCENT: usize = 99;
const MIN_ORF_LEN: usize = 30;
const MIN_COUPLED_ORF_CODONS: usize = 3;
const MIN_STEM: usize = 6;
const MAX_STEM: usize = 12;
const MIN_LOOP: usize = 3;
const MAX_LOOP: usize = 8;
const TAIL_WIN: usize = 8;
const MIN_TAIL_THYMINE_COUNT: u8 = 4;
const MIN_OPERATOR_STEM: usize = 6;
const MAX_OPERATOR_STEM: usize = 10;
const REFERENCE_TEMPERATURE_K: f64 = 310.15;
// SantaLucia & Hicks (2004), Table 4: length-dependent DNA hairpin-loop entropy
// for 3..=8 unpaired bases; ΔH = 0. DOI: 10.1146/annurev.biophys.32.110601.141800.
const HAIRPIN_LOOP_ENTROPY_CAL_PER_MOL_K: [f64; 6] =
    [-11.3, -11.3, -10.6, -12.9, -13.5, -13.9];

struct BasePlanes {
    high: Vec<u16>,
    low: Vec<u16>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RegulatoryEdge {
    regulator: usize,
    target: usize,
    motif: String,
    count: usize,
    sign: char,
}

fn exit_err(msg: &str) -> ! {
    eprintln!("Error: {}", msg);
    process::exit(1);
}

fn parse_args() -> Config {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        eprintln!("Usage: dnadec <input.fasta> [--circular] [--non-strict] [--uncertainty <report|error|omit>] [--emit <default|ir|python|verilog>]");
        process::exit(1);
    }

    let mut input_path = None;
    let mut circular = false;
    let mut non_strict = false;
    let mut emit = EmitTarget::Default;
    let mut uncertainty = UncertaintyMode::Report;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--circular" => circular = true,
            "--non-strict" => non_strict = true,
            "--uncertainty" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    exit_err("--uncertainty requires a mode (report, error, omit)");
                };
                uncertainty = match val.as_str() {
                    "report" => UncertaintyMode::Report,
                    "error" => UncertaintyMode::Error,
                    "omit" => UncertaintyMode::Omit,
                    other => exit_err(&format!("unknown --uncertainty mode '{}'", other)),
                };
            }
            "--emit" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    exit_err("--emit requires a target (default, ir, python, verilog)");
                };
                emit = match val.as_str() {
                    "default" => EmitTarget::Default,
                    "ir" => EmitTarget::Ir,
                    "python" => EmitTarget::Python,
                    "verilog" => EmitTarget::Verilog,
                    other => exit_err(&format!("unknown --emit target '{}'", other)),
                };
            }
            flag if flag.starts_with('-') => exit_err(&format!("unknown flag '{}'", flag)),
            path => input_path = Some(path.to_string()),
        }
        i += 1;
    }

    let Some(input_path) = input_path else {
        exit_err("missing input .fasta file");
    };

    Config {
        input_path,
        circular,
        non_strict,
        emit,
        uncertainty,
    }
}

fn resolve_shared_motif_edges(
    genes: &[CandidateGene],
    row_ptr: &[usize],
    switches: &[SwitchEntry],
) -> Vec<RegulatoryEdge> {
    let mut edges = Vec::new();
    for left in 0..genes.len().saturating_sub(1) {
        let right = left + 1;
        if genes[left].coupled_to.is_some()
            || genes[right].coupled_to.is_some()
            || genes[left].uncertain_evidence
            || genes[right].uncertain_evidence
            || genes[left].strand != genes[right].strand
        {
            continue;
        }
        let left_sites = &switches[row_ptr[left]..row_ptr[left + 1]];
        let right_sites = &switches[row_ptr[right]..row_ptr[right + 1]];
        for left_switch in left_sites {
            let SwitchKind::OperatorSiteCandidate(motif) = &left_switch.kind else {
                continue;
            };
            for right_switch in right_sites {
                let SwitchKind::OperatorSiteCandidate(right_motif) = &right_switch.kind else {
                    continue;
                };
                if motif != right_motif {
                    continue;
                }
                for (regulator, target, switch, linked_motif) in [
                    (left, right, right_switch, right_motif),
                    (right, left, left_switch, motif),
                ] {
                    edges.push(RegulatoryEdge {
                        regulator,
                        target,
                        motif: linked_motif.clone(),
                        count: switch.count,
                        sign: '?',
                    });
                }
            }
        }
    }
    edges.sort_by(|left, right| {
        (left.regulator, left.target, &left.motif).cmp(&(
            right.regulator,
            right.target,
            &right.motif,
        ))
    });
    edges.dedup();
    edges
}

fn read_fasta(path: &str) -> Vec<u8> {
    let content = fs::read_to_string(path)
        .unwrap_or_else(|e| exit_err(&format!("failed to read '{}': {}", path, e)));

    let mut seq = Vec::with_capacity(content.len());
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('>') {
            continue;
        }
        for b in trimmed.bytes() {
            let up = b.to_ascii_uppercase();
            if matches!(up, b'A' | b'C' | b'G' | b'T') {
                seq.push(up);
            }
        }
    }

    if seq.is_empty() {
        exit_err(&format!("no valid A/C/G/T nucleotides found in '{}'", path));
    }
    seq
}

#[inline]
fn comp_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp_base(b)).collect()
}

fn encode_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&base| base_idx(base) as u8).collect()
}

fn pack_sixmers(codes: &[u8]) -> Vec<u16> {
    if codes.len() < 6 {
        return Vec::new();
    }
    let mut packed = vec![0_u16; codes.len() - 5];
    let mut rolling = 0_u16;
    for (index, &code) in codes.iter().enumerate() {
        rolling = ((rolling << 2) | code as u16) & 0x0fff;
        if index >= 5 {
            packed[index - 5] = rolling;
        }
    }
    packed
}

const fn pack_consensus(sequence: [u8; 6]) -> u16 {
    let mut packed = 0_u16;
    let mut index = 0;
    while index < sequence.len() {
        let code = match sequence[index] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
        packed = (packed << 2) | code;
        index += 1;
    }
    packed
}

#[inline]
fn sixmer_mismatches(sequence: u16, consensus: u16) -> usize {
    let differences = sequence ^ consensus;
    ((differences | (differences >> 1)) & 0x0555).count_ones() as usize
}

fn reverse_complement_codes(codes: &[u8]) -> Vec<u8> {
    codes.iter().rev().map(|code| code ^ 3).collect()
}

fn pack_base_planes(codes: &[u8]) -> BasePlanes {
    let mut high = vec![0_u16; codes.len()];
    let mut low = vec![0_u16; codes.len()];
    if codes.len() < MAX_STEM {
        return BasePlanes { high, low };
    }

    let mut high_word = 0_u16;
    let mut low_word = 0_u16;
    for (index, &code) in codes.iter().take(MAX_STEM).enumerate() {
        high_word |= u16::from(code >> 1) << index;
        low_word |= u16::from(code & 1) << index;
    }
    for start in 0..=codes.len() - MAX_STEM {
        high[start] = high_word;
        low[start] = low_word;
        if start + MAX_STEM < codes.len() {
            high_word = (high_word >> 1)
                | (u16::from(codes[start + MAX_STEM] >> 1) << (MAX_STEM - 1));
            low_word = (low_word >> 1)
                | (u16::from(codes[start + MAX_STEM] & 1) << (MAX_STEM - 1));
        }
    }
    BasePlanes { high, low }
}

#[inline]
fn base_idx(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

// SantaLucia (1998), unified DNA nearest-neighbor ΔG°37 values at 1 M NaCl,
// indexed by the 5' dinucleotide (AA through TT). DOI: 10.1073/pnas.95.4.1460.
const NN_DG: [f64; 16] = [
    -1.00, -1.44, -1.28, -0.88, -1.45, -1.84, -2.17, -1.28, -1.30, -2.24, -1.84, -1.44, -0.58,
    -1.30, -1.45, -1.00,
];

fn duplex_dg(codes: &[u8]) -> f64 {
    codes
        .windows(2)
        .map(|pair| NN_DG[pair[0] as usize * 4 + pair[1] as usize])
        .sum()
}

fn dinucleotide_shuffle(seq: &[u8], mut seed: u64) -> Option<Vec<u8>> {
    let n = seq.len();
    if n < 2 {
        return None;
    }

    let mut out_degree = [0_usize; 4];
    for pair in seq.windows(2) {
        out_degree[base_idx(pair[0])] += 1;
    }
    let mut adjacency: [Vec<u8>; 4] =
        std::array::from_fn(|base| Vec::with_capacity(out_degree[base]));
    for pair in seq.windows(2) {
        adjacency[base_idx(pair[0])].push(pair[1]);
    }
    for edges in &mut adjacency {
        for i in (1..edges.len()).rev() {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            edges.swap(i, ((seed >> 33) as usize) % (i + 1));
        }
    }

    let mut path = Vec::with_capacity(n);
    let mut stack = Vec::with_capacity(n);
    stack.push(seq[0]);
    while let Some(&base) = stack.last() {
        if let Some(next) = adjacency[base_idx(base)].pop() {
            stack.push(next);
        } else if let Some(base) = stack.pop() {
            path.push(base);
        }
    }
    path.reverse();
    if path.len() == n && path.first() == Some(&seq[0]) && path.last() == seq.last() {
        Some(path)
    } else {
        None
    }
}

fn score_promoter_box(
    sixmers: &[u16],
    i: usize,
    box35: u16,
    box10: u16,
    spacers: std::ops::RangeInclusive<usize>,
    label: &'static str,
) -> Option<(usize, usize, &'static str)> {
    if i >= sixmers.len() {
        return None;
    }
    let matches35 = 6 - sixmer_mismatches(sixmers[i], box35);

    let mut best = None;
    for spacer in spacers {
        let pos10 = i + 6 + spacer;
        if pos10 >= sixmers.len() {
            continue;
        }
        let end10 = pos10 + 6;
        let matches10 = 6 - sixmer_mismatches(sixmers[pos10], box10);
        let matches = matches35 + matches10;
        if best.map_or(true, |(best_matches, _, _)| matches > best_matches) {
            best = Some((matches, end10, label));
        }
    }
    best
}

fn score_promoters_at(sixmers: &[u16], i: usize) -> [Option<(usize, usize, &'static str)>; 2] {
    let s70 = score_promoter_box(
        sixmers,
        i,
        pack_consensus(*b"TTGACA"),
        pack_consensus(*b"TATAAT"),
        16..=18,
        "promoter_sigma70",
    );
    let s32 = score_promoter_box(
        sixmers,
        i,
        pack_consensus(*b"CTTGAA"),
        pack_consensus(*b"CCCCAT"),
        13..=17,
        "promoter_sigma32",
    );
    [s70, s32]
}

fn score_rbs_before(sixmers: &[u16], start_codon: usize) -> Option<usize> {
    let sd = pack_consensus(*b"AGGAGG");
    let mut best: Option<usize> = None;
    for offset in 4..=15 {
        if start_codon < offset + 6 {
            continue;
        }
        let pos = start_codon - offset - 6;
        if pos >= sixmers.len() {
            continue;
        }
        let matches = 6 - sixmer_mismatches(sixmers[pos], sd);
        if best.map_or(true, |best_matches| matches > best_matches) {
            best = Some(matches);
        }
    }
    best
}

fn thymine_tail_counts(seq: &[u8]) -> Vec<u8> {
    let mut counts = vec![0_u8; seq.len() + 1];
    if seq.len() < TAIL_WIN {
        return counts;
    }
    let mut count = seq[..TAIL_WIN]
        .iter()
        .filter(|&&base| base == b'T')
        .count() as u8;
    counts[0] = count;
    for start in 1..=seq.len() - TAIL_WIN {
        count -= u8::from(seq[start - 1] == b'T');
        count += u8::from(seq[start + TAIL_WIN - 1] == b'T');
        counts[start] = count;
    }
    counts
}

fn terminator_start_candidates(tail_counts: &[u8]) -> Vec<u8> {
    let sequence_len = tail_counts.len().saturating_sub(1);
    let max_tail_start = sequence_len.saturating_sub(TAIL_WIN);
    let min_tail_offset = 2 * MIN_STEM + MIN_LOOP;
    let max_tail_offset = 2 * MAX_STEM + MAX_LOOP;
    let initial_end = max_tail_offset.min(max_tail_start);
    let mut supported_tails = 0;
    if initial_end >= min_tail_offset {
        supported_tails = (min_tail_offset..=initial_end)
            .filter(|&start| tail_counts[start] >= MIN_TAIL_THYMINE_COUNT)
            .count();
    }

    let mut candidates = vec![0; sequence_len];
    for start in 0..sequence_len {
        candidates[start] = u8::from(supported_tails > 0);

        let leaving = start + min_tail_offset;
        if leaving <= max_tail_start && tail_counts[leaving] >= MIN_TAIL_THYMINE_COUNT {
            supported_tails -= 1;
        }
        let entering = start + max_tail_offset + 1;
        if entering <= max_tail_start && tail_counts[entering] >= MIN_TAIL_THYMINE_COUNT {
            supported_tails += 1;
        }
    }
    candidates
}

#[cfg(test)]
fn count_stem_pairs(codes: &[u8], left_start: usize, right_start: usize, stem: usize) -> usize {
    let mut paired = 0;
    for k in 0..stem {
        let left = codes[left_start + k];
        let right = codes[right_start + stem - 1 - k];
        if left == (right ^ 3) {
            paired += 1;
        }
    }
    paired
}

fn count_packed_stem_pairs(
    left: &BasePlanes,
    right_reverse_complement: &BasePlanes,
    left_start: usize,
    right_reverse_start: usize,
    stem: usize,
) -> usize {
    let mask = (1_u16 << stem) - 1;
    let left_high = left.high[left_start] & mask;
    let left_low = left.low[left_start] & mask;
    let right_high = right_reverse_complement.high[right_reverse_start] & mask;
    let right_low = right_reverse_complement.low[right_reverse_start] & mask;
    let differences = (left_high ^ right_high) | (left_low ^ right_low);
    stem - (differences & mask).count_ones() as usize
}

fn hairpin_loop_initiation_dg(loop_len: usize) -> Option<f64> {
    let entropy = *HAIRPIN_LOOP_ENTROPY_CAL_PER_MOL_K.get(loop_len.checked_sub(MIN_LOOP)?)?;
    Some(-REFERENCE_TEMPERATURE_K * entropy / 1000.0)
}

fn score_hairpin_at(
    seq: &[u8],
    codes: &[u8],
    base_planes: &BasePlanes,
    reverse_planes: &BasePlanes,
    tail_counts: Option<&[u8]>,
    i: usize,
) -> Option<(f64, usize)> {
    let mut best = None;
    for stem in MIN_STEM..=MAX_STEM {
        for loop_len in MIN_LOOP..=MAX_LOOP {
            let end2 = i + 2 * stem + loop_len;
            if end2 > seq.len()
                || tail_counts.is_some() && end2 + TAIL_WIN > seq.len()
            {
                continue;
            }
            if tail_counts.map_or(false, |tails| {
                tails[end2] < MIN_TAIL_THYMINE_COUNT
            }) {
                continue;
            }
            let right_start = i + stem + loop_len;
            let reverse_start = seq.len() - right_start - stem;
            if count_packed_stem_pairs(base_planes, reverse_planes, i, reverse_start, stem) != stem
            {
                continue;
            }
            let Some(loop_dg) = hairpin_loop_initiation_dg(loop_len) else {
                continue;
            };
            let dg = duplex_dg(&codes[i..i + stem]) + loop_dg;
            if best.map_or(true, |(best_dg, _)| dg < best_dg) {
                best = Some((dg, end2));
            }
        }
    }
    best
}

fn score_terminator_at(
    seq: &[u8],
    codes: &[u8],
    base_planes: &BasePlanes,
    reverse_planes: &BasePlanes,
    tail_counts: &[u8],
    i: usize,
) -> Option<(f64, usize)> {
    let (dg, tail_start) = score_hairpin_at(
        seq,
        codes,
        base_planes,
        reverse_planes,
        Some(tail_counts),
        i,
    )?;
    Some((dg, tail_start + TAIL_WIN))
}

#[inline]
fn is_start_codon(c: &[u8]) -> bool {
    matches!(c, b"ATG" | b"GTG" | b"TTG")
}

#[inline]
fn is_stop_codon(c: &[u8]) -> bool {
    matches!(c, b"TAA" | b"TAG" | b"TGA")
}

fn find_stop_in_frame(seq: &[u8], start: usize, limit: usize) -> Option<usize> {
    let mut pos = start + 3;
    while pos + 3 <= seq.len() && pos < limit {
        if is_stop_codon(&seq[pos..pos + 3]) {
            return Some(pos + 3);
        }
        pos += 3;
    }
    None
}

fn empirical_quantile_usize(values: &mut [usize], percent: usize) -> Option<usize> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let rank = (values.len() * percent).div_ceil(100).saturating_sub(1);
    values.get(rank).copied()
}

fn empirical_quantile_f64(values: &mut [f64], percent: usize) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let rank = (values.len() * percent).div_ceil(100).saturating_sub(1);
    values.get(rank).copied()
}

fn operator_motif_at(seq: &[u8], codes: &[u8], pos: usize) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    for stem in (MIN_OPERATOR_STEM..=MAX_OPERATOR_STEM).rev() {
        let motif_len = 2 * stem;
        if pos + motif_len > seq.len() {
            continue;
        }
        let paired = (0..stem)
            .all(|k| seq[pos + k] == comp_base(seq[pos + motif_len - 1 - k]));
        if !paired {
            continue;
        }
        let dg = duplex_dg(&codes[pos..pos + stem]);
        if best.map_or(true, |(best_len, best_dg)| {
            dg < best_dg || (dg == best_dg && motif_len > best_len)
        }) {
            best = Some((motif_len, dg));
        }
    }
    best
}

fn collect_null_scores(
    seq: &[u8],
    sigma70: &mut Vec<usize>,
    sigma32: &mut Vec<usize>,
    rbs: &mut Vec<usize>,
    terminator: &mut Vec<f64>,
) {
    let codes = encode_bases(seq);
    let sixmers = pack_sixmers(&codes);
    let reverse_codes = reverse_complement_codes(&codes);
    let base_planes = pack_base_planes(&codes);
    let reverse_planes = pack_base_planes(&reverse_codes);

    for pos in 0..seq.len() {
        for (matches, _, label) in score_promoters_at(&sixmers, pos).into_iter().flatten() {
            let scores = if label == "promoter_sigma70" {
                &mut *sigma70
            } else {
                &mut *sigma32
            };
            scores.push(matches);
        }
        if pos + 3 <= seq.len() && is_start_codon(&seq[pos..pos + 3]) {
            if let Some(matches) = score_rbs_before(&sixmers, pos) {
                rbs.push(matches);
            }
        }
        if let Some((dg, _)) = score_hairpin_at(
            seq,
            &codes,
            &base_planes,
            &reverse_planes,
            None,
            pos,
        ) {
            terminator.push(dg);
        }
    }
}

fn calibrate_noise_floors(seq: &[u8]) -> NoiseFloors {
    if seq.len() < 8 {
        return NoiseFloors {
            sigma70_prom_floor: None,
            sigma32_prom_floor: None,
            rbs_floor: None,
            term_floor: None,
        };
    }

    let mut sigma70 = Vec::new();
    let mut sigma32 = Vec::new();
    let mut rbs = Vec::new();
    let mut terminator = Vec::new();
    for shuffle_index in 0..SHUFFLE_COUNT {
        let Some(decoy) = dinucleotide_shuffle(seq, SHUFFLE_SEED + shuffle_index as u64) else {
            return NoiseFloors {
                sigma70_prom_floor: None,
                sigma32_prom_floor: None,
                rbs_floor: None,
                term_floor: None,
            };
        };
        collect_null_scores(
            &decoy,
            &mut sigma70,
            &mut sigma32,
            &mut rbs,
            &mut terminator,
        );
        let reverse = reverse_complement(&decoy);
        collect_null_scores(
            &reverse,
            &mut sigma70,
            &mut sigma32,
            &mut rbs,
            &mut terminator,
        );
    }

    NoiseFloors {
        sigma70_prom_floor: empirical_quantile_usize(&mut sigma70, NULL_QUANTILE_PERCENT),
        sigma32_prom_floor: empirical_quantile_usize(&mut sigma32, NULL_QUANTILE_PERCENT),
        rbs_floor: empirical_quantile_usize(&mut rbs, NULL_QUANTILE_PERCENT),
        term_floor: empirical_quantile_f64(
            &mut terminator,
            100 - NULL_QUANTILE_PERCENT,
        ),
    }
}

fn promoter_noise_floor<'a>(label: &str, floors: &'a NoiseFloors) -> Option<usize> {
    match label {
        "promoter_sigma70" => floors.sigma70_prom_floor,
        "promoter_sigma32" => floors.sigma32_prom_floor,
        _ => None,
    }
}

/// Pass 0: Scan one strand through calibrated motif, initiation, ORF, and terminator gates.
fn scan_strand(
    seq: &[u8],
    codes: &[u8],
    sixmers: &[u16],
    base_planes: &BasePlanes,
    reverse_planes: &BasePlanes,
    strand: char,
    floors: &NoiseFloors,
    non_strict: bool,
    errors: &mut Vec<String>,
) -> Vec<CandidateGene> {
    let n = seq.len();
    let tail_counts = thymine_tail_counts(seq);
    let term_start_candidates = terminator_start_candidates(&tail_counts);
    let mut terminators: Vec<(usize, usize, f64)> = Vec::new();
    let mut i = 0;
    while i < n {
        if term_start_candidates[i] != 0 {
            if let Some((dg, end)) =
                score_terminator_at(seq, codes, base_planes, reverse_planes, &tail_counts, i)
            {
                if floors.term_floor.map_or(true, |floor| dg <= floor) {
                    terminators.push((i, end, dg));
                    i = end;
                    continue;
                }
            }
        }
        i += 1;
    }
    let mut candidates = Vec::new();
    let mut active_prom: Option<(usize, usize, usize, &'static str, bool)> = None;
    let mut term_idx = 0;
    let mut pos = 0;

    while pos < n {
        while term_idx < terminators.len() && terminators[term_idx].1 <= pos {
            term_idx += 1;
        }

        let promoter_candidates = score_promoters_at(sixmers, pos);
        let mut supported_promoters = [None, None];
        let mut supported_count = 0;
        let mut fallback_promoter = None;
        for candidate in promoter_candidates.into_iter().flatten() {
            if fallback_promoter.map_or(
                true,
                |best: (usize, usize, &'static str)| candidate.0 > best.0,
            ) {
                fallback_promoter = Some(candidate);
            }
            if promoter_noise_floor(candidate.2, floors)
                .map_or(false, |floor| candidate.0 >= floor)
            {
                supported_promoters[supported_count] = Some(candidate);
                supported_count += 1;
            }
        }
        let selected = supported_promoters
            .iter()
            .flatten()
            .copied()
            .max_by_key(|candidate| candidate.0)
            .or(fallback_promoter);
        if let Some((promoter_matches, end10, label)) = selected {
            let replaces_active = active_prom
                .as_ref()
                .map_or(true, |current| promoter_matches > current.2);
            if replaces_active {
                let uncertain_promoter = promoter_noise_floor(label, floors)
                    .map_or(true, |floor| promoter_matches < floor)
                    || supported_count > 1;
                active_prom = Some((pos, end10, promoter_matches, label, uncertain_promoter));
            }
        }

        let Some((prom_pos, prom_end10, prom_matches, prom_label, uncertain_promoter)) =
            active_prom.as_ref()
        else {
            pos += 1;
            continue;
        };

        let next_term = terminators.get(term_idx).copied();
        if let Some((t_start, _, _)) = next_term {
            if pos >= t_start {
                active_prom = None;
                pos += 1;
                continue;
            }
        }

        if pos + 3 > n || !is_start_codon(&seq[pos..pos + 3]) {
            pos += 1;
            continue;
        }

        if *prom_end10 > pos {
            pos += 1;
            continue;
        }
        let rbs_score = score_rbs_before(sixmers, pos);
        let rbs_supported = rbs_score
            .zip(floors.rbs_floor)
            .map_or(false, |(score, floor)| score >= floor);
        let has_rbs = rbs_supported;
        let leaderless_ok = pos >= *prom_end10 + 5 && pos <= *prom_end10 + 9;
        let coupled_ok = candidates.last().map_or(false, |prev: &CandidateGene| {
            prev.strand == strand && prev.end > pos && (1..=4).contains(&(prev.end - pos))
        });
        let init_mode = if has_rbs {
            "RBS"
        } else if leaderless_ok {
            "leaderless_TSS+1"
        } else if coupled_ok {
            "coupled"
        } else {
            pos += 1;
            continue;
        };

        let scan_limit = next_term.map_or(n, |(t_start, _, _)| t_start);
        let Some(stop_end) = find_stop_in_frame(seq, pos, scan_limit) else {
            pos += 1;
            continue;
        };

        let orf_len = stop_end - pos;
        let minimum_orf_len = if init_mode == "coupled" {
            MIN_COUPLED_ORF_CODONS * 3
        } else {
            MIN_ORF_LEN
        };
        if orf_len < minimum_orf_len {
            pos += 1;
            continue;
        }

        let promoter_overlaps_start = *prom_end10 > pos;
        let uncertain_rbs = init_mode == "RBS" && !rbs_supported;
        if next_term.is_none() && !non_strict {
            if !*uncertain_promoter && !uncertain_rbs && !promoter_overlaps_start {
                errors.push(format!(
                    "Error: Unterminated transcript starting at @span: {}..{} (missing downstream terminator)",
                    prom_pos, strand
                ));
            }
            active_prom = None;
            pos = stop_end;
            continue;
        }

        let term_dg = next_term.map(|(_, _, dg)| dg);
        let uncertain_terminator = next_term.map_or(false, |(_, _, dg)| {
            floors.term_floor.map_or(false, |floor| dg > floor)
        });
        let mut uncertainty_reasons = Vec::with_capacity(4);
        if *uncertain_promoter {
            uncertainty_reasons.push("WARN: UNCERTAIN_PROMOTER");
        }
        if uncertain_rbs {
            uncertainty_reasons.push("WARN: UNCERTAIN_RBS");
        }
        if uncertain_terminator {
            uncertainty_reasons.push("WARN: UNCERTAIN_TERMINATOR");
        }
        if promoter_overlaps_start {
            uncertainty_reasons.push("WARN: PROMOTER_OVERLAPS_START");
        }
        candidates.push(CandidateGene {
            start: pos,
            end: stop_end,
            strand,
            prom_pos: *prom_pos,
            prom_end: *prom_end10,
            terminator_end: next_term.map(|(_, end, _)| end),
            prom_type: prom_label.to_string(),
            init_mode: init_mode.to_string(),
            promoter_matches: *prom_matches,
            rbs_matches: rbs_score,
            terminator_dg: term_dg,
            unterminated: next_term.is_none(),
            uncertain_evidence: !uncertainty_reasons.is_empty(),
            uncertainty_reasons,
            coupled_to: None,
        });
        pos = stop_end.saturating_sub(4);
    }

    candidates
}

fn candidate_evidence_cmp(left: &CandidateGene, right: &CandidateGene) -> std::cmp::Ordering {
    let terminator_order = match (left.terminator_dg, right.terminator_dg) {
        (Some(left), Some(right)) => left.total_cmp(&right),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    };
    terminator_order
        .then_with(|| right.promoter_matches.cmp(&left.promoter_matches))
        .then_with(|| right.rbs_matches.cmp(&left.rbs_matches))
}

/// Pass 0.5: Preserve explicit TGATG overlaps and abstain on unresolved frame conflicts.
fn resolve_conflicts(
    mut candidates: Vec<CandidateGene>,
    circular: bool,
    orig_len: usize,
) -> Vec<CandidateGene> {
    candidates.sort_by(|left, right| {
        left.start
            .cmp(&right.start)
            .then_with(|| candidate_evidence_cmp(left, right))
    });

    let mut components: Vec<Vec<CandidateGene>> = Vec::new();
    for candidate in candidates {
        let overlaps_last = components.last().map_or(false, |component| {
            component
                .iter()
                .map(|member| member.end)
                .max()
                .map_or(false, |end| candidate.start < end)
        });
        if overlaps_last {
            components.last_mut().unwrap().push(candidate);
        } else {
            components.push(vec![candidate]);
        }
    }

    let mut winners = Vec::new();
    for mut component in components {
        if component.len() == 1 {
            winners.push(component.pop().unwrap());
            continue;
        }

        let is_coupled_chain = component.windows(2).all(|pair| {
            pair[0].strand == pair[1].strand
                && (1..=4).contains(&pair[0].end.saturating_sub(pair[1].start))
        });
        if !is_coupled_chain {
            if let Some(best) = component
                .into_iter()
                .min_by(candidate_evidence_cmp)
            {
                winners.push(best);
            }
            continue;
        }

        let first_index = winners.len();
        for offset in 1..component.len() {
            if component[offset].strand == '+' {
                component[offset].coupled_to = Some(first_index + offset - 1);
            } else {
                component[offset - 1].coupled_to = Some(first_index + offset);
            }
        }
        winners.extend(component);
    }

    if !circular {
        return winners;
    }

    let mut seen: HashSet<(usize, usize, char)> = HashSet::new();
    let mut deduped = Vec::new();
    for g in winners {
        if g.start >= orig_len {
            continue;
        }
        let key = (g.start % orig_len, g.end % orig_len, g.strand);
        if seen.insert(key) {
            deduped.push(g);
        }
    }
    deduped
}

fn push_dedup_switch(slice: &mut Vec<SwitchEntry>, kind: SwitchKind, position: Option<usize>) {
    if let Some(last) = slice.last_mut() {
        let within_repeat_spacing = match (&last.kind, position, last.positions.last()) {
            (SwitchKind::OperatorSiteCandidate(motif), Some(position), Some(&previous_position)) =>
                position.abs_diff(previous_position) == motif.len(),
            _ => true,
        };
        if last.kind == kind && within_repeat_spacing {
            last.count += 1;
            if let Some(position) = position {
                last.positions.push(position);
            }
            return;
        }
    }
    slice.push(SwitchEntry {
        kind,
        count: 1,
        positions: position.into_iter().collect(),
    });
}

/// Pass 1: Intergenic CSR slicer with physical-pairing operator candidates.
fn build_csr_switches(
    fwd_seq: &[u8],
    rev_seq: &[u8],
    fwd_codes: &[u8],
    rev_codes: &[u8],
    genes: &[CandidateGene],
) -> (Vec<usize>, Vec<SwitchEntry>) {
    let n = fwd_seq.len();
    let mut previous_plus = vec![None; genes.len()];
    let mut next_minus = vec![None; genes.len()];
    let mut last_plus = None;
    for (idx, gene) in genes.iter().enumerate() {
        previous_plus[idx] = last_plus;
        if gene.strand == '+' {
            last_plus = Some(idx);
        }
    }
    let mut next_negative = None;
    for idx in (0..genes.len()).rev() {
        next_minus[idx] = next_negative;
        if genes[idx].strand == '-' {
            next_negative = Some(idx);
        }
    }

    let mut row_ptr = Vec::with_capacity(genes.len() + 1);
    let mut col_switches: Vec<SwitchEntry> = Vec::new();
    row_ptr.push(0);

    for (idx, g) in genes.iter().enumerate() {
        let mut local: Vec<SwitchEntry> = Vec::new();

        if let Some(parent_idx) = g.coupled_to {
            push_dedup_switch(&mut local, SwitchKind::Coupled(parent_idx), None);
            col_switches.extend(local);
            row_ptr.push(col_switches.len());
            continue;
        }

        let (strand_seq, strand_codes, s_start) = if g.strand == '+' {
            (fwd_seq, fwd_codes, g.start)
        } else {
            (rev_seq, rev_codes, n - g.end)
        };

        let previous_idx = if g.strand == '+' {
            previous_plus[idx]
        } else {
            next_minus[idx]
        };
        let is_operon_follower = previous_idx.map_or(false, |previous_idx| {
            genes[previous_idx].prom_pos == g.prom_pos
        });
        let prev_end = previous_idx.map_or(0, |previous_idx| {
            let previous_gene_end = if g.strand == '+' {
                genes[previous_idx].end
            } else {
                n - genes[previous_idx].start
            };
            if is_operon_follower {
                previous_gene_end
            } else {
                genes[previous_idx]
                    .terminator_end
                    .unwrap_or(previous_gene_end)
            }
        });

        if !is_operon_follower {
            push_dedup_switch(&mut local, SwitchKind::Promoter(g.prom_type.clone()), None);
        }

        // Exclude the known promoter and the preceding transcription unit's terminator.
        let operator_intervals = if is_operon_follower {
            [(prev_end.min(s_start), s_start), (0, 0)]
        } else {
            [
                (prev_end.min(s_start), g.prom_pos.min(s_start)),
                (g.prom_end.max(prev_end).min(s_start), s_start),
            ]
        };
        for (region_start, region_end) in operator_intervals {
            if region_start >= region_end {
                continue;
            }
            let Some(utr) = strand_seq.get(region_start..region_end) else {
                exit_err(&format!(
                    "internal error: invalid dependency interval (strand={}, gene={}..{}, promoter_end={}, previous_end={}, region_start={}, region_end={}, gene_start={}, sequence_len={}, previous_gene={:?}, operon_follower={})",
                    g.strand,
                    g.start,
                    g.end,
                    g.prom_end,
                    prev_end,
                    region_start,
                    region_end,
                    s_start,
                    strand_seq.len(),
                    previous_idx,
                    is_operon_follower
                ));
            };
            let Some(utr_codes) = strand_codes.get(region_start..region_end) else {
                exit_err("internal error: invalid operator code interval");
            };
            let mut scan = 0;
            while scan + 2 * MIN_OPERATOR_STEM <= utr.len() {
                let Some((len, _dg)) = operator_motif_at(utr, utr_codes, scan) else {
                    scan += 1;
                    continue;
                };
                let motif = String::from_utf8_lossy(&utr[scan..scan + len]).into_owned();
                let oriented_start = region_start + scan;
                let genome_start = if g.strand == '+' {
                    oriented_start
                } else {
                    n - oriented_start - len
                };
                push_dedup_switch(
                    &mut local,
                    SwitchKind::OperatorSiteCandidate(motif),
                    Some(genome_start),
                );
                scan += len;
            }
        }

        if g.init_mode == "leaderless_TSS+1" {
            push_dedup_switch(&mut local, SwitchKind::LeaderlessTss, None);
        } else {
            push_dedup_switch(&mut local, SwitchKind::Rbs, None);
        }

        col_switches.extend(local);
        row_ptr.push(col_switches.len());
    }

    (row_ptr, col_switches)
}

fn format_switch_default(sw: &SwitchEntry) -> String {
    match &sw.kind {
        SwitchKind::Promoter(p) => format!("{}*x{}", p, sw.count),
        SwitchKind::OperatorSiteCandidate(m) => format!("!lock({})*x{}", m, sw.count),
        SwitchKind::Rbs => format!("RBS*x{}", sw.count),
        SwitchKind::LeaderlessTss => format!("leaderless_TSS+1*x{}", sw.count),
        SwitchKind::Coupled(parent) => format!("coupled(gene_{})", parent),
    }
}

fn layer_scopes_and_latches(
    gene_count: usize,
    edges: &[RegulatoryEdge],
) -> (Vec<usize>, HashSet<usize>) {
    let mut dependencies = vec![Vec::new(); gene_count];
    for edge in edges {
        if edge.regulator < gene_count && edge.target < gene_count {
            dependencies[edge.target].push(edge.regulator);
        }
    }

    let mut scopes = vec![0; gene_count];
    let mut latch_nodes = HashSet::new();
    let mut state = vec![0_u8; gene_count];
    let mut dependent_count = vec![0; gene_count];
    for node_dependencies in &dependencies {
        for &dependency in node_dependencies {
            dependent_count[dependency] += 1;
        }
    }

    let mut roots: Vec<usize> = (0..gene_count)
        .filter(|&node| dependent_count[node] == 0)
        .collect();
    roots.extend(0..gene_count);
    for root in roots {
        if state[root] != 0 {
            continue;
        }
        state[root] = 1;
        let mut stack = vec![(root, 0_usize, 0_usize)];
        while let Some(&(node, depth, next_dependency)) = stack.last() {
            if next_dependency == dependencies[node].len() {
                state[node] = 2;
                stack.pop();
                continue;
            }

            let dependency = dependencies[node][next_dependency];
            stack.last_mut().unwrap().2 += 1;
            match state[dependency] {
                0 => {
                    state[dependency] = 1;
                    scopes[dependency] = depth + 1;
                    stack.push((dependency, depth + 1, 0));
                }
                1 => {
                    if let Some(cycle_start) = stack
                        .iter()
                        .position(|(active, _, _)| *active == dependency)
                    {
                        latch_nodes.insert(stack[cycle_start].0);
                    }
                }
                _ => scopes[dependency] = scopes[dependency].min(depth + 1),
            }
        }
    }
    (scopes, latch_nodes)
}

/// Pass 2b: Emit lifted Scope/CSR IR into the requested target language
fn emit_output(
    cfg: &Config,
    orig_len: usize,
    genes: &[CandidateGene],
    row_ptr: &[usize],
    col_switches: &[SwitchEntry],
    scopes: &[usize],
    latch_nodes: &HashSet<usize>,
    regulatory_edges: &[RegulatoryEdge],
) {
    match cfg.emit {
        EmitTarget::Ir => println!("IR_MODULE {} {{", cfg.input_path),
        EmitTarget::Python => println!("# Trans-decompiled by dnadec from {}", cfg.input_path),
        EmitTarget::Verilog => println!("// Trans-decompiled by dnadec from {}", cfg.input_path),
        EmitTarget::Default => {}
    }

    for edge in regulatory_edges {
        let annotation = format!(
            "REGULATORY_EDGE gene_{} {} gene_{} [motif={}*x{}; effect={}]",
            edge.regulator,
            match edge.sign {
                '-' => "-|",
                '+' => "+|",
                _ => "?",
            },
            edge.target,
            edge.motif,
            edge.count,
            if edge.sign == '?' {
                "unknown"
            } else if edge.sign == '-' {
                "repression"
            } else {
                "activation"
            }
        );
        match cfg.emit {
            EmitTarget::Default => println!("// {}", annotation),
            EmitTarget::Ir => println!("  {}", annotation),
            EmitTarget::Python => println!("# {}", annotation),
            EmitTarget::Verilog => println!("// {}", annotation),
        }
    }
    if cfg.emit == EmitTarget::Python && !latch_nodes.is_empty() {
        println!("while circuit_active:  # scope_0_loops [LATCH]");
    }

    for i in 0..genes.len() {
        let g = &genes[i];
        let span_start = if cfg.circular {
            g.start % orig_len
        } else {
            g.start
        };
        let span_end = if cfg.circular {
            g.end % orig_len
        } else {
            g.end
        };
        let is_latch = latch_nodes.contains(&i);
        let mut warnings = Vec::new();
        if g.unterminated {
            warnings.push("WARN: UNTERMINATED_READTHROUGH");
        }
        if g.uncertain_evidence {
            warnings.push("WARN: UNCERTAIN_SEQUENCE_EVIDENCE");
            warnings.extend(g.uncertainty_reasons.iter().copied());
        }
        let warn = if warnings.is_empty() {
            String::new()
        } else {
            format!(" [{}]", warnings.join("] ["))
        };
        let switches = &col_switches[row_ptr[i]..row_ptr[i + 1]];
        if cfg.uncertainty != UncertaintyMode::Omit {
            for switch in switches {
                let SwitchKind::OperatorSiteCandidate(motif) = &switch.kind else {
                    continue;
                };
                for &position in &switch.positions {
                    let start = if cfg.circular {
                        position % orig_len
                    } else {
                        position
                    };
                    let end = if cfg.circular {
                        (position + motif.len()) % orig_len
                    } else {
                        position + motif.len()
                    };
                    let annotation = format!(
                        "CANDIDATE_OPERATOR_SITE gene_{} motif={} span={}..{} strand={}",
                        i, motif, start, end, g.strand
                    );
                    match cfg.emit {
                        EmitTarget::Default | EmitTarget::Verilog => {
                            println!("// {}", annotation)
                        }
                        EmitTarget::Ir => println!("  {}", annotation),
                        EmitTarget::Python => println!("# {}", annotation),
                    }
                }
            }
        }
        let candidate_sites: Vec<String> = if cfg.uncertainty == UncertaintyMode::Omit {
            Vec::new()
        } else {
            switches
                .iter()
                .filter_map(|sw| match &sw.kind {
                    SwitchKind::OperatorSiteCandidate(m)
                        if !regulatory_edges
                            .iter()
                            .any(|edge| {
                                edge.sign != '?'
                                    && edge.target == i
                                    && edge.motif == *m
                            }) =>
                    {
                        Some(format!("{}*x{}", m, sw.count))
                    }
                    _ => None,
                })
                .collect()
        };
        let unlinked_locks: Vec<String> = candidate_sites
            .iter()
            .enumerate()
            .map(|(index, site)| format!("!lock_{}({})", index, site))
            .collect();
        let linked_sites: Vec<&RegulatoryEdge> = regulatory_edges
            .iter()
            .filter(|edge| edge.target == i && edge.sign != '?')
            .collect();
        let candidate_sites_note = if candidate_sites.is_empty() {
            String::new()
        } else {
            format!(
                " [CANDIDATE_OPERATOR_SITES: {}]",
                candidate_sites.join(", ")
            )
        };
        let candidate_cis_edges: Vec<String> = candidate_sites
            .iter()
            .map(|site| {
                format!(
                    "operator_site({}) -> gene_{} [regulator=unknown, effect=unknown]",
                    site, i
                )
            })
            .collect();
        let candidate_cis_edges_note = if candidate_cis_edges.is_empty() {
            String::new()
        } else {
            format!(" [CANDIDATE_CIS_EDGES: {}]", candidate_cis_edges.join(", "))
        };
        let unlinked_locks_note = if unlinked_locks.is_empty() {
            String::new()
        } else {
            format!(" [UNLINKED_OPERATOR_LOCKS: {}]", unlinked_locks.join(", "))
        };
        let dependencies: Vec<&SwitchEntry> = switches
            .iter()
            .filter(|sw| match &sw.kind {
                SwitchKind::OperatorSiteCandidate(m) => {
                    linked_sites.iter().any(|edge| edge.motif == *m)
                }
                _ => true,
            })
            .collect();

        match cfg.emit {
            EmitTarget::Default => {
                let scope_tag = if is_latch {
                    format!("scope_{}_loops [LATCH]", scopes[i])
                } else {
                    format!("scope_{}", scopes[i])
                };
                let conds: Vec<String> = dependencies
                    .iter()
                    .filter_map(|sw| match &sw.kind {
                        SwitchKind::OperatorSiteCandidate(m) => linked_sites
                            .iter()
                            .find(|edge| edge.motif == *m)
                            .map(|edge| {
                                let gate = if edge.sign == '-' { "!lock" } else { "lock" };
                                format!("{}({})*x{}", gate, m, sw.count)
                            }),
                        _ => Some(format_switch_default(sw)),
                    })
                    .collect();
                println!(
                    "{}: if ({}) {{ emit(gene_{}); }} // @span: {}..{} ({}){}{}{}{}",
                    scope_tag,
                    conds.join(", "),
                    i,
                    span_start,
                    span_end,
                    g.strand,
                    warn,
                    candidate_sites_note,
                    candidate_cis_edges_note,
                    unlinked_locks_note
                );
            }
            EmitTarget::Ir => {
                let scope_tag = if is_latch {
                    format!("CALLBACK_SCOPE_{}_LOOPS [LATCH]", scopes[i])
                } else {
                    format!("CALLBACK_SCOPE_{}", scopes[i])
                };
                let deps: Vec<String> = dependencies
                    .iter()
                    .filter_map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => Some(format!("{}(x{})", p, sw.count)),
                        SwitchKind::OperatorSiteCandidate(m) => linked_sites
                            .iter()
                            .find(|edge| edge.motif == *m)
                            .map(|edge| {
                                let gate = if edge.sign == '-' { "NOT" } else { "" };
                                format!("{} lock({})(x{})", gate, m, sw.count)
                                    .trim_start()
                                    .to_string()
                            }),
                        SwitchKind::Rbs => Some(format!("RBS(x{})", sw.count)),
                        SwitchKind::LeaderlessTss => {
                            Some(format!("leaderless_TSS+1(x{})", sw.count))
                        }
                        SwitchKind::Coupled(p) => Some(format!("COUPLED(gene_{})", p)),
                    })
                    .collect();
                println!("  {}:", scope_tag);
                println!(
                    "    CALLBACK gene_{} @span({}..{}, '{}'){}",
                    i, span_start, span_end, g.strand, warn
                );
                if !candidate_sites.is_empty() {
                    println!(
                        "    CANDIDATE_OPERATOR_SITES [ {} ]",
                        candidate_sites.join(", ")
                    );
                    println!(
                        "    CANDIDATE_CIS_EDGES [ {} ]",
                        candidate_cis_edges.join(", ")
                    );
                    println!(
                        "    UNLINKED_OPERATOR_LOCKS [ {} ]",
                        unlinked_locks.join(", ")
                    );
                }
                println!("      DEPENDS_ALL [ {} ]", deps.join(", "));
                println!("      EMIT gene_{}", i);
            }
            EmitTarget::Python => {
                let conds: Vec<String> = dependencies
                    .iter()
                    .filter_map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => Some(p.clone()),
                        SwitchKind::OperatorSiteCandidate(m) => linked_sites
                            .iter()
                            .find(|edge| edge.motif == *m)
                            .map(|edge| {
                                let gate = if edge.sign == '-' { "not " } else { "" };
                                format!("{}lock_{} >= {}", gate, m, sw.count)
                            }),
                        SwitchKind::Rbs => Some("rbs".to_string()),
                        SwitchKind::LeaderlessTss => Some("leaderless_tss_1".to_string()),
                        SwitchKind::Coupled(p) => Some(format!("coupled_gene_{}", p)),
                    })
                    .collect();
                let python_indent = if latch_nodes.is_empty() { "" } else { "    " };
                println!(
                    "{}if {}:  # @span: {}..{} ({}){}",
                    python_indent,
                    conds.join(" and "),
                    span_start,
                    span_end,
                    g.strand,
                    warn
                );
                if !candidate_sites.is_empty() {
                    println!(
                        "{}# CANDIDATE_OPERATOR_SITES: {}",
                        python_indent,
                        candidate_sites.join(", ")
                    );
                    println!(
                        "{}# CANDIDATE_CIS_EDGES: {}",
                        python_indent,
                        candidate_cis_edges.join(", ")
                    );
                    println!(
                        "{}# UNLINKED_OPERATOR_LOCKS: {}",
                        python_indent,
                        unlinked_locks.join(", ")
                    );
                }
                println!("{}    gene_{} = emit(\"gene_{}\")", python_indent, i, i);
            }
            EmitTarget::Verilog => {
                let gates: Vec<String> = dependencies
                    .iter()
                    .filter_map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => Some(p.clone()),
                        SwitchKind::OperatorSiteCandidate(m) => linked_sites
                            .iter()
                            .find(|edge| edge.motif == *m)
                            .map(|edge| {
                                let gate = if edge.sign == '-' { "~" } else { "" };
                                format!("{}(lock_{}_x{})", gate, m, sw.count)
                            }),
                        SwitchKind::Rbs => Some("rbs".to_string()),
                        SwitchKind::LeaderlessTss => Some("leaderless_tss_1".to_string()),
                        SwitchKind::Coupled(p) => Some(format!("coupled_gene_{}", p)),
                    })
                    .collect();
                let latch_note = if is_latch { " // latch feedback" } else { "" };
                println!(
                    "assign gene_{} = {}; // @span: {}..{} ({}){}{}",
                    i,
                    gates.join(" & "),
                    span_start,
                    span_end,
                    g.strand,
                    latch_note,
                    warn
                );
                if !candidate_sites.is_empty() {
                    println!(
                        "// CANDIDATE_OPERATOR_SITES gene_{}: {}",
                        i,
                        candidate_sites.join(", ")
                    );
                    println!("// CANDIDATE_CIS_EDGES: {}", candidate_cis_edges.join(", "));
                    println!("// UNLINKED_OPERATOR_LOCKS: {}", unlinked_locks.join(", "));
                }
            }
        }
    }

    if cfg.emit == EmitTarget::Ir {
        println!("}}");
    }
}

fn report_calibration(floors: &NoiseFloors, sequence_len: usize) {
    eprintln!(
        "Null calibration: {} dinucleotide shuffles, seeds {}..{}; match-count q{}, ΔG q{} lower tail",
        SHUFFLE_COUNT,
        SHUFFLE_SEED,
        SHUFFLE_SEED + SHUFFLE_COUNT as u64 - 1,
        NULL_QUANTILE_PERCENT,
        100 - NULL_QUANTILE_PERCENT
    );
    let reason = if sequence_len < 8 {
        "input shorter than 8 bp"
    } else {
        "zero valid null-model hits"
    };
    for (name, threshold) in [
        ("sigma70 promoter", floors.sigma70_prom_floor.map(|v| v.to_string())),
        ("sigma32 promoter", floors.sigma32_prom_floor.map(|v| v.to_string())),
        ("RBS", floors.rbs_floor.map(|v| v.to_string())),
        ("terminator ΔG", floors.term_floor.map(|v| format!("{v:.3} kcal/mol"))),
    ] {
        if let Some(threshold) = threshold {
            eprintln!("Calibrated {name}: {threshold}");
        } else {
            eprintln!("Uncalibrated {name}: {reason}");
        }
    }
}

fn main() {
    let cfg = parse_args();
    let raw_seq = read_fasta(&cfg.input_path);
    let orig_len = raw_seq.len();

    let work_seq = if cfg.circular {
        let mut doubled = Vec::with_capacity(orig_len * 2);
        doubled.extend_from_slice(&raw_seq);
        doubled.extend_from_slice(&raw_seq);
        doubled
    } else {
        raw_seq
    };

    let rev_seq = reverse_complement(&work_seq);
    let work_codes = encode_bases(&work_seq);
    let rev_codes = reverse_complement_codes(&work_codes);
    let work_sixmers = pack_sixmers(&work_codes);
    let rev_sixmers = pack_sixmers(&rev_codes);
    let work_planes = pack_base_planes(&work_codes);
    let rev_planes = pack_base_planes(&rev_codes);
    let floors = calibrate_noise_floors(&work_seq);
    report_calibration(&floors, work_seq.len());

    let mut errors = Vec::new();
    let mut candidates = scan_strand(
        &work_seq,
        &work_codes,
        &work_sixmers,
        &work_planes,
        &rev_planes,
        '+',
        &floors,
        cfg.non_strict,
        &mut errors,
    );
    let rev_cands = scan_strand(
        &rev_seq,
        &rev_codes,
        &rev_sixmers,
        &rev_planes,
        &work_planes,
        '-',
        &floors,
        cfg.non_strict,
        &mut errors,
    );

    let total_len = work_seq.len();
    for mut g in rev_cands {
        let fwd_start = total_len - g.end;
        let fwd_end = total_len - g.start;
        g.start = fwd_start;
        g.end = fwd_end;
        candidates.push(g);
    }

    if !errors.is_empty() && !cfg.non_strict {
        for err in &errors {
            eprintln!("{}", err);
        }
        process::exit(1);
    }

    if !cfg.non_strict {
        candidates.retain(|gene| !gene.uncertain_evidence);
    }

    let winners = resolve_conflicts(candidates, cfg.circular, orig_len);
    let (row_ptr, col_switches) = build_csr_switches(
        &work_seq,
        &rev_seq,
        &work_codes,
        &rev_codes,
        &winners,
    );
    let regulatory_edges = resolve_shared_motif_edges(&winners, &row_ptr, &col_switches);
    let (scopes, latch_nodes) = layer_scopes_and_latches(winners.len(), &regulatory_edges);
    if cfg.uncertainty == UncertaintyMode::Error {
        for i in 0..winners.len() {
            if let Some(site) = col_switches[row_ptr[i]..row_ptr[i + 1]]
                .iter()
                .find_map(|sw| match &sw.kind {
                    SwitchKind::OperatorSiteCandidate(m) => Some(m),
                    _ => None,
                })
            {
                exit_err(&format!(
                    "unresolved operator-site candidate '{}' associated with gene_{}; choose --uncertainty report or omit",
                    site, i
                ));
            }
        }
    }
    emit_output(
        &cfg,
        orig_len,
        &winners,
        &row_ptr,
        &col_switches,
        &scopes,
        &latch_nodes,
        &regulatory_edges,
    );
}

#[cfg(test)]
mod tests {
    use super::{
        calibrate_noise_floors, comp_base, count_packed_stem_pairs, count_stem_pairs,
        dinucleotide_shuffle, encode_bases, empirical_quantile_usize, hairpin_loop_initiation_dg,
        layer_scopes_and_latches, operator_motif_at, pack_base_planes, pack_consensus,
        pack_sixmers, push_dedup_switch, reverse_complement_codes, score_promoters_at,
        score_terminator_at, sixmer_mismatches, terminator_start_candidates, thymine_tail_counts,
        RegulatoryEdge, SwitchEntry, SwitchKind,
    };

    #[test]
    fn sigma32_consensus_matches_its_declared_spacer() {
        let sequence = b"CCCCCTTGAAGGGGGGGGGGGGGGGCCCCAT";
        let sixmers = pack_sixmers(&encode_bases(sequence));
        let matches = score_promoters_at(&sixmers, 4);
        assert!(
            matches
                .into_iter()
                .flatten()
                .any(|(_, _, label)| label == "promoter_sigma32"),
            "{:?}",
            matches
        );
    }

    #[test]
    fn packed_sixmer_mismatch_count_matches_basewise_comparison() {
        let consensus = pack_consensus(*b"AGGAGG");
        for packed in 0..4096_u16 {
            let expected = (0..6)
                .filter(|&index| {
                    ((packed >> (10 - index * 2)) & 3) != ((consensus >> (10 - index * 2)) & 3)
                })
                .count();
            assert_eq!(sixmer_mismatches(packed, consensus), expected);
        }
    }

    fn dinucleotide_counts(seq: &[u8]) -> [usize; 16] {
        let mut counts = [0; 16];
        for pair in seq.windows(2) {
            let from = match pair[0] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            let to = match pair[1] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            counts[from * 4 + to] += 1;
        }
        counts
    }

    #[test]
    fn dinucleotide_shuffle_preserves_counts_and_is_seeded() {
        let sequence = b"ACGTTGCAACGTACGATCGTACGGTACCTAGCTAGGCTAACGT";
        let shuffled = dinucleotide_shuffle(sequence, 42).unwrap();
        assert_eq!(shuffled.len(), sequence.len());
        assert_eq!(shuffled.first(), sequence.first());
        assert_eq!(shuffled.last(), sequence.last());
        assert_eq!(
            dinucleotide_counts(&shuffled),
            dinucleotide_counts(sequence)
        );
        assert_eq!(Some(shuffled), dinucleotide_shuffle(sequence, 42));
    }

    #[test]
    fn short_or_empty_null_samples_are_reported_as_uncalibrated() {
        let floors = calibrate_noise_floors(b"ACGTACG");
        assert!(floors.sigma70_prom_floor.is_none());
        assert!(floors.sigma32_prom_floor.is_none());
        assert!(floors.rbs_floor.is_none());
        assert!(floors.term_floor.is_none());
        assert_eq!(dinucleotide_shuffle(b"A", 42), None);
    }

    #[test]
    fn empirical_quantile_uses_the_observed_null_values() {
        let mut values = [1, 2, 3, 4, 5];
        assert_eq!(empirical_quantile_usize(&mut values, 95), Some(5));
        let mut empty = [];
        assert_eq!(empirical_quantile_usize(&mut empty, 95), None);
    }

    #[test]
    fn operator_detection_requires_a_physical_palindromic_stem() {
        let motif = b"GATATATATATC";
        assert_eq!(
            operator_motif_at(motif, &encode_bases(motif), 0).map(|(len, _)| len),
            Some(motif.len())
        );
        let non_palindrome = b"TCCCTATCTCCCTATC";
        assert_eq!(
            operator_motif_at(non_palindrome, &encode_bases(non_palindrome), 0),
            None
        );
    }

    #[test]
    fn hairpin_initiation_uses_the_published_loop_length_table() {
        assert!(hairpin_loop_initiation_dg(3).unwrap() > 0.0);
        assert!(hairpin_loop_initiation_dg(8).unwrap() > 0.0);
        assert_eq!(hairpin_loop_initiation_dg(9), None);
    }

    #[test]
    fn graph_pass_marks_only_cycles_in_supplied_edges() {
        let edges = [
            RegulatoryEdge {
                regulator: 0,
                target: 1,
                motif: "AATT".to_string(),
                count: 1,
                sign: '-',
            },
            RegulatoryEdge {
                regulator: 1,
                target: 0,
                motif: "ATAT".to_string(),
                count: 1,
                sign: '-',
            },
        ];
        let (scopes, latches) = layer_scopes_and_latches(2, &edges);
        assert_eq!(scopes.len(), 2);
        assert_eq!(latches.len(), 1);
        assert!(latches.contains(&0));
        assert_eq!(scopes[1], 1);

        let (_, no_latches) = layer_scopes_and_latches(2, &[]);
        assert!(no_latches.is_empty());
    }

    #[test]
    fn packed_stem_pair_counts_match_the_base_pair_rule() {
        let bases = [b'A', b'C', b'G', b'T'];
        for &left in &bases {
            for &right in &bases {
                let mut sequence = [b'A'; 12];
                sequence[0] = left;
                sequence[11] = right;
                let codes = encode_bases(&sequence);
                let expected = usize::from(left == comp_base(right));
                assert_eq!(count_stem_pairs(&codes, 0, 11, 1), expected);
                let reverse_codes = reverse_complement_codes(&codes);
                let left_planes = pack_base_planes(&codes);
                let reverse_planes = pack_base_planes(&reverse_codes);
                assert_eq!(
                    count_packed_stem_pairs(&left_planes, &reverse_planes, 0, 0, 1),
                    expected,
                    "left={} right={}",
                    left as char,
                    right as char
                );
            }
        }
        for &left0 in &bases {
            for &left1 in &bases {
                for &right0 in &bases {
                    for &right1 in &bases {
                        let mut sequence = [b'A'; 12];
                        sequence[0] = left0;
                        sequence[1] = left1;
                        sequence[10] = right0;
                        sequence[11] = right1;
                        let codes = encode_bases(&sequence);
                        let reverse_codes = reverse_complement_codes(&codes);
                        let left_planes = pack_base_planes(&codes);
                        let reverse_planes = pack_base_planes(&reverse_codes);
                        assert_eq!(
                            count_packed_stem_pairs(&left_planes, &reverse_planes, 0, 0, 2),
                            count_stem_pairs(&codes, 0, 10, 2),
                            "left={}{} right={}{}",
                            left0 as char,
                            left1 as char,
                            right0 as char,
                            right1 as char
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn terminator_scanner_accepts_four_ts_in_an_eight_base_tail() {
        let sequence = b"GCGCGCGCGAAAGCGCGCGCTTTTTT";
        let codes = encode_bases(sequence);
        let reverse_codes = reverse_complement_codes(&codes);
        let planes = pack_base_planes(&codes);
        let reverse_planes = pack_base_planes(&reverse_codes);
        let tail_counts = thymine_tail_counts(sequence);
        assert!(
            score_terminator_at(sequence, &codes, &planes, &reverse_planes, &tail_counts, 0)
                .is_some()
        );

        let noncontiguous_tail = b"GCGCGCGCGAAAGCGCGCGCTTATTT";
        let noncontiguous_codes = encode_bases(noncontiguous_tail);
        let noncontiguous_reverse_codes = reverse_complement_codes(&noncontiguous_codes);
        let noncontiguous_planes = pack_base_planes(&noncontiguous_codes);
        let noncontiguous_reverse_planes = pack_base_planes(&noncontiguous_reverse_codes);
        let noncontiguous_counts = thymine_tail_counts(noncontiguous_tail);
        assert!(score_terminator_at(
            noncontiguous_tail,
            &noncontiguous_codes,
            &noncontiguous_planes,
            &noncontiguous_reverse_planes,
            &noncontiguous_counts,
            0
        )
        .is_some());
    }

    #[test]
    fn thymine_tail_table_matches_every_eight_base_window() {
        let sequence = b"TTTTATATTTTTAAAAT";
        let counts = thymine_tail_counts(sequence);
        for start in 0..=sequence.len() - super::TAIL_WIN {
            let expected = sequence[start..start + super::TAIL_WIN]
                .iter()
                .filter(|&&base| base == b'T')
                .count() as u8;
            assert_eq!(counts[start], expected);
        }
    }

    #[test]
    fn terminator_prefilter_matches_all_possible_hairpin_offsets() {
        let sequence = b"TTTTATATTTTTAAAATGCGCGCAAAGCGCGCTTTTTT";
        let tails = thymine_tail_counts(sequence);
        let candidates = terminator_start_candidates(&tails);
        for start in 0..sequence.len() {
            let expected = (super::MIN_STEM..=super::MAX_STEM).any(|stem| {
                (super::MIN_LOOP..=super::MAX_LOOP).any(|loop_len| {
                    let tail_start = start + 2 * stem + loop_len;
                    tail_start + super::TAIL_WIN <= sequence.len()
                        && tails[tail_start] >= super::MIN_TAIL_THYMINE_COUNT
                })
            });
            assert_eq!(candidates[start] != 0, expected);
        }
    }

    #[test]
    fn operator_repeat_compression_requires_contiguous_sites() {
        let motif = || SwitchKind::OperatorSiteCandidate("GATATATATATC".to_string());
        let mut switches: Vec<SwitchEntry> = Vec::new();
        push_dedup_switch(&mut switches, motif(), Some(100));
        push_dedup_switch(&mut switches, motif(), Some(112));
        push_dedup_switch(&mut switches, motif(), Some(150));

        assert_eq!(switches.len(), 2);
        assert_eq!(switches[0].count, 2);
        assert_eq!(switches[0].positions, [100, 112]);
        assert_eq!(switches[1].count, 1);
        assert_eq!(switches[1].positions, [150]);
    }
}
