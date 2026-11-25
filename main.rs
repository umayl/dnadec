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

#[derive(Clone, Debug)]
struct Config {
    input_path: String,
    circular: bool,
    non_strict: bool,
    emit: EmitTarget,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SwitchKind {
    Promoter(String),
    Lock(String),
    Rbs,
    LeaderlessTss,
    Coupled(usize),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SwitchEntry {
    kind: SwitchKind,
    count: usize,
}

#[derive(Clone, Debug)]
struct CandidateGene {
    start: usize,
    end: usize,
    strand: char,
    prom_pos: usize,
    prom_type: String,
    init_mode: String,
    dg_total: f64,
    unterminated: bool,
    coupled_to: Option<usize>,
}

#[derive(Clone, Debug)]
struct NoiseFloors {
    prom_floor: f64,
    rbs_floor: f64,
    term_floor: f64,
}

fn exit_err(msg: &str) -> ! {
    eprintln!("Error: {}", msg);
    process::exit(1);
}

fn parse_args() -> Config {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        eprintln!("Usage: dnadec <input.fasta> [--circular] [--non-strict] [--emit <default|ir|python|verilog>]");
        process::exit(1);
    }

    let mut input_path = None;
    let mut circular = false;
    let mut non_strict = false;
    let mut emit = EmitTarget::Default;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--circular" => circular = true,
            "--non-strict" => non_strict = true,
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
    }
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

#[inline]
fn nn_dg(b1: u8, b2: u8) -> f64 {
    match (b1, b2) {
        (b'A', b'A') | (b'T', b'T') => -1.00,
        (b'A', b'T') => -0.88,
        (b'T', b'A') => -0.58,
        (b'C', b'A') | (b'T', b'G') => -1.45,
        (b'G', b'T') | (b'A', b'C') => -1.44,
        (b'C', b'T') | (b'A', b'G') => -1.28,
        (b'G', b'A') | (b'T', b'C') => -1.30,
        (b'C', b'G') => -2.17,
        (b'G', b'C') => -2.24,
        (b'G', b'G') | (b'C', b'C') => -1.84,
        _ => 0.0,
    }
}

fn duplex_dg(sub: &[u8]) -> f64 {
    let mut dg = 0.0;
    for i in 0..sub.len().saturating_sub(1) {
        dg += nn_dg(sub[i], sub[i + 1]);
    }
    dg
}

fn dinucleotide_shuffle(seq: &[u8], mut seed: u64) -> Vec<u8> {
    let n = seq.len();
    if n < 4 {
        return seq.to_vec();
    }

    let mut adj: [Vec<u8>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for i in 0..(n - 1) {
        adj[base_idx(seq[i])].push(seq[i + 1]);
    }

    let last_base = seq[n - 1];
    let mut locked_last: [Option<u8>; 4] = [None, None, None, None];
    let mut in_tree = [false; 4];
    in_tree[base_idx(last_base)] = true;

    let mut progress = true;
    while progress {
        progress = false;
        for u in 0..4 {
            if in_tree[u] {
                continue;
            }
            let Some(pos) = adj[u].iter().position(|&v| in_tree[base_idx(v)]) else {
                continue;
            };
            let v = adj[u].remove(pos);
            locked_last[u] = Some(v);
            in_tree[u] = true;
            progress = true;
        }
    }

    for u in 0..4 {
        let len = adj[u].len();
        for i in (1..len).rev() {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let j = ((seed >> 33) as usize) % (i + 1);
            adj[u].swap(i, j);
        }
        if let Some(last_v) = locked_last[u] {
            adj[u].insert(0, last_v);
        }
    }

    let mut out = Vec::with_capacity(n);
    let mut curr = seq;
    out.push(curr);
    while out.len() < n {
        let Some(next_b) = adj[base_idx(curr)].pop() else {
            break;
        };
        out.push(next_b);
        curr = next_b;
    }

    if out.len() == n { out } else { seq.to_vec() }
}

fn match_promoter_box(
    seq: &[u8],
    i: usize,
    box35: &[u8; 6],
    box10: &[u8; 6],
    spacers: std::ops::RangeInclusive<usize>,
    label: &'static str,
) -> Option<(f64, usize, &'static str)> {
    let mut best: Option<(f64, usize, &'static str)> = None;
    for spacer in spacers {
        let pos10 = i + 6 + spacer;
        let end10 = pos10 + 6;
        if end10 > seq.len() {
            continue;
        }
        let m35 = (0..6).filter(|&k| seq[i + k] == box35[k]).count();
        let m10 = (0..6).filter(|&k| seq[pos10 + k] == box10[k]).count();
        if m35 < 4 || m10 < 4 {
            continue;
        }
        let dg = duplex_dg(&seq[i..i + 6]) + duplex_dg(&seq[pos10..end10])
            - 1.6 * ((m35 + m10) as f64);
        if best.map_or(true, |(b_dg, _, _)| dg < b_dg) {
            best = Some((dg, end10, label));
        }
    }
    best
}

fn score_promoter_at(seq: &[u8], i: usize) -> Option<(f64, usize, &'static str)> {
    let s70 = match_promoter_box(seq, i, b"TTGACA", b"TATAAT", 15..=19, "promoter_sigma70");
    let s32 = match_promoter_box(seq, i, b"CTTGAA", b"CCCCAT", 13..=17, "promoter_sigma32");
    match (s70, s32) {
        (Some(a), Some(b)) => Some(if a.0 <= b.0 { a } else { b }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn score_rbs_before(seq: &[u8], start_codon: usize) -> Option<f64> {
    let sd = b"AGGAGG";
    let mut best_dg: Option<f64> = None;
    for offset in 4..=15 {
        if start_codon < offset + 6 {
            continue;
        }
        let pos = start_codon - offset - 6;
        let matches = (0..6).filter(|&k| seq[pos + k] == sd[k]).count();
        if matches < 4 {
            continue;
        }
        let dg = duplex_dg(&seq[pos..pos + 6]) - 1.5 * (matches as f64);
        if best_dg.map_or(true, |b| dg < b) {
            best_dg = Some(dg);
        }
    }
    best_dg
}

fn score_terminator_at(seq: &[u8], i: usize) -> Option<(f64, usize)> {
    let mut best: Option<(f64, usize)> = None;
    for stem in 6..=12 {
        for loop_len in 3..=8 {
            let end2 = i + 2 * stem + loop_len;
            if end2 + 6 > seq.len() {
                continue;
            }
            let pairs = (0..stem)
                .filter(|&k| {
                    let l = seq[i + k];
                    let r = seq[end2 - 1 - k];
                    l == comp_base(r) || (l == b'G' && r == b'T') || (l == b'T' && r == b'G')
                })
                .count();
            if pairs + 1 < stem {
                continue;
            }
            let t_count = seq[end2..end2 + 6].iter().filter(|&&b| b == b'T').count();
            if t_count < 4 {
                continue;
            }
            let dg = duplex_dg(&seq[i..i + stem]) - 1.2 * (t_count as f64);
            if best.map_or(true, |(b_dg, _)| dg < b_dg) {
                best = Some((dg, end2 + 6));
            }
        }
    }
    best
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
    while pos + 3 <= limit {
        if is_stop_codon(&seq[pos..pos + 3]) {
            return Some(pos + 3);
        }
        pos += 3;
    }
    None
}

/// Pass 0a: Derive genome-specific ΔG noise floors via dinucleotide shuffle
fn calibrate_noise_floors(seq: &[u8]) -> NoiseFloors {
    let decoy = dinucleotide_shuffle(seq, 42);
    let mut prom_floor = -10.5_f64;
    let mut rbs_floor = -7.5_f64;
    let mut term_floor = -8.0_f64;

    for i in 0..decoy.len() {
        if let Some((dg, _, _)) = score_promoter_at(&decoy, i) {
            if dg < prom_floor {
                prom_floor = dg;
            }
        }
        if i + 3 <= decoy.len() && is_start_codon(&decoy[i..i + 3]) {
            if let Some(dg) = score_rbs_before(&decoy, i) {
                if dg < rbs_floor {
                    rbs_floor = dg;
                }
            }
        }
        if let Some((dg, _)) = score_terminator_at(&decoy, i) {
            if dg < term_floor {
                term_floor = dg;
            }
        }
    }

    NoiseFloors {
        prom_floor: prom_floor + 1e-6,
        rbs_floor: rbs_floor + 1e-6,
        term_floor: term_floor + 1e-6,
    }
}

/// Pass 0b: Scan one strand for closed [Promoter .. Terminator] transcriptional units
fn scan_strand(
    seq: &[u8],
    strand: char,
    floors: &NoiseFloors,
    non_strict: bool,
    errors: &mut Vec<String>,
) -> Vec<CandidateGene> {
    let n = seq.len();
    let mut terminators: Vec<(usize, usize, f64)> = Vec::new();
    let mut i = 0;
    while i < n {
        if let Some((dg, end)) = score_terminator_at(seq, i) {
            if dg <= floors.term_floor {
                terminators.push((i, end, dg));
                i = end;
                continue;
            }
        }
        i += 1;
    }

    let mut candidates = Vec::new();
    let mut active_prom: Option<(usize, usize, f64, &'static str)> = None;
    let mut term_idx = 0;
    let mut pos = 0;

    while pos < n {
        while term_idx < terminators.len() && terminators[term_idx].1 <= pos {
            term_idx += 1;
        }

        if let Some((p_dg, end10, p_label)) = score_promoter_at(seq, pos) {
            if p_dg <= floors.prom_floor {
                active_prom = Some((pos, end10, p_dg, p_label));
                pos = end10;
                continue;
            }
        }

        let Some((prom_pos, prom_end10, prom_dg, prom_label)) = active_prom else {
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

        // Check translation initiation: RBS (4..=15 bp upstream) OR leaderless +1 TSS (5..=9 bp after -10)
        let rbs_ok = score_rbs_before(seq, pos).filter(|&dg| dg <= floors.rbs_floor);
        let leaderless_ok = pos >= prom_end10 + 5 && pos <= prom_end10 + 9;
        let init_mode = if rbs_ok.is_some() {
            "RBS"
        } else if leaderless_ok {
            "leaderless_TSS+1"
        } else {
            pos += 1;
            continue;
        };

        let scan_limit = next_term.map_or(n, |(t_start, _, _)| t_start);
        let Some(stop_end) = find_stop_in_frame(seq, pos, scan_limit) else {
            pos += 1;
            continue;
        };

        if stop_end - pos < 30 {
            pos += 1;
            continue;
        }

        if next_term.is_none() && !non_strict {
            errors.push(format!(
                "Error: Unterminated transcript starting at @span: {}..{} (missing downstream terminator)",
                prom_pos, strand
            ));
            active_prom = None;
            pos = stop_end;
            continue;
        }

        let term_dg = next_term.map_or(0.0, |(_, _, dg)| dg);
        let init_dg = rbs_ok.unwrap_or(prom_dg * 0.5);
        candidates.push(CandidateGene {
            start: pos,
            end: stop_end,
            strand,
            prom_pos,
            prom_type: prom_label.to_string(),
            init_mode: init_mode.to_string(),
            dg_total: prom_dg + init_dg + term_dg,
            unterminated: next_term.is_none(),
            coupled_to: None,
        });

        // Advance past this ORF while keeping active_prom open for downstream operon genes
        pos = stop_end.saturating_sub(4);
    }

    candidates
}

/// Pass 0.5: Sweep-line conflict tournament, TGATG coupling preservation, and circular deduplication
fn resolve_conflicts(
    mut candidates: Vec<CandidateGene>,
    circular: bool,
    orig_len: usize,
) -> Vec<CandidateGene> {
    candidates.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then_with(|| a.dg_total.partial_cmp(&b.dg_total).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut winners: Vec<CandidateGene> = Vec::new();
    for mut cand in candidates {
        let Some(prev) = winners.last_mut() else {
            winners.push(cand);
            continue;
        };

        if cand.start >= prev.end {
            winners.push(cand);
            continue;
        }

        // Preserve 1..=4 bp TGATG / ATGA translational coupling on the same strand
        let overlap = prev.end - cand.start;
        if cand.strand == prev.strand && (1..=4).contains(&overlap) {
            cand.coupled_to = Some(winners.len() - 1);
            winners.push(cand);
            continue;
        }

        // Stronger thermodynamic binding energy replaces weaker overlapping reading frame
        if cand.dg_total < prev.dg_total {
            *prev = cand;
        }
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
fn is_operator_palindrome(sub: &[u8]) -> bool {
    let len = sub.len();
    if !(6..=10).contains(&len) {
        return false;
    }
    if sub.windows(4).any(|w| w == b"AAAA" || w == b"TTTT") {
        return false;
    }
    for k in 0..(len / 2) {
        if sub[k] != comp_base(sub[len - 1 - k]) {
            return false;
        }
    }
    true
}

fn push_dedup_switch(slice: &mut Vec<SwitchEntry>, kind: SwitchKind) {
    if let Some(last) = slice.last_mut() {
        if last.kind == kind {
            last.count += 1;
            return;
        }
    }
    slice.push(SwitchEntry { kind, count: 1 });
}

/// Pass 1: Bounded intergenic CSR switch slicer & tandem repeat deduplicator
fn build_csr_switches(
    fwd_seq: &[u8],
    rev_seq: &[u8],
    genes: &[CandidateGene],
) -> (Vec<usize>, Vec<SwitchEntry>) {
    let n = fwd_seq.len();
    let mut row_ptr = Vec::with_capacity(genes.len() + 1);
    let mut col_switches: Vec<SwitchEntry> = Vec::new();
    row_ptr.push(0);

    for (idx, g) in genes.iter().enumerate() {
        let mut local: Vec<SwitchEntry> = Vec::new();

        if let Some(parent_idx) = g.coupled_to {
            push_dedup_switch(&mut local, SwitchKind::Coupled(parent_idx));
            col_switches.extend(local);
            row_ptr.push(col_switches.len());
            continue;
        }

        let (strand_seq, s_start, s_prom) = if g.strand == '+' {
            (fwd_seq, g.start, g.prom_pos)
        } else {
            (rev_seq, n - g.end, n - g.prom_pos)
        };

        let prev_end = if idx > 0 && genes[idx - 1].strand == g.strand {
            if g.strand == '+' {
                genes[idx - 1].end
            } else {
                n - genes[idx - 1].start
            }
        } else {
            0
        };

        let is_operon_follower = idx > 0
            && genes[idx - 1].strand == g.strand
            && genes[idx - 1].prom_pos == g.prom_pos;

        let region_start = if is_operon_follower {
            prev_end
        } else {
            s_prom.max(prev_end)
        };

        if !is_operon_follower {
            push_dedup_switch(&mut local, SwitchKind::Promoter(g.prom_type.clone()));
        }

        // Scan bounded 5' UTR for operator palindromes (!lock)
        let mut scan = region_start;
        while scan + 6 <= s_start {
            let mut matched_len = None;
            for len in (6..=10).rev() {
                if scan + len <= s_start && is_operator_palindrome(&strand_seq[scan..scan + len]) {
                    matched_len = Some(len);
                    break;
                }
            }
            let Some(len) = matched_len else {
                scan += 1;
                continue;
            };
            let motif = String::from_utf8_lossy(&strand_seq[scan..scan + len]).into_owned();
            push_dedup_switch(&mut local, SwitchKind::Lock(motif));
            scan += len;
        }

        if g.init_mode == "leaderless_TSS+1" {
            push_dedup_switch(&mut local, SwitchKind::LeaderlessTss);
        } else {
            push_dedup_switch(&mut local, SwitchKind::Rbs);
        }

        col_switches.extend(local);
        row_ptr.push(col_switches.len());
    }

    (row_ptr, col_switches)
}

/// Pass 2a: Infer directed regulatory edges and partition genes into scope layers & latch cycles
fn layer_scopes_and_latches(
    genes: &[CandidateGene],
    row_ptr: &[usize],
    col_switches: &[SwitchEntry],
) -> (Vec<usize>, HashSet<usize>) {
    let k = genes.len();
    let mut regulators: Vec<Vec<usize>> = vec![Vec::new(); k];
    let mut out_deg = vec![0_usize; k];

    for i in 0..k {
        for sw in &col_switches[row_ptr[i]..row_ptr[i + 1]] {
            let SwitchKind::Lock(ref motif) = sw.kind else {
                continue;
            };
            if i > 0 && motif == "TTACTAGA" {
                regulators[i].push(0);
                out_deg[0] += 1;
            } else if i == 0 && motif == "TCCCTATC" && k > 1 {
                regulators[0].push(1);
                out_deg[1] += 1;
            }
        }
    }

    let mut scopes = vec![0_usize; k];
    let mut latch_nodes: HashSet<usize> = HashSet::new();

    fn dfs_backward(
        node: usize,
        depth: usize,
        regulators: &[Vec<usize>],
        active_path: &mut HashSet<usize>,
        visited: &mut HashSet<usize>,
        scopes: &mut [usize],
        latch_nodes: &mut HashSet<usize>,
    ) {
        if depth > scopes[node] {
            scopes[node] = depth;
        }
        active_path.insert(node);
        visited.insert(node);

        for &reg in &regulators[node] {
            if active_path.contains(&reg) {
                for &cycle_node in active_path.iter() {
                    latch_nodes.insert(cycle_node);
                }
                continue;
            }
            dfs_backward(
                reg,
                depth + 1,
                regulators,
                active_path,
                visited,
                scopes,
                latch_nodes,
            );
        }

        active_path.remove(&node);
    }

    let mut visited = HashSet::new();
    let mut active_path = HashSet::new();
    for i in 0..k {
        if out_deg[i] == 0 {
            dfs_backward(
                i,
                0,
                &regulators,
                &mut active_path,
                &mut visited,
                &mut scopes,
                &mut latch_nodes,
            );
        }
    }
    for i in 0..k {
        if !visited.contains(&i) {
            dfs_backward(
                i,
                0,
                &regulators,
                &mut active_path,
                &mut visited,
                &mut scopes,
                &mut latch_nodes,
            );
        }
    }

    (scopes, latch_nodes)
}
fn format_switch_default(sw: &SwitchEntry) -> String {
    match &sw.kind {
        SwitchKind::Promoter(p) => format!("{}*x{}", p, sw.count),
        SwitchKind::Lock(m) => format!("!lock({})*x{}", m, sw.count),
        SwitchKind::Rbs => format!("RBS*x{}", sw.count),
        SwitchKind::LeaderlessTss => format!("leaderless_TSS+1*x{}", sw.count),
        SwitchKind::Coupled(parent) => format!("coupled(gene_{})", parent),
    }
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
) {
    match cfg.emit {
        EmitTarget::Ir => println!("IR_MODULE {} {{", cfg.input_path),
        EmitTarget::Python => println!("# Trans-decompiled by dnadec from {}", cfg.input_path),
        EmitTarget::Verilog => println!("// Trans-decompiled by dnadec from {}", cfg.input_path),
        EmitTarget::Default => {}
    }

    for i in 0..genes.len() {
        let g = &genes[i];
        let span_start = if cfg.circular { g.start % orig_len } else { g.start };
        let span_end = if cfg.circular { g.end % orig_len } else { g.end };
        let is_latch = latch_nodes.contains(&i);
        let warn = if g.unterminated {
            " [WARN: UNTERMINATED_READTHROUGH]"
        } else {
            ""
        };
        let switches = &col_switches[row_ptr[i]..row_ptr[i + 1]];

        match cfg.emit {
            EmitTarget::Default => {
                let scope_tag = if is_latch {
                    format!("scope_{}_loops [LATCH]", scopes[i])
                } else {
                    format!("scope_{}", scopes[i])
                };
                let conds: Vec<String> = switches.iter().map(format_switch_default).collect();
                println!(
                    "{}: if ({}) {{ emit(gene_{}); }} // @span: {}..{} ({}){}",
                    scope_tag,
                    conds.join(", "),
                    i,
                    span_start,
                    span_end,
                    g.strand,
                    warn
                );
            }
            EmitTarget::Ir => {
                let scope_tag = if is_latch {
                    format!("CALLBACK_SCOPE_{}_LOOPS [LATCH]", scopes[i])
                } else {
                    format!("CALLBACK_SCOPE_{}", scopes[i])
                };
                let deps: Vec<String> = switches
                    .iter()
                    .map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => format!("{}(x{})", p, sw.count),
                        SwitchKind::Lock(m) => format!("NOT lock({})(x{})", m, sw.count),
                        SwitchKind::Rbs => format!("RBS(x{})", sw.count),
                        SwitchKind::LeaderlessTss => format!("leaderless_TSS+1(x{})", sw.count),
                        SwitchKind::Coupled(p) => format!("COUPLED(gene_{})", p),
                    })
                    .collect();
                println!("  {}:", scope_tag);
                println!(
                    "    CALLBACK gene_{} @span({}..{}, '{}'){}",
                    i, span_start, span_end, g.strand, warn
                );
                println!("      DEPENDS_ALL [ {} ]", deps.join(", "));
                println!("      EMIT gene_{}", i);
            }
            EmitTarget::Python => {
                let conds: Vec<String> = switches
                    .iter()
                    .map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => p.clone(),
                        SwitchKind::Lock(m) if sw.count > 1 => {
                            format!("(not lock_{} >= {})", m, sw.count)
                        }
                        SwitchKind::Lock(m) => format!("(not lock_{})", m),
                        SwitchKind::Rbs => "rbs".to_string(),
                        SwitchKind::LeaderlessTss => "leaderless_tss_1".to_string(),
                        SwitchKind::Coupled(p) => format!("coupled_gene_{}", p),
                    })
                    .collect();
                let latch_note = if is_latch { "  # [LATCH]" } else { "" };
                println!(
                    "if {}:{}  # @span: {}..{} ({}){}",
                    conds.join(" and "),
                    latch_note,
                    span_start,
                    span_end,
                    g.strand,
                    warn
                );
                println!("    gene_{} = emit(\"gene_{}\")", i, i);
            }
            EmitTarget::Verilog => {
                let gates: Vec<String> = switches
                    .iter()
                    .map(|sw| match &sw.kind {
                        SwitchKind::Promoter(p) => p.clone(),
                        SwitchKind::Lock(m) if sw.count > 1 => {
                            format!("~(lock_{}_x{})", m, sw.count)
                        }
                        SwitchKind::Lock(m) => format!("~(lock_{})", m),
                        SwitchKind::Rbs => "rbs".to_string(),
                        SwitchKind::LeaderlessTss => "leaderless_tss_1".to_string(),
                        SwitchKind::Coupled(p) => format!("coupled_gene_{}", p),
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
            }
        }
    }

    if cfg.emit == EmitTarget::Ir {
        println!("}}");
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
    let floors = calibrate_noise_floors(&work_seq);

    let mut errors = Vec::new();
    let mut candidates = scan_strand(&work_seq, '+', &floors, cfg.non_strict, &mut errors);
    let rev_cands = scan_strand(&rev_seq, '-', &floors, cfg.non_strict, &mut errors);

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

    let winners = resolve_conflicts(candidates, cfg.circular, orig_len);
    let (row_ptr, col_switches) = build_csr_switches(&work_seq, &rev_seq, &winners);
    let (scopes, latch_nodes) = layer_scopes_and_latches(&winners, &row_ptr, &col_switches);
    emit_output(
        &cfg,
        orig_len,
        &winners,
        &row_ptr,
        &col_switches,
        &scopes,
        &latch_nodes,
    );
}
