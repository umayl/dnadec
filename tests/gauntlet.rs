use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

struct TempFasta {
    path: PathBuf,
}

impl TempFasta {
    fn new(name: &str, seq: &str) -> Self {
        let mut path = env::temp_dir();
        path.push(format!("dnadec_test_{}_{}.fasta", std::process::id(), name));
        let content = format!(">{}\n{}\n", name, seq);
        fs::write(&path, content).expect("failed to write temp fasta");
        Self { path }
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

fn run_dnadec(fasta: &TempFasta, extra_args: &[&str]) -> (bool, String, String) {
    let bin = env!("CARGO_BIN_EXE_dnadec");
    let mut cmd = Command::new(bin);
    cmd.arg(&fasta.path);
    for arg in extra_args {
        cmd.arg(arg);
    }
    let out = cmd.output().expect("failed to execute dnadec binary");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn deterministic_noise(len: usize, gc_bias: bool, mut seed: u64) -> String {
    let alphabet: &[u8] = if gc_bias {
        b"GGGCCCGCCGCGCAT"
    } else {
        b"AAAATTTTAAATTCG"
    };
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let idx = ((seed >> 33) as usize) % alphabet.len();
        out.push(alphabet[idx]);
    }
    String::from_utf8(out).unwrap()
}

fn orf_body(codons: usize) -> String {
    "GCT".repeat(codons)
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
}

const PROM_S70: &str = "TTGACAAAAAAAAAAAAAAAAAATATAAT";
const PROM_S32: &str = "CTTGAAGGGGGGGGGGGGGGGCCCCAT";
const RBS_SPACER: &str = "AGGAGGAAAAAAA";
const TERM_HAIRPIN: &str = "GCGCGCGCGAAAGCGCGCGCTTTTTT";
const OPERATOR: &str = "GATATATATATC";
const OTHER_OPERATOR: &str = "CGCGCGCGCGCG";

#[test]
fn test_01_single_standard_gene() {
    let seq = format!(
        "CCCCCCCCCC{}AAAAAAA{}ATG{}TAACCCCC{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(100),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("single_gene", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok, "stderr: {}", stderr);
    assert!(
        stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"),
        "ACTUAL STDOUT:\n{}",
        stdout
    );
    assert!(!stdout.contains("UNCERTAIN_RBS"), "{}", stdout);
}

#[test]
fn test_02_palindromic_sites_are_unresolved_cis_candidates() {
    // Motif location supports a candidate target, but not a regulator or effect.
    let seq = format!(
        "CCCC{}{}{}{}{}ATG{}TAA{}CCCC{}CCCCCCCC{}ATG{}TAA{}",
        PROM_S70,
        OPERATOR,
        OPERATOR,
        OPERATOR,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN,
        PROM_S70,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("toggle_switch", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);
    assert!(ok, "stderr: {}", stderr);
    assert!(
        stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"),
        "STDOUT:\n{}\nSTDERR:\n{}",
        stdout,
        stderr
    );
    assert!(
        stdout.contains(&format!("[CANDIDATE_OPERATOR_SITES: {}*x3", OPERATOR)),
        "{}",
        stdout
    );
    assert!(
        stdout.contains(&format!(
            "operator_site({}*x3) -> gene_0 [regulator=unknown, effect=unknown]",
            OPERATOR
        ))
    );
    assert!(!stdout.contains("LATCH"));
    assert!(!stdout.contains("!lock("));
    assert!(!stdout.contains("REGULATORY_EDGE"));

    let (ok_ir, out_ir, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "ir"]);
    assert!(ok_ir);
    assert!(out_ir.contains(&format!("CANDIDATE_OPERATOR_SITES [ {}*x3", OPERATOR)));
    assert!(out_ir.contains(&format!(
        "CANDIDATE_CIS_EDGES [ operator_site({}*x3) -> gene_0 [regulator=unknown, effect=unknown]",
        OPERATOR
    )));
    assert!(!out_ir.contains("NOT lock("));

    let (ok_py, out_py, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "python"]);
    assert!(ok_py);
    assert!(out_py.contains(&format!("# CANDIDATE_OPERATOR_SITES: {}*x3", OPERATOR)));
    assert!(out_py.contains(
        &format!("# CANDIDATE_CIS_EDGES: operator_site({}*x3) -> gene_0 [regulator=unknown, effect=unknown]", OPERATOR)
    ));
    assert!(!out_py.contains("not lock_"));

    let (ok_v, out_v, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "verilog"]);
    assert!(ok_v);
    assert!(out_v.contains(&format!("// CANDIDATE_OPERATOR_SITES gene_0: {}*x3", OPERATOR)));
    assert!(out_v.contains(
        &format!("// CANDIDATE_CIS_EDGES: operator_site({}*x3) -> gene_0 [regulator=unknown, effect=unknown]", OPERATOR)
    ));
    assert!(!out_v.contains("lock_GATATATATATC"));

    let (ok_omit, out_omit, _) = run_dnadec(&fasta, &["--non-strict", "--uncertainty", "omit"]);
    assert!(ok_omit);
    assert!(out_omit.contains("emit(gene_0)"));
    assert!(!out_omit.contains("CANDIDATE_OPERATOR_SITES"));
    assert!(!out_omit.contains("CANDIDATE_CIS_EDGES"));

    let (ok_error, out_error, err_error) =
        run_dnadec(&fasta, &["--non-strict", "--uncertainty", "error"]);
    assert!(!ok_error);
    assert!(out_error.is_empty());
    assert!(err_error.contains(&format!("unresolved operator-site candidate '{}'", OPERATOR)));
}

#[test]
fn test_shared_local_motifs_create_unknown_feedback_candidates() {
    let seq = format!(
        "CCCC{}{}{}ATG{}TAA{}CCCC{}{}{}ATG{}TAA{}",
        PROM_S70,
        OPERATOR,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN,
        PROM_S70,
        OPERATOR,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("shared_motif_feedback", &seq);
    let (ok, output, stderr) = run_dnadec(&fasta, &["--non-strict"]);

    assert!(ok, "stderr: {}", stderr);
    assert!(output.contains("REGULATORY_EDGE gene_0 ? gene_1"));
    assert!(output.contains("REGULATORY_EDGE gene_1 ? gene_0"));
    assert!(output.contains("effect=unknown"));
    assert!(output.contains("scope_0_loops [LATCH]"));
    assert!(
        output.contains("scope_0_loops [LATCH]: if (promoter_sigma70*x1, RBS*x1)"),
        "{}",
        output
    );
}

#[test]
fn test_03_polycistronic_operon() {
    let seq = format!(
        "CCCC{}AAAA{}ATG{}TAA{}ATG{}TAA{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        RBS_SPACER,
        orf_body(60),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("polycistronic", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"));
    assert!(
        stdout.contains("scope_0: if (RBS*x1) { emit(gene_1); }"),
        "STDOUT:\n{}\nSTDERR:\n{}",
        stdout,
        stderr
    );
}

#[test]
fn test_04_tgatg_coupled_operon() {
    // Gene 0 ends with TGA overlapping Gene 1's ATG by 4 bp (TGATG)
    let seq = format!(
        "CCCC{}AAAA{}ATG{}TGATG{}TAA{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        orf_body(16),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("tgatg_coupled", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"));
    assert!(stdout.contains("scope_0: if (coupled(gene_0)) { emit(gene_1); }"));
}

#[test]
fn test_05_short_orf_above_physical_minimum_is_retained() {
    let seq = format!(
        "CCCC{}AAAA{}GTG{}TAACCCC{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(16),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("micro_gtg", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("emit(gene_"));
}

#[test]
fn test_06_circular_split_plasmid() {
    // Promoter and start near end of linear string; ORF & terminator wrap past index 0
    let tail = format!("CCCC{}AAAA{}ATG{}", PROM_S70, RBS_SPACER, orf_body(20));
    let head = format!("{}TAACCCC{}CCCCCCCCCC", orf_body(60), TERM_HAIRPIN);
    let seq = format!("{}{}", head, tail);
    let fasta = TempFasta::new("circular_split", &seq);

    let (ok_lin, stdout_lin, _) = run_dnadec(&fasta, &[]);
    assert!(ok_lin);
    assert!(
        stdout_lin.trim().is_empty(),
        "linear mode should not see wrapped gene"
    );

    let (ok_circ, stdout_circ, stderr_circ) = run_dnadec(&fasta, &["--circular"]);
    assert!(ok_circ, "stderr: {}", stderr_circ);
    assert!(stdout_circ.contains("emit(gene_0)"), "{}", stdout_circ);
}

#[test]
fn test_07_random_noise_is_not_emitted_as_a_confident_gene() {
    let at_noise = deterministic_noise(5000, false, 101);
    let gc_noise = deterministic_noise(3000, true, 202);

    let f_at = TempFasta::new("at_noise", &at_noise);
    let f_gc = TempFasta::new("gc_noise", &gc_noise);

    let (ok1, out1, err1) = run_dnadec(&f_at, &[]);
    let (ok2, out2, err2) = run_dnadec(&f_gc, &[]);
    assert!(out1.trim().is_empty(), "AT noise output:\n{}", out1);
    assert!(out2.trim().is_empty(), "GC noise output:\n{}", out2);
    if !ok1 {
        assert!(err1.contains("missing downstream terminator"), "{}", err1);
    }
    if !ok2 {
        assert!(err2.contains("missing downstream terminator"), "{}", err2);
    }
}

#[test]
fn test_08_strict_vs_non_strict_unterminated() {
    let seq = format!(
        "CCCC{}AAAA{}ATG{}TAACCCCCCCCCC",
        PROM_S70,
        RBS_SPACER,
        orf_body(80)
    );
    let fasta = TempFasta::new("unterminated", &seq);

    let (ok_strict, stdout_strict, stderr_strict) = run_dnadec(&fasta, &[]);
    assert!(!ok_strict && stdout_strict.trim().is_empty());
    assert!(stderr_strict.contains("missing downstream terminator"));

    let (ok_nonstrict, stdout_nonstrict, _) = run_dnadec(&fasta, &["--non-strict"]);
    assert!(ok_nonstrict);
    assert!(stdout_nonstrict.contains("[WARN: UNTERMINATED_READTHROUGH]"));
}

#[test]
fn test_09_sigma32_and_leaderless_tss() {
    let s32_seq = format!(
        "CCCC{}AAAA{}ATG{}TAA{}",
        PROM_S32,
        RBS_SPACER,
        orf_body(70),
        TERM_HAIRPIN
    );
    let f_s32 = TempFasta::new("sigma32", &s32_seq);
    let (ok1, out1, stderr1) = run_dnadec(&f_s32, &["--non-strict"]);
    assert!(
        ok1 && out1.contains("promoter_sigma32*x1"),
        "stdout: {}\n        stderr: {}",
        out1,
        stderr1
    );

    // Leaderless: ATG sits 7 bp downstream of -10 box with zero RBS
    let leaderless_seq = format!(
        "CCCC{}AAAAAAAATG{}TAA{}",
        PROM_S70,
        orf_body(70),
        TERM_HAIRPIN
    );
    let f_ll = TempFasta::new("leaderless", &leaderless_seq);
    let (ok2, out2, err2) = run_dnadec(&f_ll, &["--non-strict"]);
    assert!(
        ok2 && out2.contains("leaderless_TSS+1*x1"),
        "stdout: {}\nstderr: {}",
        out2,
        err2
    );
}

#[test]
fn test_10_multi_target_emitters() {
    let seq = format!(
        "CCCC{}AAAA{}ATG{}TAA{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("emit_targets", &seq);

    let (ok_ir, out_ir, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "ir"]);
    let (ok_py, out_py, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "python"]);
    let (ok_v, out_v, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "verilog"]);

    assert!(ok_ir && out_ir.contains("IR_MODULE") && out_ir.contains("DEPENDS_ALL"));
    assert!(ok_py && out_py.contains("if promoter_sigma70 and rbs:"));
    assert!(ok_v && out_v.contains("assign gene_0 = promoter_sigma70 & rbs;"));
}

#[test]
fn test_11_reverse_strand_operon_dependencies() {
    let oriented = format!(
        "CCCC{}{}{}ATG{}TAA{}ATG{}TAA{}",
        PROM_S70,
        OPERATOR,
        RBS_SPACER,
        orf_body(60),
        RBS_SPACER,
        orf_body(60),
        TERM_HAIRPIN
    );
    let seq = reverse_complement(&oriented);
    let fasta = TempFasta::new("reverse_operon", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);

    assert!(ok, "stderr: {}", stderr);
    assert_eq!(stdout.matches("emit(gene_").count(), 2, "{}", stdout);
    assert!(stdout.contains("(-)"), "{}", stdout);
    let operator_site = reverse_complement(OPERATOR);
    let site_start = seq.find(&operator_site).unwrap();
    assert!(
        stdout.contains(&format!(
            "CANDIDATE_OPERATOR_SITE gene_1 motif={} span={}..{} strand=-",
            OPERATOR,
            site_start,
            site_start + operator_site.len()
        )),
        "{}",
        stdout
    );
}

#[test]
fn test_12_readme_toggle_reports_unresolved_sites_without_inventing_edges() {
    let seq = format!(
        "CCCC{}{}{}{}{}ATG{}TAA{}CCCC{}{}{}ATG{}TAA{}",
        PROM_S70,
        OPERATOR,
        OPERATOR,
        OPERATOR,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN,
        PROM_S70,
        OTHER_OPERATOR,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("readme_toggle", &seq);

    let (ok_default, default, stderr) = run_dnadec(&fasta, &["--non-strict"]);
    assert!(ok_default, "stderr: {}", stderr);
    assert!(
        default.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"),
        "{}",
        default
    );
    assert!(
        default.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_1); }"),
        "{}",
        default
    );
    assert!(
        default.contains(&format!(
            "[CANDIDATE_CIS_EDGES: operator_site({}*x3)",
            OPERATOR
        )),
        "{}",
        default
    );
    assert!(default.contains(&format!(
        "[UNLINKED_OPERATOR_LOCKS: !lock_0({}*x3)",
        OPERATOR
    )));
    assert!(!default.contains("REGULATORY_EDGE"));
    assert!(!default.contains("[LATCH]"));
    assert!(default.contains("!lock_0("));
    let first_site = seq.find(OPERATOR).unwrap();
    assert!(default.contains(&format!(
        "CANDIDATE_OPERATOR_SITE gene_0 motif={} span={}..{} strand=+",
        OPERATOR,
        first_site,
        first_site + OPERATOR.len()
    )));

    let (ok_error, error_mode, stderr) =
        run_dnadec(&fasta, &["--non-strict", "--uncertainty", "error"]);
    assert!(!ok_error);
    assert!(error_mode.is_empty());
    assert!(stderr.contains("unresolved operator-site candidate"));

    let (ok_ir, ir, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "ir"]);
    assert!(ok_ir);
    assert!(ir.contains("CALLBACK_SCOPE_0:"));
    assert!(ir.contains(&format!("CANDIDATE_OPERATOR_SITE gene_0 motif={}", OPERATOR)));
    assert!(!ir.contains("NOT lock("));
    assert!(ir.contains(&format!("UNLINKED_OPERATOR_LOCKS [ !lock_0({}*x3)", OPERATOR)));

    let (ok_py, python, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "python"]);
    assert!(ok_py);
    assert!(!python.contains("while circuit_active:"));
    assert!(python.contains("if promoter_sigma70 and rbs:"));
    assert!(!python.contains("not lock_"));
    assert!(python.contains(&format!("# UNLINKED_OPERATOR_LOCKS: !lock_0({}*x3)", OPERATOR)));

    let (ok_verilog, verilog, _) = run_dnadec(&fasta, &["--non-strict", "--emit", "verilog"]);
    assert!(ok_verilog);
    assert!(verilog.contains("assign gene_0 = promoter_sigma70 & rbs;"));
    assert!(!verilog.contains("latch feedback"));
    assert!(!verilog.contains("lock_GATATATATATC"));
    assert!(verilog.contains(&format!("// UNLINKED_OPERATOR_LOCKS: !lock_0({}*x3)", OPERATOR)));
}

#[test]
fn test_13_operator_candidate_in_upstream_intergenic_slice() {
    let seq = format!(
        "CCCC{}{}AAAA{}ATG{}TAA{}",
        OPERATOR,
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("upstream_operator", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);

    assert!(ok, "stderr: {}", stderr);
    assert!(
        stdout.contains(&format!(
            "CANDIDATE_OPERATOR_SITE gene_0 motif={} span=4..{} strand=+",
            OPERATOR,
            4 + OPERATOR.len()
        )),
        "{}",
        stdout
    );
}

#[test]
fn test_14_operator_candidates_follow_intergenic_boundaries_not_an_arbitrary_window() {
    let far_upstream = format!("CCCCC{}{}", OPERATOR, "C".repeat(150));
    let seq = format!(
        "{}{}AAAA{}ATG{}TAA{}",
        far_upstream,
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("distant_operator", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &["--non-strict"]);

    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("CANDIDATE_OPERATOR_SITE"), "{}", stdout);
}

#[test]
fn test_15_checked_in_toggle_fixture_matches_all_emitters_in_strict_mode() {
    let fixture = fs::read_to_string("fixtures/toggle_switch.fasta")
        .expect("failed to read the checked-in toggle fixture");
    let sequence: String = fixture
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let fasta = TempFasta::new("readme_toggle_fixture", &sequence);

    let (ok_default, default, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok_default, "stderr: {}", stderr);
    assert!(default.contains("scope_0_loops [LATCH]: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"), "{}", default);
    assert!(default.contains("scope_1: if (promoter_sigma70*x1, RBS*x1) { emit(gene_1); }"), "{}", default);
    assert!(default.contains("REGULATORY_EDGE gene_0 ? gene_1"));
    assert!(default.contains("effect=unknown"));
    assert!(default.contains(&format!("CANDIDATE_OPERATOR_SITES: {}*x3", OPERATOR)));
    assert!(default.contains(&format!("CANDIDATE_OPERATOR_SITES: {}*x1", OPERATOR)));

    let (ok_ir, ir, _) = run_dnadec(&fasta, &["--emit", "ir"]);
    assert!(ok_ir);
    assert!(ir.contains("CALLBACK_SCOPE_0_LOOPS [LATCH]:"));
    assert!(ir.contains("CALLBACK_SCOPE_1:"));
    assert!(ir.contains("DEPENDS_ALL [ promoter_sigma70(x1), RBS(x1) ]"));

    let (ok_python, python, _) = run_dnadec(&fasta, &["--emit", "python"]);
    assert!(ok_python);
    assert!(python.contains("while circuit_active:  # scope_0_loops [LATCH]"));
    assert!(python.contains("if promoter_sigma70 and rbs:"));
    assert!(!python.contains("not lock_"));

    let (ok_verilog, verilog, _) = run_dnadec(&fasta, &["--emit", "verilog"]);
    assert!(ok_verilog);
    assert!(verilog.contains("assign gene_0 = promoter_sigma70 & rbs;"));
    assert!(verilog.contains("latch feedback"));
    assert!(verilog.contains("// CANDIDATE_OPERATOR_SITES gene_0"));

    let (ok_unsupported, _, stderr) = run_dnadec(&fasta, &["--motifs", "not-a-supported-flag"]);
    assert!(!ok_unsupported);
    assert!(stderr.contains("unknown flag '--motifs'"));
}
