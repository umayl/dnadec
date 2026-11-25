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

const PROM_S70: &str = "TTGACAAAAAAAAAAAAAAAAAATATAAT";
const PROM_S32: &str = "CTTGAAGGGGGGGGGGGGGGGCCCCAT";
const RBS_SPACER: &str = "AGGAGGAAAAAAA";
const TERM_HAIRPIN: &str = "GCGCGCGCGAAAACGCGCGCTCTTTTTT";

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
    assert!(stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"));
}

#[test]
fn test_02_gardner_toggle_switch_latch() {
    // Gene 0 has 3x tandem TCCCTATC operator; Gene 1 has TTACTAGA operator
    let seq = format!(
        "CCCC{}TCCCTATCTCCCTATCTCCCTATC{}ATG{}TAA{}CCCC{}TTACTAGA{}ATG{}TAA{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN,
        PROM_S70,
        RBS_SPACER,
        orf_body(80),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("toggle_switch", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("scope_0_loops [LATCH]: if (promoter_sigma70*x1, !lock(TCCCTATC)*x3, RBS*x1) { emit(gene_0); }"));
    assert!(stdout.contains("scope_1: if (promoter_sigma70*x1, !lock(TTACTAGA)*x1, RBS*x1) { emit(gene_1); }"));
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
    assert!(stdout.contains("scope_0: if (RBS*x1) { emit(gene_1); }"));
}

#[test]
fn test_04_tgatg_coupled_operon() {
    // Gene 0 ends with TGA overlapping Gene 1's ATG by 4 bp (TGATG)
    let seq = format!(
        "CCCC{}AAAA{}ATG{}TGATG{}TAA{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(60),
        orf_body(60),
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("tgatg_coupled", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("scope_0: if (promoter_sigma70*x1, RBS*x1) { emit(gene_0); }"));
    assert!(stdout.contains("scope_0: if (coupled(gene_0)) { emit(gene_1); }"));
}

#[test]
fn test_05_micro_peptide_and_gtg_start() {
    let seq = format!(
        "CCCC{}AAAA{}GTG{}TAACCCC{}",
        PROM_S70,
        RBS_SPACER,
        orf_body(16), // 54 bp total ORF
        TERM_HAIRPIN
    );
    let fasta = TempFasta::new("micro_gtg", &seq);
    let (ok, stdout, stderr) = run_dnadec(&fasta, &[]);
    assert!(ok, "stderr: {}", stderr);
    assert!(stdout.contains("emit(gene_0)"));
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
    assert!(stdout_lin.trim().is_empty(), "linear mode should not see wrapped gene");

    let (ok_circ, stdout_circ, stderr_circ) = run_dnadec(&fasta, &["--circular"]);
    assert!(ok_circ, "stderr: {}", stderr_circ);
    assert!(stdout_circ.contains("emit(gene_0)"));
}

#[test]
fn test_07_random_noise_zero_false_positives() {
    let at_noise = deterministic_noise(5000, false, 101);
    let gc_noise = deterministic_noise(3000, true, 202);

    let f_at = TempFasta::new("at_noise", &at_noise);
    let f_gc = TempFasta::new("gc_noise", &gc_noise);

    let (ok1, out1, _) = run_dnadec(&f_at, &[]);
    let (ok2, out2, _) = run_dnadec(&f_gc, &[]);
    assert!(ok1 && out1.trim().is_empty(), "false positive in AT noise: {}", out1);
    assert!(ok2 && out2.trim().is_empty(), "false positive in GC noise: {}", out2);
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

    let (ok_strict, _, stderr_strict) = run_dnadec(&fasta, &[]);
    assert!(!ok_strict, "strict mode must fail on unterminated transcript");
    assert!(stderr_strict.contains("Unterminated transcript"));

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
    let (ok1, out1, _) = run_dnadec(&f_s32, &[]);
    assert!(ok1 && out1.contains("promoter_sigma32*x1"));

    // Leaderless: ATG sits 7 bp downstream of -10 box with zero RBS
    let leaderless_seq = format!(
        "CCCC{}AAAAAAAATG{}TAA{}",
        PROM_S70,
        orf_body(70),
        TERM_HAIRPIN
    );
    let f_ll = TempFasta::new("leaderless", &leaderless_seq);
    let (ok2, out2, _) = run_dnadec(&f_ll, &[]);
    assert!(ok2 && out2.contains("leaderless_TSS+1*x1"));
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

    let (ok_ir, out_ir, _) = run_dnadec(&fasta, &["--emit", "ir"]);
    let (ok_py, out_py, _) = run_dnadec(&fasta, &["--emit", "python"]);
    let (ok_v, out_v, _) = run_dnadec(&fasta, &["--emit", "verilog"]);

    assert!(ok_ir && out_ir.contains("IR_MODULE") && out_ir.contains("DEPENDS_ALL"));
    assert!(ok_py && out_py.contains("if promoter_sigma70 and rbs:"));
    assert!(ok_v && out_v.contains("assign gene_0 = promoter_sigma70 & rbs;"));
}
