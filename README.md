**Author:** Oladapo Oyelaja  
**Date:** 25 November 2026  

# `dnadec`: A Deterministic DNA-to-IR Prototype

Look at a raw `.fasta` file and you see an unbroken wall of four letters: `A`, `T`, `C`, and `G`. Many established gene-prediction pipelines use statistical models to rank candidate features. `dnadec` explores an alternative deterministic sequence-rule approach—reducing a DNA strand to ordered spatial intervals and lifting sequence gates into an intermediate representation (IR) without trained weights or species databases.

The compiler pass-through is a chain of physical and structural sequence gates: promoter evidence opens a transcriptional unit; an initiation signal admits a start codon; an in-frame stop defines a callback; a terminator closes the unit. Later passes consume only evidence that survived earlier gates. `dnadec` is an experimental biophysical static analyzer, not a validated whole-genome annotator.

---

# I. Architecture and Primitives

`dnadec` operates in pure Rust `std` with zero external dependencies across linear sequence passes (`O(N + K log K)` time, where `N` is sequence length and `K ≪ N` is candidate frame count; `O(N)` space):

- **Multi-Shuffle Null Calibration (Altschul–Erickson Dinucleotide Shuffle):** Pass 0 constructs a 4-node Eulerian graph (`A, C, G, T`) with a locked final-edge spanning tree and runs `10` deterministic dinucleotide-preserving shuffles (`seed = 42`) to establish genome-specific background floors for promoter match counts, RBS match counts, and hairpin stem free energies. Inputs too short to form a valid shuffle or yielding zero null candidates are explicitly reported as uncalibrated (`None`) rather than falling back to hardcoded floors.
- **Promoter and Initiation Gates:** σ70 (`TTGACA` / `TATAAT`, spacer `16..=18 bp`) and σ32 (`CTTGAA` / `CCCCAT`, spacer `13..=17 bp`) boxes require at least `4/6` matches per box. Shine–Dalgarno RBS scoring requires at least `4/6` matches to `AGGAGG` at `4..=15 bp` upstream of a start codon (`ATG`/`GTG`/`TTG`). Callbacks require an in-frame `STOP` (`TAA`/`TAG`/`TGA`) and a minimum coding span of `30 bp` (10 codons), or `1..=4 bp` translational overlap (`coupled(gene_i)`), or leaderless spacing (`5..=9 bp` downstream of the `-10` box).
- **Transcriptional-Unit Closure:** Candidate rho-independent terminators evaluate inverted-repeat stems (`6..=12 bp`) and loops (`3..=8 bp`) using nearest-neighbor Watson–Crick base-pair stacking energies (`ΔG`) paired with a downstream poly-T tract (`≥ 4` Ts in an `8 bp` window). Default strict mode halts on unsupported read-through spans; `--non-strict` emits them with an explicit `[WARN: UNTERMINATED_READTHROUGH]` annotation.
- **Unwrapped Circular Topology (`--circular`):** Circular plasmids and chromosomes are unrolled once (`seq + seq`), keeping all interval overlap math strictly linear across passes and deferring `% orig_len` until final `@span` emission.
- **CSR Dependency Slicing & Scope Latches:** Pass 1 scans bounded promoter/5'-UTR neighborhoods (`-120..+25 bp`) for inverted-repeat operator candidates (`6..=10 bp` stems) and compresses identical adjacent hits into multiplicity terms (`*xN`) only when start coordinates lie within `35 bp`. Pass 2 traces shared local operator links (emitted with unknown regulatory sign unless resolved) and quarantines feedback cycles into `scope_N_loops [LATCH]`. Unlinked operator sites remain isolated as `!lock_N` guards.

---

# II. Case Study: Synthetic Toggle Fixture

Passes 0 through 2 lower sequence coordinates into a structured Intermediate Representation of callback dependencies and feedback cycles, supporting structured text (`default`), `ir`, `python`, and `verilog` emitters.

```bash
dnadec fixtures/toggle_switch.fasta
```

### 1. Lifted Scope / CSR Intermediate Representation (`--emit ir`)

```text
IR_MODULE fixtures/toggle_switch.fasta {
  REGULATORY_EDGE gene_0 -- gene_1 [motif=TTACTAGA*x1]
  REGULATORY_EDGE gene_1 -- gene_0 [motif=TCCCTATC*x3]
  CANDIDATE_OPERATOR_SITE gene_0 motif=TCCCTATC span=84..92 strand=+
  CANDIDATE_OPERATOR_SITE gene_0 motif=TCCCTATC span=92..100 strand=+
  CANDIDATE_OPERATOR_SITE gene_0 motif=TCCCTATC span=100..108 strand=+
  CALLBACK_SCOPE_0_LOOPS [LATCH]:
    CALLBACK gene_0 @span(170..872, '+')
      DEPENDS_ALL [ promoter_sigma70(x1), NOT lock(TCCCTATC)(x3), RBS(x1) ]
      EMIT gene_0
  CANDIDATE_OPERATOR_SITE gene_1 motif=TTACTAGA span=948..956 strand=+
  CALLBACK_SCOPE_1:
    CALLBACK gene_1 @span(1032..1734, '+')
      DEPENDS_ALL [ promoter_sigma70(x1), NOT lock(TTACTAGA)(x1), RBS(x1) ]
      EMIT gene_1
}
```

### 2. Multi-Target Emitter Examples

#### Default Structured Output (`--emit default`)
```rust
scope_0_loops [LATCH]: if (promoter_sigma70*x1, !lock(TCCCTATC)*x3, RBS*x1) { emit(gene_0); } // @span: 170..872 (+)
scope_1: if (promoter_sigma70*x1, !lock(TTACTAGA)*x1, RBS*x1) { emit(gene_1); } // @span: 1032..1734 (+)
```

#### Python-Style Output (`--emit python`)
```python
# Trans-decompiled by dnadec from fixtures/toggle_switch.fasta
while circuit_active:  # scope_0_loops [LATCH]
    if promoter_sigma70 and not lock_TCCCTATC >= 3 and rbs:  # @span: 170..872 (+)
        gene_0 = emit("gene_0")
    if promoter_sigma70 and not lock_TTACTAGA >= 1 and rbs:  # @span: 1032..1734 (+)
        gene_1 = emit("gene_1")
```

#### Verilog-Style Output (`--emit verilog`)
```verilog
// Trans-decompiled by dnadec from fixtures/toggle_switch.fasta
assign gene_0 = promoter_sigma70 & ~(lock_TCCCTATC_x3) & rbs; // @span: 170..872 (+) // latch feedback
assign gene_1 = promoter_sigma70 & ~(lock_TTACTAGA_x1) & rbs; // @span: 1032..1734 (+)
```
---

# III. Emergent Circuit Logic

- **Implicit Boolean Conjunction (`AND`) and Steric Negation (`NOT`):** Dependencies surviving the sequence gates within the same 5' CSR slice form a conjunction (`AND`), while candidate operator palindromes emit as steric `!lock(...)` guards.
- **Spatial Tandem Repeat Compression:** Identical adjacent operator hits within `35 bp` collapse into a single multiplicity term (`!lock(TCCCTATC)*x3`), while separated sites remain distinct entries.
- **Feedback Cycle Quarantine (`[LATCH]`):** During Pass 2's backward depth-first traversal, an `active_path` recursion stack detects bidirectional feedback cycles and isolates them into `scope_0_loops [LATCH]`.
- **Operon & Coupling Propagation:** Downstream callbacks inside polycistronic transcripts sharing one upstream promoter and one downstream terminator emit clean `if (RBS*x1)` or `if (coupled(gene_0))` guards without requiring internal promoters.

---

# IV. Benchmarks, Scaling, and Usage

Built in pure Rust `std` with zero external dependencies, `dnadec` is evaluated against both synthetic constructs and the full *Escherichia coli* K-12 reference genome (`NC_000913.3`, `4,641,652 bp`):

- **Wild-genome benchmark (`benchmarks/regulondb_eval.py`):** On the `4,641,652 bp` *E. coli* K-12 chromosome, the release build achieved a log-log linear scaling exponent of **`0.981`** with a median runtime of **`25.19 s`** in non-strict mode (including `10` full-genome dinucleotide shuffles). Strict mode exited with status `1` due to unclosed read-through spans. In non-strict mode (5,594 predicted ORFs vs. 4,435 RegulonDB references):
  - **Gene detection (`IoU ≥ 0.50`):** **26.86% F1** (`24.08%` precision, `30.37%` recall; `1,347` matched genes)
  - **Gene detection (`IoU ≥ 0.90`):** **19.36% F1** (`17.36%` precision, `21.89%` recall; `971` matched genes)
  - **Gene detection (`IoU = 1.00`):** **17.29% F1** (`15.50%` precision, `19.55%` recall; `867` exact matches)
  - **Operator-site matching (`within 10 bp`):** **2.50% F1** (`4.61%` precision, `1.71%` recall; `1,020` predicted sites vs. `2,742` confirmed/strong references)
  - **Signed regulatory edges:** **0** signed edges predicted against `4,057` confirmed/strong references (unlinked operators remain abstained `!lock_N` guards).
- **Synthetic benchmark (`benchmarks/synthetic_eval.py`):** Across `80` held-out synthetic constructs (`40` positive, `40` negative), `dnadec` achieved **59.79% gene F1** (`50.88%` precision, `72.50%` recall) and **72.44% operator-site F1** (`95.83%` precision, `58.23%` recall).

### Quick Start

```bash
git clone https://github.com/umayl/dnadec.git
cd dnadec
cargo build --release

# Decompile a linear circuit (Strict Mode by default)
./target/release/dnadec fixtures/toggle_switch.fasta

# Decompile a circular plasmid (handles genes split across index 0)
./target/release/dnadec plasmid.fasta --circular

# Emit IR, Python, or Verilog HDL
./target/release/dnadec fixtures/toggle_switch.fasta --emit ir
./target/release/dnadec fixtures/toggle_switch.fasta --emit python
./target/release/dnadec fixtures/toggle_switch.fasta --emit verilog

# Run synthetic and RegulonDB genome benchmarks
python benchmarks/synthetic_eval.py
python benchmarks/regulondb_eval.py

# Permit Rho-terminated or partial contigs with inline warnings
./target/release/dnadec wild_contig.fasta --non-strict
```
---

# V. Biological Boundaries and Future Work

Because `dnadec` relies strictly on 1D sequence-intrinsic rules without external proteome context, clear physical boundaries define its current scope:

- **Protein-Assisted & Non-Canonical Regulation:** Wild bacterial operons relying on external activator proteins (`CRP`, `AraC`) to rescue weak `-35` boxes, non-Shine–Dalgarno initiation, or **Rho helicase** termination without a hairpin cannot be resolved from 1D sequence geometry alone.
- **Trans-Acting Repressor Identity:** While `dnadec` identifies inverted-repeat operator candidates (`!lock(...)`), linking an unannotated palindrome to a trans-acting repressor encoded elsewhere on the chromosome—or determining activation vs. repression sign—cannot be established from a 1D DNA palindrome alone; unknown sites remain isolated as unlinked `!lock_N` guards.
- **Analytical Null Calibration:** Replacing the `10`-sample Monte Carlo dinucleotide shuffle with an exact first-order Markov transition expectation (`4×4` dinucleotide matrix) is a natural next step to eliminate sampling variance and cut full-genome runtime to under a second.
- **Prokaryotic Scope:** Eukaryotic genomes—characterized by intron splicing, chromatin accessibility, and distal enhancers—fall outside linear prokaryotic coordinate scanning.

---

# VI. Conclusion

`dnadec` demonstrates how deterministic 1D sequence intervals can be lifted into a structured CSR/scope intermediate representation and lowered into multi-target circuit code without trained statistical models—while empirically mapping the boundary where 1D sequence physics stops and 3D cellular/proteome context begins.

---

# Acknowledgements & Conceptual Independence

The four-pass biophysical decompilation pipeline, transcriptional-unit closure rules, and CSR-to-scope control-flow mapping presented in `dnadec` were conceptualized and engineered independently from first principles in zero-dependency standard-library Rust.