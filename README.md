**Author:** Oladapo Oyelaja  
**Date:** 24 November 2026  

# A Deterministic Biophysical Trans-Decompiler for Prokaryotic Genomes (`dnadec`)

Look at a raw `.fasta` file and you see an unbroken wall of four letters: `A`, `T`, `C`, and `G`. For decades, computational biology has treated that string as an evolutionary puzzle best solved through probability—training Hidden Markov Models, codon-bias matrices, and deep neural networks to infer where genes sit. Those tools output flat coordinate tables, telling you *where* open reading frames are located while leaving the underlying program completely opaque.

Strip away the Latin nomenclature and reduce a DNA strand to its bare physical events along a coordinate axis: `x1, x2_start, x2_end, x3_start, x4_end`. When you examine a prokaryotic genome purely as ordered spatial intervals governed by nearest-neighbor base-pairing thermodynamics (ΔG), a structured hierarchy emerges directly from the sequence. From those flat coordinates alone, you can deterministically derive the **unknowns** (upstream promoters, Shine–Dalgarno sites, and operator palindromes), isolate the **callbacks** (expressed reading frames), trace their **dependencies** (bounded intergenic switch slices), and reconstruct the **broader control flow** of the circuit—nested conditional guards, polycistronic operons, and local or motif-linked feedback latches.

Once higher-level control flow falls straight out of raw coordinate intervals without trained statistical models or arbitrary energy thresholds, the abstraction becomes clear: a prokaryotic genome behaves as **thermodynamic bytecode**—and lifting its static logic does not require machine learning. It requires a decompiler.

---

# I. Architecture and Primitives

`dnadec` lifts raw nucleotide sequences into a structured Intermediate Representation (IR) across four linear-scan passes (`O(N + K log K)` time, where `N` is sequence length and `K ≪ N` is candidate frame count; `O(N)` space), deriving energy baselines directly from sequence physics:

- **Self-Calibrating Noise Floor (Altschul–Erickson Dinucleotide Shuffle):** Rather than relying on hardcoded ΔG cutoffs that break across AT-rich and GC-rich organisms, Pass 0 constructs a 4-node Eulerian graph (`A, C, G, T`) with a locked final-edge spanning tree and deterministically shuffles the sequence (`seed = 42`). Scoring the physical motif scanners against this chemically identical decoy derives the genome's own ΔG background ceiling automatically.
- **Transcriptional Units (`[Promoter .. Terminator]`):** A promoter (σ70 or σ32) opens a continuous transcriptional span that remains active across downstream stop codons until closed by a rho-independent hairpin terminator (`stem: 6..=12 bp`, `loop: 3..=8 bp`, trailing `TTTT` tract).
  - In **Default (Strict) Mode**, any promoter lacking a closing downstream terminator fails compilation as an unterminated read-through error.
  - Under **`--non-strict`**, open or Rho-terminated transcripts compile with an explicit `[WARN: UNTERMINATED_READTHROUGH]` annotation.
- **Callbacks & Translational Coupling:** Inside an open transcriptional unit, `START` (`ATG`/`GTG`/`TTG`) → `STOP` (`TAA`/`TAG`/`TGA`) callbacks initiate via either a Shine–Dalgarno `RBS` (`4..=15 bp` upstream), a leaderless `+1` transcription start (`5..=9 bp` downstream of the `-10` box), or a `1..=4 bp` `TGATG` translational overlap (`coupled(gene_i)`).
- **Unwrapped Circular Topology (`--circular`):** Circular plasmids and chromosomes are unrolled once (`seq + seq`), keeping all interval overlap math strictly linear across passes and deferring `% orig_len` until final `@span` emission.
- **CSR Dependency Slicing & Scope Latches:** Pass 1 extracts the upstream unknowns (`promoter`, `!lock_N` operator palindromes, and `RBS`) bounding each winning callback into a flat Compressed Sparse Row (`CSR`) dependency array (`callback_dependencies`) and collapses tandem repeats into multiplicity factors (`*xN`). Pass 2 walks backward from terminal callbacks (`callback_scope_0`) to layer known regulatory edges (`callback_scope_1..N`) and quarantines active-path feedback cycles into `callback_scope_N_loops [LATCH]`.

---

# II. Case Study: Gardner Bistable Toggle Switch (Trans-Decompilation)

Because Passes 0 through 2 lower raw nucleotide coordinates into a **domain-agnostic Intermediate Representation (IR)**—consisting of flat CSR dependency slices, topological callback scopes, and quarantined feedback latches—`dnadec` operates as a **trans-decompiler**. Once the `.fasta` input is lifted into the Scope/CSR IR, the final emission pass can target any imperative, hardware-description, or formal logic language without re-scanning the genome.

Given a raw `.fasta` sequence encoding a mutual-repression toggle switch (`gene_0` repressing `gene_1`, and `gene_1` repressing `gene_0` via a 3x tandem operator array):

```bash
dnadec fixtures/toggle_switch.fasta
```

### 1. Lifted Scope / CSR Intermediate Representation (`--emit ir`)

```text
IR_MODULE toggle_switch {
  CALLBACK_SCOPE_0_LOOPS [LATCH]:
    CALLBACK gene_0 @span(170..872, '+')
      DEPENDS_ALL [ promoter_sigma70(x1), NOT lock(TCCCTATC)(x3), RBS(x1) ]
      EMIT gene_0

  CALLBACK_SCOPE_1:
    CALLBACK gene_1 @span(1032..1734, '+')
      DEPENDS_ALL [ promoter_sigma70(x1), NOT lock(TTACTAGA)(x1), RBS(x1) ]
      EMIT gene_1
}
```

### 2. Multi-Target Code Emission (Trans-Decompilation)

From that single lifted IR, backend emitters lower the genetic circuit directly into arbitrary target languages:

#### Default Structured Output (`--emit default`)
```rust
scope_0_loops [LATCH]: if (promoter_sigma70*x1, !lock(TCCCTATC)*x3, RBS*x1) { emit(gene_0); } // @span: 170..872 (+)
scope_1: if (promoter_sigma70*x1, !lock(TTACTAGA)*x1, RBS*x1) { emit(gene_1); } // @span: 1032..1734 (+)
```

#### Target A: Executable Python Simulation (`--emit python`)
```python
# Trans-decompiled from toggle_switch.fasta
while circuit_active:  # scope_0_loops [LATCH]
    if promoter_sigma70 and (not lock_TCCCTATC >= 3) and rbs:
        gene_0 = emit("gene_0")  # @span: 170..872 (+)
    if promoter_sigma70 and (not lock_TTACTAGA) and rbs:
        gene_1 = emit("gene_1")  # @span: 1032..1734 (+)
```

#### Target B: Hardware Description (`--emit verilog`)
```verilog
// Trans-decompiled hardware gate equivalent
assign gene_0 = promoter_sigma70 & ~(lock_TCCCTATC_x3) & rbs; // latch feedback
assign gene_1 = promoter_sigma70 & ~(lock_TTACTAGA)    & rbs;
```
---

# III. Emergent Circuit Logic

Several emergent compiler behaviors appear directly in the decompiled output above:

- **Implicit Boolean Conjunction (`AND`) and Steric Negation (`NOT`):** Every unknown residing within the same 5' UTR CSR slice must be satisfied simultaneously for RNA polymerase and the ribosome to trigger the downstream callback, naturally forming an `AND` guard where bound repressor palindromes act as steric `NOT` gates (`!lock(...)`).
- **Tandem Repeat Compression:** Adjacent identical palindromes upstream of `gene_0` are collapsed during Pass 1's sliding-window deduplication into a single cooperative multiplicity term (`!lock(TCCCTATC)*x3`), reducing CSR memory footprint without losing binding stoichiometry.
- **Feedback Cycle Quarantine (`[LATCH]`):** During Pass 2's backward depth-first traversal, maintaining an `active_path` recursion stack detects bidirectional repression cycles and isolates them into `scope_0_loops [LATCH]` rather than emitting an invalid acyclic tree.
- **Operon & Coupling Propagation:** In polycistronic transcripts sharing one upstream promoter and one downstream terminator, downstream callbacks emit clean `if (RBS*x1)` or `if (coupled(gene_0))` guards without requiring synthetic internal promoters.

---

# IV. Complexity, Scaling, and Usage

Built in pure Rust `std` with zero external dependencies, `dnadec` bridges human-readable circuit inspection and high-throughput static analysis:

- **Linear Sequence Scaling (`O(N)` Space, `O(N + K log K)` Time):** Every pass executes via linear scans, bounded physical windows, sweep-line interval sorting over `K` candidates, and contiguous CSR arrays (~80–110 bytes/bp), allowing full bacterial chromosomes to compile in tens of milliseconds on commodity hardware.
- **Zero Trained Weights or Species Databases:** Because thermodynamic stability (ΔG) and dinucleotide self-calibration derive directly from the input sequence, `dnadec` analyzes novel bacteriophages, unannotated metagenomic contigs, and synthetic plasmids without pre-trained HMMs.
- **Compile-Time Circuit Linting:** Treating unclosed transcripts as hard errors by default brings static verification discipline to synthetic plasmid design prior to DNA synthesis.

### Quick Start

```bash
git clone https://github.com/umayl/dnadec.git
cd dnadec
cargo build --release

# Decompile a linear circuit (Strict Mode by default)
./target/release/dnadec fixtures/toggle_switch.fasta

# Decompile a circular plasmid (handles genes split across index 0)
./target/release/dnadec plasmid.fasta --circular

# Trans-decompile to IR, Python, or Verilog HDL
./target/release/dnadec fixtures/toggle_switch.fasta --emit ir
./target/release/dnadec fixtures/toggle_switch.fasta --emit python
./target/release/dnadec fixtures/toggle_switch.fasta --emit verilog

# Permit Rho-terminated or partial contigs with inline warnings
./target/release/dnadec wild_contig.fasta --non-strict
```

---

# V. Biological Boundaries and Future Work

Because `dnadec` relies strictly on sequence-intrinsic physical thermodynamics without external proteome context, several biological boundaries define its current scope:

- **Protein-Assisted Regulation:** Wild bacterial operons that rely on external activator proteins (`CRP`, `AraC`) to rescue weak `-35` promoter boxes, or **Rho helicase** to terminate unstructured transcripts without a hairpin, cannot be resolved from DNA folding energy alone without `--non-strict` mode.
- **Trans-Acting Repressor Resolution:** While `dnadec` deterministically identifies operator palindromes (`!lock(...)`), linking an unannotated palindrome to a specific trans-acting repressor protein encoded elsewhere on the chromosome requires known motif dictionaries or structural protein–DNA docking; unknown operators remain isolated as unlinked steric guards (`!lock_N`).
- **Prokaryotic Scope:** Eukaryotic genomes—characterized by intron splicing, chromatin accessibility, and distal enhancers tens of kilobases away—fall outside linear prokaryotic coordinate scanning.

---

# VI. Conclusion

Treating DNA solely through opaque statistical models obscures the clean, deterministic control flow governing prokaryotic genetic circuits. By reducing raw nucleotide strings to physical coordinate intervals and deriving unknowns, callbacks, dependencies, and feedback latches through nearest-neighbor thermodynamics and self-calibrating dinucleotide backgrounds, `dnadec` lifts prokaryotic genomes into a universal intermediate representation ready for multi-target trans-decompilation in linear time.

---

# Acknowledgements & Conceptual Independence

The four-pass biophysical decompilation pipeline, transcriptional-unit closure rules, and CSR-to-scope control-flow mapping presented in `dnadec` were conceptualized and engineered independently from first principles. Rather than adapting existing statistical gene-prediction heuristics or machine-learning classifiers, the system was built from scratch in standard-library Rust by deriving control-flow structure directly from physical sequence intervals.
