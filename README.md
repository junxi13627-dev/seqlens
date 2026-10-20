# SeqLens 🧬

![Build](https://github.com/junxi13627-dev/seqlens/actions/workflows/build.yml/badge.svg)
![Java](https://img.shields.io/badge/Java-17%2B-orange)
![License](https://img.shields.io/badge/license-MIT-blue)

**SeqLens** is a lightweight, dependency-free command-line toolkit for DNA and protein sequence analysis, written in pure Java.

It reads standard FASTA files and covers the everyday tasks of a molecular biology workflow: sequence statistics, reverse complements, transcription, translation, open reading frame (ORF) discovery and k-mer profiling.

## Features

| Command      | What it does |
|--------------|--------------|
| `stats`      | Length, A/C/G/T/N counts, GC content and estimated melting temperature (Tm) |
| `revcomp`    | Reverse complement, with full IUPAC ambiguity-code support |
| `transcribe` | DNA → mRNA |
| `translate`  | DNA → protein in any of the three forward reading frames (standard genetic code) |
| `orfs`       | Finds ORFs (ATG → stop) in all **six** reading frames, reported in forward-strand coordinates |
| `kmers`      | Counts overlapping k-mers and shows the most frequent ones |

- Multi-record, multi-line FASTA parsing (comments and blank lines handled)
- Input from a file, from standard input (`-`), or straight from the command line (`--seq`)
- No runtime dependencies — just a JDK
- Unit-tested with JUnit 5 and built automatically with GitHub Actions

## Quick start

Requirements: **Java 17+** and **Maven 3.8+**.

```bash
git clone https://github.com/junxi13627-dev/seqlens.git
cd seqlens
mvn package            # compiles, runs the tests, builds target/seqlens.jar
java -jar target/seqlens.jar stats examples/sample.fasta
```

## Usage examples

### Sequence statistics

```text
$ java -jar target/seqlens.jar stats examples/sample.fasta
id                       length        A        C        G        T        N      GC%     Tm(C)
contig_1                    681      181      154      169      177        0   47.43%      83.4
gc_rich_region              400       57      136      151       56        0   71.75%      92.6
fragment_with_gap           256       54       49       65       78       10   46.34%      81.2
```

### Find open reading frames

```text
$ java -jar target/seqlens.jar orfs --min-aa 45 examples/sample.fasta
id                   strand frame    start      end   aa_len  protein
contig_1                  +     1      121      381       86  MVEVEKHNCSSQNSLTLTVVTLDQSRIEVIHFARLTIPDMITAIWSTMGIRLELVYS...
contig_1                  -     1      472      621       49  MDRKTSWDSDSTQCYSCIITVSRASCSLYSNVPCYSYQYTSSIYLQYIN
contig_1                  +     2       20      166       48  MHQCGRVEQTPLETGLTIRSSVMKVTVRATFDSEWWKLKNITVLHRIL
gc_rich_region            -     3       15      197       60  MQRMRISGRRPSAGITPRRYSAAFDPDLPRVSRLGVGPNVHQGVRRVPRAPPTPSQVLGA

4 ORF(s) found with >= 45 amino acids.
```

`start` and `end` are 1-based, inclusive and always given on the forward strand (the stop codon is included), so ORFs on the `-` strand can be located directly in the original sequence.

### Translate a sequence

```text
$ java -jar target/seqlens.jar translate --seq ATGGCCATTGTAATGGGCCGCTGA
>input_frame1
MAIVMGR*
```

### Other commands

```bash
java -jar target/seqlens.jar revcomp examples/sample.fasta
java -jar target/seqlens.jar transcribe --seq ATGCTTAGC
java -jar target/seqlens.jar kmers --k 4 --top 5 examples/sample.fasta
cat my_genes.fasta | java -jar target/seqlens.jar stats -
```

Run `java -jar target/seqlens.jar --help` for all options.

## How it works

- **GC content** is calculated over unambiguous bases only (A/C/G/T), so `N` gaps don't skew the result.
- **Melting temperature** uses the Wallace rule (`2·(A+T) + 4·(G+C)`) for oligos under 14 nt and `64.9 + 41·(G+C−16.4)/N` for longer sequences. It's an estimate meant for quick primer checks, not a substitute for nearest-neighbor thermodynamics.
- **Translation** uses NCBI translation table 1 (the standard code). Stops appear as `*`, and codons containing ambiguous bases translate to `X`.
- **ORF finding** scans all six frames. Each ATG not already inside an ORF in the same frame starts a new ORF that runs to the first in-frame stop; ORFs with no stop codon before the sequence ends are skipped. Results are sorted longest first.

## Project structure

```
seqlens/
├── src/main/java/io/github/seqlens/
│   ├── Main.java          # command-line interface
│   ├── FastaReader.java   # FASTA parser
│   ├── FastaRecord.java   # one FASTA entry
│   ├── DnaTools.java      # GC content, reverse complement, Tm, k-mers
│   ├── Translator.java    # genetic code and translation
│   ├── OrfFinder.java     # six-frame ORF search
│   └── Orf.java           # ORF result record
├── src/test/java/...      # JUnit 5 tests
├── examples/sample.fasta  # synthetic demo sequences
└── .github/workflows/     # CI: build + test on every push
```

## Roadmap

- [ ] Alternative genetic codes (e.g. vertebrate mitochondrial)
- [ ] Protein statistics (molecular weight, isoelectric point)
- [ ] Pairwise alignment (Needleman–Wunsch / Smith–Waterman)
- [ ] FASTQ support with quality-score summaries

## License

[MIT](LICENSE) © 2026 Junxi Song

*The sequences in `examples/sample.fasta` are synthetic and generated for demonstration only.*
