# Changelog

## Unreleased

- The k-mer search and the banded alignment of `build_matchlist` moved to
  the new package NucleotideSearch.jl (`search_self`), now a dependency.
  Match lists are identical bit for bit (verified on the example data with
  nine settings) and equally fast. The alignment unit tests moved with the
  code. Scores and penalties above 10,000 in absolute value are now
  rejected.
- `read_fasta` now calls `NucleotideSearch.read_fasta` (FASTX.jl): faster on
  large or line-wrapped files, all whitespace removed from sequences, and an
  `ArgumentError` for malformed files, including an empty record directly
  followed by another header.

## 0.1.0 (2026-10-19)

First version.

- `lulu`: port of the R function `lulu::lulu` 0.1.0, with identical results
  and log files (verified against R 4.5.0 on 20 test cases and on the full
  example data set with four settings).
- `build_matchlist`: pure Julia, multi-threaded construction of the match
  list (k-mer ranking and banded global alignment, VSEARCH-like options).
- I/O: `read_otutable`, `write_otutable`, `read_matchlist`, `write_matchlist`,
  `read_fasta`, `write_otu_map`.
- DataFrames.jl package extension.
- `lulu` restructured into small type-stable functions; the log keyword is
  `logfile` (not `log`, which would shadow `Base.log`); `NaN` thresholds are
  rejected.
- Sparse OTU tables: `OTUTable` keeps a `SparseMatrixCSC`;
  `read_otutable(...; storage = :sparse)`; `sparse(table)`/`Matrix(table)`;
  sparse kernels for co-occurrence and abundance ratios; `lulu(...; storage)`.
- `lulu` without a log searches parents on all threads (identical results).
- `build_matchlist`: vectorised (SIMD) anti-diagonal alignment with 16-bit
  scores (32-bit when needed); each pair aligned once and shared by both
  directions, so identities are symmetric; search split into small
  functions.
