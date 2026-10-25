# LULU.jl

LULU.jl is a pure Julia implementation of **LULU**, the post-clustering
curation algorithm for amplicon (metabarcoding) data by Frøslev et al. (2017).
It ports the R package [`lulu`](https://github.com/tobiasgf/lulu) 0.1.0.

LULU removes erroneous OTUs (or ASVs) from an OTU table without discarding
rare but real ones. A "daughter" OTU is merged into a "parent" OTU when the
two sequences are similar, the parent is at least as widespread, the
daughter (almost) always co-occurs with it, and the daughter is less abundant.
No reference database is needed.

> Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K.,
> Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation
> of DNA amplicon data yields reliable biodiversity estimates.
> *Nature Communications*, 8, 1188. <https://doi.org/10.1038/s41467-017-01312-x>

Main features:

- `lulu`: the curation algorithm. It gives the same results as the R function
  (curated table, OTU map, retained and discarded OTUs) and writes the same
  log file. On the tutorial data set it takes 0.01–0.02 s instead of 70 s
  (0.2 s when it also writes the 766,354-line log). It accepts sparse tables
  and searches parents on several threads.
- `build_matchlist`: builds the match list in pure Julia, as the
  BLASTN/VSEARCH step of the LULU workflow does. It uses a vectorised (SIMD)
  alignment and runs on several threads: 2–3 s for the tutorial data set,
  where BLASTN takes 61 s.
- Readers and writers for OTU tables, match lists, FASTA files and OTU maps.
  [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) is supported
  through a package extension.

Dependencies: the standard libraries `Dates`, `Printf` and `SparseArrays`,
and [NucleotideSearch.jl](../NucleotideSearch.jl) (pure Julia k-mer search and
banded alignment, used by `build_matchlist`). Julia ≥ 1.10.

## Installation

The package is not registered, and neither is NucleotideSearch.jl. Install
both from local copies, NucleotideSearch.jl first:

```julia
using Pkg
Pkg.develop(path = "path/to/NucleotideSearch.jl")
Pkg.develop(path = "path/to/LULU.jl")
```

With Julia ≥ 1.11, `Project.toml` also points to `../NucleotideSearch.jl`
(`[sources]`), so `julia --project=.` in a checkout next to
NucleotideSearch.jl works without the first step.

## Quick start

```julia
using LULU

otutab    = read_otutable("otutable.tsv")     # OTUs × samples, ids in the first column
matchlist = read_matchlist("match_list.txt")  # query id, hit id, % similarity

result = lulu(otutab, matchlist)
# Same defaults as R: minimum_ratio_type = "min", minimum_ratio = 1,
# minimum_match = 84, minimum_relative_cooccurence = 0.95

result.curated_table      # OTUTable: daughters summed into their parents
result.curated_count      # number of retained OTUs
result.otu_map            # per-OTU summary (total, spread, parent_id, curated, rank)
write_otutable("curated_table.tsv", result.curated_table)
```

Without an external match list, build it from the OTU sequences:

```julia
matchlist = build_matchlist("otu_sequences.fasta")   # start Julia with -t auto
result = lulu(otutab, matchlist)
```

`examples/quickstart.jl` runs the whole workflow on a subset of the example
data of the R package (`julia --project=. examples/quickstart.jl`).

### Large tables

OTU tables are mostly zeros: most OTUs occur in few samples. For large tables,
read the counts into a sparse matrix:

```julia
otutab = read_otutable("otutable.tsv"; storage = :sparse)
result = lulu(otutab, matchlist)   # start Julia with -t auto
```

A sparse table uses memory only for non-zero counts, and co-occurrence and
abundance ratios are computed over the samples where an OTU is present. The
curated table is sparse too. `sparse(table)` and `Matrix(table)` convert
between the two storages (`using SparseArrays`). The results are identical.

On a synthetic table with 30,000 OTUs, 1,000 samples (4.5% non-zero counts)
and 837,000 match-list entries, `lulu` takes 0.25 s with a sparse table on 12
threads (0.46 s on 1 thread). The same table stored densely takes 0.6–0.8 s
and allocates 5–6 times more memory.

Without a log, the search for parents runs on all threads of the Julia
session. Whether a candidate is accepted does not depend on the parents
assigned before, so the first accepted candidate of every OTU is found in
parallel; the chains of errors are then resolved in processing order. The
result does not depend on the number of threads. With a log the search is
sequential, so that the log follows the processing order as in R.

### Inputs

**OTU table.** OTUs are rows and samples are columns. The values are read
counts, or any other non-negative numbers. `read_otutable` reads a
tab-separated file with sample names in the header and OTU ids in the first
column. The header may have an empty first field, as in the LULU example
data. `OTUTable(ids, samples, counts)` builds a table in memory, and
`OTUTable(df)` builds one from a `DataFrame` (ids in the first column).

**Match list.** It has one line per pair of OTUs: the query id (a potential
daughter), the hit id (a potential parent) and the similarity in percent.
`read_matchlist` reads whitespace-separated files, such as the BLASTN output
`-outfmt '6 qseqid sseqid pident'` or the VSEARCH output
`--userfields query+target+id`. Hits equal to `*`, self-hits and ids missing
from the table are ignored.

The README of the R package suggests these external commands:

```sh
makeblastdb -in OTU_sequences.fasta -parse_seqids -dbtype nucl
blastn -db OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt \
       -qcov_hsp_perc 80 -perc_identity 84 -query OTU_sequences.fasta

vsearch --usearch_global OTU_sequences.fasta --db OTU_sequences.fasta --self --id .84 \
        --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 \
        --query_cov .9 --maxhits 10
```

## The algorithm

This is what `lulu` does. Each step follows the R code.

1. Compute, for each OTU, its total read count and its *spread* (the number of
   samples in which it is present). Sort the OTUs by decreasing spread, then
   by decreasing total. Ties keep the input order.
2. Process the OTUs from top to bottom. Each one is a potential daughter.
3. Its *hits* are the match-list entries where it is the query and the
   similarity is **strictly greater** than `minimum_match`.
4. Its *potential parents* are the hits whose spread is at least as large as
   the daughter's. They are tested in processing order.
5. A potential parent is accepted if:
   - the relative co-occurrence is `>= minimum_relative_cooccurence`. This is
     the fraction of the daughter's samples in which the parent is also
     present;
   - the parent/daughter read ratio over the daughter's samples is
     **strictly greater** than `minimum_ratio`. With `minimum_ratio_type =
     "min"` the minimum ratio is used, with `"avg"` the mean ratio.
6. If the accepted parent was processed earlier, the daughter takes the
   parent's own `parent_id`, so chains of errors end at the ultimate parent.
   Otherwise the daughter takes the parent's id. This happens when the parent
   has the same spread but a smaller total.
7. If no parent is accepted, the OTU is its own parent.
8. Rows with the same `parent_id` are summed into the curated table.

Step 6 has a consequence that comes from the R code. A daughter can be
assigned to a parent with the same spread that is itself merged later. In
that case the daughter's reads stay under that id: the id appears in
`curated_otus` and in `curated_table`, but its `otu_map` entry says `"merged"`.
LULU.jl keeps this behaviour.

### Arguments and result

```julia
lulu(otutable, matchlist; minimum_ratio_type = "min", minimum_ratio = 1,
     minimum_match = 84, minimum_relative_cooccurence = 0.95,
     logfile = nothing, progress = false)
```

| Argument | Meaning |
|---|---|
| `minimum_ratio_type` | `"min"` or `"avg"` (strings or symbols) |
| `minimum_ratio` | minimum parent/daughter abundance ratio, compared with `>` |
| `minimum_match` | minimum similarity in **percent**, compared with `>` |
| `minimum_relative_cooccurence` | minimum co-occurrence rate, compared with `>=` |
| `logfile` | `nothing`; `true` writes `lulu.log_YYYYmmdd_HHMMSS` in the working directory, as R does; a file path; or an `IO` |
| `progress` | print `progress: N%` for each OTU, as R does |
| `storage` | reads used during the search: `:auto` (sparse if the table is sparse or less than a quarter of its counts are non-zero), `:dense` or `:sparse`; affects speed and memory only |

Other methods: `lulu(counts::AbstractMatrix, ids, matchlist; samples)`,
`lulu("otutable.tsv", "match_list.txt")`, and data frames for either input
(`lulu(df, matchlist_df; id_col = 1)`).

The result is a `LuluResult`. It has the same components as the list
returned by R:

| Field | Content |
|---|---|
| `curated_table` | `OTUTable` with the curated OTUs, rows sorted by id |
| `curated_count`, `curated_otus` | retained OTUs (sorted) |
| `discarded_count`, `discarded_otus` | merged OTUs (in processing order) |
| `runtime` | elapsed time in seconds |
| `minimum_match`, `minimum_relative_cooccurence`, `minimum_ratio_type`, `minimum_ratio` | settings used |
| `otu_map` | `OTUMap` with `ids`, `total`, `spread`, `parent_id`, `curated` (`"parent"`/`"merged"`), `rank` |
| `original_table` | the input table, rows in processing order |

`DataFrame(result.curated_table)` and `DataFrame(result.otu_map)` convert the
results to data frames.

### Differences from the R function

The results and logs are the same as in R. LULU.jl even emulates R's
extended-precision `rowSums`, `sum` and `mean`, and R's printing of numbers in
the log. The only differences are:

- No log file and no progress output unless requested. R always writes both.
- An OTU with no reads that has potential parents is kept as its own parent.
  R stops with `missing value where TRUE/FALSE needed`.
- An invalid `minimum_ratio_type` or a `NaN` threshold raises an error. R
  silently uses `"min"` for any other ratio type.
- Ids are sorted by Unicode code point, as dplyr does. R's `table()` uses the
  locale collation for `curated_otus`. The two orders differ only for ids that
  mix upper and lower case or contain punctuation.
- `read_otutable` keeps sample names as written. R's `read.csv` applies
  `make.names` (e.g. `1A` becomes `X1A`).
- If `minimum_match` is between 0 and 1, LULU.jl warns that it looks like a
  fraction instead of a percentage.

## Building the match list in Julia

```julia
build_matchlist(ids, sequences; min_identity = 84, iddef = 1, query_cov = 0.9,
                maxhits = 10, maxaccepts = 0, maxrejects = 32, wordlength = 8,
                minwordmatches = 12, match = 2, mismatch = -4, gap_open = 20,
                gap_extend = 2, terminal_gap_open = 2, terminal_gap_extend = 1)
build_matchlist("sequences.fasta"; kwargs...)
build_matchlist(["id1" => "ACGT…", …]; kwargs...)
```

The defaults reproduce the settings of the VSEARCH command in the LULU README.
The search follows the same strategy as `vsearch --usearch_global --self`:

- Targets are ranked by the number of unique 8-mers they share with the
  query.
- They are aligned in that order with a global alignment. It uses affine gap
  penalties, with separate penalties for terminal gaps.
- The search stops after `maxaccepts` accepted hits (0: no limit) or after
  `maxrejects` rejected candidates (0: no limit).
- The best `maxhits` hits are kept (0: all).

A hit is accepted when its identity is `>= min_identity` (percent) and the
fraction of query bases aligned to target bases is `>= query_cov`.

`iddef` selects the identity definition, as in VSEARCH:

| `iddef` | Identity |
|---|---|
| `0` | matches / length of the shorter sequence |
| `1`, `4` | matches / alignment columns |
| `2` | matches / alignment columns excluding terminal gaps |
| `3` | 1 − (mismatches + gap openings) / length of the longer sequence |

The alignment is restricted to a diagonal band. The band is wide enough for
every alignment that can be accepted under `iddef`, `min_identity` and
`query_cov`. An accepted
alignment that touches the band edge is recomputed without the band. The
output does not depend on the number of threads (verified with 1 and 12
threads). `maxrejects = 0, minwordmatches = 0` gives an exhaustive search.

Implementation notes (the search is `search_self` of NucleotideSearch.jl):

- The dynamic programming matrix is filled by anti-diagonals. The cells of an
  anti-diagonal are independent, so the inner loop is vectorised (SIMD) with
  16-bit scores, or 32-bit scores when 16 bits could overflow (very long
  sequences). This is about 5 times faster than a row-by-row loop.
- Each pair of sequences is aligned once, in a fixed orientation, and the
  result is used for both directions (query → target and target → query).
  This saves about 38% of the alignments. The identity of a pair is
  therefore the same in both directions.
- Most aligned candidates are accepted hits (77% on the example data), whose
  identity has to be computed anyway, so a stricter k-mer filter would save
  little; exact filters (base composition, q-gram lemma) would remove 0.5%
  of the rejected candidates at 84% identity.

The implementation is independent of VSEARCH and BLASTN, so the hits and
identities are close to theirs but not identical. On the example data of the
R package (2,425 OTUs, 130 samples, 12 threads), compared with the BLASTN
match list of the tutorial (LULU retains 459 OTUs):

| `build_matchlist` settings | Time | Entries | OTUs retained by LULU | Same parent/merged status as with BLASTN |
|---|---|---|---|---|
| defaults (README VSEARCH command) | 1.6–2.7 s | 20,977 | 573 | 2,293 of 2,425 |
| `maxhits = 0` | 1.6–2.9 s | 227,820 | 526 | 2,334 of 2,425 |
| `iddef = 2, query_cov = 0.8, maxhits = 0` (closest to BLASTN) | 2.3–4.2 s | 234,726 | 468 | 2,392 of 2,425 |

BLASTN with 8 threads took 61 s on the same machine. With `iddef = 1`,
terminal gaps count against identity, so length variants of ITS2 sequences
get lower identities and fewer OTUs are merged. `maxhits = 10` also removes
potential parents.

## Validation

- **`test/`** (`julia --project=. -e 'using Pkg; Pkg.test()'`, 1,249 tests,
  no network, about 1 minute). It compares LULU.jl with the original R
  function on 20 cases, generated with `tools/make_test_data.jl` and
  `tools/make_reference.R`:
  - 3 synthetic tables × 5 settings;
  - a table of relative abundances (floats) with an OTU without reads;
  - a subset of the example data (every 8th OTU, BLASTN matches).

  The synthetic cases include ties, chains of errors, parents with the same
  spread, `*` hits, self-hits, unknown ids and duplicated entries. The OTU
  map, the curated table, the OTU lists and the **full log text** are
  identical. Every case is also run with the parallel search, with the
  sparse kernels and with a sparse table (with and without log). Unit tests
  check the documented behaviour and the I/O. The banded aligner is checked
  against an unbanded Gotoh alignment on 300 random pairs, with 16- and
  32-bit scores.
- **`validation/example_data.jl`** (`julia -t auto --project=.
  validation/example_data.jl`) uses the full example data set: 2,425 OTUs and
  a BLASTN match list with 320,803 entries. It needs the files created by
  `tools/make_validation_data.sh` (BLAST+, R with dplyr). Results for four
  settings (default; `avg`; strict; loose) are identical to R, including the
  logs (766,354 lines for the default settings), also with the parallel
  search and with a sparse table. With the default settings,
  459 OTUs are retained and 1,966 merged, as in the tutorial of the R
  package.

## License and credits

The algorithm and the original R code are by Tobias Guldberg Frøslev
(<https://github.com/tobiasgf/lulu>). LULU.jl is a derivative work released
under the GNU Lesser General Public License v3 (see `LICENSE`, `COPYING`
and `NOTICE.md`).
