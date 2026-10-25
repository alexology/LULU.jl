# Notice

LULU.jl is a Julia port of the R package **lulu** 0.1.0, "Post Clustering
Curation of Amplicon Data for Reliable Biodiversity Metrics", by Tobias
Guldberg Frøslev: <https://github.com/tobiasgf/lulu> (commit cloned on
2026-10-19).

- `src/curation.jl` follows `R/Functions.R` (function `lulu`) line by line,
  including the text of the log messages.
- `test/data/example_subset_*` and the full validation use the example data
  (`Example_data/otutable_test.txt`, `Example_data/centroids_test.txt`) of
  the R package. They come from Frøslev et al. (2017, *Nature Communications*
  8:1188). The match list was computed with BLASTN 2.17.0+ using the command
  of the LULU tutorial.
- The reference results in `test/reference/` were produced by running the
  original R code (R 4.5.0, dplyr 1.1.4) with `tools/make_reference.R`.

Licence: the `DESCRIPTION` file of the R package states "GNU General Public
License v3 (GPL-3)", while its `LICENSE` file contains the text of the GNU
Lesser General Public License v3. LULU.jl follows the `LICENSE` file and is
distributed under the GNU LGPL v3 (`LICENSE`). The LGPL v3 supplements the
GNU GPL v3, whose text is in `COPYING`.

`build_matchlist` (`src/matchlist.jl`, `src/alignment.jl`) is new code. It
takes its search strategy and option names from VSEARCH
(<https://github.com/torognes/vsearch>), but contains no VSEARCH code.
