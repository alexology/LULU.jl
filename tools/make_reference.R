# Runs the original R function lulu() (lulu 0.1.0, R/Functions.R) on every case
# listed in test/data/cases.tsv and writes the results and log files to
# test/reference/. Run from the package root:
#
#     Rscript tools/make_reference.R [path/to/lulu/R/Functions.R] [data dir] [output dir]
#
# Default source: reference/lulu/R/Functions.R (a clone of
# https://github.com/tobiasgf/lulu, see tools/make_validation_data.sh).
# Requires dplyr. Doubles are written with 17 significant digits so that
# Julia can compare them exactly.

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args) >= 1) args[1] else "reference/lulu/R/Functions.R"
data_dir <- normalizePath(if (length(args) >= 2) args[2] else "test/data")
out_dir <- if (length(args) >= 3) args[3] else "test/reference"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir)
source(src)
suppressPackageStartupMessages(library(dplyr))

fmt <- function(df) {
  for (k in seq_along(df)) if (is.double(df[[k]])) df[[k]] <- sprintf("%.17g", df[[k]])
  df
}

cases <- read.delim(file.path(data_dir, "cases.tsv"), colClasses = "character")
work <- tempfile("lulu_ref_")
dir.create(work)
old <- setwd(work)
for (k in seq_len(nrow(cases))) {
  cs <- cases[k, ]
  otutab <- read.csv(file.path(data_dir, cs$otutable), sep = "\t", header = TRUE,
                     as.is = TRUE, row.names = 1)
  matchlist <- read.table(file.path(data_dir, cs$matchlist), header = FALSE, as.is = TRUE,
                          stringsAsFactors = FALSE)
  unlink(list.files(pattern = "^lulu[.]log_"))
  invisible(capture.output(res <- lulu(otutab, matchlist,
    minimum_ratio_type = cs$minimum_ratio_type,
    minimum_ratio = as.numeric(cs$minimum_ratio),
    minimum_match = as.numeric(cs$minimum_match),
    minimum_relative_cooccurence = as.numeric(cs$minimum_relative_cooccurence))))
  prefix <- file.path(out_dir, cs$name)
  m <- res$otu_map
  write.table(fmt(data.frame(OTUid = row.names(m), m)), paste0(prefix, "_otu_map.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ct <- res$curated_table
  write.table(fmt(data.frame(OTUid = row.names(ct), ct, check.names = FALSE)),
              paste0(prefix, "_curated_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(res$curated_otus, paste0(prefix, "_curated_otus.txt"))
  writeLines(res$discarded_otus, paste0(prefix, "_discarded_otus.txt"))
  log <- list.files(pattern = "^lulu[.]log_")
  stopifnot(length(log) == 1)
  txt <- gsub("\r\n", "\n", readChar(log, file.size(log), useBytes = TRUE), fixed = TRUE)
  writeBin(charToRaw(txt), paste0(prefix, ".log"))
  cat(sprintf("%-28s %4d retained, %4d merged\n", cs$name, res$curated_count, res$discarded_count))
}
unlink(list.files(pattern = "^lulu[.]log_"))
setwd(old)
unlink(work, recursive = TRUE)
