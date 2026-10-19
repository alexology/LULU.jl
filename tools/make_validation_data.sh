#!/usr/bin/env sh
# Prepares reference/work/ for validation/example_data.jl and for the example
# subset in tools/make_test_data.jl:
#   - clones the original R package (https://github.com/tobiasgf/lulu) into
#     reference/lulu;
#   - builds the BLASTN match list of the example data as in the LULU tutorial;
#   - runs the R function lulu() with four settings and keeps results and logs.
# Requires git, BLAST+ (makeblastdb, blastn), R with dplyr. Run from the
# package root:  sh tools/make_validation_data.sh [path/to/Rscript]
set -e
RSCRIPT=${1:-Rscript}
mkdir -p reference/work
[ -d reference/lulu ] || git clone --depth 1 https://github.com/tobiasgf/lulu.git reference/lulu
cd reference/work
cp ../lulu/Example_data/otutable_test.txt ../lulu/Example_data/centroids_test.txt .
makeblastdb -in centroids_test.txt -parse_seqids -dbtype nucl > /dev/null
blastn -db centroids_test.txt -outfmt '6 qseqid sseqid pident' -out match_list.txt \
       -qcov_hsp_perc 80 -perc_identity 84 -query centroids_test.txt -num_threads 4
cat > cases.tsv <<'CASES'
name	otutable	matchlist	minimum_ratio_type	minimum_ratio	minimum_match	minimum_relative_cooccurence
blast_default	otutable_test.txt	match_list.txt	min	1	84	0.95
blast_avg	otutable_test.txt	match_list.txt	avg	1	84	0.95
blast_strict	otutable_test.txt	match_list.txt	min	10	90	1
blast_loose	otutable_test.txt	match_list.txt	avg	0.5	80	0.6
CASES
cd ../..
"$RSCRIPT" tools/make_reference.R reference/lulu/R/Functions.R reference/work reference/work
