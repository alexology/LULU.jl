# Full-size validation on the example data of the R package lulu (2,425 OTUs,
# 130 samples; Frøslev et al. 2017).
#
# Requires the files produced by `tools/make_validation_data.sh`:
# reference/work/{otutable_test.txt, centroids_test.txt, match_list.txt}
# (BLASTN match list) and the R results `blast_*` with their logs.
#
#     julia -t auto --project=. validation/example_data.jl

using LULU
using SparseArrays

const WORK = joinpath(@__DIR__, "..", "reference", "work")
isfile(joinpath(WORK, "match_list.txt")) ||
    error("missing $(joinpath(WORK, "match_list.txt")): run tools/make_validation_data.sh first")

function read_r_map(path)
    lines = readlines(path)[2:end]
    f = split.(lines, '\t')
    (ids = String.(getindex.(f, 1)), total = parse.(Float64, getindex.(f, 2)),
     spread = parse.(Int, getindex.(f, 3)), parent_id = String.(getindex.(f, 4)),
     curated = String.(getindex.(f, 5)), rank = parse.(Int, getindex.(f, 6)))
end

function compare(label, res, prefix; logfile = nothing)
    ok = true
    check(name, cond) = (cond || (ok = false; println("  MISMATCH: ", name)))
    m = read_r_map(joinpath(WORK, prefix * "_otu_map.tsv"))
    om = res.otu_map
    check("otu_map ids (order)", om.ids == m.ids)
    check("otu_map total", om.total == m.total)
    check("otu_map spread", om.spread == m.spread)
    check("otu_map parent_id", om.parent_id == m.parent_id)
    check("otu_map curated", om.curated == m.curated)
    check("otu_map rank", om.rank == m.rank)
    ct = read_otutable(joinpath(WORK, prefix * "_curated_table.tsv"))
    check("curated_table", res.curated_table.ids == ct.ids && res.curated_table.counts == ct.counts &&
                           res.curated_table.samples == ct.samples)
    check("curated_otus", res.curated_otus == readlines(joinpath(WORK, prefix * "_curated_otus.txt")))
    check("discarded_otus", res.discarded_otus == readlines(joinpath(WORK, prefix * "_discarded_otus.txt")))
    if logfile !== nothing
        rlog = replace(read(joinpath(WORK, prefix * ".log"), String), "\r\n" => "\n")
        check("log text", rlog == read(logfile, String))
    end
    println(rpad(label, 36), ok ? "identical to R" : "DIFFERENT from R",
            "  (", res.curated_count, " retained, ", res.discarded_count, " merged, ",
            round(res.runtime; digits = 2), " s)")
    ok
end

table = read_otutable(joinpath(WORK, "otutable_test.txt"))
blast = read_matchlist(joinpath(WORK, "match_list.txt"))
println("OTU table: ", table, "; BLASTN match list: ", blast)

settings = [
    ("blast_default", (;)),
    ("blast_avg", (minimum_ratio_type = "avg",)),
    ("blast_strict", (minimum_ratio = 10, minimum_match = 90, minimum_relative_cooccurence = 1.0)),
    ("blast_loose", (minimum_ratio_type = "avg", minimum_ratio = 0.5, minimum_match = 80,
                     minimum_relative_cooccurence = 0.6)),
]
allok = true
sparse_table = sparse(table)
for (prefix, kw) in settings
    isfile(joinpath(WORK, prefix * "_otu_map.tsv")) || (println("skipped $prefix (no R result)"); continue)
    logfile = tempname()
    res = lulu(table, blast; logfile, kw...)
    global allok &= compare(prefix, res, prefix;
                            logfile = isfile(joinpath(WORK, prefix * ".log")) ? logfile : nothing)
    rm(logfile; force = true)
    # Parallel search (no log) and sparse table.
    for (label, tbl) in (("parallel, no log", table), ("sparse table", sparse_table))
        global allok &= compare(prefix * ", " * label, lulu(tbl, blast; kw...), prefix)
    end
end

# Pure Julia match list instead of BLASTN.
ids, seqs = read_fasta(joinpath(WORK, "centroids_test.txt"))
for (label, kw) in [("README VSEARCH settings", (;)),
                    ("no maxhits limit", (maxhits = 0,)),
                    ("BLAST-like: iddef 2, query_cov 0.8", (iddef = 2, query_cov = 0.8, maxhits = 0))]
    t = @elapsed ml = build_matchlist(ids, seqs; kw...)
    isempty(kw) && write_matchlist(joinpath(WORK, "julia_matchlist_t$(Threads.nthreads()).tsv"), ml; digits = nothing)
    res = lulu(table, ml)
    rb = lulu(table, blast)
    agree = count(res.otu_map.curated .== rb.otu_map.curated[indexin(res.otu_map.ids, rb.otu_map.ids)])
    println("build_matchlist (", label, "): ", length(ml), " entries in ", round(t; digits = 1),
            " s with ", Threads.nthreads(), " threads; LULU retains ", res.curated_count,
            " OTUs (BLASTN: ", rb.curated_count, "); parent/merged status agrees for ",
            agree, " of ", notus(table), " OTUs")
end
println(allok ? "All comparisons with R: identical." : "Some comparisons with R differ.")
