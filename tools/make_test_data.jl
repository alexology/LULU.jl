# Generates the inputs of the R parity tests in test/data/ (synthetic tables
# and a subset of the example data of the R package). Run from the package
# root, then run tools/make_reference.R:
#
#     julia --project=. tools/make_test_data.jl
#     Rscript tools/make_reference.R
#
# The subset of the example data requires reference/work/ (see
# tools/make_validation_data.sh).

using LULU, Random

const DATA = joinpath(@__DIR__, "..", "test", "data")

# Synthetic community: parents with random profiles; daughters are scaled,
# noisy copies of a parent (sometimes of another daughter), plus unrelated
# OTUs and duplicated rows (ties in spread and total).
function synthetic(rng; nparents = 15, ndaughters = 50, nother = 10, nsamples = 14)
    rows = Vector{Vector{Int}}()
    ids = String[]
    for p in 1:nparents
        push!(rows, [rand(rng) < 0.6 ? rand(rng, 50:5000) : 0 for _ in 1:nsamples])
        push!(ids, "par$(lpad(p, 2, '0'))")
    end
    parent_of = Dict{String,String}()
    for d in 1:ndaughters
        src = rand(rng) < 0.25 && length(rows) > nparents ? rand(rng, nparents+1:length(rows)) :
                                                           rand(rng, 1:nparents)
        f = rand(rng, (0.002, 0.01, 0.05, 0.2, 0.8, 1.2))
        r = [x > 0 ? round(Int, x * f * (0.5 + rand(rng))) : 0 for x in rows[src]]
        for j in eachindex(r)          # dropouts and spurious presences
            rand(rng) < 0.1 && (r[j] = 0)
            rand(rng) < 0.04 && (r[j] = rand(rng, 1:20))
        end
        # R fails on an OTU without reads that has potential parents
        all(==(0), r) && (r[rand(rng, eachindex(r))] = 1)
        push!(rows, r)
        push!(ids, "dau$(lpad(d, 2, '0'))")
        parent_of[ids[end]] = ids[src]
    end
    for o in 1:nother
        push!(rows, [rand(rng) < 0.3 ? rand(rng, 1:300) : 0 for _ in 1:nsamples])
        push!(ids, "oth$(lpad(o, 2, '0'))")
    end
    for k in 1:4                       # exact duplicates: ties
        src = rand(rng, 1:length(rows))
        push!(rows, copy(rows[src])); push!(ids, "dup$(k)_" * ids[src])
        parent_of[ids[end]] = ids[src]
    end
    perm = randperm(rng, length(ids))
    counts = permutedims(reduce(hcat, rows))[perm, :]
    ids = ids[perm]
    # Match list.
    ml = Tuple{String,String,Float64}[]
    for (d, p) in parent_of
        push!(ml, (d, p, round(rand(rng, 78:0.001:100); digits = 3)))
        rand(rng) < 0.5 && push!(ml, (p, d, round(rand(rng, 78:0.001:100); digits = 3)))
    end
    for _ in 1:150
        a, b = rand(rng, ids), rand(rng, ids)
        push!(ml, (a, b, round(rand(rng, 75:0.001:100); digits = 3)))
    end
    push!(ml, (ids[1], ids[1], 100.0))            # self-hit
    push!(ml, (ids[2], "*", 0.0))                 # VSEARCH no-hit
    push!(ml, (ids[3], "unknown_otu", 99.0))      # id absent from the table
    append!(ml, ml[1:5])                          # duplicated entries
    shuffle!(rng, ml)
    OTUTable(ids, ["s$(lpad(j, 2, '0'))" for j in 1:nsamples], counts), MatchList(ml)
end

cases = Tuple{String,String,String,Vector{String}}[]  # name, table, matchlist, settings
settings = Dict(
    "default" => ["min", "1", "84", "0.95"],
    "avg" => ["avg", "1", "84", "0.95"],
    "strict" => ["min", "10", "90", "1"],
    "loose" => ["avg", "0.5", "80", "0.6"],
    "lowcooc" => ["min", "0.9", "84", "0.5"])

for seed in 1:3
    rng = MersenneTwister(seed)
    t, m = synthetic(rng)
    write_otutable(joinpath(DATA, "synthetic$(seed)_table.tsv"), t)
    write_matchlist(joinpath(DATA, "synthetic$(seed)_matchlist.tsv"), m)
    for s in ("default", "avg", "strict", "loose", "lowcooc")
        push!(cases, ("synthetic$(seed)_$s", "synthetic$(seed)_table.tsv",
                      "synthetic$(seed)_matchlist.tsv", settings[s]))
    end
end

# Relative abundances (floating-point counts) and an OTU without reads and
# without hits.
let rng = MersenneTwister(11)
    t, m = synthetic(rng)
    rel = t.counts ./ max.(sum(t.counts; dims = 1), 1)
    rel = vcat(rel, zeros(1, size(rel, 2)))
    tf = OTUTable([t.ids; "zero_otu"], t.samples, rel)
    write_otutable(joinpath(DATA, "relative_table.tsv"), tf)
    write_matchlist(joinpath(DATA, "relative_matchlist.tsv"), m)
    for s in ("default", "avg", "loose")
        push!(cases, ("relative_$s", "relative_table.tsv", "relative_matchlist.tsv", settings[s]))
    end
end

# Every 8th OTU of the example data with the BLASTN matches among them.
work = joinpath(@__DIR__, "..", "reference", "work")
if isfile(joinpath(work, "match_list.txt"))
    full = read_otutable(joinpath(work, "otutable_test.txt"))
    keep = 1:8:notus(full)
    sub = OTUTable(full.ids[keep], full.samples, full.counts[keep, :])
    ids = Set(sub.ids)
    bl = read_matchlist(joinpath(work, "match_list.txt"))
    k = [i for i in 1:length(bl) if bl.query[i] in ids && bl.hit[i] in ids]
    write_otutable(joinpath(DATA, "example_subset_table.tsv"), sub)
    write_matchlist(joinpath(DATA, "example_subset_matchlist.tsv"),
                    MatchList(bl.query[k], bl.hit[k], bl.match[k]))
    for s in ("default", "avg")
        push!(cases, ("example_subset_$s", "example_subset_table.tsv",
                      "example_subset_matchlist.tsv", settings[s]))
    end
else
    @warn "reference/work/match_list.txt not found: example subset not regenerated"
    for s in ("default", "avg")
        push!(cases, ("example_subset_$s", "example_subset_table.tsv",
                      "example_subset_matchlist.tsv", settings[s]))
    end
end

open(joinpath(DATA, "cases.tsv"), "w") do io
    println(io, join(["name", "otutable", "matchlist", "minimum_ratio_type", "minimum_ratio",
                      "minimum_match", "minimum_relative_cooccurence"], '\t'))
    for (n, t, m, s) in cases
        println(io, join([n, t, m, s...], '\t'))
    end
end
println("Wrote $(length(cases)) cases to ", joinpath(DATA, "cases.tsv"))
