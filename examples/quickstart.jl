# Quick start with the files shipped in test/data (a subset of the example
# data of the R package lulu: 304 OTUs, 130 samples, BLASTN match list).
#
#     julia --project=. examples/quickstart.jl

using LULU

data = joinpath(@__DIR__, "..", "test", "data")
otutab = read_otutable(joinpath(data, "example_subset_table.tsv"))
matchlist = read_matchlist(joinpath(data, "example_subset_matchlist.tsv"))
println(otutab)
println(matchlist)

result = lulu(otutab, matchlist)   # minimum_ratio_type = "min", minimum_ratio = 1,
                                   # minimum_match = 84, minimum_relative_cooccurence = 0.95
println(result)
println("First retained OTUs: ", first(result.curated_otus, 3))
println("First merged OTUs:   ", first(result.discarded_otus, 3))

m = result.otu_map
println("\nOTU map (first 5 rows):")
println(join(["OTUid", "total", "spread", "parent_id", "curated", "rank"], '\t'))
for i in 1:5
    println(join((first(m.ids[i], 10) * "…", Int(m.total[i]), m.spread[i],
                  first(m.parent_id[i], 10) * "…", m.curated[i], m.rank[i]), '\t'))
end

# Write the curated table and the OTU map.
out = mktempdir()
write_otutable(joinpath(out, "curated_table.tsv"), result.curated_table)
write_otu_map(joinpath(out, "otu_map.tsv"), result.otu_map)
println("\nWritten to ", out)
