# Comparison with the results and log files of the original R function
# (test/reference/, produced by tools/make_reference.R).

function read_r_map(path)
    f = split.(readlines(path)[2:end], '\t')
    (ids = String.(getindex.(f, 1)), total = parse.(Float64, getindex.(f, 2)),
     spread = parse.(Int, getindex.(f, 3)), parent_id = String.(getindex.(f, 4)),
     curated = String.(getindex.(f, 5)), rank = parse.(Int, getindex.(f, 6)))
end

@testset "Parity with R lulu 0.1.0" begin
    cases = split.(readlines(joinpath(DATA, "cases.tsv"))[2:end], '\t')
    @test length(cases) == 20
    for c in cases
        name, tfile, mfile, rtype, ratio, mm, cooc = c
        @testset "$name" begin
            table = read_otutable(joinpath(DATA, tfile))
            ml = read_matchlist(joinpath(DATA, mfile))
            io = IOBuffer()
            res = lulu(table, ml; minimum_ratio_type = rtype, minimum_ratio = parse(Float64, ratio),
                       minimum_match = parse(Float64, mm),
                       minimum_relative_cooccurence = parse(Float64, cooc), logfile = io)
            m = read_r_map(joinpath(REF, name * "_otu_map.tsv"))
            om = res.otu_map
            @test om.ids == m.ids
            @test om.total == m.total
            @test om.spread == m.spread
            @test om.parent_id == m.parent_id
            @test om.curated == m.curated
            @test om.rank == m.rank
            ct = read_otutable(joinpath(REF, name * "_curated_table.tsv"))
            @test res.curated_table.ids == ct.ids
            @test res.curated_table.samples == ct.samples
            @test res.curated_table.counts == ct.counts
            @test res.curated_otus == readlines(joinpath(REF, name * "_curated_otus.txt"))
            @test res.discarded_otus == readlines(joinpath(REF, name * "_discarded_otus.txt"))
            @test res.curated_count == length(res.curated_otus)
            @test res.discarded_count == length(res.discarded_otus)
            @test String(take!(io)) == read(joinpath(REF, name * ".log"), String)
            # Original table in processing order, and reads are conserved.
            @test res.original_table.ids == m.ids
            @test sum(res.curated_table.counts; dims = 1) ≈ sum(table.counts; dims = 1)
            # Same results with the parallel search (no log), the sparse kernels
            # and a sparse table; same log with a sparse table.
            settings = (minimum_ratio_type = rtype, minimum_ratio = parse(Float64, ratio),
                        minimum_match = parse(Float64, mm),
                        minimum_relative_cooccurence = parse(Float64, cooc))
            sio = IOBuffer()
            for (tbl, kw) in ((table, (storage = :dense,)), (table, (storage = :sparse,)),
                              (sparse(table), (;)), (sparse(table), (logfile = sio,)))
                r2 = lulu(tbl, ml; settings..., kw...)
                @test r2.otu_map.ids == m.ids
                @test r2.otu_map.parent_id == m.parent_id
                @test r2.otu_map.total == m.total && r2.otu_map.spread == m.spread
                @test r2.otu_map.rank == m.rank
                @test r2.curated_table.ids == ct.ids
                @test Matrix(r2.curated_table.counts) == ct.counts
                @test issparse(r2.curated_table.counts) == issparse(tbl.counts)
                @test issparse(r2.original_table.counts) == issparse(tbl.counts)
                @test r2.discarded_otus == res.discarded_otus
            end
            @test String(take!(sio)) == read(joinpath(REF, name * ".log"), String)
        end
    end
end
