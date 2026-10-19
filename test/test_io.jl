@testset "Types and input/output" begin
    @test_throws DimensionMismatch OTUTable(["a"], ["s1", "s2"], [1 2 3])
    @test_throws ArgumentError OTUTable(["a", "a"], ["s1"], [1; 2;;])
    @test_throws ArgumentError OTUTable(["a"], ["s1"], [-1;;])
    @test_throws ArgumentError OTUTable(["a"], ["s1"], [NaN;;])
    @test_throws DimensionMismatch MatchList(["a"], ["b", "c"], [1.0])
    t = OTUTable(["o1", "o2"], ["s1", "s2", "s3"], [1 0 3; 4 5 0])
    @test notus(t) == 2 && nsamples(t) == 3 && size(t) == (2, 3)
    @test occursin("2 OTUs and 3 samples", sprint(show, t))

    mktempdir() do dir
        p = joinpath(dir, "t.tsv")
        write_otutable(p, t)
        @test read_otutable(p) == t
        @test readline(p) == "OTUid\ts1\ts2\ts3"
        # Header without a name for the id column (R write.table style) and
        # with an empty first field (the LULU example file), CRLF, quotes.
        write(p, "s1\ts2\n\"o1\"\t1\t2\r\no2\t3\t4\n\n")
        r = read_otutable(p)
        @test r.ids == ["o1", "o2"] && r.samples == ["s1", "s2"] && r.counts == [1 2; 3 4]
        write(p, "\ts1\ts2\no1\t1.5\t2\n")
        r = read_otutable(p)
        @test r.counts isa Matrix{Float64} && r.counts == [1.5 2.0]
        write(p, "\ts1\ts2\no1\t1\tx\n")
        @test_throws ArgumentError read_otutable(p)
        write(p, "\ts1\ts2\no1\t1\t2\no2\t1\n")
        @test_throws ArgumentError read_otutable(p)
        write(p, "a\tb\tc\td\no1\t1\n")
        @test_throws ArgumentError read_otutable(p)

        # Sparse reading and writing.
        write(p, "\ts1\ts2\ts3\no1\t0\t3\t0\no2\t0\t0\t0\no3\t5\t0\t1\n")
        rs = read_otutable(p; storage = :sparse)
        rd = read_otutable(p)
        @test issparse(rs.counts) && nnz(rs.counts) == 3 && rs == rd
        @test rs.counts isa SparseMatrixCSC{Int}
        write(p, "\ts1\ts2\no1\t0\t2.5\no2\t1\t0\n")
        rs = read_otutable(p; storage = :sparse)
        @test rs.counts isa SparseMatrixCSC{Float64} && nnz(rs.counts) == 2
        @test Matrix(rs.counts) == [0.0 2.5; 1.0 0.0]
        @test_throws ArgumentError read_otutable(p; storage = :compressed)
        big = OTUTable(["a", "b", "c"], ["s1", "s2", "s3", "s4"], [0 0 7 0; 1 0 0 2; 0 0 0 0])
        p2 = joinpath(dir, "t2.tsv")
        write_otutable(p, big); write_otutable(p2, sparse(big))
        @test read(p, String) == read(p2, String)
        @test read_otutable(p2; storage = :sparse) == big

        m = MatchList([("o1", "o2", 90.12345), ("o2", "o1", 100.0)])
        q = joinpath(dir, "m.tsv")
        write_matchlist(q, m)
        @test readlines(q) == ["o1\to2\t90.123", "o2\to1\t100.000"]
        write_matchlist(q, m; digits = nothing)
        @test read_matchlist(q) == m
        # Whitespace separators, comments, extra columns (as read.table).
        write(q, "# comment\no1  o2 95.5 extra\n\no2\to1\t88 # note\n")
        @test read_matchlist(q) == MatchList([("o1", "o2", 95.5), ("o2", "o1", 88.0)])
        write(q, "o1 o2\n")
        @test_throws ArgumentError read_matchlist(q)
        write(q, "o1 o2 high\n")
        @test_throws ArgumentError read_matchlist(q)

        f = joinpath(dir, "s.fasta")
        write(f, ">a desc\nACGT\nacg\n\n>b\nTTTT\n")
        @test read_fasta(f) == (["a", "b"], ["ACGTACG", "TTTT"])
        write(f, "ACGT\n>a\nA\n")
        @test_throws ArgumentError read_fasta(f)

        r = lulu(t, MatchList([("o1", "o2", 99.0)]))
        mp = joinpath(dir, "map.tsv")
        write_otu_map(mp, r.otu_map)
        @test readlines(mp)[1] == "OTUid\ttotal\tspread\tparent_id\tcurated\trank"
        @test length(readlines(mp)) == 3
    end
end
