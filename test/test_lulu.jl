@testset "lulu: small examples" begin
    # P is the most widespread OTU; D and Z are its daughters.
    t = OTUTable(["D", "P", "Z"], ["a", "b", "c"], [1 2 0; 10 20 5; 0 0 3])
    ml = MatchList([("D", "P", 90.0), ("Z", "P", 90.0)])
    r = lulu(t, ml)
    @test r.otu_map.ids == ["P", "D", "Z"]
    @test r.otu_map.parent_id == ["P", "P", "P"]
    @test r.otu_map.curated == ["parent", "merged", "merged"]
    @test r.otu_map.total == [35.0, 3.0, 3.0]
    @test r.otu_map.spread == [3, 2, 1]
    @test r.otu_map.rank == [1, 2, 3]
    @test r.curated_otus == ["P"]
    @test r.discarded_otus == ["D", "Z"]
    @test r.curated_table.counts == [11 22 8]
    @test r.curated_table isa OTUTable{Int}
    @test r.minimum_match == 84 && r.minimum_relative_cooccurence == 0.95
    @test r.minimum_ratio_type == "min" && r.minimum_ratio == 1
    @test r.runtime >= 0
    @test occursin("1 OTUs retained, 2 merged", sprint(show, r))

    # Strict comparisons: similarity must be > minimum_match, ratio > minimum_ratio.
    @test lulu(t, MatchList([("D", "P", 84.0)])).otu_map.parent_id[2] == "D"
    @test lulu(t, MatchList([("D", "P", 84.01)])).otu_map.parent_id[2] == "P"
    @test lulu(t, ml; minimum_ratio = 10).otu_map.parent_id[2] == "D"   # min ratio = 10
    @test lulu(t, ml; minimum_ratio = 9.99).otu_map.parent_id[2] == "P"
    # Only query -> hit is used: the reverse direction has no effect.
    @test lulu(t, MatchList([("P", "D", 99.0)])).curated_count == 3
    # Self-hits, "*" and unknown ids are ignored.
    @test lulu(t, MatchList([("D", "D", 100.0), ("D", "*", 0.0), ("D", "X", 99.0)])).curated_count == 3

    # Co-occurrence >= threshold: Z occurs in 1 sample out of 1 with P.
    t2 = OTUTable(["P", "D"], ["a", "b", "c", "d"], [10 10 10 10; 1 1 1 0])
    t3 = OTUTable(["P", "D"], ["a", "b", "c", "d"], [10 10 10 0; 1 1 1 1])
    m2 = MatchList([("D", "P", 95.0)])
    @test lulu(t2, m2).otu_map.parent_id == ["P", "P"]
    # D present in 4 samples, P in 3: same spread rule excludes nothing here
    # (spread 4 > 3, so P cannot be a potential parent of D).
    @test lulu(t3, m2).otu_map.parent_id == ["D", "P"]

    # Relative co-occurrence below threshold.
    t4 = OTUTable(["P", "D"], ["a", "b", "c", "d", "e"], [10 10 10 10 0; 1 1 1 0 1])
    @test lulu(t4, m2).otu_map.parent_id[2] == "D"                        # 3/4 < 0.95
    @test lulu(t4, m2; minimum_relative_cooccurence = 0.75).otu_map.parent_id[2] == "D"  # min ratio 0
    @test lulu(t4, m2; minimum_relative_cooccurence = 0.75, minimum_ratio_type = :avg,
               minimum_ratio = 1).otu_map.parent_id[2] == "P"            # mean (10+10+10+0)/4 = 7.5
end

@testset "lulu: chains and same-spread parents" begin
    # A is processed before B (same spread, higher total) and gets B's id as
    # parent; B is then merged into P. As in R, A keeps parent_id "B", so B
    # appears among the curated OTUs although it is "merged".
    t = OTUTable(["A", "B", "P"], ["a", "b", "c"], [100 1 0; 50 50 0; 1000 1000 1000])
    ml = MatchList([("A", "B", 95.0), ("B", "P", 95.0)])
    r = lulu(t, ml; minimum_ratio_type = "avg")
    @test r.otu_map.ids == ["P", "A", "B"]
    @test r.otu_map.parent_id == ["P", "B", "P"]
    @test r.otu_map.curated == ["parent", "merged", "merged"]
    @test r.curated_otus == ["B", "P"]
    @test r.discarded_otus == ["A"]
    @test r.curated_table.ids == ["B", "P"]
    @test r.curated_table.counts == [100 1 0; 1050 1050 1000]

    # Chain: C -> B (processed earlier) -> P gives C the parent of B.
    t = OTUTable(["P", "B", "C"], ["a", "b", "c"], [1000 1000 1000; 100 100 0; 10 10 0])
    ml = MatchList([("B", "P", 95.0), ("C", "B", 95.0)])
    r = lulu(t, ml)
    @test r.otu_map.parent_id == ["P", "P", "P"]
    @test r.curated_table.counts == [1110 1110 1000]

    # Potential parents are tested in processing order; later ones are only logged.
    t = OTUTable(["P1", "P2", "D"], ["a", "b"], [100 100; 50 50; 1 1])
    ml = MatchList([("D", "P2", 99.0), ("D", "P1", 90.0)])
    io = IOBuffer()
    r = lulu(t, ml; logfile = io)
    @test r.otu_map.parent_id == ["P1", "P2", "P1"]
    s = String(take!(io))
    @test occursin("\n---hits: P2 \n---hits: P1\n---potential parent: P1 \n---potential parent: P2\n" *
                   "------checking: P1\n------relative cooccurence: 1 which is sufficient!\n" *
                   "------min avg abundance: 100 which is OK!\nSETTING D to be an ERROR of P1\n\n" *
                   "------checking: P2", s)
    @test endswith(s, "------checking: P2")
    @test occursin("####processing: P1 #####\n---hits: \n---potential parent: \nNo parent found!\n", s)

    # Ties in spread and total keep the input order.
    t = OTUTable(["x", "y", "z"], ["a", "b"], [1 1; 1 1; 2 2])
    @test lulu(t, MatchList(String[], String[], Float64[])).otu_map.ids == ["z", "x", "y"]
    @test lulu(t, MatchList(String[], String[], Float64[])).otu_map.rank == [1, 2, 3]
end

@testset "lulu: documented differences from R" begin
    # An OTU without reads that has potential parents: R stops with an error,
    # here it is kept.
    t = OTUTable(["P", "Z"], ["a", "b"], [5 5; 0 0])
    r = lulu(t, MatchList([("Z", "P", 99.0)]))
    @test r.otu_map.parent_id == ["P", "Z"]
    @test r.otu_map.spread == [2, 0]
    @test_throws ArgumentError lulu(t, MatchList([("Z", "P", 99.0)]); minimum_ratio_type = "max")
    @test_throws ArgumentError lulu(t, MatchList([("Z", "P", 99.0)]); minimum_ratio = NaN)
    @test_throws ArgumentError lulu(t, MatchList([("Z", "P", 99.0)]); minimum_match = NaN)
    @test lulu(t, MatchList([("Z", "P", 99.0)]); minimum_ratio_type = :avg).minimum_ratio_type == "avg"
    @test lulu(t, MatchList([("Z", "P", 99.0)]); logfile = false).curated_count == 2
    # Empty table.
    e = lulu(OTUTable(String[], ["a"], zeros(Int, 0, 1)), MatchList(String[], String[], Float64[]))
    @test e.curated_count == 0 && e.discarded_count == 0 && size(e.curated_table) == (0, 1)
    @test_logs (:warn, r"looks like a fraction") lulu(t, MatchList([("Z", "P", 99.0)]); minimum_match = 0.84)
end

@testset "lulu: float tables, matrices, files, log file" begin
    c = [0.5 0.25; 0.01 0.02; 0.1 0.0]
    t = OTUTable(c; ids = ["a", "b", "c"])
    @test t.samples == ["S1", "S2"]
    r = lulu(t, MatchList([("b", "a", 99.0)]))
    @test r.curated_table isa OTUTable{Float64}
    @test r.curated_table.counts ≈ [0.51 0.27; 0.1 0.0]
    @test lulu(c, ["a", "b", "c"], MatchList([("b", "a", 99.0)])).otu_map.parent_id ==
          r.otu_map.parent_id
    mktempdir() do dir
        cd(dir) do
            r1 = lulu(t, MatchList([("b", "a", 99.0)]); logfile = true)
            logs = filter(startswith("lulu.log_"), readdir())
            @test length(logs) == 1
            @test occursin(r"^lulu\.log_\d{8}_\d{6}$", logs[1])
            @test occursin("SETTING b to be an ERROR of a", read(logs[1], String))
            write_otutable("t.tsv", t)
            write_matchlist("m.tsv", MatchList([("b", "a", 99.0)]))
            r2 = lulu("t.tsv", "m.tsv"; logfile = "my.log")
            @test r2.otu_map.parent_id == r1.otu_map.parent_id
            @test isfile("my.log")
        end
    end
    redirect_stdout(devnull) do
        @test lulu(t, MatchList([("b", "a", 99.0)]); progress = true).curated_count == 2
    end
end

@testset "lulu: sparse tables and parallel search" begin
    t = OTUTable(["D", "P", "Z", "Q"], ["a", "b", "c", "d"],
                 [1 2 0 0; 10 20 5 0; 0 0 3 0; 0 7 0 9])
    ml = MatchList([("D", "P", 90.0), ("Z", "P", 90.0), ("Q", "P", 95.0)])
    st = sparse(t)
    @test issparse(st) && !issparse(t) && st == t
    @test Matrix(st) == t && !issparse(Matrix(st))
    @test occursin("sparse, 8 stored counts", sprint(show, st))
    for kw in ((;), (storage = :sparse,), (storage = :dense,), (minimum_ratio_type = :avg,))
        a = lulu(t, ml; kw...); b = lulu(st, ml; kw...)
        @test a.otu_map.parent_id == b.otu_map.parent_id
        @test a.curated_table == b.curated_table
        @test issparse(b.curated_table.counts)
    end
    @test_throws ArgumentError lulu(t, ml; storage = :compressed)

    # Explicitly stored zeros are not presences.
    c = SparseMatrixCSC(4, 4, [1, 3, 6, 8, 9], [1, 2, 1, 2, 3, 2, 3, 4],
                        [1, 10, 2, 20, 0, 5, 3, 9])
    sz = OTUTable(["D", "P", "Z", "Q"], ["a", "b", "c", "d"], c)
    @test nnz(sz.counts) == 8 && Matrix(sz.counts) == [1 2 0 0; 10 20 5 0; 0 0 3 0; 0 0 0 9]
    r1 = lulu(sz, ml); r2 = lulu(Matrix(sz), ml)
    @test r1.otu_map.spread == r2.otu_map.spread
    @test r1.otu_map.parent_id == r2.otu_map.parent_id
    @test r1.curated_table == r2.curated_table
    @test_throws ArgumentError OTUTable(["a"], ["s"], sparse([-1;;]))

    # Ratios over the daughter's samples, iterable twice.
    reads = sparse([1 10; 0 20; 2 0; 0 0])   # samples × OTUs: daughter 1, parent 2
    r = LULU._abundance_ratios(reads, 1, 2)
    @test collect(r) == [10.0, 0.0] == collect(r)
    @test LULU._rmean(r) == 5.0
    @test LULU._relative_cooccurrence(reads, 1, 2, 2) == 0.5
    @test LULU._relative_cooccurrence(Matrix(reads), 1, 2, 2) == 0.5
    @test collect(LULU._abundance_ratios(Matrix(reads), 1, 2)) == [10.0, 0.0]

    # Float sparse table: group sums as R's long double sums.
    f = sparse([0.1 0.0; 0.2 0.5; 0.3 0.0])
    tf = OTUTable(["a", "b", "c"], ["s1", "s2"], f)
    rf = lulu(tf, MatchList([("a", "b", 99.0), ("c", "b", 99.0)]))
    @test rf.curated_table.counts isa SparseMatrixCSC{Float64}
    @test rf.curated_table.counts == lulu(Matrix(tf), MatchList([("a", "b", 99.0), ("c", "b", 99.0)])).curated_table.counts

    # The result does not depend on the parallel decomposition.
    @test LULU._resolve([1, 1, 0], 3, 2) == 1
    @test LULU._resolve([1, 0, 0], 2, 3) == 3
    @test LULU._resolve([1, 0, 0], 2, 0) == 2
end

@testset "R numerics" begin
    # as.character() of doubles in R 4.5.
    for (x, s) in [(0.68, "0.68"), (1.0, "1"), (27.88888888888889, "27.8888888888889"),
                   (1e5, "1e+05"), (123456.0, "123456"), (1e-4, "1e-04"), (0.001, "0.001"),
                   (1234567890123456.0, "1234567890123456"), (3.6864099744707345e15, "3686409974470735"), (0.1 + 0.2, "0.3"),
                   (1 / 3, "0.333333333333333"), (2 / 3, "0.666666666666667"), (100.0, "100"),
                   (0.0, "0"), (NaN, "NaN"), (Inf, "Inf"), (-2.5, "-2.5"), (1e15, "1e+15"),
                   (1e300, "1e+300"), (0.000123, "0.000123")]
        @test LULU._rformat(x) == s
    end
    @test LULU._rmean([1.0, 2.0, 4.0]) == 7 / 3
    @test LULU._rmean([0.1, 0.2, 0.3]) == 0.2
    @test LULU._rsum([0.1, 0.2, 0.3]) == 0.6   # double accumulation would give 0.6000000000000001
end
