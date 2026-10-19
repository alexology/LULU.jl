using Random

function mutate(rng, s; subs = 0, indels = 0)
    v = collect(s)
    for _ in 1:subs
        k = rand(rng, eachindex(v)); v[k] = rand(rng, setdiff("ACGT", v[k]))
    end
    for _ in 1:indels
        (isempty(v) || rand(rng, Bool)) ? insert!(v, rand(rng, 1:length(v)+1), rand(rng, "ACGT")) :
                          deleteat!(v, rand(rng, eachindex(v)))
    end
    String(v)
end

@testset "build_matchlist" begin
    rng = MersenneTwister(7)
    base = String(rand(rng, collect("ACGT"), 200))
    seqs = ["base" => base,
            "snp1" => mutate(rng, base; subs = 1),
            "snp10" => mutate(rng, base; subs = 10),
            "short" => base[1:170],
            "other" => String(rand(rng, collect("ACGT"), 200)),
            "indel" => base[1:100] * base[104:end]]
    ml = build_matchlist(seqs)
    pairs = Dict((ml.query[k], ml.hit[k]) => ml.match[k] for k in 1:length(ml))
    @test pairs[("base", "snp1")] == 99.5
    @test pairs[("snp1", "base")] == 99.5
    @test pairs[("base", "snp10")] == 95
    @test pairs[("short", "base")] == 85          # 170 / 200 columns (iddef 1)
    @test !haskey(pairs, ("base", "short"))       # query coverage 0.85 < 0.9
    @test pairs[("base", "indel")] == 197 / 200 * 100
    @test !any(==("other"), ml.query) && !any(==("other"), ml.hit)
    @test all(ml.query .!= ml.hit)
    @test all(ml.match .>= 84)
    # Hits of each query sorted by decreasing identity.
    for q in unique(ml.query)
        @test issorted(ml.match[ml.query.==q]; rev = true)
    end

    ml2 = build_matchlist(seqs; iddef = 2)
    p2 = Dict((ml2.query[k], ml2.hit[k]) => ml2.match[k] for k in 1:length(ml2))
    @test p2[("short", "base")] == 100            # terminal gaps excluded
    ml0 = build_matchlist(seqs; query_cov = 0.8)
    @test any(k -> ml0.query[k] == "base" && ml0.hit[k] == "short", 1:length(ml0))
    ml1 = build_matchlist(seqs; maxhits = 1)
    @test all(q -> count(==(q), ml1.query) == 1, unique(ml1.query))
    @test length(build_matchlist(seqs; min_identity = 99.9)) == 0
    @test build_matchlist(seqs; maxrejects = 0, minwordmatches = 0, maxhits = 0) ==
          build_matchlist(seqs; maxhits = 0)
    @test build_matchlist(seqs; maxaccepts = 1).query |> unique |>
          qs -> all(q -> count(==(q), build_matchlist(seqs; maxaccepts = 1).query) == 1, qs)

    mktempdir() do dir
        f = joinpath(dir, "s.fasta")
        open(f, "w") do io
            for (id, s) in seqs
                println(io, ">", id, "\n", s)
            end
        end
        @test build_matchlist(f) == ml
    end
    @test build_matchlist(first.(seqs), last.(seqs)) == ml

    @test_throws ArgumentError build_matchlist(["a", "a"], ["ACGT", "ACGT"])
    @test_throws DimensionMismatch build_matchlist(["a"], ["ACGT", "ACGT"])
    @test_throws ArgumentError build_matchlist(seqs; iddef = 5)
    @test_throws ArgumentError build_matchlist(seqs; min_identity = 140)
    @test_throws ArgumentError build_matchlist(seqs; query_cov = 2)
    @test_throws ArgumentError build_matchlist(seqs; wordlength = 20)
    @test_throws ArgumentError build_matchlist(seqs; gap_open = -1)
    @test_logs (:warn, r"percentage") build_matchlist(seqs; min_identity = 0.9)

    # End-to-end: the match list feeds lulu.
    t = OTUTable(first.(seqs), ["s1", "s2", "s3"],
                 [1000 800 900; 10 8 9; 20 0 5; 3 2 1; 50 60 70; 5 5 5])
    r = lulu(t, ml)
    @test r.curated_otus == ["base", "other"]
end
