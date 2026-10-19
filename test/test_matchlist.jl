using Random

# Reference: unbanded Gotoh alignment score with the same gap rules.
function brute_score(a, b, sc)
    la, lb = length(a), length(b)
    NEG = typemin(Int) ÷ 4
    H = fill(NEG, la + 1, lb + 1); E = copy(H); F = copy(H)
    H[1, 1] = 0
    for i in 0:la, j in 0:lb
        (i == 0 && j == 0) && continue
        hterm = i == 0 || i == la
        vterm = j == 0 || j == lb
        e = j > 0 ? max(H[i+1, j] - (hterm ? sc.terminal_gap_open + sc.terminal_gap_extend :
                                              sc.gap_open + sc.gap_extend),
                        E[i+1, j] - (hterm ? sc.terminal_gap_extend : sc.gap_extend)) : NEG
        f = i > 0 ? max(H[i, j+1] - (vterm ? sc.terminal_gap_open + sc.terminal_gap_extend :
                                              sc.gap_open + sc.gap_extend),
                        F[i, j+1] - (vterm ? sc.terminal_gap_extend : sc.gap_extend)) : NEG
        d = (i > 0 && j > 0) ?
            H[i, j] + ((a[i] == b[j] && a[i] < 4) ? sc.match : sc.mismatch) : NEG
        E[i+1, j+1] = e; F[i+1, j+1] = f; H[i+1, j+1] = max(d, e, f)
    end
    H[la+1, lb+1]
end

# Score of the path left in the workspace by `_align` (operations in reverse).
function path_score(ops, a, b, sc)
    la, lb = length(a), length(b)
    i = j = 0
    s = 0
    prev = 0x00
    for o in Iterators.reverse(ops)
        if o == UInt8('M')
            i += 1; j += 1
            s += (a[i] == b[j] && a[i] < 4) ? sc.match : sc.mismatch
        elseif o == UInt8('D')
            term = i == 0 || i == la
            s -= (term ? sc.terminal_gap_extend : sc.gap_extend) +
                 (prev == o ? 0 : (term ? sc.terminal_gap_open : sc.gap_open))
            j += 1
        else
            term = j == 0 || j == lb
            s -= (term ? sc.terminal_gap_extend : sc.gap_extend) +
                 (prev == o ? 0 : (term ? sc.terminal_gap_open : sc.gap_open))
            i += 1
        end
        prev = o
    end
    @assert i == la && j == lb
    s
end

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

@testset "Banded alignment" begin
    rng = MersenneTwister(42)
    sc = LULU._Scoring(2, -4, 20, 2, 2, 1)
    ws = LULU._Workspace()
    for trial in 1:300
        a0 = String(rand(rng, collect("ACGT"), rand(rng, 1:70)))
        b0 = rand(rng) < 0.3 ? String(rand(rng, collect("ACGTN"), rand(rng, 1:70))) :
             mutate(rng, a0; subs = rand(rng, 0:8), indels = rand(rng, 0:4))
        isempty(b0) && (b0 = "A")
        a = LULU._encode(a0); b = LULU._encode(b0)
        st = LULU._align(ws, a, b, -length(a), length(b), sc)
        @test path_score(ws.ops, a, b, sc) == brute_score(a, b, sc)
        @test st.columns == length(ws.ops)
        @test st.matches + st.mismatches <= min(length(a), length(b))
        @test !st.touches_band_edge
        # 16- and 32-bit scores give the same alignment.
        @test LULU._align(ws, a, b, -length(a), length(b), sc, Int16) ==
              LULU._align(ws, a, b, -length(a), length(b), sc, Int32)
    end
    @test LULU._score_type(300, 300, sc) === Int16
    @test LULU._score_type(5_000, 5_000, sc) === Int32
    # Long sequences (32-bit scores): same identity as computed by the band.
    long = String(rand(rng, collect("ACGT"), 3_000))
    la = LULU._encode(long); lb = LULU._encode(mutate(rng, long; subs = 30, indels = 3))
    st32 = LULU._align(ws, la, lb, -length(la), length(lb), sc)
    @test st32.matches >= 2_900
    @test path_score(ws.ops, la, lb, sc) >= 2 * 2_900 - 1_000
    a = LULU._encode("ACGTACGTAC"); b = LULU._encode("ACGTTCGTAC")
    st = LULU._align(ws, a, b, -2, 2, sc)
    @test (st.matches, st.mismatches, st.columns, st.terminal_gap_columns, st.gap_openings) ==
          (9, 1, 10, 0, 0)
    @test LULU._identity(st, 10, 10, 1) == 90
    # Terminal gaps: b is a with 5 extra bases at the end.
    a = LULU._encode("ACGTACGTACGTACGTACGT"); b = LULU._encode("ACGTACGTACGTACGTACGTGGGGG")
    st = LULU._align(ws, a, b, -20, 25, sc)
    @test (st.matches, st.columns, st.terminal_gap_columns, st.gap_openings) == (20, 25, 5, 1)
    @test LULU._identity(st, 20, 25, 0) == 100
    @test LULU._identity(st, 20, 25, 1) == 80
    @test LULU._identity(st, 20, 25, 2) == 100
    @test LULU._identity(st, 20, 25, 3) ≈ 96
    @test LULU._identity(st, 20, 25, 4) == 80
    # N never matches.
    st = LULU._align(ws, LULU._encode("ANGT"), LULU._encode("ANGT"), -4, 4, sc)
    @test st.matches == 3 && st.mismatches == 1
    @test LULU._encode("acgu") == UInt8[0, 1, 2, 3]
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
