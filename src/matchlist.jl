# Pure Julia construction of the LULU match list (all OTU sequences searched
# against each other), replacing the external BLASTN/VSEARCH step. The search
# follows the strategy of `vsearch --usearch_global --self`: candidate targets
# are ranked by the number of shared unique k-mers and aligned in that order
# until `maxaccepts` hits are accepted or `maxrejects` candidates are rejected.

function _unique_kmers(s::Vector{UInt8}, k::Int)
    out = UInt32[]
    mask = (UInt32(1) << (2k)) - UInt32(1)
    w = UInt32(0); valid = 0
    for x in s
        if x > 0x03
            valid = 0; w = UInt32(0)
            continue
        end
        w = ((w << 2) | x) & mask
        valid += 1
        valid >= k && push!(out, w)
    end
    unique!(sort!(out))
end

# Inverted index in CSR form: targets containing each k-mer.
function _kmer_index(kmers::Vector{Vector{UInt32}}, k::Int)
    nk = 1 << (2k)
    offsets = zeros(Int, nk + 1)
    for ks in kmers, w in ks
        offsets[w+2] += 1
    end
    cumsum!(offsets, offsets)
    postings = Vector{Int32}(undef, offsets[end])
    fillpos = offsets[1:nk] .+ 1
    for (t, ks) in enumerate(kmers), w in ks
        postings[fillpos[w+1]] = t
        fillpos[w+1] += 1
    end
    offsets, postings
end

function _identity(st::_AlignStats, la::Int, lb::Int, iddef::Int)
    if iddef == 0
        return 100 * st.matches / min(la, lb)
    elseif iddef == 1 || iddef == 4
        return 100 * st.matches / st.columns
    elseif iddef == 2
        den = st.columns - st.terminal_gap_columns
        return den == 0 ? 0.0 : 100 * st.matches / den
    else
        return 100 * (1 - (st.mismatches + st.gap_openings) / max(la, lb))
    end
end

# Settings shared by all queries of a search.
struct _SearchSettings
    min_identity::Float64   # fraction
    id_threshold::Float64   # percent, with a tolerance for rounding
    iddef::Int
    query_cov::Float64
    maxhits::Int
    maxaccepts::Int
    maxrejects::Int
    minwordmatches::Int
    scoring::_Scoring
end

# Diagonal band for sequences of lengths `la` and `lb`. Under iddef 1/4 an
# alignment with identity >= t has at most `slack` gap columns, so it stays
# within (slack -/+ (lb - la)) / 2 of the main diagonal. Other definitions get a
# wider band.
function _band(la::Int, lb::Int, s::_SearchSettings)
    t = s.min_identity
    slack = t > 0 ? floor(Int, (1 - t) / t * min(la, lb)) : la + lb
    delta = lb - la
    if s.iddef == 1 || s.iddef == 4
        max(-la, min(0, delta, -fld(slack - delta, 2))), min(lb, max(0, delta, fld(slack + delta, 2)))
    else
        max(-la, min(0, delta) - slack), min(lb, max(0, delta) + slack)
    end
end

# Exact test on the lengths alone: false if no alignment can be accepted.
_length_compatible(la::Int, lb::Int, s::_SearchSettings) =
    la > 0 && lb > 0 && min(la, lb) >= s.query_cov * la &&
    !((s.iddef == 1 || s.iddef == 4) && 100 * min(la, lb) / max(la, lb) < s.id_threshold)

# Alignment of sequences `x` < `y`. A pair is always aligned in this canonical
# orientation, so (q, t) and (t, q) share one alignment (identity and number of
# aligned pairs are symmetric). An accepted alignment that touches the band edge
# may not be optimal: it is recomputed without band.
function _pair_stats(ws::_Workspace, seqs, x::Int, y::Int, s::_SearchSettings)
    a = seqs[x]; b = seqs[y]
    la = length(a); lb = length(b)
    klo, khi = _band(la, lb, s)
    st = _align(ws, a, b, klo, khi, s.scoring)
    if st.touches_band_edge && _identity(st, la, lb, s.iddef) >= s.id_threshold
        st = _align(ws, a, b, -la, lb, s.scoring)
    end
    st
end

# Alignments shared between the two directions of a pair, in lock-striped
# dictionaries. An entry is removed when it is used the second time.
struct _PairCache
    tables::Vector{Dict{Int,_AlignStats}}
    locks::Vector{ReentrantLock}
end

_PairCache(nstripes::Int) =
    _PairCache([Dict{Int,_AlignStats}() for _ in 1:nstripes], [ReentrantLock() for _ in 1:nstripes])

function _get_or_align!(compute, cache::_PairCache, key::Int)
    stripe = mod1(key, length(cache.tables))
    table = cache.tables[stripe]
    lk = cache.locks[stripe]
    cached = lock(() -> pop!(table, key, nothing), lk)
    cached === nothing || return cached
    st = compute()   # outside the lock; a concurrent duplicate gives the same result
    lock(lk) do
        haskey(table, key) ? delete!(table, key) : (table[key] = st)
    end
    st
end

# Per-task memory of the search.
struct _SearchState
    ws::_Workspace
    shared::Vector{Int32}    # shared unique k-mers with the current query
    touched::Vector{Int}
    candidates::Vector{Int}
end

_SearchState(n::Int) = _SearchState(_Workspace(), zeros(Int32, n), Int[], Int[])

# Targets sharing enough unique k-mers with query `q`, by decreasing number of
# shared k-mers (ties: index).
function _rank_candidates!(st::_SearchState, q::Int, kmers, offsets, postings,
                           s::_SearchSettings)
    shared = st.shared; touched = empty!(st.touched); candidates = empty!(st.candidates)
    for w in kmers[q], p in (offsets[w+1]+1):offsets[w+2]
        tg = postings[p]
        shared[tg] == 0 && push!(touched, tg)
        shared[tg] += Int32(1)
    end
    threshold = min(s.minwordmatches, length(kmers[q]))
    if threshold == 0
        append!(candidates, eachindex(shared))
    else
        append!(candidates, Iterators.filter(tg -> shared[tg] >= threshold, touched))
    end
    sort!(candidates; alg = MergeSort, by = tg -> (-shared[tg], tg))
    for tg in touched
        shared[tg] = 0
    end
    candidates
end

# Accepted hits (target, identity) of query `q`.
function _search_query(st::_SearchState, q::Int, seqs, kmers, offsets, postings,
                       cache::_PairCache, s::_SearchSettings)
    n = length(seqs)
    la = length(seqs[q])
    accepted = Tuple{Int,Float64}[]
    rejected = 0
    for tg in _rank_candidates!(st, q, kmers, offsets, postings, s)
        tg == q && continue
        lb = length(seqs[tg])
        ok = false
        if _length_compatible(la, lb, s)
            x, y = minmax(q, tg)
            stats = _get_or_align!(cache, (x - 1) * n + y) do
                _pair_stats(st.ws, seqs, x, y, s)
            end
            id = _identity(stats, la, lb, s.iddef)
            ok = id >= s.id_threshold && stats.matches + stats.mismatches >= s.query_cov * la
            ok && push!(accepted, (tg, id))
        end
        if ok
            s.maxaccepts > 0 && length(accepted) >= s.maxaccepts && break
        else
            rejected += 1
            s.maxrejects > 0 && rejected >= s.maxrejects && break
        end
    end
    sort!(accepted; alg = MergeSort, by = hit -> -hit[2])
    s.maxhits > 0 && length(accepted) > s.maxhits && resize!(accepted, s.maxhits)
    accepted
end

"""
    build_matchlist(ids, sequences; min_identity = 84, iddef = 1, query_cov = 0.9,
                    maxhits = 10, maxaccepts = 0, maxrejects = 32, wordlength = 8,
                    minwordmatches = 12, match = 2, mismatch = -4, gap_open = 20,
                    gap_extend = 2, terminal_gap_open = 2, terminal_gap_extend = 1)
        -> MatchList
    build_matchlist(fasta_path; kwargs...)

Build the match list required by [`lulu`](@ref) by searching every sequence
(query) against all the others (targets), in pure Julia. It replaces the
external command suggested in the LULU README,

    vsearch --usearch_global OTU_sequences.fasta --db OTU_sequences.fasta --self
            --id .84 --iddef 1 --userout match_list.txt
            -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

whose settings are the defaults here. The search strategy is the same as in
VSEARCH, but the implementation is independent: hits and identities are
close to, not identical with, those of VSEARCH or BLASTN.

For each query:

1. targets are ranked by the number of unique k-mers (`wordlength`) shared
   with the query, decreasing; targets sharing fewer than `minwordmatches`
   (or all of the query's unique k-mers, if fewer) are ignored;
2. candidates are aligned in that order with a banded global alignment
   (affine gaps; `match`, `mismatch`, gap penalties as positive numbers; a gap
   of length `n` costs `open + n * extend`; gaps at the ends of either
   sequence use the terminal penalties);
3. a candidate is accepted if its identity (percent) is `>= min_identity` and
   the fraction of query bases aligned to target bases is `>= query_cov`;
   otherwise it is rejected. The search stops after `maxaccepts` accepts (0:
   no limit) or `maxrejects` rejects (0: no limit);
4. the accepted hits are sorted by decreasing identity and the first
   `maxhits` (0: all) are kept.

With `maxrejects = 0` and `minwordmatches = 0` every pair is aligned
(exhaustive search; slow for large data sets).

Identity definitions (`iddef`, as in VSEARCH):
`0`: matches / length of the shorter sequence; `1` (default) and `4`:
matches / alignment columns; `2`: matches / alignment columns excluding
terminal gaps; `3`: `1 - (mismatches + gap openings) / length of the longer
sequence`. Bases other than A, C, G, T/U never count as matches.

The alignment is restricted to a diagonal band wide enough for every
alignment with identity `>= min_identity` under `iddef` 1. Queries run in
parallel on `Threads.nthreads()` threads; the output does not depend on the
number of threads. Self-hits are never reported.
"""
function build_matchlist(ids::AbstractVector{<:AbstractString},
                         sequences::AbstractVector{<:AbstractString};
                         min_identity::Real = 84, iddef::Integer = 1, query_cov::Real = 0.9,
                         maxhits::Integer = 10, maxaccepts::Integer = 0,
                         maxrejects::Integer = 32, wordlength::Integer = 8,
                         minwordmatches::Integer = 12, match::Integer = 2,
                         mismatch::Integer = -4, gap_open::Integer = 20,
                         gap_extend::Integer = 2, terminal_gap_open::Integer = 2,
                         terminal_gap_extend::Integer = 1)
    length(ids) == length(sequences) ||
        throw(DimensionMismatch("ids and sequences must have the same length"))
    allunique(ids) || throw(ArgumentError("sequence ids must be unique"))
    0 <= min_identity <= 100 ||
        throw(ArgumentError("min_identity must be a percentage between 0 and 100"))
    0 < min_identity <= 1 &&
        @warn "min_identity = $min_identity is a percentage: did you mean $(100 * min_identity)?"
    iddef in 0:4 || throw(ArgumentError("iddef must be 0, 1, 2, 3 or 4"))
    0 <= query_cov <= 1 || throw(ArgumentError("query_cov must be between 0 and 1"))
    3 <= wordlength <= 12 || throw(ArgumentError("wordlength must be between 3 and 12"))
    min(maxhits, maxaccepts, maxrejects, minwordmatches) >= 0 ||
        throw(ArgumentError("maxhits, maxaccepts, maxrejects and minwordmatches must be >= 0"))
    min(gap_open, gap_extend, terminal_gap_open, terminal_gap_extend) >= 0 ||
        throw(ArgumentError("gap penalties must be given as non-negative numbers"))

    settings = _SearchSettings(min_identity / 100, Float64(min_identity) - 1e-9, iddef,
                               query_cov, maxhits, maxaccepts, maxrejects, minwordmatches,
                               _Scoring(match, mismatch, gap_open, gap_extend,
                                        terminal_gap_open, terminal_gap_extend))
    k = Int(wordlength)
    seqs = [_encode(s) for s in sequences]
    n = length(seqs)
    kmers = [_unique_kmers(s, k) for s in seqs]
    offsets, postings = _kmer_index(kmers, k)
    cache = _PairCache(64)
    results = Vector{Vector{Tuple{Int,Float64}}}(undef, n)
    chunks = Iterators.partition(1:n, max(1, cld(n, 8 * Threads.nthreads())))
    tasks = map(chunks) do range
        Threads.@spawn begin
            state = _SearchState(n)
            for q in range
                results[q] = _search_query(state, q, seqs, kmers, offsets, postings, cache,
                                           settings)
            end
        end
    end
    foreach(wait, tasks)

    query = String[]; hit = String[]; ident = Float64[]
    for q in 1:n, (tg, id) in results[q]
        push!(query, ids[q]); push!(hit, ids[tg]); push!(ident, id)
    end
    MatchList(query, hit, ident)
end

function build_matchlist(fasta_path::AbstractString; kwargs...)
    ids, seqs = read_fasta(fasta_path)
    build_matchlist(ids, seqs; kwargs...)
end

"""
    build_matchlist(records::AbstractVector{<:Pair}; kwargs...)

Match list from `id => sequence` pairs.
"""
build_matchlist(records::AbstractVector{<:Pair}; kwargs...) =
    build_matchlist([string(first(r)) for r in records], [string(last(r)) for r in records];
                    kwargs...)
