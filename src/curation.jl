# The LULU curation algorithm: port of `lulu()` from R/Functions.R of the R
# package lulu 0.1.0 (Tobias Guldberg Frøslev). The control flow, the
# comparisons (`>` or `>=`) and the log messages follow the R code exactly.
#
# Internally OTUs are identified by their position in processing order, and
# the reads are stored as a samples × OTUs matrix (dense or sparse) so that
# each OTU is a contiguous column.

# ---------------------------------------------------------------------------
# Settings

# Validated thresholds of a LULU run.
struct _Criteria
    minimum_ratio_type::Symbol   # :min or :avg
    minimum_ratio::Float64
    minimum_match::Float64
    minimum_relative_cooccurence::Float64
end

function _criteria(minimum_ratio_type, minimum_ratio, minimum_match,
                   minimum_relative_cooccurence)
    type = Symbol(minimum_ratio_type)
    type in (:min, :avg) ||
        throw(ArgumentError("minimum_ratio_type must be \"min\" or \"avg\", got \"$type\""))
    for (name, value) in (("minimum_ratio", minimum_ratio), ("minimum_match", minimum_match),
                          ("minimum_relative_cooccurence", minimum_relative_cooccurence))
        isnan(value) && throw(ArgumentError("$name must be a number, got NaN"))
    end
    _Criteria(type, minimum_ratio, minimum_match, minimum_relative_cooccurence)
end

# Statistic of the parent/daughter read ratios: R `min()` or `mean()`.
_ratio_statistic(c::_Criteria) = c.minimum_ratio_type === :avg ? _rmean : minimum
_ratio_label(c::_Criteria) =
    c.minimum_ratio_type === :avg ? "\n------mean avg abundance: " : "\n------min avg abundance: "

function _warn_if_fraction(minimum_match::Float64, matchlist::MatchList)
    if 0 < minimum_match <= 1 && any(>(1), matchlist.match)
        @warn "minimum_match = $minimum_match looks like a fraction, but LULU expects a " *
              "percentage (e.g. 84) like the similarities in the match list"
    end
    nothing
end

# ---------------------------------------------------------------------------
# Log: every function has a no-op method for `nothing` (no log requested), so
# the algorithm calls them unconditionally and nothing is formatted when no
# log is written.

_log(::Nothing, args...) = nothing
_log(io::IO, args...) = (print(io, args...); nothing)

# A number printed as R's `paste0` does.
_log_number(::Nothing, prefix, x) = nothing
_log_number(io::IO, prefix, x) = (print(io, prefix, _rformat(x)); nothing)

# R's `cat(paste0(prefix, xs))`: one element per entry, separated by a space;
# an empty vector gives the prefix alone.
_log_vector(::Nothing, prefix, xs) = nothing
function _log_vector(io::IO, prefix, xs)
    isempty(xs) && return _log(io, prefix)
    for (k, x) in enumerate(xs)
        k > 1 && print(io, ' ')
        print(io, prefix, x)
    end
    nothing
end

# Runs `f(io)` with the log destination: `nothing`/`false` (no log), an `IO`,
# a file path (appended to) or `true` (file named as in R).
_with_log(f, ::Nothing, start) = f(nothing)
_with_log(f, io::IO, start) = f(io)
_with_log(f, path::AbstractString, start) = open(f, path, "a")
_with_log(f, flag::Bool, start) =
    flag ? _with_log(f, "lulu.log_" * Dates.format(start, "yyyymmdd_HHMMSS"), start) :
           f(nothing)

# ---------------------------------------------------------------------------
# Statistics and processing order

# Total read count of an OTU, computed as R's `rowSums` (zeros do not change
# the sum, so the stored values of a sparse column suffice).
_total_reads(x::AbstractVector{<:Integer}) = _rsum_exact(x)
_total_reads(x::AbstractVector) = _rsum(x)

# Total read count and spread (number of samples with reads) of every OTU of
# an OTUs × samples table. For a sparse table the stored values of each OTU are
# read from the transposed matrix (one column per OTU).
function _otu_statistics(counts::AbstractMatrix)
    total = [_total_reads(row) for row in eachrow(counts)]
    spread = [count(>(0), row) for row in eachrow(counts)]
    total, spread
end

function _otu_statistics(counts::SparseMatrixCSC)
    by_otu = _otu_columns(counts)
    stored(i) = view(nonzeros(by_otu), nzrange(by_otu, i))
    total = [_total_reads(stored(i)) for i in axes(by_otu, 2)]
    spread = [count(>(0), stored(i)) for i in axes(by_otu, 2)]
    total, spread
end

# Processing order: decreasing spread, then decreasing total; ties keep the
# input order (R `order(spread, total, decreasing = TRUE)`, a stable sort).
_processing_order(spread, total) =
    sortperm(collect(zip(spread, total)); rev = true, alg = MergeSort)

# Rank by decreasing total, ties in processing order
# (R `rank(-total, ties.method = "first")`).
_total_rank(total) = invperm(sortperm(total; rev = true, alg = MergeSort))

# Match-list entries used by LULU: no "*" hits, no self-hits, similarity
# strictly greater than `minimum_match`.
_usable_hit(query, hit, match, minimum_match) =
    hit != "*" && hit != query && match > minimum_match

# Hit ids of each OTU (in processing order), in match-list order, for the log
# (which lists them all, including ids missing from the table).
function _hits_by_otu(matchlist::MatchList, position::Dict{String,Int}, minimum_match)
    hits = [String[] for _ in 1:length(position)]
    for (query, hit, match) in zip(matchlist.query, matchlist.hit, matchlist.match)
        _usable_hit(query, hit, match, minimum_match) || continue
        i = get(position, query, 0)
        i > 0 && push!(hits[i], hit)
    end
    hits
end

# Positions of the hits of each OTU (0 for ids missing from the table), in
# match-list order. The string lookups run in parallel.
function _hit_positions(matchlist::MatchList, position::Dict{String,Int}, minimum_match)
    m = length(matchlist)
    query_pos = zeros(Int, m)
    hit_pos = zeros(Int, m)
    _foreach_chunk(1:m) do range
        for k in range
            query = matchlist.query[k]; hit = matchlist.hit[k]
            _usable_hit(query, hit, matchlist.match[k], minimum_match) || continue
            query_pos[k] = get(position, query, 0)
            query_pos[k] > 0 && (hit_pos[k] = get(position, hit, 0))
        end
    end
    hits = [Int[] for _ in 1:length(position)]
    for k in 1:m
        query_pos[k] > 0 && push!(hits[query_pos[k]], hit_pos[k])
    end
    hits
end

# Runs `f(subrange)` on chunks of `range`, on all threads when there are more
# than one. `f` must only write to disjoint locations.
function _foreach_chunk(f, range::AbstractUnitRange; chunks_per_thread::Int = 16)
    nchunks = min(length(range), chunks_per_thread * Threads.nthreads())
    if Threads.nthreads() == 1 || nchunks <= 1
        f(range)
    else
        tasks = [Threads.@spawn f(r) for r in Iterators.partition(range, cld(length(range), nchunks))]
        foreach(wait, tasks)
    end
    nothing
end

# Positions of the hits that are at most `last` (OTUs with a spread at least
# as large as the daughter's), sorted and unique (R `which(... %in% hits)`).
# `hit_positions` holds 0 for ids missing from the table.
function _potential_parents!(out::Vector{Int}, hit_positions::Vector{Int}, last::Int)
    empty!(out)
    for p in hit_positions
        1 <= p <= last && push!(out, p)
    end
    unique!(sort!(out))
end

# Number of OTUs with a spread >= that of OTU `line` (spread is non-increasing).
_last_relevant(spread, line::Int) = searchsortedlast(spread, spread[line]; rev = true)

# ---------------------------------------------------------------------------
# Kernels on the reads (samples × OTUs; one column per OTU, processing order)

# Fraction of the samples of daughter `d` in which candidate `p` is present.
function _relative_cooccurrence(reads::AbstractMatrix, d::Int, p::Int, spread::Int)
    shared = 0
    @inbounds for s in axes(reads, 1)
        shared += (reads[s, d] > 0) & (reads[s, p] > 0)
    end
    shared / spread
end

# Sparse columns: merge of the two sorted lists of samples.
function _relative_cooccurrence(reads::SparseMatrixCSC, d::Int, p::Int, spread::Int)
    rows = rowvals(reads); vals = nonzeros(reads)
    i, iend = first(nzrange(reads, d)), last(nzrange(reads, d))
    j, jend = first(nzrange(reads, p)), last(nzrange(reads, p))
    shared = 0
    @inbounds while i <= iend && j <= jend
        ri = rows[i]; rj = rows[j]
        if ri == rj
            shared += (vals[i] > 0) & (vals[j] > 0)
            i += 1; j += 1
        elseif ri < rj
            i += 1
        else
            j += 1
        end
    end
    shared / spread
end

# Parent/daughter read ratios over the samples where the daughter is present,
# in sample order. The result can be iterated twice (needed by `_rmean`).
_abundance_ratios(reads::AbstractMatrix, d::Int, p::Int) =
    (Float64(reads[s, p]) / Float64(reads[s, d]) for s in axes(reads, 1) if reads[s, d] > 0)
_abundance_ratios(reads::SparseMatrixCSC, d::Int, p::Int) = _SparseRatios(reads, d, p)

struct _SparseRatios{T,Ti}
    reads::SparseMatrixCSC{T,Ti}
    d::Int
    p::Int
end

Base.IteratorSize(::Type{<:_SparseRatios}) = Base.SizeUnknown()
Base.eltype(::Type{<:_SparseRatios}) = Float64

function Base.iterate(r::_SparseRatios,
                      state = (first(nzrange(r.reads, r.d)), first(nzrange(r.reads, r.p))))
    i, j = state
    rows = rowvals(r.reads); vals = nonzeros(r.reads)
    iend = last(nzrange(r.reads, r.d)); jend = last(nzrange(r.reads, r.p))
    @inbounds while i <= iend
        if vals[i] > 0                       # stored zeros are not presences
            s = rows[i]
            while j <= jend && rows[j] < s
                j += 1
            end
            parent = (j <= jend && rows[j] == s) ? vals[j] : zero(eltype(vals))
            return Float64(parent) / Float64(vals[i]), (i + 1, j)
        end
        i += 1
    end
    nothing
end

# Reads used by the kernels: samples × OTUs (one column per OTU), built from
# the OTUs × samples `table` in the requested storage.
function _kernel_reads(table::AbstractMatrix, storage::Symbol)
    if storage === :auto   # sparse unless at least a quarter of the counts are non-zero
        storage = count(!iszero, table) < length(table) / 4 ? :sparse : :dense
    end
    storage === :sparse ? _otu_columns(sparse(table)) : permutedims(table)
end
_kernel_reads(table::SparseMatrixCSC, storage::Symbol) =
    storage === :dense ? permutedims(Matrix(table)) : _otu_columns(table)

# ---------------------------------------------------------------------------
# Parent search

# Whether `candidate` satisfies both criteria for daughter `line`, with log
# messages.
function _is_parent(logio, reads, line::Int, candidate::Int, spread::Int, c::_Criteria,
                    statistic::F) where {F}
    cooccurrence = _relative_cooccurrence(reads, line, candidate, spread)
    _log_number(logio, "\n------relative cooccurence: ", cooccurrence)
    cooccurrence >= c.minimum_relative_cooccurence || return false
    _log(logio, " which is sufficient!")
    ratio = statistic(_abundance_ratios(reads, line, candidate))
    _log_number(logio, _ratio_label(c), ratio)
    ratio > c.minimum_ratio || return false
    _log(logio, " which is OK!")
    true
end

# Parent of `line` given the first accepted candidate `found` (0: none). A
# parent processed earlier passes on its own parent (chains end at the ultimate
# parent); a parent processed later (same spread, smaller total) is used as it
# is, as in R.
_resolve(parents::Vector{Int}, line::Int, found::Int) =
    found == 0 ? line : found < line ? parents[found] : found

# First candidate accepted as parent of `line` (0: none), without log.
function _first_accepted(reads, line::Int, candidates, spread, c::_Criteria,
                         statistic::F) where {F}
    for candidate in candidates
        _is_parent(nothing, reads, line, candidate, spread[line], c, statistic) &&
            return candidate
    end
    0
end

# Parallel search: whether a candidate is accepted does not depend on the
# parents assigned before, so the first accepted candidate of every OTU is
# found independently; the chains are then resolved in processing order.
function _assign_parents_parallel(reads, spread, hit_positions, c::_Criteria,
                                  statistic::F) where {F}
    n = length(spread)
    found = zeros(Int, n)
    _foreach_chunk(1:n) do range
        candidates = Int[]
        for line in range
            _potential_parents!(candidates, hit_positions[line], _last_relevant(spread, line))
            found[line] = _first_accepted(reads, line, candidates, spread, c, statistic)
        end
    end
    parents = Vector{Int}(undef, n)
    for line in 1:n
        parents[line] = _resolve(parents, line, found[line])
    end
    parents
end

# Sequential search with the log written by R. Candidates after the accepted
# one are only logged, as in R.
function _assign_parents_logged(logio::IO, reads, ids, spread, hits, hit_positions,
                                c::_Criteria, statistic::F) where {F}
    n = length(ids)
    # hits: hit ids of each OTU as strings, listed in the log.
    parents = Vector{Int}(undef, n)
    candidates = Int[]
    for line in 1:n
        _log(logio, "\n####processing: ", ids[line], " #####")
        _log_vector(logio, "\n---hits: ", hits[line])
        _potential_parents!(candidates, hit_positions[line], _last_relevant(spread, line))
        _log_vector(logio, "\n---potential parent: ", view(ids, candidates))
        found = 0
        for candidate in candidates
            _log(logio, "\n------checking: ", ids[candidate])
            found == 0 || continue
            _is_parent(logio, reads, line, candidate, spread[line], c, statistic) || continue
            found = candidate
            _log(logio, "\nSETTING ", ids[line], " to be an ERROR of ",
                 ids[_resolve(parents, line, found)], "\n")
        end
        found == 0 && _log(logio, "\nNo parent found!\n")
        parents[line] = _resolve(parents, line, found)
    end
    parents
end

_assign_parents(::Nothing, reads, ids, spread, matchlist, position, hit_positions, c,
                statistic) =
    _assign_parents_parallel(reads, spread, hit_positions, c, statistic)
_assign_parents(logio::IO, reads, ids, spread, matchlist, position, hit_positions, c,
                statistic) =
    _assign_parents_logged(logio, reads, ids, spread,
                           _hits_by_otu(matchlist, position, c.minimum_match),
                           hit_positions, c, statistic)

# ---------------------------------------------------------------------------
# Result

# Column sums of the rows `rows` of `c`, as dplyr `summarise(sum)`.
_group_sum(c::AbstractMatrix{T}, rows, j) where {T<:Integer} =
    T(sum(i -> Int128(c[i, j]), rows; init = Int128(0)))
_group_sum(c::AbstractMatrix{T}, rows, j) where {T<:AbstractFloat} =
    T(_rsum(c[i, j] for i in rows))
_group_sum(c::AbstractMatrix{T}, rows, j) where {T<:Real} =
    sum(c[i, j] for i in rows; init = zero(T))

# Rows of `table` (OTUs × samples, processing order) summed by group; groups
# are the sorted parent ids.
function _merge_by_parent(table::AbstractMatrix, group::Vector{Int}, ngroups::Int)
    members = [Int[] for _ in 1:ngroups]
    for (i, g) in enumerate(group)
        push!(members[g], i)
    end
    [_group_sum(table, rows, j) for rows in members, j in axes(table, 2)]
end

# Sparse table: per sample, the stored values are added to their group in row
# (processing) order; zeros do not change the sums.
_accumulator(::Type{<:Integer}) = Int128
_accumulator(::Type{<:AbstractFloat}) = BigFloat   # R's long double (see rnumeric.jl)
_accumulator(::Type{T}) where {T} = T

function _merge_by_parent(table::SparseMatrixCSC{T,Ti}, group::Vector{Int},
                          ngroups::Int) where {T,Ti}
    setprecision(BigFloat, _LDOUBLE_PRECISION) do
        A = _accumulator(T)
        acc = [zero(A) for _ in 1:ngroups]
        used = falses(ngroups)
        touched = Int[]
        I = Ti[]; J = Ti[]; V = T[]
        rows = rowvals(table); vals = nonzeros(table)
        for j in axes(table, 2)
            for k in nzrange(table, j)
                g = group[rows[k]]
                used[g] || (used[g] = true; push!(touched, g))
                acc[g] += vals[k]
            end
            for g in sort!(touched)
                push!(I, g); push!(J, j); push!(V, T(acc[g]))
                acc[g] = zero(A); used[g] = false
            end
            empty!(touched)
        end
        sparse(I, J, V, ngroups, size(table, 2))
    end
end

function _lulu_result(table::AbstractMatrix{T}, ids, samples, total, spread, parents,
                      c::_Criteria, runtime) where {T}
    parent_ids = ids[parents]
    curated = [p == i ? "parent" : "merged" for (i, p) in enumerate(parents)]
    curated_otus = sort!(unique(parent_ids))
    group_of = Dict(id => k for (k, id) in enumerate(curated_otus))
    group = [group_of[id] for id in parent_ids]
    curated_counts = _merge_by_parent(table, group, length(curated_otus))
    kept = Set(curated_otus)
    discarded_otus = filter(!in(kept), ids)
    LuluResult(OTUTable(copy(curated_otus), copy(samples), curated_counts),
               length(curated_otus), curated_otus,
               length(discarded_otus), discarded_otus,
               runtime, c.minimum_match, c.minimum_relative_cooccurence,
               String(c.minimum_ratio_type), c.minimum_ratio,
               OTUMap(ids, total, spread, parent_ids, curated, _total_rank(total)),
               OTUTable(ids, copy(samples), table))
end

# ---------------------------------------------------------------------------
# Public interface

"""
    lulu(otutable::OTUTable, matchlist::MatchList;
         minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84,
         minimum_relative_cooccurence = 0.95, logfile = nothing, progress = false,
         storage = :auto)
        -> LuluResult

Post-clustering curation of an OTU table with the LULU algorithm
(Frøslev et al. 2017). Port of the R function `lulu::lulu`.

OTUs are processed in order of decreasing spread (number of samples where the
OTU is present), then decreasing total read count (ties keep the input order).
Each OTU is treated as a potential daughter. Its potential parents are the OTUs
listed as its hits in `matchlist` with similarity `> minimum_match` (percent)
and with a spread at least as large as its own. They are tested in processing
order; the first one that satisfies both criteria becomes the parent:

1. relative co-occurrence: the fraction of the daughter's samples in which the
   parent is also present is `>= minimum_relative_cooccurence`;
2. abundance ratio: parent/daughter read ratio over the daughter's samples,
   its minimum (`minimum_ratio_type = "min"`) or mean (`"avg"`), is
   `> minimum_ratio`.

If the parent was processed before the daughter, the daughter receives the
parent's own `parent_id` (so chains of errors end at the ultimate parent);
otherwise (a parent with the same spread but lower in the processing order) it
receives the parent's id. Rows with the same `parent_id` are then summed into
the curated table.

# Arguments
- `otutable`: OTU table (see [`OTUTable`](@ref), [`read_otutable`](@ref)),
  dense or sparse.
- `matchlist`: pairwise similarities (see [`MatchList`](@ref),
  [`read_matchlist`](@ref), [`build_matchlist`](@ref)).
- `minimum_ratio_type`: `"min"` (default) or `"avg"` (also as `Symbol`).
- `minimum_ratio`: minimum parent/daughter abundance ratio (default 1).
- `minimum_match`: minimum similarity in percent (default 84).
- `minimum_relative_cooccurence`: minimum relative co-occurrence
  (default 0.95).
- `logfile`: `nothing` (default, no log); `true` for a file named
  `lulu.log_YYYYmmdd_HHMMSS` in the working directory, as the R package does;
  a file path (appended to); or an `IO`. The log text is identical to the
  one written by R.
- `progress`: print `progress: N%` for each OTU, as R does (default `false`).
- `storage`: representation of the reads during the parent search: `:auto`
  (default; sparse if the table is sparse or less than a quarter of its counts
  are non-zero), `:dense` or `:sparse`. It affects speed and memory only.

Without a log, the parent search runs on `Threads.nthreads()` threads; the
result does not depend on the number of threads. With a log it is sequential,
so that the log follows the processing order. The curated table has the same
storage (dense or sparse) as `otutable`.

# Differences from R
- No log file and no progress messages unless requested.
- An OTU with no reads that has potential parents is kept as its own parent
  (R stops with "missing value where TRUE/FALSE needed").
- `minimum_ratio_type` other than `"min"`/`"avg"`, or a `NaN` threshold, is an
  error (R silently uses `"min"`).
- Ids are sorted by Unicode code point (`curated_otus`, rows of
  `curated_table`), as dplyr does; `table()` in R uses the locale collation,
  which differs only for ids that mix upper and lower case or use
  punctuation.

See [`LuluResult`](@ref) for the returned components.
"""
function lulu(otutable::OTUTable, matchlist::MatchList;
              minimum_ratio_type::Union{AbstractString,Symbol} = "min",
              minimum_ratio::Real = 1, minimum_match::Real = 84,
              minimum_relative_cooccurence::Real = 0.95,
              logfile::Union{Nothing,Bool,AbstractString,IO} = nothing,
              progress::Bool = false, storage::Symbol = :auto)
    t0 = time()
    start = Dates.now()
    storage in (:auto, :dense, :sparse) ||
        throw(ArgumentError("storage must be :auto, :dense or :sparse, got :$storage"))
    criteria = _criteria(minimum_ratio_type, minimum_ratio, minimum_match,
                         minimum_relative_cooccurence)
    _warn_if_fraction(criteria.minimum_match, matchlist)

    total0, spread0 = _otu_statistics(otutable.counts)
    order = _processing_order(spread0, total0)
    ids = otutable.ids[order]
    total = total0[order]
    spread = spread0[order]
    table = otutable.counts[order, :]        # OTUs × samples, processing order
    reads = _kernel_reads(table, storage)    # samples × OTUs

    position = Dict(id => i for (i, id) in enumerate(ids))
    hit_positions = _hit_positions(matchlist, position, criteria.minimum_match)
    statistic = _ratio_statistic(criteria)
    progress && _print_progress(length(ids))
    parents = _with_log(logfile, start) do logio
        _assign_parents(logio, reads, ids, spread, matchlist, position, hit_positions,
                        criteria, statistic)
    end
    _lulu_result(table, ids, otutable.samples, total, spread, parents, criteria, time() - t0)
end

# The progress lines printed by R, one per OTU.
function _print_progress(n::Int)
    for line in 1:n
        println("progress: ", round(Int, line / n * 100), "%")
    end
end

"""
    lulu(counts::AbstractMatrix, ids, matchlist; samples, kwargs...)

Run [`lulu`](@ref) on a count matrix (OTUs × samples) with OTU ids `ids`.
"""
lulu(counts::AbstractMatrix, ids::AbstractVector, matchlist::MatchList;
     samples = ["S$j" for j in axes(counts, 2)], kwargs...) =
    lulu(OTUTable(ids, samples, counts), matchlist; kwargs...)

"""
    lulu(otutable_path::AbstractString, matchlist_path::AbstractString; kwargs...)

Read the OTU table and the match list from files (see [`read_otutable`](@ref)
and [`read_matchlist`](@ref)) and run [`lulu`](@ref).
"""
lulu(otutable_path::AbstractString, matchlist_path::AbstractString; kwargs...) =
    lulu(read_otutable(otutable_path), read_matchlist(matchlist_path); kwargs...)
