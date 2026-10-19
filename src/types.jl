# Data structures: OTU table, match list, OTU map and the result of `lulu`.

"""
    OTUTable(ids, samples, counts)

An OTU table: OTUs as rows, samples as columns, read counts as values
(the R `otutable` data frame with OTU ids as row names).

- `ids`: unique OTU ids (row names);
- `samples`: sample names (column names);
- `counts`: `length(ids) × length(samples)` matrix of non-negative real numbers
  (integers or floats). A `SparseMatrixCSC` is kept sparse, which saves memory
  and time for large tables where most counts are zero; any other matrix is
  converted to a dense `Matrix`.

Fields: `ids::Vector{String}`, `samples::Vector{String}`, `counts::M` (a
`Matrix{T}` or a `SparseMatrixCSC{T}`). `sparse(table)` and `Matrix(table)`
convert between the two storages.
"""
struct OTUTable{T<:Real,M<:AbstractMatrix{T}}
    ids::Vector{String}
    samples::Vector{String}
    counts::M
    function OTUTable{T,M}(ids::Vector{String}, samples::Vector{String},
                           counts::M) where {T<:Real,M<:AbstractMatrix{T}}
        size(counts) == (length(ids), length(samples)) ||
            throw(DimensionMismatch("counts has size $(size(counts)), expected " *
                                    "($(length(ids)), $(length(samples))) from ids and samples"))
        if !allunique(ids)
            seen = Set{String}()
            dup = first(x for x in ids if !(x in seen || (push!(seen, x); false)))
            throw(ArgumentError("OTU ids must be unique; duplicated id: \"$dup\""))
        end
        for x in _stored_values(counts)
            (isnan(x) || x < 0) &&
                throw(ArgumentError("counts must be non-negative numbers; found $x"))
        end
        new{T,M}(ids, samples, counts)
    end
end

# Values that need checking: all entries, or only the stored ones.
_stored_values(c::AbstractMatrix) = c
_stored_values(c::SparseMatrixCSC) = nonzeros(c)

_count_storage(c::SparseMatrixCSC) = c
_count_storage(c::AbstractMatrix) = Matrix(c)

function OTUTable(ids::AbstractVector, samples::AbstractVector,
                  counts::AbstractMatrix{T}) where {T<:Real}
    c = _count_storage(counts)
    OTUTable{T,typeof(c)}(String[string(x) for x in ids], String[string(x) for x in samples], c)
end

"""
    OTUTable(counts::AbstractMatrix; ids, samples)

Build a table from a matrix; default ids are `"OTU1"`, `"OTU2"`, … and default
sample names `"S1"`, `"S2"`, ….
"""
function OTUTable(counts::AbstractMatrix;
                  ids = ["OTU$i" for i in axes(counts, 1)],
                  samples = ["S$j" for j in axes(counts, 2)])
    OTUTable(ids, samples, counts)
end

SparseArrays.sparse(t::OTUTable) = OTUTable(t.ids, t.samples, sparse(t.counts))
Base.Matrix(t::OTUTable) = OTUTable(t.ids, t.samples, Matrix(t.counts))
SparseArrays.issparse(t::OTUTable) = issparse(t.counts)

"""
    notus(t::OTUTable)

Number of OTUs (rows).
"""
notus(t::OTUTable) = length(t.ids)

"""
    nsamples(t::OTUTable)

Number of samples (columns).
"""
nsamples(t::OTUTable) = length(t.samples)

Base.size(t::OTUTable) = size(t.counts)
Base.:(==)(a::OTUTable, b::OTUTable) =
    a.ids == b.ids && a.samples == b.samples && a.counts == b.counts

function Base.show(io::IO, t::OTUTable{T}) where {T}
    kind = issparse(t.counts) ? " (sparse, $(nnz(t.counts)) stored counts)" : ""
    print(io, "OTUTable{$T} with $(notus(t)) OTUs and $(nsamples(t)) samples", kind)
end

"""
    MatchList(query, hit, match)

Pairwise sequence similarities between OTUs, as produced by BLASTN, VSEARCH or
[`build_matchlist`](@ref): one entry per (query OTU, hit OTU) pair, with the
similarity `match` in percent (e.g. `96.5`).

In `lulu` the query is the potential daughter and the hit a potential parent.
Hits equal to `"*"` (VSEARCH "no hit") and self-hits are ignored.

Fields: `query::Vector{String}`, `hit::Vector{String}`, `match::Vector{Float64}`.
"""
struct MatchList
    query::Vector{String}
    hit::Vector{String}
    match::Vector{Float64}
    function MatchList(query::Vector{String}, hit::Vector{String}, match::Vector{Float64})
        length(query) == length(hit) == length(match) ||
            throw(DimensionMismatch("query, hit and match must have the same length"))
        new(query, hit, match)
    end
end

MatchList(query::AbstractVector, hit::AbstractVector, match::AbstractVector{<:Real}) =
    MatchList(String[string(x) for x in query], String[string(x) for x in hit],
              Float64[Float64(x) for x in match])

"""
    MatchList(entries)

Build a match list from an iterable of `(query, hit, match)` tuples.
"""
function MatchList(entries)
    q = String[]; h = String[]; m = Float64[]
    for (a, b, c) in entries
        push!(q, string(a)); push!(h, string(b)); push!(m, Float64(c))
    end
    MatchList(q, h, m)
end

Base.length(m::MatchList) = length(m.query)
Base.:(==)(a::MatchList, b::MatchList) =
    a.query == b.query && a.hit == b.hit && a.match == b.match

function Base.show(io::IO, m::MatchList)
    print(io, "MatchList with $(length(m)) entries")
end

"""
    OTUMap

Per-OTU curation summary (the R `otu_map` data frame), with one entry per OTU
in processing order (decreasing spread, then decreasing total read count).

Fields (all vectors of the same length):

- `ids`: OTU id;
- `total`: total read count;
- `spread`: number of samples where the OTU is present (count > 0);
- `parent_id`: id of the OTU it was merged into, or its own id;
- `curated`: `"parent"` (kept) or `"merged"`;
- `rank`: rank by decreasing total read count (ties: processing order).
"""
struct OTUMap
    ids::Vector{String}
    total::Vector{Float64}
    spread::Vector{Int}
    parent_id::Vector{String}
    curated::Vector{String}
    rank::Vector{Int}
end

Base.length(m::OTUMap) = length(m.ids)

function Base.show(io::IO, m::OTUMap)
    print(io, "OTUMap with $(length(m)) OTUs ($(count(==("parent"), m.curated)) parents)")
end

"""
    LuluResult

Result of [`lulu`](@ref), with the same components as the list returned by the
R function:

- `curated_table::OTUTable`: daughters summed into their parents; rows sorted
  by id;
- `curated_count::Int`, `curated_otus::Vector{String}` (sorted): retained OTUs;
- `discarded_count::Int`, `discarded_otus::Vector{String}` (processing order):
  OTUs merged into a parent;
- `runtime::Float64`: elapsed time in seconds;
- `minimum_match`, `minimum_relative_cooccurence`: settings used;
- `otu_map::OTUMap`: per-OTU summary;
- `original_table::OTUTable`: input table, rows in processing order.

The other settings are stored in `minimum_ratio_type` (`"min"` or `"avg"`)
and `minimum_ratio`.
"""
struct LuluResult{T<:Real,M<:AbstractMatrix{T}}
    curated_table::OTUTable{T,M}
    curated_count::Int
    curated_otus::Vector{String}
    discarded_count::Int
    discarded_otus::Vector{String}
    runtime::Float64
    minimum_match::Float64
    minimum_relative_cooccurence::Float64
    minimum_ratio_type::String
    minimum_ratio::Float64
    otu_map::OTUMap
    original_table::OTUTable{T,M}
end

function Base.show(io::IO, r::LuluResult)
    print(io, "LuluResult: $(r.curated_count) OTUs retained, $(r.discarded_count) merged ",
          "(minimum_match = $(r.minimum_match), minimum_relative_cooccurence = ",
          "$(r.minimum_relative_cooccurence), minimum_ratio_type = \"$(r.minimum_ratio_type)\", ",
          "minimum_ratio = $(r.minimum_ratio))")
end
