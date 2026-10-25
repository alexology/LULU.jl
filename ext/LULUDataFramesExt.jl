# DataFrames.jl integration, loaded automatically when DataFrames is loaded.

module LULUDataFramesExt

using LULU
using DataFrames: DataFrames, AbstractDataFrame, DataFrame, ncol, nrow

"""
    OTUTable(df::AbstractDataFrame; id_col = 1)

OTU table from a data frame with one row per OTU: the OTU ids in column
`id_col` (index or name) and the counts in all the other columns.
"""
function LULU.OTUTable(df::AbstractDataFrame; id_col = 1)
    idname = id_col isa Integer ? names(df)[id_col] : string(id_col)
    idname in names(df) || throw(ArgumentError("column $id_col not found"))
    samples = [n for n in names(df) if n != idname]
    for s in samples
        eltype(df[!, s]) <: Real ||
            throw(ArgumentError("sample column \"$s\" is not numeric (eltype $(eltype(df[!, s])))"))
    end
    counts = isempty(samples) ? zeros(Int, nrow(df), 0) :
             reduce(hcat, [collect(df[!, s]) for s in samples])
    LULU.OTUTable(string.(df[!, idname]), samples, counts)
end

"""
    MatchList(df::AbstractDataFrame)

Match list from the first three columns of a data frame (query id, hit id,
percent similarity), whatever their names.
"""
function LULU.MatchList(df::AbstractDataFrame)
    ncol(df) >= 3 || throw(ArgumentError("a match list needs at least 3 columns"))
    LULU.MatchList(string.(df[!, 1]), string.(df[!, 2]), Float64.(df[!, 3]))
end

"""
    DataFrame(t::OTUTable; id_col = :OTUid)

Data frame with a column of OTU ids followed by one column per sample.
Sample names must be unique and differ from `id_col`; choose another
`id_col` when a sample is named `OTUid`.
"""
function DataFrames.DataFrame(t::LULU.OTUTable; id_col = :OTUid)
    length(unique(t.samples)) == length(t.samples) ||
        throw(ArgumentError("sample names must be unique to convert to a DataFrame"))
    string(id_col) in t.samples &&
        throw(ArgumentError("id_col '$id_col' is also a sample name; choose another id_col"))
    df = DataFrame(string(id_col) => copy(t.ids))
    for (j, s) in enumerate(t.samples)
        df[!, s] = collect(t.counts[:, j])   # dense column also for sparse tables
    end
    df
end

"""
    DataFrame(m::OTUMap)

Data frame with columns `OTUid`, `total`, `spread`, `parent_id`, `curated`,
`rank` (the R `otu_map`).
"""
DataFrames.DataFrame(m::LULU.OTUMap) =
    DataFrame(OTUid = m.ids, total = m.total, spread = m.spread, parent_id = m.parent_id,
              curated = m.curated, rank = m.rank)

"""
    DataFrame(m::MatchList)

Data frame with columns `query`, `hit`, `match`.
"""
DataFrames.DataFrame(m::LULU.MatchList) = DataFrame(query = m.query, hit = m.hit, match = m.match)

LULU.lulu(otutable::AbstractDataFrame, matchlist::Union{LULU.MatchList,AbstractDataFrame};
          id_col = 1, kwargs...) =
    LULU.lulu(LULU.OTUTable(otutable; id_col), _matchlist(matchlist); kwargs...)

LULU.lulu(otutable::LULU.OTUTable, matchlist::AbstractDataFrame; kwargs...) =
    LULU.lulu(otutable, LULU.MatchList(matchlist); kwargs...)

_matchlist(m::LULU.MatchList) = m
_matchlist(m::AbstractDataFrame) = LULU.MatchList(m)

end # module
