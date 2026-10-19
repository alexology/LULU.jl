# Reading and writing OTU tables, match lists, FASTA files and OTU maps.

_unquote(s::AbstractString) =
    (length(s) >= 2 && s[1] == '"' && s[end] == '"') ? String(s[2:end-1]) : String(s)

# Counts parsed line by line: integers while every value is an integer, then
# floating point. With `keep_zeros = false` only non-zero values are stored
# (with their row and column), for a sparse matrix.
mutable struct _CountParser
    keep_zeros::Bool
    all_int::Bool
    ints::Vector{Int}
    floats::Vector{Float64}
    rows::Vector{Int}
    cols::Vector{Int}
end

_CountParser(keep_zeros::Bool) = _CountParser(keep_zeros, true, Int[], Float64[], Int[], Int[])

function _push_count!(p::_CountParser, field::AbstractString, row::Int, col::Int, path)
    if p.all_int
        v = tryparse(Int, field)
        if v !== nothing
            (p.keep_zeros || v != 0) || return p
            push!(p.ints, v)
            p.keep_zeros || (push!(p.rows, row); push!(p.cols, col))
            return p
        end
        p.floats = Float64.(p.ints)
        p.all_int = false
    end
    x = tryparse(Float64, field)
    x === nothing && throw(ArgumentError(
        "$path: cannot parse \"$field\" as a number (row $(row + 1), column $(col + 1))"))
    (p.keep_zeros || x != 0) || return p
    push!(p.floats, x)
    p.keep_zeros || (push!(p.rows, row); push!(p.cols, col))
    p
end

_values(p::_CountParser) = p.all_int ? p.ints : p.floats

function _count_matrix(p::_CountParser, nrow::Int, ncol::Int)
    v = _values(p)
    p.keep_zeros ? permutedims(reshape(v, ncol, nrow)) : sparse(p.rows, p.cols, v, nrow, ncol)
end

"""
    read_otutable(path; delim = '\t', storage = :dense) -> OTUTable

Read an OTU table: a header line with sample names, then one line per OTU with
the OTU id in the first column and the counts in the following ones (the
format read in the LULU tutorial with
`read.csv(path, sep = "\t", header = TRUE, row.names = 1)`).

The header may or may not contain a name for the id column (e.g. an empty first
field or `"OTUid"`). Counts are read as `Int` when all values are integers,
otherwise as `Float64`. Double quotes around fields are removed; empty lines
are skipped. Sample names are kept as written (R would apply `make.names`).

`storage = :sparse` stores the counts in a `SparseMatrixCSC` (only non-zero
values are kept while reading): recommended for large tables in which most
OTUs occur in few samples.
"""
function read_otutable(path::AbstractString; delim::Union{Char,AbstractString} = '\t',
                       storage::Symbol = :dense)
    storage in (:dense, :sparse) ||
        throw(ArgumentError("storage must be :dense or :sparse, got :$storage"))
    lines = filter(!isempty, map(l -> rstrip(l, '\r'), readlines(path)))
    isempty(lines) && throw(ArgumentError("$path: empty file"))
    header = [_unquote(x) for x in split(lines[1], delim)]
    ids = String[]
    parser = _CountParser(storage === :dense)
    ncol = -1
    for (row, line) in enumerate(@view lines[2:end])
        fields = split(line, delim)
        ncol < 0 && (ncol = length(fields) - 1)
        length(fields) == ncol + 1 ||
            throw(ArgumentError("$path: line $(row + 1) has $(length(fields)) fields, " *
                                "expected $(ncol + 1)"))
        push!(ids, _unquote(fields[1]))
        for col in 1:ncol
            _push_count!(parser, _unquote(fields[col+1]), row, col, path)
        end
    end
    ncol < 0 && (ncol = length(header) - 1)
    samples = if length(header) == ncol + 1
        header[2:end]
    elseif length(header) == ncol
        header
    else
        throw(ArgumentError("$path: the header has $(length(header)) fields, " *
                            "but data lines have $(ncol + 1)"))
    end
    OTUTable(ids, samples, _count_matrix(parser, length(ids), ncol))
end

"""
    write_otutable(path, table::OTUTable; delim = '\\t', id_header = "OTUid")

Write an OTU table in the format read by [`read_otutable`](@ref).
"""
function write_otutable(path::AbstractString, t::OTUTable; delim = '\t',
                        id_header::AbstractString = "OTUid")
    by_otu = _otu_columns(t.counts)   # samples × OTUs: one column per output line
    open(path, "w") do io
        println(io, join([id_header; t.samples], delim))
        for i in 1:notus(t)
            print(io, t.ids[i])
            _write_counts(io, by_otu, i, delim)
            println(io)
        end
    end
    path
end

# Counts of samples × OTUs matrices, one column per OTU.
_otu_columns(c::AbstractMatrix) = permutedims(c)
_otu_columns(c::SparseMatrixCSC) = sparse(transpose(c))

_write_counts(io::IO, c::AbstractMatrix, i::Int, delim) =
    foreach(x -> print(io, delim, _format_count(x)), view(c, :, i))

function _write_counts(io::IO, c::SparseMatrixCSC{T}, i::Int, delim) where {T}
    rows = rowvals(c); vals = nonzeros(c)
    next = 1                                  # next sample to write
    for k in nzrange(c, i)
        for _ in next:rows[k]-1
            print(io, delim, _format_count(zero(T)))
        end
        print(io, delim, _format_count(vals[k]))
        next = rows[k] + 1
    end
    for _ in next:size(c, 1)
        print(io, delim, _format_count(zero(T)))
    end
end

_format_count(x::Integer) = string(x)
_format_count(x::Real) = isinteger(x) && abs(x) < 1e15 ? string(Int(x)) : string(x)

"""
    read_matchlist(path) -> MatchList

Read a match list: whitespace-separated lines with query OTU id, hit OTU id and
percent similarity, e.g. BLASTN `-outfmt '6 qseqid sseqid pident'` or VSEARCH
`--userfields query+target+id`. Like R's `read.table(path, header = FALSE)`:
any whitespace separates fields, text after `#` is ignored, columns after the
third are ignored.
"""
function read_matchlist(path::AbstractString)
    q = String[]; h = String[]; m = Float64[]
    for (k, raw) in enumerate(readlines(path))
        line = raw
        c = findfirst('#', line)
        c === nothing || (line = line[1:prevind(line, c)])
        f = split(line)
        isempty(f) && continue
        length(f) >= 3 ||
            throw(ArgumentError("$path: line $k has $(length(f)) fields, expected at least 3"))
        v = tryparse(Float64, f[3])
        v === nothing && throw(ArgumentError("$path: line $k: cannot parse \"$(f[3])\" as a number"))
        push!(q, _unquote(f[1])); push!(h, _unquote(f[2])); push!(m, v)
    end
    MatchList(q, h, m)
end

"""
    write_matchlist(path, matchlist::MatchList; digits = 3)

Write a tab-separated match list without header (query, hit, percent
similarity rounded to `digits` decimals, as in BLASTN output). Use
`digits = nothing` to write full precision.
"""
function write_matchlist(path::AbstractString, m::MatchList; digits::Union{Int,Nothing} = 3)
    open(path, "w") do io
        for k in 1:length(m)
            v = digits === nothing ? string(m.match[k]) : _printf_fixed(digits, m.match[k])
            println(io, m.query[k], '\t', m.hit[k], '\t', v)
        end
    end
    path
end

"""
    read_fasta(path) -> (ids::Vector{String}, sequences::Vector{String})

Read a FASTA file. The id is the first whitespace-delimited word of the header
line; sequences may span several lines and are converted to upper case.
"""
function read_fasta(path::AbstractString)
    ids = String[]; seqs = String[]
    buf = IOBuffer()
    for raw in readlines(path)   # readlines closes the file even if parsing fails
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, '>')
            isempty(ids) || push!(seqs, uppercase(String(take!(buf))))
            w = split(line[2:end])
            push!(ids, isempty(w) ? "" : String(w[1]))
        else
            isempty(ids) && throw(ArgumentError("$path: sequence data before the first header"))
            print(buf, line)
        end
    end
    isempty(ids) || push!(seqs, uppercase(String(take!(buf))))
    ids, seqs
end

"""
    write_otu_map(path, otu_map::OTUMap; delim = '\\t')

Write the OTU map as a table with columns `OTUid`, `total`, `spread`,
`parent_id`, `curated`, `rank`.
"""
function write_otu_map(path::AbstractString, m::OTUMap; delim = '\t')
    open(path, "w") do io
        println(io, join(("OTUid", "total", "spread", "parent_id", "curated", "rank"), delim))
        for i in 1:length(m)
            println(io, join((m.ids[i], _format_count(m.total[i]), m.spread[i],
                              m.parent_id[i], m.curated[i], m.rank[i]), delim))
        end
    end
    path
end
