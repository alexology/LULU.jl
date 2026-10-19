# Pure Julia construction of the LULU match list (all OTU sequences searched
# against each other), replacing the external BLASTN/VSEARCH step. The search
# itself (k-mer prefilter and banded alignment) is `search_self` of
# NucleotideSearch.jl, with the strategy of `vsearch --usearch_global --self`.

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
close to, not identical with, those of VSEARCH or BLASTN. The search is
`search_self` of NucleotideSearch.jl.

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
    3 <= wordlength <= 12 || throw(ArgumentError("wordlength must be between 3 and 12"))
    # The other arguments are checked by NucleotideSearch.
    scoring = NucleotideSearch.Scoring(; match, mismatch, gap_open, gap_extend,
                                       terminal_gap_open, terminal_gap_extend)
    db = NucleotideSearch.SequenceDB(ids, sequences; k = wordlength)
    results = NucleotideSearch.search_self(db; min_identity, iddef, query_cov, maxhits,
                                           maxaccepts, maxrejects, minwordmatches, scoring)
    query = String[]; hit = String[]; ident = Float64[]
    for hits in results, h in hits
        push!(query, ids[h.query]); push!(hit, ids[h.target]); push!(ident, h.identity)
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
