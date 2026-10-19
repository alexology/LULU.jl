"""
    LULU

Pure Julia implementation of LULU, the post-clustering curation algorithm for
amplicon data by Frøslev et al. (2017, *Nature Communications* 8:1188),
ported from the R package `lulu` 0.1.0 (<https://github.com/tobiasgf/lulu>).

LULU merges "daughter" OTUs (presumed artefacts) into sequence-similar, more
widespread "parent" OTUs with which they consistently co-occur.

Main entry points:

- [`lulu`](@ref): the curation algorithm (port of the R function `lulu()`);
- [`build_matchlist`](@ref): pure Julia replacement for the external
  BLASTN/VSEARCH step that produces the match list;
- [`read_otutable`](@ref), [`read_matchlist`](@ref), [`read_fasta`](@ref),
  [`write_otutable`](@ref), [`write_matchlist`](@ref), [`write_otu_map`](@ref).
"""
module LULU

using Dates: Dates
using Printf: Printf, @sprintf
using SparseArrays: SparseArrays, SparseMatrixCSC, issparse, nnz, nonzeros, nzrange,
                    rowvals, sparse

export OTUTable, MatchList, OTUMap, LuluResult
export lulu, build_matchlist
export read_otutable, write_otutable, read_matchlist, write_matchlist,
       read_fasta, write_otu_map
export notus, nsamples

include("types.jl")
include("rnumeric.jl")
include("io.jl")
include("curation.jl")
include("alignment.jl")
include("matchlist.jl")

end # module
