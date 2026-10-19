using LULU
using Test
using DataFrames
using SparseArrays

const DATA = joinpath(@__DIR__, "data")
const REF = joinpath(@__DIR__, "reference")

include("test_parity.jl")
include("test_lulu.jl")
include("test_io.jl")
include("test_matchlist.jl")
include("test_dataframes.jl")
