using DataAxesFormats
using Documenter
using MetacellsGraphs
using NestedTests
using Random
using SomeGraphs
using Test

test_prefixes(ARGS)
abort_on_first_failure(true)

Random.seed!(123456)

nested_test("doctests") do
    DocMeta.setdocmeta!(MetacellsGraphs, :DocTestSetup, :(using MetacellsGraphs); recursive = true)
    return doctest(MetacellsGraphs; manual = false)
end

include("scatter_graphs.jl")
include("heatmap_graphs.jl")
