using Documenter
using MetacellsGraphs
using Random
using Test

Random.seed!(123456)

@testset "doctests" begin
    DocMeta.setdocmeta!(MetacellsGraphs, :DocTestSetup, :(using MetacellsGraphs); recursive = true)
    return doctest(MetacellsGraphs; manual = false)
end

@testset "hello" begin
    @test hello_metacells_graphs() == "Hello from MetacellsGraphs!"
end
