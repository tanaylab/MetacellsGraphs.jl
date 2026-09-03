"""
Draw the figures worth looking at of a metacells repository.

This sits between [`SomeGraphs`](https://github.com/tanaylab/SomeGraphs.jl), which draws heatmaps, scatters, bars and
distributions from arrays and knows nothing of metacells, and [`Metacells`](https://github.com/tanaylab/Metacells.jl),
which computes what a metacells repository holds and draws nothing. A figure here reads a repository through
`DataAxesFormats`, names its properties as `Metacells` does, and draws with `SomeGraphs`.

It does the arithmetic a figure needs and keeps the result in the figure: means, quantiles, normalizations, a fold
factor from a median, the order a clustered heatmap is drawn in. It also gathers across repositories, which is how a
figure compares the rounds of sharpening. What it does not do is produce data the repository should have held - if a
number is worth storing, it belongs in a computation in `Metacells`, and nothing here writes to a repository.

The figures are meant to be asked for from a notebook: a call takes a repository and a few names and returns a figure,
writes no files, and needs no output directory.
"""
module MetacellsGraphs

export hello_metacells_graphs

using DataAxesFormats
using Metacells
using SomeGraphs
using TanayLabUtilities

"""
    hello_metacells_graphs()::AbstractString

Say hello.

A placeholder, so that the package has something in it: the gates have nothing to look at otherwise. Delete it when the
first real figure arrives.

```jldoctest
println(hello_metacells_graphs())

# output

Hello from MetacellsGraphs!
```
"""
function hello_metacells_graphs()::AbstractString
    return "Hello from MetacellsGraphs!"
end

end  # module
