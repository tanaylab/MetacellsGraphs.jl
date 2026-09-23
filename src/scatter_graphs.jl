"""
Scatter graphs of a metacells repository.

A graph function fills in the data the graph is about, and nothing else. In particular it does not color the points.
The caller says how to color them in a second call, so that any property can be the color:

    graph = umap_graph(daf)
    fill_type!(points_colors_vector_fields(graph), daf)

The second call must name the same `axis` and `entries` as the first, or the vectors will not match. The same goes for
the sizes of the points, and for anything else the graph leaves empty.
"""
module ScatterGraphs

export gene_base_delta_correlations_graph
export gene_gene_graph
export umap_graph

using DataAxesFormats
using SomeGraphs
using TanayLabUtilities

using ..DataSources

"""
    gene_gene_graph(
        daf::DafReader;
        axis::AbstractString = "metacell",
        x_gene::AbstractString,
        y_gene::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = $(GENE_FRACTION_REGULARIZATION_FOR_GRAPHS),
    )::PointsGraph

The expression of `x_gene` against `y_gene`, a point per each of the `entries` of the `axis`, on log scale.

The `axis` can be any axis which has a `linear_fraction` per gene, typically `metacell` or `block`. A block is coarser,
so the same pair of genes shows fewer points and less scatter.

The hover of each point names its entry and gives both fractions.
"""
function gene_gene_graph(
    daf::DafReader;
    axis::AbstractString = "metacell",
    x_gene::AbstractString,
    y_gene::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
)::PointsGraph
    graph = points_graph()
    fill_gene_expression!(x_axis_vector_fields(graph), daf; gene = x_gene, axis, entries, gene_fraction_regularization)  # NOJET
    fill_gene_expression!(y_axis_vector_fields(graph), daf; gene = y_gene, axis, entries, gene_fraction_regularization)
    return graph
end

"""
    umap_graph(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::PointsGraph

The 2D UMAP embedding, a point per each of the `entries` of the `axis`.

The `axis` can be any axis which has UMAP coordinates, typically `metacell` or `block`. A block is placed at the mean
position of its metacells, so the same embedding shows fewer points closer to its center.

The coordinates are the arbitrary output of the UMAP projection, so the axes are drawn without ticks or a grid. Only
which points are near which other points means anything.

The hover of each point names its entry and gives its coordinates.
"""
function umap_graph(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::PointsGraph
    graph = points_graph()
    fill_umap!(x_axis_vector_fields(graph), daf; coordinate = "x", axis, entries)
    fill_umap!(y_axis_vector_fields(graph), daf; coordinate = "y", axis, entries)
    return graph
end

"""
    gene_base_delta_correlations_graph(;
        daf::DafReader,
        base_daf::DafReader,
        gene::AbstractString,
        axis::AbstractString = "base_block",
        base_axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = $(GENE_FRACTION_REGULARIZATION_FOR_GRAPHS),
    )::PointsGraph

What the metacells did to one gene, a point per each of the `entries`: how much its correlation with the cells changed,
against how much of the gene there is to say anything about.

Both repositories hold the correlations along the `axis`. The `base_daf` holds the expression along its `base_axis`
instead, because what is a base block there is also a block. The two must list the same names in the same order.

The `x` axis is the gene's correlation between the cells of the entry's neighborhood and their punctuated metacells,
minus what it is in the `base_daf`. A point right of zero is an entry whose metacells describe the gene better than the
base does. The `y` axis is how much of the gene the entry expresses, on log scale, since moving the correlation of a
gene which is barely there says less than moving one which is everywhere.

Only the entries whose base correlation is not zero are shown. A zero there is what a gene saying nothing about an
entry looks like, rather than a correlation which happens to be zero. A gene which says nothing about any entry has
nothing to show and is an error rather than an empty graph. The hidden entries stay in the data, masked, and do not
affect the ranges of the axes.

The hover of each point names its entry and gives both correlations, the change between them, and the fraction of the
gene.
"""
function gene_base_delta_correlations_graph(;
    daf::DafReader,
    base_daf::DafReader,
    gene::AbstractString,
    axis::AbstractString = "base_block",
    base_axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
)::PointsGraph
    @assert axis_vector(base_daf, "gene") == axis_vector(daf, "gene")
    @assert axis_vector(base_daf, base_axis) == axis_vector(daf, axis)

    graph = points_graph()
    graph.configuration.x_axis.scale.include_hidden = false
    graph.configuration.y_axis.scale.include_hidden = false

    entities = points_entities(graph)
    fill_gene_correlation!(entities, base_daf; gene, axis, entries, title = "base correlation")
    fill_gene_correlation!(entities, daf; gene, axis, entries)
    fill_gene_correlation_change!(x_axis_vector_fields(graph), daf, base_daf; gene, axis, entries)
    fill_gene_expression!(
        y_axis_vector_fields(graph),
        base_daf;
        gene,
        axis = base_axis,
        entries,
        gene_fraction_regularization,
    )

    # A zero base correlation is the gene saying nothing about the entry, so such entries are hidden.
    is_shown_per_entry = get_gene_correlation_vector(base_daf; gene, axis, entries) .!= 0
    if !any(is_shown_per_entry)
        error("no base block correlates the gene: $(gene)\nof the base daf data: $(base_daf.name)")
    end
    put_vector_mask_data!(entities, is_shown_per_entry)

    return graph
end

end  # module
