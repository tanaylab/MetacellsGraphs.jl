"""
Scatter graphs of a metacells repository.
"""
module ScatterGraphs

export blocks_gene_gene_graph
export blocks_umap_graph
export gene_base_delta_correlations_graph
export metacells_gene_gene_graph
export metacells_umap_graph

using DataAxesFormats
using Metacells
using SomeGraphs
using TanayLabUtilities

using ..Utilities

# Needed because of JET:
import Metacells.Contracts.base_block_axis
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.type_axis
import Metacells.Contracts.vector_of_color_per_type
import Metacells.Contracts.vector_of_type_per_block
import Metacells.Contracts.vector_of_type_per_metacell
import Metacells.Contracts.vector_of_umap_x_per_block
import Metacells.Contracts.vector_of_umap_x_per_metacell
import Metacells.Contracts.vector_of_umap_y_per_block
import Metacells.Contracts.vector_of_umap_y_per_metacell

"""
The regularization added to a fraction before taking its log, so that a gene seen in no cell of an entity has a value
to be drawn at instead of falling off the bottom of the axis. It is deliberately small: a fraction this size is one UMI
in a hundred thousand, which is below anything these graphs are read for.
"""
GENE_FRACTION_REGULARIZATION_FOR_GRAPHS::Float64 = 1e-5

"""
    metacells_gene_gene_graph(
        daf::DafReader;
        x_gene::AbstractString,
        y_gene::AbstractString,
        gene_fraction_regularization::Real = $(DEFAULT.gene_fraction_regularization),
    )::PointsGraph

The expression of one gene against another, a point per metacell, on log scale.

If the metacells have a type, the points are colored by it, using the colors of the type axis, and the types are shown
in a legend. If they do not, they are all the same color; nothing else about the graph changes.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_type_per_metacell(OptionalInput),
        vector_of_color_per_type(OptionalInput),
    ],
) function metacells_gene_gene_graph(
    daf::DafReader;
    x_gene::AbstractString,
    y_gene::AbstractString,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
)::PointsGraph
    return gene_gene_graph(daf; axis = "metacell", x_gene, y_gene, gene_fraction_regularization)
end

"""
    blocks_gene_gene_graph(
        daf::DafReader;
        x_gene::AbstractString,
        y_gene::AbstractString,
        gene_fraction_regularization::Real = $(DEFAULT.gene_fraction_regularization),
    )::PointsGraph

The expression of one gene against another, a point per block, on log scale.

This is [`metacells_gene_gene_graph`](@ref) of the blocks the metacells were grouped into, and reads the same way; a
block is coarser, so the same pair of genes shows fewer points and less scatter.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_type_per_block(OptionalInput),
        vector_of_color_per_type(OptionalInput),
    ],
) function blocks_gene_gene_graph(
    daf::DafReader;
    x_gene::AbstractString,
    y_gene::AbstractString,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
)::PointsGraph
    return gene_gene_graph(daf; axis = "block", x_gene, y_gene, gene_fraction_regularization)
end

# One gene against another, a point per entry of some axis. The two graphs differ only in which axis the points are,
# and every axis which has a `linear_fraction` per gene also has an optional `type`, so they are one function.
function gene_gene_graph(
    daf::DafReader;
    axis::AbstractString,
    x_gene::AbstractString,
    y_gene::AbstractString,
    gene_fraction_regularization::Real,
)::PointsGraph
    @assert gene_fraction_regularization >= 0

    points_xs = get_matrix(daf, "gene", axis, "linear_fraction")[x_gene, :].array
    points_ys = get_matrix(daf, "gene", axis, "linear_fraction")[y_gene, :].array
    points_colors, colors_configuration = points_type_colors(daf, axis)

    return points_graph(;
        x_axis_title = "$(x_gene) fraction",
        y_axis_title = "$(y_gene) fraction",
        points_colors_title = "type",
        points_xs,
        points_ys,
        points_colors,
        points_hovers = axis_vector(daf, axis),
        configuration = PointsGraphConfiguration(;
            x_axis = AxisConfiguration(; log_scale = Log2Scale, log_regularization = gene_fraction_regularization),
            y_axis = AxisConfiguration(; log_scale = Log2Scale, log_regularization = gene_fraction_regularization),
            points = ScattersConfiguration(; colors = colors_configuration),
        ),
    )
end

"""
    metacells_umap_graph(daf::DafReader)::PointsGraph

The 2D UMAP embedding of the metacells, a point per metacell.

If the metacells have a type, the points are colored by it, using the colors of the type axis, and the types are shown
in a legend. If they do not, they are all the same color; nothing else about the graph changes.

The coordinates are the arbitrary output of the UMAP projection, so the axes are drawn without ticks or a grid; only
which points are near which other points means anything.

$(CONTRACT)
"""
@computation Contract(;
    axes = [metacell_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        vector_of_umap_x_per_metacell(RequiredInput),
        vector_of_umap_y_per_metacell(RequiredInput),
        vector_of_type_per_metacell(OptionalInput),
        vector_of_color_per_type(OptionalInput),
    ],
) function metacells_umap_graph(daf::DafReader)::PointsGraph
    return umap_graph(daf; axis = "metacell")
end

"""
    blocks_umap_graph(daf::DafReader)::PointsGraph

The 2D UMAP embedding of the blocks, a point per block.

This is [`metacells_umap_graph`](@ref) of the blocks the metacells were grouped into, and reads the same way; a block
is placed at the mean position of its metacells, so the same embedding shows fewer points closer to its center.

$(CONTRACT)
"""
@computation Contract(;
    axes = [block_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        vector_of_umap_x_per_block(RequiredInput),
        vector_of_umap_y_per_block(RequiredInput),
        vector_of_type_per_block(OptionalInput),
        vector_of_color_per_type(OptionalInput),
    ],
) function blocks_umap_graph(daf::DafReader)::PointsGraph
    return umap_graph(daf; axis = "block")
end

# The UMAP embedding, a point per entry of some axis. As with `gene_gene_graph`, the two graphs differ only in which
# axis the points are.
function umap_graph(daf::DafReader; axis::AbstractString)::PointsGraph
    points_colors, colors_configuration = points_type_colors(daf, axis)

    return points_graph(;
        x_axis_title = "UMAP x",
        y_axis_title = "UMAP y",
        points_colors_title = "type",
        points_xs = get_vector(daf, axis, "umap_x").array,
        points_ys = get_vector(daf, axis, "umap_y").array,
        points_colors,
        points_hovers = axis_vector(daf, axis),
        configuration = PointsGraphConfiguration(;
            x_axis = AxisConfiguration(; show_ticks = false, show_grid = false),
            y_axis = AxisConfiguration(; show_ticks = false, show_grid = false),
            points = ScattersConfiguration(; colors = colors_configuration),
        ),
    )
end

"""
    gene_base_delta_correlations_graph(;
        daf::DafReader,
        base_daf::DafReader,
        gene::AbstractString,
        gene_fraction_regularization::Real = $(DEFAULT.gene_fraction_regularization),
    )::PointsGraph

What the metacells did to one gene, a point per base block: how much its correlation with the cells changed, against
how much of the gene there is in the block to say anything about.

The `x` axis is the gene's correlation between the cells of the base block's neighborhood and their punctuated
metacells, minus what it is in the `base_daf` - so a point right of zero is a block whose metacells describe the gene
better than the base does. The `y` axis is how much of the gene the base block expresses, on log scale, since moving
the correlation of a gene which is barely there says less than moving one which is everywhere.

Only the base blocks whose base correlation is not zero are shown. A zero there is what a gene saying nothing about a
block looks like, rather than a correlation which happens to be zero; a gene which says nothing about any block has
nothing to show and is an error rather than an empty graph.

If the blocks have a type, the points are colored by it, using the colors of the type axis. Each point's hover names
its block and type and gives both correlations along with the change between them.

# Daf

$(CONTRACT1)

# Base

$(CONTRACT2)
"""
@computation Contract(;
    name = "daf",
    link = Metacells,
    axes = [gene_axis(RequiredInput), base_block_axis(RequiredInput)],
    data = [
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
    ],
) Contract(;
    name = "base_daf",
    link = Metacells,
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_type_per_block(OptionalInput),
        vector_of_color_per_type(OptionalInput),
    ],
) function gene_base_delta_correlations_graph(;
    daf::DafReader,
    base_daf::DafReader,
    gene::AbstractString,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
)::PointsGraph
    @assert gene_fraction_regularization >= 0
    @assert axis_vector(base_daf, "gene") == axis_vector(daf, "gene")
    @assert axis_vector(base_daf, "block") == axis_vector(daf, "base_block")

    base_correlation_per_base_block = correlation_per_base_block(base_daf, gene)
    delta_correlation_per_base_block = correlation_per_base_block(daf, gene) .- base_correlation_per_base_block

    block_indices = findall(base_correlation_per_base_block .!= 0)
    if isempty(block_indices)
        error("no base block correlates the gene: $(gene)\nof the base daf data: $(base_daf.name)")
    end

    type_per_block, colors_configuration = points_type_colors(base_daf, "block")
    linear_fraction_per_block = get_matrix(base_daf, "gene", "block", "linear_fraction")[gene, :].array

    return points_graph(;
        x_axis_title = "$(gene) correlation change",
        y_axis_title = "$(gene) fraction",
        points_colors_title = "type",
        points_xs = delta_correlation_per_base_block[block_indices],
        points_ys = linear_fraction_per_block[block_indices],
        points_colors = type_per_block === nothing ? nothing : type_per_block[block_indices],
        points_hovers = entries_hovers(
            "block" => axis_vector(base_daf, "block")[block_indices],
            "type" => type_per_block === nothing ? nothing : type_per_block[block_indices],
            "base correlation" => rounded(base_correlation_per_base_block[block_indices]),
            "correlation" => rounded(
                base_correlation_per_base_block[block_indices] .+ delta_correlation_per_base_block[block_indices],
            ),
            "change" => rounded(delta_correlation_per_base_block[block_indices]),
        ),
        configuration = PointsGraphConfiguration(;
            y_axis = AxisConfiguration(; log_scale = Log2Scale, log_regularization = gene_fraction_regularization),
            points = ScattersConfiguration(; colors = colors_configuration),
        ),
    )
end

# The correlation of one gene with the cells of each base block's neighborhood.
function correlation_per_base_block(daf::DafReader, gene::AbstractString)::AbstractVector{<:AbstractFloat}
    return get_matrix(daf, "gene", "base_block", "correlation_between_base_neighborhood_cells_and_punctuated_metacells")[
        gene,
        :,
    ].array
end

# Correlations as a reader reads them, rather than as many digits as a `Float32` prints.
function rounded(values::AbstractVector{<:AbstractFloat})::Vector{Float32}
    return round.(values; digits = 3)
end

# The color of each point of some axis, and the configuration for drawing it. A type which has no color to be drawn in
# is not shown at all, so such a graph is drawn as if it had no types.
function points_type_colors(
    daf::DafReader,
    axis::AbstractString,
)::Tuple{Maybe{AbstractVector{<:AbstractString}}, ColorsConfiguration}
    type_per_point, colors_configuration = type_colors(daf, axis; show_legend = true)
    if colors_configuration === nothing
        return (nothing, ColorsConfiguration())
    else
        return (type_per_point, colors_configuration)
    end
end

end  # module
