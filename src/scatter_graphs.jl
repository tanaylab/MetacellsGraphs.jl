"""
Scatter graphs of a metacells repository.
"""
module ScatterGraphs

export blocks_gene_gene_graph
export blocks_umap_graph
export metacells_gene_gene_graph
export metacells_umap_graph

using DataAxesFormats
using Metacells
using SomeGraphs
using TanayLabUtilities

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
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
    points_colors, colors_configuration = type_colors(daf, axis)

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
    points_colors, colors_configuration = type_colors(daf, axis)

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

# The color of each point of some axis, and the configuration for drawing it. The type is optional, and so is coloring
# by it. Without the colors of the type axis there is nothing to map a type name to, so a repository which has types
# but no colors for them is drawn as if it had no types at all.
function type_colors(
    daf::DafReader,
    axis::AbstractString,
)::Tuple{Maybe{AbstractVector{<:AbstractString}}, ColorsConfiguration}
    type_per_point = get_vector(daf, axis, "type"; default = nothing)
    if type_per_point === nothing || !has_axis(daf, "type")
        return (nothing, ColorsConfiguration())
    else
        return (
            type_per_point.array,
            ColorsConfiguration(; palette = get_vector(daf, "type", "color"), show_legend = true),
        )
    end
end

end  # module
