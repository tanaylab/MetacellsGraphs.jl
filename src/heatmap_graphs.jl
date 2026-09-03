"""
Heatmap graphs of a metacells repository.
"""
module HeatmapGraphs

export markers_blocks_heatmap_graph
export markers_metacells_heatmap_graph
export skeletons_blocks_heatmap_graph
export skeletons_metacells_heatmap_graph

using DataAxesFormats
using Metacells
using SomeGraphs
using Statistics
using TanayLabUtilities

using ..Utilities

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.type_axis
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_color_per_type
import Metacells.Contracts.vector_of_global_flow_order_per_type
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_marker_rank_per_gene
import Metacells.Contracts.vector_of_type_per_block
import Metacells.Contracts.vector_of_type_per_metacell

"""
The extreme of the color scale of the fold factors. A fold factor is a log2, so a gene is drawn in full blue where it
is 3 below its median (an eighth of it), in white where it is at it, and in full red where it is 3 above it (eight
times it); anything more extreme is drawn in the same color as this. All the heatmaps use the same extreme, so that
they can be compared to each other.
"""
MAX_FOLD_FOR_GRAPHS::Float64 = 3

"""
The fraction of a grouped axis given over to the gaps between its groups. A gap is measured in entries, so a gap of a
fixed number of them is invisible in an axis holding thousands and is most of an axis holding a few dozen. Sizing the
gaps so that together they take this fraction of the axis makes them visible whatever the axis holds, while an axis
holding many small groups keeps the single entry it takes to separate them at all instead of becoming mostly gaps.
"""
TOTAL_GROUPS_GAP_FRACTION::Float64 = 1 / 20

"""
    markers_metacells_heatmap_graph(
        daf::DafReader;
        markers_count::Integer = $(DEFAULT.markers_count),
        group_by_type::Bool = $(DEFAULT.group_by_type),
        group_by_block::Bool = $(DEFAULT.group_by_block),
        use_global_flow_order::Bool = $(DEFAULT.use_global_flow_order),
    )::HeatmapGraph

The expression of the `markers_count` best marker genes in each metacell, as the fold factor of the gene from its
median across all the metacells.

The genes are the rows and the metacells are the columns; both are clustered. There are too many metacells to name
them, so each is identified in its hover, together with the block and the type it belongs to. If the metacells have a
type, it is also shown as a color annotation of the columns, using the colors of the type axis.

Given `group_by_type` and/or `group_by_block`, the columns are grouped by the type and/or by the block of each
metacell, which both constrains the clustering to the groups and puts a gap between them. Given both, the grouping is
hierarchical: the blocks are grouped inside the types, which is how they nest. Given `use_global_flow_order`, the type
groups run in the global flow order of the types, so that every graph lays them out the same way; otherwise the
clustering places them, as it always places the blocks.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_marker_rank_per_gene(RequiredInput),
        vector_of_block_per_metacell(OptionalInput),
        vector_of_type_per_metacell(OptionalInput),
        vector_of_color_per_type(OptionalInput),
        vector_of_global_flow_order_per_type(OptionalInput),
    ],
) function markers_metacells_heatmap_graph(
    daf::DafReader;
    markers_count::Integer = 100,
    group_by_type::Bool = false,
    group_by_block::Bool = false,
    use_global_flow_order::Bool = false,
)::HeatmapGraph
    return genes_heatmap_graph(
        daf;
        axis = "metacell",
        gene_indices = top_marker_gene_indices(daf, markers_count),
        rows_axis_title = "Marker genes",
        columns_axis_title = "Metacells",
        group_by_type,
        group_by_block,
        use_global_flow_order,
    )
end

"""
    skeletons_metacells_heatmap_graph(
        daf::DafReader;
        group_by_type::Bool = $(DEFAULT.group_by_type),
        group_by_block::Bool = $(DEFAULT.group_by_block),
        use_global_flow_order::Bool = $(DEFAULT.use_global_flow_order),
    )::HeatmapGraph

The expression of the skeleton genes in each metacell, as the fold factor of the gene from its median across all the
metacells.

This is [`markers_metacells_heatmap_graph`](@ref) of the skeleton genes, and reads the same way. The skeletons are the
genes the blocks were computed from, so this shows what the blocks were told apart by, rather than everything which
tells the metacells apart.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_block_per_metacell(OptionalInput),
        vector_of_type_per_metacell(OptionalInput),
        vector_of_color_per_type(OptionalInput),
        vector_of_global_flow_order_per_type(OptionalInput),
    ],
) function skeletons_metacells_heatmap_graph(
    daf::DafReader;
    group_by_type::Bool = false,
    group_by_block::Bool = false,
    use_global_flow_order::Bool = false,
)::HeatmapGraph
    return genes_heatmap_graph(
        daf;
        axis = "metacell",
        gene_indices = skeleton_gene_indices(daf),
        rows_axis_title = "Skeleton genes",
        columns_axis_title = "Metacells",
        group_by_type,
        group_by_block,
        use_global_flow_order,
    )
end

"""
    markers_blocks_heatmap_graph(
        daf::DafReader;
        markers_count::Integer = $(DEFAULT.markers_count),
        group_by_type::Bool = $(DEFAULT.group_by_type),
        use_global_flow_order::Bool = $(DEFAULT.use_global_flow_order),
    )::HeatmapGraph

The expression of the `markers_count` best marker genes in each block, as the fold factor of the gene from its median
across all the blocks.

This is [`markers_metacells_heatmap_graph`](@ref) of the blocks the metacells were grouped into, and reads the same
way, except that a block belongs to nothing below its type, so there is no `group_by_block` and the hover names just
the block and its type.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_marker_rank_per_gene(RequiredInput),
        vector_of_type_per_block(OptionalInput),
        vector_of_color_per_type(OptionalInput),
        vector_of_global_flow_order_per_type(OptionalInput),
    ],
) function markers_blocks_heatmap_graph(
    daf::DafReader;
    markers_count::Integer = 100,
    group_by_type::Bool = false,
    use_global_flow_order::Bool = false,
)::HeatmapGraph
    return genes_heatmap_graph(
        daf;
        axis = "block",
        gene_indices = top_marker_gene_indices(daf, markers_count),
        rows_axis_title = "Marker genes",
        columns_axis_title = "Blocks",
        group_by_type,
        group_by_block = false,
        use_global_flow_order,
    )
end

"""
    skeletons_blocks_heatmap_graph(
        daf::DafReader;
        group_by_type::Bool = $(DEFAULT.group_by_type),
        use_global_flow_order::Bool = $(DEFAULT.use_global_flow_order),
    )::HeatmapGraph

The expression of the skeleton genes in each block, as the fold factor of the gene from its median across all the
blocks.

This is [`markers_blocks_heatmap_graph`](@ref) of the skeleton genes, and reads the same way. The skeletons are the
genes the blocks were computed from, so this shows what the blocks were told apart by, rather than everything which
tells them apart.

$(CONTRACT)
"""
@computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput), type_axis(OptionalInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_type_per_block(OptionalInput),
        vector_of_color_per_type(OptionalInput),
        vector_of_global_flow_order_per_type(OptionalInput),
    ],
) function skeletons_blocks_heatmap_graph(
    daf::DafReader;
    group_by_type::Bool = false,
    use_global_flow_order::Bool = false,
)::HeatmapGraph
    return genes_heatmap_graph(
        daf;
        axis = "block",
        gene_indices = skeleton_gene_indices(daf),
        rows_axis_title = "Skeleton genes",
        columns_axis_title = "Blocks",
        group_by_type,
        group_by_block = false,
        use_global_flow_order,
    )
end

# The best `markers_count` marker genes. A rank of 1 is the gene which best distinguishes between the cell states, and
# a non-marker gene is ranked above any count that can be asked for, so this is simply the genes ranked that high; a
# repository with fewer markers than that shows all of them.
function top_marker_gene_indices(daf::DafReader, markers_count::Integer)::Vector{Int}
    @assert markers_count > 0
    return findall(get_vector(daf, "gene", "marker_rank").array .<= markers_count)
end

# The skeleton genes, all of them; unlike the markers they are not ranked, and there are few enough of them to show.
function skeleton_gene_indices(daf::DafReader)::Vector{Int}
    return findall(get_vector(daf, "gene", "is_skeleton").array)
end

# The expression of some genes over the entries of some axis. The graphs differ only in which genes are the rows and
# which axis is the columns, and every axis which has a `log_linear_fraction` per gene also has an optional `type`, so
# they are one function.
function genes_heatmap_graph(
    daf::DafReader;
    axis::AbstractString,
    gene_indices::AbstractVector{<:Integer},
    rows_axis_title::AbstractString,
    columns_axis_title::AbstractString,
    group_by_type::Bool,
    group_by_block::Bool,
    use_global_flow_order::Bool,
)::HeatmapGraph
    @assert group_by_type || !use_global_flow_order "use_global_flow_order without group_by_type"

    log_linear_fraction_per_gene_per_entry = get_matrix(daf, "gene", axis, "log_linear_fraction").array
    @views log_linear_fraction_per_shown_gene_per_entry = log_linear_fraction_per_gene_per_entry[gene_indices, :]
    fold_per_shown_gene_per_entry = Float32.(  # NOJET
        log_linear_fraction_per_shown_gene_per_entry .- median(log_linear_fraction_per_shown_gene_per_entry; dims = 2),
    )

    name_per_entry = axis_vector(daf, axis)
    type_per_entry, colors_configuration = type_colors(daf, axis)

    # Only the metacells belong to a block; a block belongs to nothing below its type.
    block_per_entry = axis == "block" ? nothing : get_vector(daf, axis, "block"; default = nothing)

    columns_groups, columns_subgroups = columns_grouping(
        daf;
        axis,
        type_per_entry,
        block_per_entry,
        group_by_type,
        group_by_block,
        use_global_flow_order,
    )

    return heatmap_graph(;
        x_axis_title = columns_axis_title,
        y_axis_title = rows_axis_title,
        entries_colors_title = "log2 fold from median",
        entries_values = fold_per_shown_gene_per_entry,
        rows_names = axis_vector(daf, "gene")[gene_indices],
        # There are too many entries to name them; they are still identified in the hover, together with everything
        # they belong to, so the nesting the graph is grouped by can be read off any of its columns.
        columns_hovers = entries_hovers(axis => name_per_entry, "block" => block_per_entry, "type" => type_per_entry),
        columns_annotations = type_annotations(type_per_entry, colors_configuration),
        columns_groups,
        columns_subgroups,
        configuration = HeatmapGraphConfiguration(;
            figure = FigureConfiguration(; margins = MarginsConfiguration(; left = 100, bottom = 100)),
            rows_reorder = length(gene_indices) > 1 ? OptimalHclust : nothing,
            columns_reorder = length(name_per_entry) > 1 ? OptimalHclust : nothing,
            columns_groups_gap = groups_gap(length(name_per_entry), columns_groups),
            origin = HeatmapTopLeft,
            entries_colors = ColorsConfiguration(;
                palette = "BuWtRd",
                axis = AxisConfiguration(; minimum = -MAX_FOLD_FOR_GRAPHS, maximum = MAX_FOLD_FOR_GRAPHS),
                show_legend = true,
            ),
        ),
    )
end

# The groups and the subgroups of the columns. Numbered groups are laid out in the order of their numbers, which is how
# the global flow order pins the types to the same places in every graph; named groups are laid out by the clustering,
# as the blocks always are.
function columns_grouping(
    daf::DafReader;
    axis::AbstractString,
    type_per_entry::Maybe{AbstractVector{<:AbstractString}},
    block_per_entry::Maybe{AbstractVector{<:AbstractString}},
    group_by_type::Bool,
    group_by_block::Bool,
    use_global_flow_order::Bool,
)::Tuple{Maybe{AbstractVector}, Maybe{AbstractVector}}
    @assert !group_by_type || type_per_entry !== nothing "group_by_type without a type per $(axis)"
    @assert !group_by_block || block_per_entry !== nothing "group_by_block without a block per $(axis)"

    if group_by_type
        type_groups = use_global_flow_order ? global_flow_groups(daf, type_per_entry) : type_per_entry  # NOJET
        return (type_groups, group_by_block ? block_per_entry : nothing)
    elseif group_by_block
        return (block_per_entry, nothing)
    else
        return (nothing, nothing)
    end
end

# The global flow order of the type of each entry, compacted to cover `1:N` over just the types which appear, as the
# clustering requires the numbers of the groups to.
function global_flow_groups(daf::DafReader, type_per_entry::AbstractVector{<:AbstractString})::Vector{UInt32}
    @assert has_vector(daf, "type", "global_flow_order") "use_global_flow_order without a global_flow_order per type"
    order_per_type_name =
        Dict{AbstractString, UInt32}(zip(axis_vector(daf, "type"), get_vector(daf, "type", "global_flow_order").array))
    order_per_entry = UInt32[order_per_type_name[type] for type in type_per_entry]
    group_per_order =
        Dict{UInt32, UInt32}(order => UInt32(group) for (group, order) in enumerate(sort(unique(order_per_entry))))
    return [group_per_order[order] for order in order_per_entry]
end

# The gap between the groups of an axis, in entries. An ungrouped axis has no gaps, and is given the smallest gap there
# is rather than a special case of its own.
function groups_gap(n_entries::Integer, group_per_entry::Maybe{AbstractVector})::Int
    if group_per_entry === nothing
        return 1
    end
    n_gaps = max(length(unique(group_per_entry)) - 1, 1)
    return max(1, round(Int, n_entries * TOTAL_GROUPS_GAP_FRACTION / n_gaps))
end

# The hover of each entry of an axis, naming the entry and everything it belongs to, skipping whatever the repository
# does not have.
function entries_hovers(fields::Pair{<:AbstractString, <:Maybe{AbstractVector}}...)::Vector{AbstractString}
    known_fields = Pair{AbstractString, AbstractVector}[field for field in fields if field[2] !== nothing]
    return AbstractString[
        join(String["$(label): $(value_per_entry[entry_index])" for (label, value_per_entry) in known_fields], "<br>")
        for entry_index in eachindex(known_fields[1][2])
    ]
end

# The color annotation of the type of each entry, or no annotation at all when there are no types to show or no colors
# to show them in.
function type_annotations(
    type_per_entry::Maybe{AbstractVector{<:AbstractString}},
    colors_configuration::Maybe{ColorsConfiguration},
)::Vector{AnnotationData}
    if type_per_entry === nothing || colors_configuration === nothing
        return AnnotationData[]
    else
        return [
            AnnotationData(;
                title = "type",
                values = type_per_entry,
                hovers = type_per_entry,
                colors = colors_configuration,
            ),
        ]
    end
end

end  # module
