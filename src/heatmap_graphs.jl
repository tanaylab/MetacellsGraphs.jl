"""
Heatmap graphs of a metacells repository.

A graph function fills in the data the graph is about, and nothing else. In particular it does not group the columns
and does not annotate them. The caller says how, in a second call, so that any property can group or annotate:

    graph = genes_heatmap_graph(daf; genes = get_skeleton_gene_indices(daf))
    fill_type!(columns_annotations_colors_vector_fields(graph, add_columns_annotation!(graph)), daf)
    fill_global_flow_order!(columns_groups_vector_data_fields(graph), daf)
    fill_block!(columns_subgroups_vector_data_fields(graph), daf)

Grouping a level constrains the clustering to it and puts a gap between the groups. Naming both levels nests them: the
blocks are grouped inside the types, which is how they nest. A level given numbers is laid out in the order of these
numbers, and a level given names is laid out by the clustering.

The second call must name the same `axis` as the first. It needn't name the `entries`, since it fills the entries the
graph was named after.
"""
module HeatmapGraphs

export genes_heatmap_graph
export put_genes_heatmap_configuration!

using DataAxesFormats
using SomeGraphs
using TanayLabUtilities

using ..DataSources

"""
    genes_heatmap_graph(
        daf::DafReader;
        axis::AbstractString = "metacell",
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        rows_axis_title::Maybe{AbstractString} = "Genes",
        columns_axis_title::Maybe{AbstractString} = "\$(uppercasefirst(axis))s",
        max_fold::Real = $(MAX_FOLD_FOR_GRAPHS),
    )::HeatmapGraph

The expression of the `genes` in each of the `entries` of the `axis`, as the fold factor of the gene from its median
across the `entries`.

The genes are the rows and the entries are the columns; both are clustered. Pick the genes with a data source such as
[`get_top_marker_gene_indices`](@ref MetacellsGraphs.DataSources.get_top_marker_gene_indices) or
[`get_skeleton_gene_indices`](@ref MetacellsGraphs.DataSources.get_skeleton_gene_indices), and name them in the
`rows_axis_title`.

There are too many columns to label, so they are named in the hover of each cell instead. The rows are labelled.
"""
function genes_heatmap_graph(
    daf::DafReader;
    axis::AbstractString = "metacell",
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    rows_axis_title::Maybe{AbstractString} = "Genes",
    columns_axis_title::Maybe{AbstractString} = "$(uppercasefirst(axis))s",
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
)::HeatmapGraph
    graph = heatmap_graph()
    fill_genes_fold_matrix!(entries_matrix_fields(graph), daf; genes, axis, entries, max_fold)
    n_genes, n_entries = size(graph.data.entries.matrix)
    put_genes_heatmap_configuration!(graph.configuration; n_genes, n_entries, rows_axis_title, columns_axis_title)
    return graph
end

"""
    put_genes_heatmap_configuration!(
        configuration::HeatmapGraphConfiguration;
        n_genes::Integer,
        n_entries::Integer,
        rows_axis_title::Maybe{AbstractString} = "Genes",
        columns_axis_title::Maybe{AbstractString} = "Metacells",
    )::Nothing

Lay out a heatmap of `n_genes` genes (the rows) by `n_entries` entries (the columns), as
[`genes_heatmap_graph`](@ref) shows it. The rows are labelled and named by the `rows_axis_title`. There are too many
columns to label, so they are named only by the `columns_axis_title` and in the hover of each cell. Both sides are
clustered, except a side with a single entry, which has nothing to be ordered against.
"""
function put_genes_heatmap_configuration!(
    configuration::HeatmapGraphConfiguration;
    n_genes::Integer,
    n_entries::Integer,
    rows_axis_title::Maybe{AbstractString} = "Genes",
    columns_axis_title::Maybe{AbstractString} = "Metacells",
)::Nothing
    configuration.figure.margins.left = 100
    configuration.origin = HeatmapTopLeft
    configuration.rows.title = rows_axis_title
    configuration.rows.order_source = n_genes > 1 ? OptimalTreeReorder : nothing
    configuration.columns.title = columns_axis_title
    configuration.columns.show_ticks = false
    configuration.columns.order_source = n_entries > 1 ? OptimalTreeReorder : nothing
    return nothing
end

end  # module
