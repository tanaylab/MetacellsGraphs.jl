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

The second call must name the same `axis` and `entries` as the first, or the vectors will not match.
"""
module HeatmapGraphs

export genes_heatmap_graph

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
    graph = heatmap_graph(;
        configuration = HeatmapGraphConfiguration(;
            figure = FigureConfiguration(; margins = MarginsConfiguration(; left = 100, bottom = 100)),
            origin = HeatmapTopLeft,
            rows = HeatmapAxisConfiguration(; title = rows_axis_title),
            columns = HeatmapAxisConfiguration(; title = columns_axis_title, show_ticks = false),
        ),
    )

    fill_genes_fold_matrix!(entries_matrix_fields(graph), daf; genes, axis, entries, max_fold)  # NOJET

    # A single entry has nothing to be ordered against, so there is nothing to cluster.
    n_genes, n_entries = size(graph.data.entries.matrix)
    graph.configuration.rows.reorder = n_genes > 1 ? OptimalHclust : nothing
    graph.configuration.columns.reorder = n_entries > 1 ? OptimalHclust : nothing

    return graph
end

end  # module
