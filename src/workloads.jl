# Compile the uses of the data sources and the graphs which a notebook makes, when the package is precompiled, so that
# their first use in a session is fast.
module Workloads

using DataAxesFormats
using PrecompileTools
using SomeGraphs
using SparseArrays

using ..DataSources
using ..HeatmapGraphs
using ..ScatterGraphs

import SomeGraphs.Common.graph_to_json

# Write a small chain of two repositories into the `directory`, holding the properties the data sources read, with the
# same types as a metacells repository holds them. The `top` repository's base is the `base` repository next to it.
function write_workload_daf(directory::AbstractString)::Nothing
    genes = ["A", "B", "C", "D", "E"]
    metacells = ["M1", "M2", "M3", "M4"]

    base = FilesDaf(joinpath(directory, "base"), "w"; name = "base")
    add_axis!(base, "gene", genes)
    add_axis!(base, "metacell", metacells)
    add_axis!(base, "type", ["T1", "T2"])
    set_vector!(base, "type", "color", ["red", "blue"])
    set_vector!(base, "type", "global_flow_order", UInt32[2, 1])
    set_vector!(base, "gene", "is_lateral", sparsevec([2], [true], 5))  # NOJET
    set_vector!(base, "gene", "is_marker", [true, true, true, false, false])
    set_vector!(base, "gene", "marker_rank", UInt32[1, 2, 3, typemax(UInt32), typemax(UInt32)])
    set_vector!(base, "gene", "is_skeleton", [true, true, false, false, false])
    set_vector!(base, "metacell", "type", ["T1", "T1", "T2", "T2"])
    set_vector!(base, "metacell", "n_cells", UInt32[10, 20, 30, 40])
    set_vector!(base, "metacell", "total_UMIs", UInt32[100, 200, 300, 400])
    set_vector!(base, "metacell", "is_doublet", [false, true, false, false])
    set_vector!(base, "metacell", "umap_x", Float32[0, 1, 2, 3])
    set_vector!(base, "metacell", "umap_y", Float32[3, 2, 1, 0])
    fractions = Float32[
        0.1 0.2 0.3 0.4
        0.4 0.3 0.2 0.1
        0.0 0.1 0.0 0.1
        0.2 0.2 0.2 0.2
        0.1 0.1 0.1 0.1
    ]
    set_matrix!(base, "gene", "metacell", "linear_fraction", fractions)
    set_matrix!(base, "gene", "metacell", "log_linear_fraction", log2.(fractions .+ Float32(1e-5)))

    top = FilesDaf(joinpath(directory, "top"), "w"; name = "top")
    add_axis!(top, "gene", genes)
    add_axis!(top, "metacell", metacells)
    add_axis!(top, "block", ["B1", "B2"])  # NOJET
    set_base_daf_links!(top, base)
    set_vector!(top, "metacell", "block", ["B1", "B1", "B2", "B2"])
    set_vector!(top, "block", "type", ["T1", "T2"])
    set_vector!(top, "block", "n_metacells", UInt32[2, 2])
    set_vector!(top, "block", "n_cells", UInt32[30, 70])
    set_vector!(top, "block", "total_UMIs", UInt32[300, 700])
    set_matrix!(top, "gene", "block", "linear_fraction", fractions[:, [1, 3]])
    set_matrix!(top, "gene", "block", "log_linear_fraction", log2.(fractions[:, [1, 3]] .+ Float32(1e-5)))  # NOJET
    set_matrix!(top, "block", "block", "is_in_neighborhood", [true false; false true])
    is_in_environment = sparse([true false; true false; false true; false true])
    set_matrix!(top, "metacell", "block", "is_in_environment", is_in_environment)  # NOJET
    return nothing
end

# Make the graphs a notebook makes, from a chain of repositories opened as a notebook opens them, filling them by each
# data source a notebook offers, and convert each to JSON (which is what showing a graph does). This runs when the
# package is precompiled.
function run_workload()::Nothing
    mktempdir() do directory
        write_workload_daf(directory)
        daf = complete_daf(joinpath(directory, "top"))  # NOJET
        for axis in ("metacell", "block")
            points = gene_gene_graph(daf; axis, x_gene = "A", y_gene = "B")
            fill_type!(points_colors_vector_fields(points), daf; axis)  # NOLINT
            fill_n_cells!(points_sizes_vector_fields(points), daf; axis)  # NOLINT
            fill_total_UMIs!(points_entities(points), daf; axis)  # NOLINT
            fill_mean_total_UMIs_per_cell!(points_entities(points), daf; axis)
            fill_gene_expression!(points_entities(points), daf; axis, gene = "C")
            fill_global_flow_order!(points_entities(points), daf; axis)
            graph_to_json(points)  # NOJET

            for heatmap in (
                genes_heatmap_graph(daf; axis),
                genes_heatmap_graph(
                    daf;
                    axis,
                    genes = get_top_marker_gene_indices(daf; markers_count = 3),
                    entries = get_is_in_neighborhood_vector(daf; axis, block = "B1") .|
                              get_is_in_categories_vector(daf; axis, property = "type", categories = ["T2"]),
                ),
            )
                add_genes_heatmap_default_rows_annotations!(heatmap, daf)
                add_genes_heatmap_default_columns_annotations!(heatmap, daf; axis)
                annotation = rows_annotations_colors_vector_fields(heatmap, add_rows_annotation!(heatmap))
                fill_is_marker!(annotation, daf)  # NOLINT
                graph_to_json(heatmap)  # NOJET
            end
        end

        points = gene_gene_graph(daf; axis = "metacell", x_gene = "A", y_gene = "B")
        fill_block!(points_colors_vector_fields(points), daf; axis = "metacell")  # NOLINT
        fill_boolean!(points_entities(points), daf; axis = "metacell", property = "is_doublet")
        fill_umap!(points_entities(points), daf; axis = "metacell", coordinate = "x")
        get_is_in_environment_vector(daf; axis = "metacell", block = "B1")
        get_skeleton_gene_indices(daf)
        graph_to_json(points)

        points = gene_gene_graph(daf; axis = "block", x_gene = "A", y_gene = "B")
        fill_n_metacells!(points_sizes_vector_fields(points), daf; axis = "block")  # NOLINT
        fill_mean_cells_per_metacell!(points_entities(points), daf; axis = "block")
        fill_mean_total_UMIs_per_metacell!(points_entities(points), daf; axis = "block")
        graph_to_json(points)
        return nothing
    end
    return nothing
end

@compile_workload begin
    run_workload()
end

end  # module
