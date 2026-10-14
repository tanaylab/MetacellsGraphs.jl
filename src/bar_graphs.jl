"""
Bar graphs of a metacells repository.

Unlike the other graphs, a graph here picks its own bars. It reads a report computed by `Metacells`, keeps the genes
worth showing, and draws the report's values about them. As with the other graphs, annotations are a second call, e.g.
[`add_changed_genes_default_annotations!`](@ref).
"""
module BarGraphs

export add_changed_genes_default_annotations!
export degraded_genes_graph
export improved_genes_graph
export put_changed_genes_configuration!

using DataAxesFormats
using DataFrames
using Metacells
using SomeGraphs
using TanayLabUtilities

using ..DataSources

"""
    improved_genes_graph(;
        daf::DafReader,
        base_daf::DafReader,
        genes_count::Integer = 30,
        regulators_count::Integer = 5,
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells improved in the most of the base
neighborhoods, as a butterfly of the percentage of neighborhoods each degraded in (to the left) and improved in (to the
right).

Everything drawn here is a column of
[`compute_gene_report`](@extref Metacells Metacells.AnalyzeGenes.compute_gene_report), which is what says when a gene
is improved or degraded in a base block and which regulators it moves with; the graph picks the genes worth showing
and draws them. A gene which is not a marker is not in the report and is therefore not shown.

The hover of each bar gives the percentage of base neighborhoods it shows. It then says how often the gene is in no
module at all in that side's base blocks, and which `regulators_count` regulators it is most often in a module with
there - that is, what the sharpening moved the gene with, where the bars say only how often it moved it. Each wing asks
about its own base blocks, so the two hovers of a gene name different regulators.

Lateral and regulator genes are shown like any other. The graph doesn't annotate them; that is a second call (see
[`add_changed_genes_default_annotations!`](@ref)). A lateral gene high in this graph is worth looking at: it is a gene
the analysis was told to ignore, which the metacells nevertheless describe better than the base does.
"""
function improved_genes_graph(;
    daf::DafReader,
    base_daf::DafReader,
    genes_count::Integer = 30,
    regulators_count::Integer = 5,
)::SeriesBarsGraph
    return changed_genes_graph(; daf, base_daf, genes_count, regulators_count, by_improved = true)
end

"""
    degraded_genes_graph(;
        daf::DafReader,
        base_daf::DafReader,
        genes_count::Integer = 30,
        regulators_count::Integer = 5,
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells degraded in the most of the base
neighborhoods.

This is [`improved_genes_graph`](@ref) picking its genes by the other side, and reads the same way; the two graphs show
the same two series and differ only in which of them decides what is worth showing.
"""
function degraded_genes_graph(;
    daf::DafReader,
    base_daf::DafReader,
    genes_count::Integer = 30,
    regulators_count::Integer = 5,
)::SeriesBarsGraph
    return changed_genes_graph(; daf, base_daf, genes_count, regulators_count, by_improved = false)
end

# The genes whose correlation changed the most, by either of the two directions. The graphs differ only in which
# direction picks them, and both show both, so they are one function.
function changed_genes_graph(;
    daf::DafReader,
    base_daf::DafReader,
    genes_count::Integer,
    regulators_count::Integer,
    by_improved::Bool,
)::SeriesBarsGraph
    @assert genes_count > 0
    @assert regulators_count >= 0

    report = compute_gene_report(; daf, base_daf, regulators_count)  # NOJET
    shown_report = get_top_changed_genes_report(report; genes_count, by_improved)

    graph = series_bars_graph()
    put_changed_genes_configuration!(graph.configuration)

    fill_column_names_data!(bars_entities(graph), shown_report; column = "gene")

    # The series are the two sides of the butterfly, so their colors say which side a bar is on. The hover of a bar
    # gives its value as a percent, as the axis does, then what the sharpening moved the gene with.
    for (prefix, side_name, color) in (("deg", "degraded", "darkred"), ("imp", "improved", "darkblue"))
        series = add_series!(graph, SeriesData(; name = side_name, color))
        fill_column_vector_data!(  # NOJET
            series.values,
            shown_report;
            column = "$(prefix)_f",
            title = "Base neighborhoods",
            hover_format = percent_format(; digits = 1),
        )
        fill_module_regulators_hovers!(series.entities, shown_report; prefix, side_name, regulators_count)
    end

    return graph
end

"""
    put_changed_genes_configuration!(configuration::SeriesBarsGraphConfiguration)::Nothing

Lay out a butterfly of genes, as [`improved_genes_graph`](@ref) and [`degraded_genes_graph`](@ref) show it. The genes
are horizontal bars, and the two series are mirrored around them. The values are percentages of the base
neighborhoods.
"""
function put_changed_genes_configuration!(configuration::SeriesBarsGraphConfiguration)::Nothing
    configuration.values_orientation = HorizontalValues
    configuration.mirrored = true
    configuration.value_axis.title = "Base neighborhoods"
    configuration.value_axis.scale.percent = true
    configuration.bar_axis.title = "Genes"
    return nothing
end

"""
    add_changed_genes_default_annotations!(
        graph::SeriesBarsGraph;
        daf::DafReader,
        base_daf::DafReader,
    )::Nothing

Add the default annotations of the bars (the genes) of an [`improved_genes_graph`](@ref) or a
[`degraded_genes_graph`](@ref) of the `daf` and the `base_daf`: whether each gene is lateral, a regulator, and a marker
(see [`fill_is_lateral!`](@ref MetacellsGraphs.DataSources.fill_is_lateral!) and its siblings). The graph measures the
metacells against the base, so the annotations are of the `base_daf`. An annotation whose data the `base_daf` doesn't
have is left out.
"""
function add_changed_genes_default_annotations!(
    graph::SeriesBarsGraph;
    daf::DafReader,  # NOLINT
    base_daf::DafReader,
)::Nothing
    for (property, fill_property!) in (
        ("is_lateral", fill_is_lateral!),  # NOLINT
        ("is_regulator", fill_is_regulator!),  # NOLINT
        ("is_marker", fill_is_marker!),  # NOLINT
    )
        if has_vector(base_daf, "gene", property)
            fill_property!(annotations_colors_vector_fields(graph, add_annotation!(graph)), base_daf)
        end
    end
    return nothing
end

# The `genes_count` genes of the `report` which changed in the most base blocks, on the improved or the degraded side.
# They are in ascending order, because the bars of a horizontal graph run from the bottom up, and the gene the graph is
# named for belongs at its top. A stable sort of a report which arrives in marker rank order means genes which moved in
# equally many base blocks are ranked by how good a marker they are, rather than by whatever the sort happens to do.
function get_top_changed_genes_report(report::DataFrame; genes_count::Integer, by_improved::Bool)::DataFrame
    sorted_report = sort(report, by_improved ? "imp_f" : "deg_f"; alg = MergeSort)
    return last(sorted_report, genes_count)
end

end  # module
