"""
Bar graphs of a metacells repository.

Unlike the other graphs, a graph here picks its own bars, so there is nothing for a caller to fill in afterwards. It
reads a report computed by `Metacells`, keeps the genes worth showing, and draws everything about them.
"""
module BarGraphs

export degraded_genes_graph
export improved_genes_graph

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

The hover of each bar says how often the gene is in no module at all in that side's base blocks, and which
`regulators_count` regulators it is most often in a module with there - that is, what the sharpening moved the gene
with, where the bars say only how often it moved it. Each wing asks about its own base blocks, so the two hovers of a
gene name different regulators.

Lateral and regulator genes are shown like any other, and are marked as such in the annotations between the two sides.
A lateral gene high in this graph is worth looking at: it is a gene the analysis was told to ignore, which the
metacells nevertheless describe better than the base does.
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

    # Ascending, because the bars of a horizontal graph run from the bottom up, and the gene the graph is named for
    # belongs at its top. A stable sort of a report which arrives in marker rank order means genes which moved in
    # equally many base blocks are ranked by how good a marker they are, rather than by whatever the sort happens to do.
    sort!(report, by_improved ? "imp_f" : "deg_f"; alg = MergeSort)
    shown_report = last(report, genes_count)

    graph = series_bars_graph(;
        configuration = SeriesBarsGraphConfiguration(;
            values_orientation = HorizontalValues,
            value_axis = AxisConfiguration(;
                title = "Base neighborhoods",
                scale = ScaleConfiguration(; percent = true),
            ),
            bar_axis = AxisConfiguration(; title = "Genes"),
            mirrored = true,
        ),
    )

    fill_column_names_data!(bars_entities(graph), shown_report; column = "gene")

    # The series are the two sides of the butterfly, so their colors say which side a bar is on. The bar shows the
    # value, so its hover says only what the sharpening moved the gene with.
    for (prefix, side_name, color) in (("deg", "degraded", "darkred"), ("imp", "improved", "darkblue"))
        series = SeriesData(; name = side_name, color)
        add_series!(graph, series)
        fill_column_vector_data!(series.values, shown_report; column = "$(prefix)_f")
        fill_module_regulators_hovers!(series.bars, shown_report; prefix, side_name, regulators_count)
    end

    for (title, column) in (("is lateral", "lat?"), ("is regulator", "reg?"))
        index = add_annotation!(graph)
        fill_column_boolean_annotation!(annotations_colors_vector_fields(graph, index), shown_report; column, title)
    end

    return graph
end

end  # module
