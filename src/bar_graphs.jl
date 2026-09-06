"""
Bar graphs of a metacells repository.
"""
module BarGraphs

export degraded_genes_graph
export improved_genes_graph

using DataAxesFormats
using DataFrames
using Metacells
using SomeGraphs
using TanayLabUtilities

# Needed because of JET:
import Metacells.Contracts.base_block_axis
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_correlation_between_markers_per_gene_per_gene
import Metacells.Contracts.matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_regulator_per_gene
import Metacells.Contracts.matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_regulator_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_regulator_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_is_transcription_factor_per_gene
import Metacells.Contracts.vector_of_marker_rank_per_gene
import Metacells.Contracts.vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_gene
import Metacells.Contracts.vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_gene

"""
    improved_genes_graph(;
        daf::DafReader,
        base_daf::DafReader,
        genes_count::Integer = $(DEFAULT.genes_count),
        regulators_count::Integer = $(DEFAULT.regulators_count),
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells improved in the most of the base
neighborhoods, as a butterfly of the percentage of neighborhoods each degraded in (to the left) and improved in (to the
right).

Everything drawn here is a column of
[`compute_gene_report`](@extref Metacells Metacells.AnalyzeGenes.compute_gene_report), which is what says when a gene
is improved or degraded in a base block and which regulators it moves with; the graph picks the genes worth showing
and draws them. A gene which is not a marker is not in the report and is therefore not shown.

The hover of each bar says how often the gene is in no module at all in that side's base blocks, and which
`regulators_count` regulators it most often shares one with there - that is, what the sharpening moved the gene with,
where the bars say only how often it moved it. Each wing asks about its own base blocks, so the two hovers of a gene
name different regulators.

Lateral and regulator genes are shown like any other, and are marked as such in the annotations between the two sides.
A lateral gene high in this graph is worth looking at: it is a gene the analysis was told to ignore, which the
metacells nevertheless describe better than the base does.

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
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_marker_rank_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_transcription_factor_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
        vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_gene(RequiredInput),
        vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_gene(RequiredInput),
        matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_regulator_per_gene(
            RequiredInput,
        ),
        matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_regulator_per_gene(
            RequiredInput,
        ),
    ],
) Contract(;
    name = "base_daf",
    link = Metacells,
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
    ],
) function improved_genes_graph(;
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
        genes_count::Integer = $(DEFAULT.genes_count),
        regulators_count::Integer = $(DEFAULT.regulators_count),
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells degraded in the most of the base
neighborhoods.

This is [`improved_genes_graph`](@ref) picking its genes by the other side, and reads the same way; the two graphs show
the same two series and differ only in which of them decides what is worth showing.

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
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_marker_rank_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_transcription_factor_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
        vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_gene(RequiredInput),
        vector_of_mean_no_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_gene(RequiredInput),
        matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_improved_base_blocks_per_regulator_per_gene(
            RequiredInput,
        ),
        matrix_of_mean_shared_module_fraction_in_base_neighborhood_cells_at_degraded_base_blocks_per_regulator_per_gene(
            RequiredInput,
        ),
    ],
) Contract(;
    name = "base_daf",
    link = Metacells,
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
    ],
) function degraded_genes_graph(;
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
        bar_axis_title = "Genes",
        value_axis_title = "Base neighborhoods",
        series_bars_values = [shown_report[!, "deg_f"], shown_report[!, "imp_f"]],
        bars_names = shown_report[!, "gene"],
        series_bars_hovers = [
            module_sharing_hovers(shown_report, "deg", "degraded", regulators_count),
            module_sharing_hovers(shown_report, "imp", "improved", regulators_count),
        ],
        series_names = ["degraded", "improved"],
        series_colors = ["darkred", "darkblue"],
        configuration = SeriesBarsGraphConfiguration(;
            values_orientation = HorizontalValues,
            value_axis = AxisConfiguration(; percent = true),
            mirrored = true,
        ),
    )

    # The annotations are not part of what the constructor takes, so they are attached to the data.
    graph.data.bars_annotations = [
        mask_annotation("is lateral", shown_report[!, "lat?"]),
        mask_annotation("is regulator", shown_report[!, "reg?"]),
    ]

    return graph
end

# The hover of each shown gene in one of the two series, saying how often the gene is in no module at all in the base
# blocks of that side, and which regulators it most often shares one with there. Each side is asked about its own base
# blocks, so the two wings of a gene say different things.
function module_sharing_hovers(
    shown_report::DataFrame,
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::Vector{AbstractString}
    return AbstractString[
        join(gene_hover_lines(shown_report, row_index, prefix, side_name, regulators_count), "<br>") for
        row_index in 1:nrow(shown_report)
    ]
end

# The lines of one gene's hover: the gene, how often it is in no module, and then a line per regulator it shares one
# with. The report pads a side which gives the gene fewer regulators than asked for with empty names, which are the
# regulators there are nothing to say about.
function gene_hover_lines(
    shown_report::DataFrame,
    row_index::Integer,
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::Vector{String}
    lines = String[
        shown_report[row_index, "gene"],
        "$(side_name): in no module in $(percent(shown_report[row_index, "$(prefix)_no_mod_f"])) of the cells",
    ]
    for rank in 1:regulators_count
        regulator_name = shown_report[row_index, "$(prefix)_reg$(rank)"]
        if regulator_name != ""
            push!(lines, "- $(regulator_name): $(percent(shown_report[row_index, "$(prefix)_reg$(rank)_f"]))")
        end
    end
    return lines
end

# Fractions in the units the value axis shows, so that the hover and the bars are read the same way.
function percent(fraction::AbstractFloat)::String
    return "$(round(100 * fraction; digits = 1))%"
end

# One mask of the genes, shown between the two sides. A mask is drawn as a category rather than as a number, so that it
# is read as what it says rather than as a quantity.
function mask_annotation(title::AbstractString, mask_per_gene::AbstractVector{Bool})::AnnotationData
    values = [mask ? "yes" : "no" for mask in mask_per_gene]
    return AnnotationData(;
        title,
        values,
        hovers = values,
        colors = ColorsConfiguration(; palette = Dict("yes" => "black", "no" => "lightgrey")),
    )
end

end  # module
