"""
Bar graphs of a metacells repository.
"""
module BarGraphs

export declined_genes_graph
export improved_genes_graph

using DataAxesFormats
using Metacells
using SomeGraphs
using TanayLabUtilities

# Needed because of JET:
import Metacells.Contracts.base_block_axis
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_regulator_per_gene

"""
How much of a difference in correlation is a difference at all. A gene whose correlation in a base block moved by less
than this between the two repositories is neither improved nor declined there. This is the threshold the pipeline
reports its own counts by, so it is a constant rather than a parameter: a graph drawn by a different threshold is not
comparable to those.
"""
MIN_CHANGED_CORRELATION_FOR_GRAPHS::Float64 = 0.05

"""
    improved_genes_graph(;
        daf::DafReader,
        base_daf::DafReader,
        genes_count::Integer = $(DEFAULT.genes_count),
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells improved in the most of the base
neighborhoods, as a butterfly of the percentage of neighborhoods each declined in (to the left) and improved in (to the
right).

A gene is improved in a base block when its correlation between the cells of that block's neighborhood and their
punctuated metacells is at least $(MIN_CHANGED_CORRELATION_FOR_GRAPHS) above what it is in the `base_daf`, and declined
when it is at least that far below. Only the base blocks the gene is correlated in at all - the ones its correlation in
the `base_daf` is not zero in - are counted, so a gene which is a marker in a few neighborhoods is judged on those
rather than diluted by the ones it says nothing about. A gene correlated in no base block is not shown.

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
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
    ],
) Contract(;
    name = "base_daf",
    link = Metacells,
    # The gene masks are read from the base, as they are for every computation comparing two repositories.
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            RequiredInput,
        ),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
    ],
) function improved_genes_graph(; daf::DafReader, base_daf::DafReader, genes_count::Integer = 30)::SeriesBarsGraph
    return changed_genes_graph(; daf, base_daf, genes_count, by_improved = true)
end

"""
    declined_genes_graph(;
        daf::DafReader,
        base_daf::DafReader,
        genes_count::Integer = $(DEFAULT.genes_count),
    )::SeriesBarsGraph

The `genes_count` genes whose correlation with their cells the metacells declined in the most of the base
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
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
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
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
    ],
) function declined_genes_graph(; daf::DafReader, base_daf::DafReader, genes_count::Integer = 30)::SeriesBarsGraph
    return changed_genes_graph(; daf, base_daf, genes_count, by_improved = false)
end

# The genes whose correlation changed the most, by either of the two directions. The graphs differ only in which
# direction picks them, and both show both, so they are one function.
function changed_genes_graph(;
    daf::DafReader,
    base_daf::DafReader,
    genes_count::Integer,
    by_improved::Bool,
)::SeriesBarsGraph
    @assert genes_count > 0
    @assert axis_vector(base_daf, "gene") == axis_vector(daf, "gene")
    @assert axis_vector(base_daf, "block") == axis_vector(daf, "base_block")

    improved_fraction_per_gene, declined_fraction_per_gene = changed_fractions_per_gene(daf, base_daf)

    # Ascending, because the bars of a horizontal graph run from the bottom up, and the gene the graph is named for
    # belongs at its top.
    gene_indices =
        shown_gene_indices(by_improved ? improved_fraction_per_gene : declined_fraction_per_gene, genes_count)

    is_lateral_per_gene = get_vector(base_daf, "gene", "is_lateral").array
    is_regulator_per_gene = get_vector(base_daf, "gene", "is_regulator").array

    graph = series_bars_graph(;
        bar_axis_title = "Genes",
        value_axis_title = "Base neighborhoods",
        series_bars_values = [declined_fraction_per_gene[gene_indices], improved_fraction_per_gene[gene_indices]],
        bars_names = axis_vector(base_daf, "gene")[gene_indices],
        series_names = ["declined", "improved"],
        series_colors = ["darkred", "darkblue"],
        configuration = SeriesBarsGraphConfiguration(;
            values_orientation = HorizontalValues,
            value_axis = AxisConfiguration(; percent = true),
            mirrored = true,
        ),
    )

    # The annotations are not part of what the constructor takes, so they are attached to the data.
    graph.data.bars_annotations = [
        mask_annotation("is lateral", is_lateral_per_gene[gene_indices]),
        mask_annotation("is regulator", is_regulator_per_gene[gene_indices]),
    ]

    return graph
end

# The fraction of the base blocks each gene is correlated in that it improved in, and the fraction it declined in. A
# gene correlated in no base block is `NaN` in both, which keeps it out of the genes shown.
function changed_fractions_per_gene(daf::DafReader, base_daf::DafReader)::Tuple{Vector{Float32}, Vector{Float32}}
    base_correlation_per_gene_per_base_block = get_matrix(
        base_daf,
        "gene",
        "base_block",
        "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
    ).array
    correlation_per_gene_per_base_block =
        get_matrix(daf, "gene", "base_block", "correlation_between_base_neighborhood_cells_and_punctuated_metacells").array

    # A gene which was not correlated in a base block has a zero there, which is what the pipeline takes as "this gene
    # says nothing here" - so it is neither improved nor declined there, and does not count towards either fraction.
    is_correlated_per_gene_per_base_block = base_correlation_per_gene_per_base_block .!= 0
    n_correlated_per_gene = vec(sum(is_correlated_per_gene_per_base_block; dims = 2))

    difference_per_gene_per_base_block = correlation_per_gene_per_base_block .- base_correlation_per_gene_per_base_block
    n_improved_per_gene = vec(
        sum(
            is_correlated_per_gene_per_base_block .&
            (difference_per_gene_per_base_block .>= MIN_CHANGED_CORRELATION_FOR_GRAPHS);
            dims = 2,
        ),
    )
    n_declined_per_gene = vec(
        sum(
            is_correlated_per_gene_per_base_block .&
            (difference_per_gene_per_base_block .<= -MIN_CHANGED_CORRELATION_FOR_GRAPHS);
            dims = 2,
        ),
    )

    return (
        Float32.(n_improved_per_gene ./ n_correlated_per_gene),
        Float32.(n_declined_per_gene ./ n_correlated_per_gene),
    )
end

# The genes with the highest fraction, worst first, skipping the genes which have no fraction at all.
function shown_gene_indices(fraction_per_gene::AbstractVector{<:AbstractFloat}, genes_count::Integer)::Vector{Int}
    gene_indices = findall(.!isnan.(fraction_per_gene))
    gene_indices = gene_indices[sortperm(fraction_per_gene[gene_indices])]
    return gene_indices[max(1, length(gene_indices) - genes_count + 1):end]
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
