"""
Data sources for the graphs of a metacells repository.

A data source is a function extracting data from `Daf` metacells repositories and inserting it into a `SomeGraphs`
graph. This uses the data source views (e.g. `VectorDataFields`, `MatrixDataFields`) provided for this purpose. These
allow using the data for any graph role (coordinates, bar sizes, colors, annotations, etc.).

A data source view has two halves: the `data` (what to show) and the `configuration` (how to show it). Each half is
written the same way, by a `get_` which fetches the data from `daf` and a `put_` which writes the data into a `graph`,
with a `fill_` function that composes them:

| function                               | takes                        | does                                           |
|:-------------------------------------- |:---------------------------- |:---------------------------------------------- |
| `get_<something>_vector`, `..._matrix` | `daf` and what to query      | returns the data                               |
| `put_<something>_data!`                | sinks, the data              | writes the data fields to the `graph`          |
| `fill_<something>_data!`               | sinks, `daf`, what to query  | `get_` then `put_`                             |
| `get_<something>_colors`, ...          | `daf` and what to query      | returns what the configuration needs           |
| `put_<something>_configuration!`       | sinks, how to show it        | writes the configuration fields to the `graph` |
| `fill_<something>_configuration!`      | sinks, `daf`, how to show it | `get_` (if needed) then `put_`                 |
| `fill_<something>!`                    | sinks, `daf`, everything     | composes whichever of the above exist          |

In most cases there's no need for getting data from `daf` for the configuration. In this case there's only a put
function; there's no get and no need for a fill function.

Generic functions (e.g., `axis_names`, `axis_vector`, `axes_matrix`) only deal with data so have no configuration
variants and no overall fill wrapper. Specific data (e.g., `gene_expression`, `type`) use these under the hood
but add specific configuration for how to display the data in the graph.

A `get_` and a `fill_` are named after where the data comes from, so a generic one says which `daf` axes it queries
(`axis_vector` for one, `axes_matrix` for two). A `put_` is named after where the data goes, so a generic one says
which shape of sink it writes:

| function                 | writes                                |
|:------------------------ |:------------------------------------- |
| `put_vector_data!`       | values and a hover, into vector sinks |
| `put_matrix_data!`       | values and a hover, into matrix sinks |
| `put_vector_names_data!` | names, into vector sinks              |
| `put_matrix_names_data!` | names of both axes, into matrix sinks |

A specific `put_` is named after its source instead, since there is only one shape it can take. Nothing about a
generic `put_` depends on `daf`, so the same one serves a `DataFrame` or any other source of a vector.

The fill and put functions take [`Sinks`](@extref SomeGraphs SomeGraphs.Sources.Sinks): a single graph struct, or a
tuple or vector of them. This allows easily reusing the same data in multiple places in the graph (e.g., values,
labels, hovers). A data function writes only the sinks which are a
[`DataSink`](@extref SomeGraphs SomeGraphs.Sources.DataSink), and a configuration function only those which are a
[`ConfigurationSink`](@extref SomeGraphs SomeGraphs.Sources.ConfigurationSink), so a mixed collection is fine and
either may match nothing at all.

In general the data and the configuration functions take different sets of parameters (e.g., axis for data and gene
regularization for configuration of log scale). The `title` is a notable exception that gets passed to both sides, as it
is used for computing hover values (data) and as an axis/legend title (configuration).
"""
module DataSources

export EMPTY_TYPE_COLOR
export fill_axes_matrix_data!
export fill_axes_names_data!
export fill_axis_names_data!
export fill_axis_vector_data!
export fill_block!
export fill_boolean_annotation!
export fill_column_boolean_annotation!
export fill_column_names_data!
export fill_column_vector_data!
export fill_gene_correlation!
export fill_gene_correlation_change!
export fill_gene_expression!
export fill_genes_expression_matrix!
export fill_genes_fold_matrix!
export fill_global_flow_order!
export fill_mean_cells_per_metacell!
export fill_mean_total_UMIs_per_cell!
export fill_mean_total_UMIs_per_metacell!
export fill_module_regulators_hovers!
export fill_n_cells!
export fill_n_metacells!
export fill_total_UMIs!
export fill_type!
export fill_umap!
export GENE_FRACTION_REGULARIZATION_FOR_GRAPHS
export get_axes_matrix
export get_axis_entries_vector
export get_axis_vector
export get_block_vector
export get_boolean_annotation_vector
export get_column_vector
export get_gene_correlation_change_vector
export get_gene_correlation_vector
export get_gene_expression_vector
export get_genes_expression_matrix
export get_genes_fold_matrix
export get_global_flow_order_vector
export get_matrix_query
export get_mean_cells_per_metacell_vector
export get_mean_total_UMIs_per_cell_vector
export get_mean_total_UMIs_per_metacell_vector
export get_module_regulators_hovers
export get_n_cells_vector
export get_n_metacells_vector
export get_skeleton_gene_indices
export get_top_marker_gene_indices
export get_total_UMIs_vector
export get_type_colors
export get_type_vector
export get_umap_vector
export get_vector_query
export MAX_FOLD_FOR_GRAPHS
export put_boolean_annotation_configuration!
export put_count_configuration!
export put_gene_correlation_change_configuration!
export put_genes_expression_configuration!
export put_genes_fold_configuration!
export put_matrix_data!
export put_matrix_names_data!
export put_type_configuration!
export put_umap_configuration!
export put_umap_data!
export put_vector_data!
export put_vector_mask_data!
export put_vector_names_data!

using DataAxesFormats
using DataFrames
using SomeGraphs
using Statistics
using TanayLabUtilities

"""
The regularization added to a gene fraction before taking its log.
"""
GENE_FRACTION_REGULARIZATION_FOR_GRAPHS::Float64 = 1e-5

"""
The extreme of the color scale of the log-base-2 fold factors between gene expression and median gene expression.
"""
MAX_FOLD_FOR_GRAPHS::Float64 = 3

"""
The color to give to entities without any type annotation (empty string type). This should be different from the color
given to any explicit type. The `"magenta"` color was chosen as it isn't often used and visibly stands out.
"""
EMPTY_TYPE_COLOR::AbstractString = "magenta"

"""
    get_total_UMIs_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
    )::AbstractVector{<:Real}

Get the total UMIs of each of the `entries` of the `daf` `axis`. By default this is the direct property of the axis
entries. However, by specifying `via`, it is possible to access indirect counts (e.g., the total UMIs of the block of
the metacell of each cell).

An entry whose `via` chain breaks gets the `empty_value`, which is a count like any other, so the log scale has
something to work with.
"""
function get_total_UMIs_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
)::AbstractVector{<:Real}
    return get_axis_vector(daf; axis, query_suffix = ": total_UMIs", via, empty_value, entries)
end

"""
    put_count_configuration!(sinks::Sinks; title::Maybe{AbstractString} = nothing)::Nothing

Show a count on a log-base-2 scale, in the `YlOrRd` color scale, named by the `title`. This is how every count-like
quantity is shown: how many cells, how many metacells, how many UMIs, and the means of these. They all span orders of
magnitude, which is what the log scale is for.

The color scale runs from yellow through orange to red, so the low end is still visible; one which starts at white
would lose it against the background.

There is no regularization, because anything real has at least one of whatever is being counted.
"""
function put_count_configuration!(
    sinks::Sinks;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_count_configuration!(sink; title, show_legend)
    end
    return nothing
end

function put_count_configuration!(
    axis::AxisConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,  # NOLINT
)::Nothing
    put_count_configuration!(axis.scale)
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_count_configuration!(
    colors::ColorsConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    colors.palette = "YlOrRd"
    colors.show_legend = show_legend
    put_count_configuration!(colors.scale)
    if title !== nothing
        colors.title = title
    end
    return nothing
end

function put_count_configuration!(
    scale::ScaleConfiguration;
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = true,  # NOLINT
)::Nothing
    scale.log_base = Log2Base
    return nothing
end

"""
    fill_total_UMIs!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
        title::Maybe{AbstractString} = "total UMIs",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the total UMIs of each of the `entries` of the `daf` `axis`, shown in log base 2 in the `YlOrRd`
color scale, and name the entities after the `entries`.
"""
function fill_total_UMIs!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
    title::Maybe{AbstractString} = "total UMIs",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_total_UMIs_vector(daf; axis, entries, via, empty_value); title)  # NOJET
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_n_cells_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
    )::AbstractVector{<:Real}

Get the number of cells of each of the `entries` of the `daf` `axis`. Both metacells and blocks have one.
"""
function get_n_cells_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
)::AbstractVector{<:Real}
    return get_axis_vector(daf; axis, query_suffix = ": n_cells", via, empty_value, entries)
end

"""
    fill_n_cells!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
        title::Maybe{AbstractString} = "cells",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the number of cells of each of the `entries` of the `daf` `axis`, shown as a count, and name the
entities after the `entries`.
"""
function fill_n_cells!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
    title::Maybe{AbstractString} = "cells",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_n_cells_vector(daf; axis, entries, via, empty_value); title)
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_n_metacells_vector(
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
    )::AbstractVector{<:Real}

Get the number of metacells of each of the `entries` of the `daf` `axis`. Only blocks have one, so unlike most sources
the default `axis` here is the block.
"""
function get_n_metacells_vector(
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
)::AbstractVector{<:Real}
    return get_axis_vector(daf; axis, query_suffix = ": n_metacells", via, empty_value, entries)
end

"""
    fill_n_metacells!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
        title::Maybe{AbstractString} = "metacells",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the number of metacells of each of the `entries` of the `daf` `axis`, shown as a count, and name
the entities after the `entries`.
"""
function fill_n_metacells!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
    title::Maybe{AbstractString} = "metacells",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_n_metacells_vector(daf; axis, entries, via, empty_value); title)
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_mean_cells_per_metacell_vector(
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the mean number of cells per metacell of each of the `entries` of the `daf` `axis`, which is its number of cells
divided by its number of metacells.
"""
function get_mean_cells_per_metacell_vector(
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_n_cells_vector(daf; axis, entries, via) ./ get_n_metacells_vector(daf; axis, entries, via)
end

"""
    fill_mean_cells_per_metacell!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean cells per metacell",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the mean number of cells per metacell of each of the `entries` of the `daf` `axis`, shown as a
count, and name the entities after the `entries`.
"""
function fill_mean_cells_per_metacell!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean cells per metacell",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_mean_cells_per_metacell_vector(daf; axis, entries, via); title)
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_mean_total_UMIs_per_metacell_vector(
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the mean total UMIs per metacell of each of the `entries` of the `daf` `axis`, which is its total UMIs divided by
its number of metacells.
"""
function get_mean_total_UMIs_per_metacell_vector(
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_total_UMIs_vector(daf; axis, entries, via) ./ get_n_metacells_vector(daf; axis, entries, via)
end

"""
    fill_mean_total_UMIs_per_metacell!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean UMIs per metacell",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the mean total UMIs per metacell of each of the `entries` of the `daf` `axis`, shown as a count,
and name the entities after the `entries`.
"""
function fill_mean_total_UMIs_per_metacell!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean UMIs per metacell",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_mean_total_UMIs_per_metacell_vector(daf; axis, entries, via); title)
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_mean_total_UMIs_per_cell_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the mean total UMIs per cell of each of the `entries` of the `daf` `axis`, which is its total UMIs divided by its
number of cells. Both metacells and blocks have both.
"""
function get_mean_total_UMIs_per_cell_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_total_UMIs_vector(daf; axis, entries, via) ./ get_n_cells_vector(daf; axis, entries, via)
end

"""
    fill_mean_total_UMIs_per_cell!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean UMIs per cell",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the mean total UMIs per cell of each of the `entries` of the `daf` `axis`, shown as a count, and
name the entities after the `entries`.
"""
function fill_mean_total_UMIs_per_cell!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean UMIs per cell",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_mean_total_UMIs_per_cell_vector(daf; axis, entries, via); title)
    put_count_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_top_marker_gene_indices(daf::DafReader; markers_count::Integer)::Vector{Int}

Get the indices of the `markers_count` best marker genes of a `daf` repository, to pass as the `genes` of another data
source.

A rank of 1 is the gene which best distinguishes between the cell states, and a non-marker gene is ranked above any
count that can be asked for, so this is simply the genes ranked that high; a repository with fewer markers than that
gives all of them.
"""
function get_top_marker_gene_indices(daf::DafReader; markers_count::Integer)::Vector{Int}
    @assert markers_count > 0
    return findall(get_axis_vector(daf; axis = "gene", query_suffix = ": marker_rank") .<= markers_count)
end

"""
    get_skeleton_gene_indices(daf::DafReader)::Vector{Int}

Get the indices of the skeleton genes of a `daf` repository, to pass as the `genes` of another data source.

Unlike the markers these are not ranked, and there are few enough of them to show them all. They are the genes the
blocks were computed from, so they show what the blocks were told apart by.
"""
function get_skeleton_gene_indices(daf::DafReader)::Vector{Int}
    return findall(get_axis_vector(daf; axis = "gene", query_suffix = ": is_skeleton"))
end

"""
    get_block_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::AbstractString = "",
    )::AbstractVector{<:AbstractString}

Get the block of each of the `entries` of the `daf` `axis`. By default this is the direct property of the axis
entries. However, by specifying `via`, it is possible to access indirect blocks (e.g., the block of the metacell of
each cell). An entry whose chain breaks gets the `empty_value`, which names no block.
"""
function get_block_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::AbstractString = "",
)::AbstractVector{<:AbstractString}
    return get_axis_vector(daf; axis, query_suffix = ": block", via, empty_value, entries)
end

"""
    fill_block!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::AbstractString = "",
        title::Maybe{AbstractString} = "block",
    )::Nothing

Fill the `sinks` with the block of each of the `entries` of the `daf` `axis`, and name the entities after the
`entries`. Blocks are shown as a grouping of an axis, as a hover, or both; grouping by them puts a gap between the
blocks, and the hover says which block an entry is in.

Aim this at a whole data source view to get the values and a hover line, or at its `values` alone to get only the
values.
"""
function fill_block!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::AbstractString = "",
    title::Maybe{AbstractString} = "block",
)::Nothing
    put_vector_data!(sinks, get_block_vector(daf; axis, entries, via, empty_value); title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_global_flow_order_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        type_property::AbstractString = "type",
        empty_value::Real = 0,
    )::AbstractVector{<:Real}

Get the global flow order of the type of each of the `entries` of the `daf` `axis`. An entry with no type gets the
`empty_value`, so it forms a group of its own at that end of the order.

The order is a property of the type, so all the entries of a type get the same number. Grouping an axis by it
therefore groups the entries by their type, and lays the groups out in the flow order, so that every graph shows the
types in the same order.
"""
function get_global_flow_order_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    type_property::AbstractString = "type",
    empty_value::Real = 0,
)::AbstractVector{<:Real}
    return get_axis_vector(
        daf;
        axis,
        query_suffix = ": global_flow_order",
        via = (type_property,),
        empty_value,
        entries,
    )
end

"""
    fill_global_flow_order!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        type_property::AbstractString = "type",
        empty_value::Real = 0,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the global flow order of the type of each of the `entries` of the `daf` `axis`, and name the
entities after the `entries`.
"""
function fill_global_flow_order!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    type_property::AbstractString = "type",
    empty_value::Real = 0,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    put_vector_data!(sinks, get_global_flow_order_vector(daf; axis, entries, type_property, empty_value); title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_boolean_annotation_vector(
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractString}

Get the Boolean mask `property` of each of the `entries` of the `daf` `axis`, as the *strings* `"true"` and `"false"`.
Coloring by a Boolean is coloring by two categories, and a categorical palette is keyed by strings.
"""
function get_boolean_annotation_vector(
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractString}
    return string.(get_axis_vector(daf; axis, query_suffix = ": $(escape_value(property))", entries))
end

"""
    put_boolean_annotation_configuration!(
        sinks::Sinks;
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = false,
    )::Nothing

Show Boolean *string* annotations: `true` in black and `false` in light grey, named by the `title`. The values
themselves are put in by [`put_vector_data!`](@ref).

Black against light grey has only two values and is named by its title, so `show_legend` is off.
"""
function put_boolean_annotation_configuration!(
    sinks::Sinks;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = false,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_boolean_annotation_configuration!(sink; title, show_legend)
    end
    return nothing
end

function put_boolean_annotation_configuration!(
    colors::ColorsConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = false,
)::Nothing
    colors.palette = Dict("true" => "black", "false" => "lightgrey")
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

"""
    fill_boolean_annotation!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = property,
        show_legend::Bool = false,
    )::Nothing

Fill the `sinks` with the Boolean mask `property` of each of the `entries` of the `daf` `axis`, shown as `true` in
black and `false` in light grey, and name the entities after the `entries`.
"""
function fill_boolean_annotation!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = property,
    show_legend::Bool = false,
)::Nothing
    put_vector_data!(sinks, get_boolean_annotation_vector(daf; axis, property, entries); title)
    put_boolean_annotation_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    fill_gene_expression!(
        sinks::Sinks,
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
        title::Maybe{AbstractString} = "\$(gene) fraction",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the `gene` expression level (linear fraction) per each of the `entries` of some `daf` `axis`,
shown in log base 2 using the `gene_fraction_regularization`, and name the entities after the `entries`.
"""
function fill_gene_expression!(
    sinks::Sinks,
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "$(gene) fraction",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_gene_expression_vector(daf; gene, axis, entries); title)
    put_genes_expression_configuration!(sinks; gene_fraction_regularization, title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_gene_expression_vector(
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the vector of the `gene` expression level (linear fraction) per each of the `entries` of some `daf` `axis`.
"""
function get_gene_expression_vector(
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_axis_vector(daf; axis, query_suffix = ":: linear_fraction @ gene = $(escape_value(gene))", entries)
end

"""
    get_gene_correlation_vector(
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "base_block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the correlation of the `gene` between the cells of the neighborhood of each of the `entries` of some `daf` `axis`
and their punctuated metacells. A zero is the gene saying nothing about that entry, rather than a correlation which
happens to be zero.
"""
function get_gene_correlation_vector(
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_axis_vector(
        daf;
        axis,
        query_suffix = ":: correlation_between_base_neighborhood_cells_and_punctuated_metacells " *
                       "@ gene = $(escape_value(gene))",
        entries,
    )
end

"""
    fill_gene_correlation!(
        sinks::Sinks,
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "base_block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = "correlation",
    )::Nothing

Fill the `sinks` with the correlation of the `gene` per each of the `entries` of some `daf` `axis`, and name the
entities after the `entries`.
"""
function fill_gene_correlation!(
    sinks::Sinks,
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "correlation",
)::Nothing
    put_vector_data!(sinks, get_gene_correlation_vector(daf; gene, axis, entries); title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_gene_correlation_change_vector(
        daf::DafReader,
        base_daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "base_block",
        base_axis::AbstractString = "base_block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get how much the correlation of the `gene` changed between the `base_daf` and the `daf`, per each of the `entries`. A
value above zero is an entry whose metacells describe the gene better than the base does.

The same entries are the `axis` of the `daf` and the `base_axis` of the `base_daf`. They must therefore list the same
names in the same order.
"""
function get_gene_correlation_change_vector(
    daf::DafReader,
    base_daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    base_axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    @assert axis_vector(daf, axis) == axis_vector(base_daf, base_axis)
    return get_gene_correlation_vector(daf; gene, axis, entries) .-
           get_gene_correlation_vector(base_daf; gene, axis = base_axis, entries)
end

"""
    put_gene_correlation_change_configuration!(
        sinks::Sinks;
        title::Maybe{AbstractString} = "correlation change",
        show_legend::Bool = true,
    )::Nothing

Name the correlation change in the configuration, and include it in the legend if `show_legend` where it is the colors.
There is nothing else to say about how to show it: it is a difference of two correlations, so it is already on the
scale it is read in.
"""
function put_gene_correlation_change_configuration!(
    sinks::Sinks;
    title::Maybe{AbstractString} = "correlation change",
    show_legend::Bool = true,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_gene_correlation_change_configuration!(sink; title, show_legend)
    end
    return nothing
end

function put_gene_correlation_change_configuration!(
    axis::AxisConfiguration;
    title::Maybe{AbstractString} = "correlation change",
    show_legend::Bool = true,  # NOLINT
)::Nothing
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_gene_correlation_change_configuration!(
    colors::ColorsConfiguration;
    title::Maybe{AbstractString} = "correlation change",
    show_legend::Bool = true,
)::Nothing
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

"""
    fill_gene_correlation_change!(
        sinks::Sinks,
        daf::DafReader,
        base_daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "base_block",
        base_axis::AbstractString = "base_block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = "correlation change",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with how much the correlation of the `gene` changed between the `base_daf` and the `daf`, per each of
the `entries`, and name the entities after them.
"""
function fill_gene_correlation_change!(
    sinks::Sinks,
    daf::DafReader,
    base_daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    base_axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "correlation change",
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, get_gene_correlation_change_vector(daf, base_daf; gene, axis, base_axis, entries); title)
    put_gene_correlation_change_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    fill_genes_expression_matrix!(
        sinks::Sinks,
        daf::DafReader;
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
        title::Maybe{AbstractString} = "fraction",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the `genes` expression level (linear fraction) per each of the `entries` of some `daf` `axis`,
shown in log base 2 using the `gene_fraction_regularization`, and name the rows and the columns. The `genes` are the
rows and the `entries` are the columns.
"""
function fill_genes_expression_matrix!(
    sinks::Sinks,
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
)::Nothing
    put_matrix_data!(sinks, get_genes_expression_matrix(daf; genes, axis, entries); title)
    put_genes_expression_configuration!(sinks; gene_fraction_regularization, title, show_legend)
    fill_axes_names_data!(
        sinks,
        daf;
        rows_axis = "gene",
        columns_axis = axis,
        row_entries = genes,
        column_entries = entries,
    )
    return nothing
end

"""
    get_genes_expression_matrix(
        daf::DafReader;
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractMatrix{<:AbstractFloat}

Get the matrix of the `genes` expression level (linear fraction) per each of the `entries` of some `daf` `axis`. The
`genes` are the rows and the `entries` are the columns.
"""
function get_genes_expression_matrix(
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractMatrix{<:AbstractFloat}
    return get_axes_matrix(
        daf;
        rows_axis = "gene",
        columns_axis = axis,
        query_suffix = ":: linear_fraction",
        row_entries = genes,
        column_entries = entries,
    )
end

"""
    put_genes_expression_configuration!(
        sinks::Sinks;
        gene_fraction_regularization::Real = $(GENE_FRACTION_REGULARIZATION_FOR_GRAPHS),
        title::Maybe{AbstractString} = "fraction",
        show_legend::Bool = true,
    )::Nothing

Show gene expression (linear fraction) on a log-base-2 scale, regularized by `gene_fraction_regularization`, named by
the `title`. The fractions themselves are put in by [`put_vector_data!`](@ref) or
[`put_matrix_data!`](@ref).

The `show_legend` applies where the fractions are the colors. An axis is read off its own ticks and ignores it.
"""
function put_genes_expression_configuration!(
    sinks::Sinks;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
)::Nothing
    @assert gene_fraction_regularization >= 0
    visit_configuration_sinks(sinks) do sink
        return put_genes_expression_configuration!(sink; gene_fraction_regularization, title, show_legend)
    end
    return nothing
end

function put_genes_expression_configuration!(
    axis::AxisConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,  # NOLINT
)::Nothing
    put_genes_expression_configuration!(axis.scale; gene_fraction_regularization)
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_genes_expression_configuration!(
    colors::ColorsConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
)::Nothing
    put_genes_expression_configuration!(colors.scale; gene_fraction_regularization)
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

function put_genes_expression_configuration!(
    scale::ScaleConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = true,  # NOLINT
)::Nothing
    scale.log_base = Log2Base
    scale.log_regularization = gene_fraction_regularization
    return nothing
end

"""
    put_umap_data!(
        sinks::Sinks,
        coordinate_per_entry::AbstractVector{<:AbstractFloat};
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Put a UMAP `coordinate_per_entry` into the `sinks`, as the values of a role.
"""
function put_umap_data!(
    sinks::Sinks,
    coordinate_per_entry::AbstractVector{<:AbstractFloat};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    visit_data_sinks(sinks) do sink
        return put_umap_data!(sink, coordinate_per_entry; title)
    end
    return nothing
end

function put_umap_data!(
    values::VectorValuesData,
    coordinate_per_entry::AbstractVector{<:AbstractFloat};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    values.vector = coordinate_per_entry
    return nothing
end

# A UMAP coordinate is the arbitrary output of the projection, so there is nothing worth saying about one in a hover.
function put_umap_data!(
    ::VectorEntitiesData,
    ::AbstractVector{<:AbstractFloat};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    return nothing
end

"""
    fill_umap!(
        sinks::Sinks,
        daf::DafReader;
        coordinate::AbstractString,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = "UMAP \$(uppercase(coordinate))",
    )::Nothing

Fill the `sinks` with the UMAP `coordinate` per each of the `entries` of some `daf` `axis`, shown without ticks or a
grid, and name the entities after the `entries`.
"""
function fill_umap!(
    sinks::Sinks,
    daf::DafReader;
    coordinate::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "UMAP $(uppercase(coordinate))",
)::Nothing
    put_umap_data!(sinks, get_umap_vector(daf; coordinate, axis, entries); title)
    put_umap_configuration!(sinks; title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_umap_vector(
        daf::DafReader;
        coordinate::AbstractString,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractFloat}

Get the vector of the UMAP `coordinate` per each of the `entries` of some `daf` `axis`.
"""
function get_umap_vector(
    daf::DafReader;
    coordinate::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractFloat}
    return get_axis_vector(daf; axis, query_suffix = ": $(escape_value("umap_$(coordinate)"))", entries)
end

"""
    put_umap_configuration!(sinks::Sinks; title::Maybe{AbstractString} = nothing)::Nothing

Show UMAP coordinates without ticks or a grid, named by the `title`. Only which points are near which other points
means anything, so the values themselves are not worth labelling.
"""
function put_umap_configuration!(sinks::Sinks; title::Maybe{AbstractString} = nothing)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_umap_configuration!(sink; title)
    end
    return nothing
end

function put_umap_configuration!(axis::AxisConfiguration; title::Maybe{AbstractString} = nothing)::Nothing
    if title !== nothing
        axis.title = title
    end
    axis.show_ticks = false
    axis.show_grid = false
    return nothing
end

"""
    fill_genes_fold_matrix!(
        sinks::Sinks,
        daf::DafReader;
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        max_fold::Real = MAX_FOLD_FOR_GRAPHS,
        title::Maybe{AbstractString} = "fold from median",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the `genes` fold (log base-2 minus the median) per each of the `entries` of some `daf` `axis`,
shown in the range -`max_fold` (blue) to 0 (white) to +`max_fold` (red), and name the rows and the columns. The `genes`
are the rows and the `entries` are the columns. This is based on the `log_linear_fraction` per gene per `axis`, which is
already a log, so there's no regularization to apply here.
"""
function fill_genes_fold_matrix!(
    sinks::Sinks,
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fold from median",
    show_legend::Bool = true,
)::Nothing
    put_matrix_data!(sinks, get_genes_fold_matrix(daf; genes, axis, entries); title)  # NOJET
    put_genes_fold_configuration!(sinks; max_fold, title, show_legend)
    fill_axes_names_data!(
        sinks,
        daf;
        rows_axis = "gene",
        columns_axis = axis,
        row_entries = genes,
        column_entries = entries,
    )
    return nothing
end

"""
    get_genes_fold_matrix(
        daf::DafReader;
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractMatrix{<:AbstractFloat}

Get the matrix of the `genes` fold (log base-2 minus the median) per each of the `entries` of some `daf` `axis`. The
`genes` are the rows and the `entries` are the columns. This is based on the `log_linear_fraction` per gene per `axis`.
The median of each gene is taken across the `entries`.
"""
function get_genes_fold_matrix(
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractMatrix{<:AbstractFloat}
    log_linear_fraction_per_gene_per_entry = get_axes_matrix(
        daf;
        rows_axis = "gene",
        columns_axis = axis,
        query_suffix = ":: log_linear_fraction",
        row_entries = genes,
        column_entries = entries,
    )
    return log_linear_fraction_per_gene_per_entry .- median(log_linear_fraction_per_gene_per_entry; dims = 2)
end

"""
    put_genes_fold_configuration!(
        sinks::Sinks;
        max_fold::Real = MAX_FOLD_FOR_GRAPHS,
        title::Maybe{AbstractString} = "fold from median",
        show_legend::Bool = true,
    )::Nothing

Show the genes fold (log base-2 minus the median) in the range -`max_fold` (blue) to 0 (white) to +`max_fold` (red),
named by the `title`, and include it in the legend if `show_legend`. The folds themselves are put in by
[`put_matrix_data!`](@ref).
"""
function put_genes_fold_configuration!(
    sinks::Sinks;
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fold from median",
    show_legend::Bool = true,
)::Nothing
    @assert max_fold > 0
    visit_configuration_sinks(sinks) do sink
        return put_genes_fold_configuration!(sink; max_fold, title, show_legend)
    end
    return nothing
end

function put_genes_fold_configuration!(
    colors::ColorsConfiguration;
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fold from median",
    show_legend::Bool = true,
)::Nothing
    colors.palette = "BuWtRd"
    colors.show_legend = show_legend
    colors.scale.minimum = -max_fold
    colors.scale.maximum = max_fold
    if title !== nothing
        colors.title = title
    end
    return nothing
end

"""
    fill_axis_names_data!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::Nothing

Name the entities of the `sinks` after each of the `entries` of the `daf` `axis`.

This fetches the names with [`get_axis_entries_vector`](@ref), then hands them to [`put_vector_names_data!`](@ref).
"""
function fill_axis_names_data!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::Nothing
    put_vector_names_data!(sinks, get_axis_entries_vector(daf; axis, entries))
    return nothing
end

"""
    put_vector_names_data!(sinks::Sinks, name_per_entry::AbstractVector{<:AbstractString})::Nothing

Name the entities of the `sinks` after the `name_per_entry`. Where the graph has room to label the entities, these
become their tick labels. They are also the first line of the hover of each entity.
"""
function put_vector_names_data!(sinks::Sinks, name_per_entry::AbstractVector{<:AbstractString})::Nothing
    visit_data_sinks(sinks) do sink
        return put_vector_names_data!(sink, name_per_entry)
    end
    return nothing
end

function put_vector_names_data!(entities::VectorEntitiesData, name_per_entry::AbstractVector{<:AbstractString})::Nothing
    entities.names = name_per_entry
    return nothing
end

# A name identifies an entity, so it belongs to the entities rather than to any one role's values.
function put_vector_names_data!(::VectorValuesData, ::AbstractVector{<:AbstractString})::Nothing
    return nothing
end

"""
    put_vector_mask_data!(
        sinks::Sinks,
        is_shown_per_entry::Union{AbstractVector{Bool}, BitVector},
    )::Nothing

Hide the entities of the `sinks` which are not shown by the `is_shown_per_entry` mask. Hidden entities are still part
of the data, so they take part in whatever is computed from it, unless the relevant configuration says otherwise.
"""
function put_vector_mask_data!(sinks::Sinks, is_shown_per_entry::Union{AbstractVector{Bool}, BitVector})::Nothing
    visit_data_sinks(sinks) do sink
        return put_vector_mask_data!(sink, is_shown_per_entry)
    end
    return nothing
end

function put_vector_mask_data!(
    entities::VectorEntitiesData,
    is_shown_per_entry::Union{AbstractVector{Bool}, BitVector},
)::Nothing
    entities.mask = is_shown_per_entry
    return nothing
end

# A mask hides an entity, so it belongs to the entities rather than to any one role's values.
function put_vector_mask_data!(::VectorValuesData, ::Union{AbstractVector{Bool}, BitVector})::Nothing
    return nothing
end

"""
    get_axis_vector(
        daf::DafReader;
        axis::AbstractString,
        query_suffix::AbstractString,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Maybe{StorageScalarBase} = nothing,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:StorageScalarBase}

Get the result of the query `@ axis via query_suffix` for each of the `entries` of the `daf` `axis`.

The `via` is a chain of axes to walk before asking anything, so that a property of something an entry belongs to can be
asked for (e.g., the block of the metacell of each cell). An `empty_value` is what an entry whose chain breaks gets; a
`nothing` leaves the query without a default, so such an entry is an error.

The `query_suffix` is whatever a vector query says at the end of that walk, e.g. `": block"` for a property or
`":: linear_fraction @ gene = FOXA1"` for a slice of a matrix. It is written by the caller, so it is used as it is;
everything else is escaped.
"""
function get_axis_vector(
    daf::DafReader;
    axis::AbstractString,
    query_suffix::AbstractString,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Maybe{StorageScalarBase} = nothing,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:StorageScalarBase}
    return get_vector_query(
        daf;
        query = "@ $(escape_value(axis))$(via_query(via, empty_value)) $(query_suffix)",
        indices = entries_indices(daf, axis, entries),
    )
end

# The part of a query which walks `via` a chain of axes, giving an entry whose chain breaks the `empty_value`.
function via_query(
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}},
    empty_value::Maybe{StorageScalarBase},
)::AbstractString
    if via === nothing
        return ""
    end
    return join(String[" : $(escape_value(via_axis))$(empty_query(empty_value))" for via_axis in via])
end

function empty_query(::Nothing)::AbstractString
    return ""
end

function empty_query(empty_value::AbstractString)::AbstractString
    return " ?? $(escape_value(empty_value))"
end

function empty_query(empty_value::StorageScalarBase)::AbstractString
    return " ?? $(empty_value)"
end

"""
    put_vector_data!(
        sinks::Sinks,
        value_per_entry::AbstractVector{<:StorageScalarBase};
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Put a `value_per_entry` of an axis into the `sinks`: as the values of a role, and as a hover line on the entities. Use
this when you have the data already; [`fill_axis_vector_data!`](@ref) fetches it from a `daf` repository for you.

The `title` prefixes the hover line, as `title: value`. Nothing is configured here, because how to show a value depends
on which property it is; a specific data source says that in its own `put_..._configuration!`.
"""
function put_vector_data!(
    sinks::Sinks,
    value_per_entry::AbstractVector{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    visit_data_sinks(sinks) do sink
        return put_vector_data!(sink, value_per_entry; title)
    end
    return nothing
end

function put_vector_data!(
    values::VectorValuesData,
    value_per_entry::AbstractVector{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    values.vector = value_per_entry
    return nothing
end

function put_vector_data!(
    entities::VectorEntitiesData,
    value_per_entry::AbstractVector{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    add_hovers!(entities, hover_strings(value_per_entry); title)
    return nothing
end

"""
    put_matrix_data!(
        sinks::Sinks,
        value_per_row_per_column::AbstractMatrix{<:StorageScalarBase};
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Put a `value_per_row_per_column` of two axes into the `sinks`: as the values of the entries, and as a hover line on
each entry. The matrix twin of [`put_vector_data!`](@ref).

The row and the column entities are not written here. They belong to their axes and are sized by one axis each, so they
are named separately; see [`put_matrix_names_data!`](@ref).
"""
function put_matrix_data!(
    sinks::Sinks,
    value_per_row_per_column::AbstractMatrix{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    visit_data_sinks(sinks) do sink
        return put_matrix_data!(sink, value_per_row_per_column; title)
    end
    return nothing
end

function put_matrix_data!(
    values::MatrixValuesData,
    value_per_row_per_column::AbstractMatrix{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    values.matrix = value_per_row_per_column
    return nothing
end

function put_matrix_data!(
    entities::MatrixEntitiesData,
    value_per_row_per_column::AbstractMatrix{<:StorageScalarBase};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    add_hovers!(entities, hover_strings(value_per_row_per_column); title)
    return nothing
end

"""
    put_matrix_names_data!(
        sinks::Sinks,
        name_per_row::AbstractVector{<:AbstractString},
        name_per_column::AbstractVector{<:AbstractString},
    )::Nothing

Name the rows and the columns of the `sinks` after the `name_per_row` and the `name_per_column`. The matrix twin of
[`put_vector_names_data!`](@ref).

The two axes are named together here rather than through
[`visit_data_sinks`](@extref SomeGraphs SomeGraphs.Sources.visit_data_sinks), which doesn't walk into them because
they are sized by one axis each while the entries are sized by both.
"""
function put_matrix_names_data!(
    sinks::Sinks,
    name_per_row::AbstractVector{<:AbstractString},
    name_per_column::AbstractVector{<:AbstractString},
)::Nothing
    for sink in matrix_data_sinks(sinks)
        put_vector_names_data!(sink.rows_entities, name_per_row)
        put_vector_names_data!(sink.columns_entities, name_per_column)
    end
    return nothing
end

"""
    fill_axes_names_data!(
        sinks::Sinks,
        daf::DafReader;
        rows_axis::AbstractString,
        columns_axis::AbstractString,
        row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::Nothing

Name the rows of the `sinks` after the `row_entries` of the `daf` `rows_axis`, and their columns after the
`column_entries` of its `columns_axis`.

This fetches the names with [`get_axis_entries_vector`](@ref), then hands them to [`put_matrix_names_data!`](@ref).
"""
function fill_axes_names_data!(
    sinks::Sinks,
    daf::DafReader;
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::Nothing
    put_matrix_names_data!(
        sinks,
        get_axis_entries_vector(daf; axis = rows_axis, entries = row_entries),
        get_axis_entries_vector(daf; axis = columns_axis, entries = column_entries),
    )
    return nothing
end

"""
    get_axes_matrix(
        daf::DafReader;
        rows_axis::AbstractString,
        columns_axis::AbstractString,
        query_suffix::AbstractString,
        row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractMatrix{<:StorageScalarBase}

Get the result of the query `@ rows_axis @ columns_axis query_suffix` for each of the `row_entries` and
`column_entries` of the two `daf` axes. The matrix twin of [`get_axis_vector`](@ref).
"""
function get_axes_matrix(
    daf::DafReader;
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    query_suffix::AbstractString,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractMatrix{<:StorageScalarBase}
    return get_matrix_query(
        daf;
        query = "@ $(escape_value(rows_axis)) @ $(escape_value(columns_axis)) $(query_suffix)",
        row_indices = entries_indices(daf, rows_axis, row_entries),
        column_indices = entries_indices(daf, columns_axis, column_entries),
    )
end

"""
    fill_axes_matrix_data!(
        sinks::Sinks,
        daf::DafReader;
        rows_axis::AbstractString,
        columns_axis::AbstractString,
        query_suffix::AbstractString,
        row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the result of the query `@ rows_axis @ columns_axis query_suffix`, and name the rows and the
columns. The matrix twin of [`fill_axis_vector_data!`](@ref).
"""
function fill_axes_matrix_data!(
    sinks::Sinks,
    daf::DafReader;
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    query_suffix::AbstractString,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    put_matrix_data!(
        sinks,
        get_axes_matrix(daf; rows_axis, columns_axis, query_suffix, row_entries, column_entries);
        title,
    )
    fill_axes_names_data!(sinks, daf; rows_axis, columns_axis, row_entries, column_entries)
    return nothing
end

# The `MatrixDataFields` of the `sinks`, which are the only ones which have row and column entities to name.
function matrix_data_sinks(sinks::Sinks)::Vector{MatrixDataFields}
    matrix_sinks = MatrixDataFields[]
    for sink in (sinks isa Union{Tuple, AbstractVector} ? sinks : (sinks,))
        if sink isa MatrixFields
            push!(matrix_sinks, sink.data)
        elseif sink isa MatrixDataFields
            push!(matrix_sinks, sink)
        end
    end
    return matrix_sinks
end

# The values as the strings a hover shows. Only strings can be a hover, so anything else is converted.
function hover_strings(value_per_entry::AbstractArray{<:AbstractString})::AbstractArray{<:AbstractString}
    return value_per_entry
end

function hover_strings(value_per_entry::AbstractArray{<:StorageScalarBase})::AbstractArray{<:AbstractString}
    return string.(value_per_entry)
end

"""
    fill_axis_vector_data!(
        sinks::Sinks,
        daf::DafReader;
        axis::AbstractString,
        query_suffix::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the result of the query `@ axis query_suffix` for each of the `entries` of the `daf` `axis`, and
name the entities after the `entries`.

This fetches the data with [`get_axis_vector`](@ref), then hands it to [`put_vector_data!`](@ref).
"""
function fill_axis_vector_data!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString,
    query_suffix::AbstractString,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Maybe{StorageScalarBase} = nothing,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    put_vector_data!(sinks, get_axis_vector(daf; axis, query_suffix, via, empty_value, entries); title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    get_column_vector(frame::DataFrame; column::AbstractString)::AbstractVector{<:StorageScalarBase}

Get the `column` of a `frame`, one value per row. This is the `DataFrame` twin of [`get_axis_vector`](@ref), for the
graphs whose data is a report computed by `Metacells` rather than a property of a `daf` repository.
"""
function get_column_vector(frame::DataFrame; column::AbstractString)::AbstractVector{<:StorageScalarBase}
    return frame[!, column]
end

"""
    fill_column_vector_data!(
        sinks::Sinks,
        frame::DataFrame;
        column::AbstractString,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the `column` of a `frame`, one value per row.

The entities are not named here, unlike [`fill_axis_vector_data!`](@ref); a frame has no one column which identifies
its rows, so say which one does with [`fill_column_names_data!`](@ref).
"""
function fill_column_vector_data!(
    sinks::Sinks,
    frame::DataFrame;
    column::AbstractString,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    put_vector_data!(sinks, get_column_vector(frame; column); title)
    return nothing
end

"""
    fill_column_names_data!(sinks::Sinks, frame::DataFrame; column::AbstractString)::Nothing

Name the entities of the `sinks` after the `column` of a `frame`, one name per row.
"""
function fill_column_names_data!(sinks::Sinks, frame::DataFrame; column::AbstractString)::Nothing
    put_vector_names_data!(sinks, string.(get_column_vector(frame; column)))
    return nothing
end

"""
    fill_column_boolean_annotation!(
        sinks::Sinks,
        frame::DataFrame;
        column::AbstractString,
        title::Maybe{AbstractString} = column,
        show_legend::Bool = false,
    )::Nothing

Fill the `sinks` with the Boolean `column` of a `frame`, shown as `true` in black and `false` in light grey. The
`DataFrame` twin of [`fill_boolean_annotation!`](@ref).
"""
function fill_column_boolean_annotation!(
    sinks::Sinks,
    frame::DataFrame;
    column::AbstractString,
    title::Maybe{AbstractString} = column,
    show_legend::Bool = false,
)::Nothing
    put_vector_data!(sinks, string.(get_column_vector(frame; column)); title)
    put_boolean_annotation_configuration!(sinks; title, show_legend)
    return nothing
end

"""
    get_module_regulators_hovers(
        frame::DataFrame;
        prefix::AbstractString,
        side_name::AbstractString,
        regulators_count::Integer,
    )::AbstractVector{<:AbstractString}

Get a hover per row of a gene report `frame`, saying how often the gene is in no module at all in the base blocks of
one side, and which regulators it is most often in a module with there. The `prefix` picks that side's columns (`imp`
or `deg`), and the `side_name` names it in the text.

The report pads a side which gives a gene fewer regulators than the `regulators_count` asked for with empty names,
which are the regulators there is nothing to say about, and which are left out.
"""
function get_module_regulators_hovers(
    frame::DataFrame;
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::AbstractVector{<:AbstractString}
    @assert regulators_count >= 0
    return AbstractString[
        join(module_regulators_lines(frame, row_index; prefix, side_name, regulators_count), "<br>") for
        row_index in 1:nrow(frame)
    ]
end

# The lines of one gene's hover: how often it is in no module, and then a line per regulator it is in a module with.
function module_regulators_lines(
    frame::DataFrame,
    row_index::Integer;
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::Vector{String}
    lines = String["$(side_name): in no module in $(percent_text(frame[row_index, "$(prefix)_no_mod_f"])) of the cells"]
    for rank in 1:regulators_count
        regulator_name = frame[row_index, "$(prefix)_reg$(rank)"]
        if regulator_name != ""
            push!(lines, "- $(regulator_name): $(percent_text(frame[row_index, "$(prefix)_reg$(rank)_f"]))")
        end
    end
    return lines
end

# A fraction in the units the value axis shows, so that a hover and the bars are read the same way. This is not
# `TanayLabUtilities.percent`, which rounds to whole percentages and says `<1%` and `>99%` at the extremes.
function percent_text(fraction::AbstractFloat)::String
    return "$(round(100 * fraction; digits = 1))%"
end

"""
    fill_module_regulators_hovers!(
        sinks::Sinks,
        frame::DataFrame;
        prefix::AbstractString,
        side_name::AbstractString,
        regulators_count::Integer,
    )::Nothing

Add a hover per row of a gene report `frame` to the `sinks`, saying which regulators the gene is in a module with in
the base blocks of one side. The lines carry their own labels, so nothing is prefixed to them.
"""
function fill_module_regulators_hovers!(
    sinks::Sinks,
    frame::DataFrame;
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::Nothing
    put_vector_data!(sinks, get_module_regulators_hovers(frame; prefix, side_name, regulators_count))
    return nothing
end

"""
    fill_type!(
        sinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        type_property::AbstractString = "type",
        type_axis::AbstractString = "type",
        via::Maybe{Union{Tuple{<:AbstractString}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = type_axis,
        show_legend::Bool = true,
        empty_type_color::AbstractString = EMPTY_TYPE_COLOR,
    )::Nothing

Fill the `sinks` with the `type_property` of each of the `entries` of the `daf` `axis`, and name the entities after the
`entries`. By default the type is the direct property of the axis entries. However, by specifying `via`, it is possible
to access indirect types (e.g., the type of the block of the metacell of each cell). Types are shown using the `color`
of each entry of the `type_axis`, with an additional `empty_type_color` for `entries` w/ no type (empty string).

This fetches the data with [`get_type_vector`](@ref) and [`get_type_colors`](@ref), puts the types in with
[`put_vector_data!`](@ref), and says how to show them with [`put_type_configuration!`](@ref).
"""
function fill_type!(
    sinks::Sinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    type_property::AbstractString = "type",
    type_axis::AbstractString = "type",
    via::Maybe{Union{Tuple{<:AbstractString}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = type_axis,
    show_legend::Bool = true,
    empty_type_color::Maybe{AbstractString} = EMPTY_TYPE_COLOR,
)::Nothing
    put_vector_data!(sinks, get_type_vector(daf; axis, entries, type_property, via); title)
    put_type_configuration!(sinks, get_type_colors(daf; type_axis, empty_type_color); title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

"""
    put_type_configuration!(
        sinks::Sinks,
        color_per_type::CategoricalColors;
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = true,
    )::Nothing

Show the types of the `sinks` using the `color_per_type` palette, named by the `title`. The types themselves are put in
by [`put_vector_data!`](@ref).

A color says nothing without a legend to read it by, so `show_legend` is on. Turn it off for an annotation of a graph
whose colors already name the same types.
"""
function put_type_configuration!(
    sinks::Sinks,
    color_per_type::CategoricalColors;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_type_configuration!(sink, color_per_type; title, show_legend)
    end
    return nothing
end

function put_type_configuration!(
    colors::ColorsConfiguration,
    color_per_type::CategoricalColors;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    colors.palette = color_per_type
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

"""
    get_type_vector(
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        type_property::AbstractString = "type",
        via::Maybe{Union{Tuple{<:AbstractString}, AbstractVector{<:AbstractString}}} = nothing,
    )::AbstractVector{<:AbstractString}

Get the `type_property` for each of the `entries` of the `daf` `axis`. By default this is the direct property of the
axis entries. However, by specifying `via`, it is possible to access indirect types (e.g., the type of the block of the
metacell of each cell).
"""
function get_type_vector(
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    type_property::AbstractString = "type",
    via::Maybe{Union{Tuple{<:AbstractString}, AbstractVector{<:AbstractString}}} = nothing,
)::AbstractVector{<:AbstractString}
    return get_axis_vector(daf; axis, query_suffix = ": $(escape_value(type_property))", via, empty_value = "", entries)
end

"""
    get_type_colors(
        daf::DafReader;
        type_axis::AbstractString = "type",
        empty_type_color::Maybe{AbstractString} = EMPTY_TYPE_COLOR,
    )::Dict{AbstractString, AbstractString}

Get the palette mapping each entry of the `daf` `type_axis` to its `color`, with an additional `empty_type_color` for
entries w/ no type (empty string). The `type_axis` must exist and must have a `color`; asking to color by a type is the
caller's decision, so a repository which can't answer is an error rather than an empty palette.

A `nothing` `empty_type_color` leaves the empty type out of the palette, which says you expect every entry to have a
type. An entry which doesn't is then rejected, because its type isn't a key of the palette.
"""
function get_type_colors(
    daf::DafReader;
    type_axis::AbstractString = "type",
    empty_type_color::Maybe{AbstractString} = EMPTY_TYPE_COLOR,
)::Dict{AbstractString, AbstractString}
    color_per_type = get_vector(daf, type_axis, "color")
    color_palette = Dict{AbstractString, AbstractString}(zip(names(color_per_type)[1], color_per_type))
    if empty_type_color !== nothing
        color_palette[""] = empty_type_color
    end
    return color_palette
end

"""
    get_axis_entries_vector(
        daf::DafReader;
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::AbstractVector{<:AbstractString}

Get the name for each of the `entries` of the `daf` `axis`. If no `entries` are specified, returns the names of all the
entries.
"""
function get_axis_entries_vector(
    daf::DafReader;
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::AbstractVector{<:AbstractString}
    names = entries_names(daf, axis, entries)
    if names === nothing
        names = axis_vector(daf, axis)
    end
    return names
end

"""
    get_vector_query(
        daf::DafReader;
        query::AbstractString,
        indices::Maybe{AbstractVector{<:Integer}} = nothing,
    )::AbstractVector{<:StorageScalarBase}

Execute a vector `query` on a `daf` repository and optionally access specific `indices` in it. This always returns a
dense vector. If indices are specified this returns a copy of the data, not a view.
"""
function get_vector_query(
    daf::DafReader;
    query::AbstractString,
    indices::Maybe{AbstractVector{<:Integer}} = nothing,
)::AbstractVector{<:StorageScalarBase}
    data = daf[query].array
    if indices !== nothing
        data = data[indices]
    else
        data = densify(data)
    end
    return data
end

"""
    get_matrix_query(
        daf::DafReader;
        query::AbstractString,
        row_indices::Maybe{AbstractVector{<:Integer}} = nothing,
        column_indices::Maybe{AbstractVector{<:Integer}} = nothing,
    )::AbstractMatrix{<:StorageScalarBase}

Execute a matrix `query` on a `daf` repository and optionally access specific `row_indices` and/or `column_indices` in
it. This always returns a dense matrix. If indices are specified this returns a copy of the data, not a view.
"""
function get_matrix_query(
    daf::DafReader;
    query::AbstractString,
    row_indices::Maybe{AbstractVector{<:Integer}} = nothing,
    column_indices::Maybe{AbstractVector{<:Integer}} = nothing,
)::AbstractMatrix{<:StorageScalarBase}
    data = daf[query].array

    is_view = false
    if column_indices !== nothing
        @views data = data[:, column_indices]
        is_view = true
    end
    if row_indices !== nothing
        @views data = data[row_indices, :]
        is_view = true
    end

    if is_view
        data = Matrix(data)
    else
        data = densify(data)
    end

    return data
end

# The names of the `entries` of an `axis`, given either their names or their indices. A `nothing` means all of them,
# which is left to the caller to resolve.
function entries_names(::DafReader, ::AbstractString, ::Nothing)::Nothing
    return nothing
end

function entries_names(
    daf::DafReader,
    axis::AbstractString,
    entries::AbstractVector{<:Integer},
)::AbstractVector{<:AbstractString}
    return axis_vector(daf, axis)[entries]
end

function entries_names(
    ::DafReader,
    ::AbstractString,
    entries::AbstractVector{<:AbstractString},
)::AbstractVector{<:AbstractString}
    return entries
end

# The indices of the `entries` of an `axis`, given either their names or their indices. A `nothing` means all of them,
# which is left to the caller to resolve.
function entries_indices(::DafReader, ::AbstractString, ::Nothing)::Nothing
    return nothing
end

function entries_indices(::DafReader, ::AbstractString, entries::AbstractVector{<:Integer})::AbstractVector{<:Integer}
    return entries
end

function entries_indices(
    daf::DafReader,
    axis::AbstractString,
    entries::AbstractVector{<:AbstractString},
)::AbstractVector{<:Integer}
    return axis_indices(daf, axis, entries)
end

end  # module
