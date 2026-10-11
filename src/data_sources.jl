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
(`axis_vector` for one, `axes_matrix` for two). The generic `put_` functions depend on nothing but the graph, so they
are in [`SomeGraphs.Sources`](@extref SomeGraphs SomeGraphs.Sources) (e.g., `put_vector_data!`,
`put_matrix_names_data!`). They are named after where the data goes, so they say which shape of sink they write, and
the same one serves a `DataFrame` or any other source of a vector. A specific `put_` here is named after its source
instead, since there is only one shape it can take.

The fill and put functions take [`Sinks`](@extref SomeGraphs SomeGraphs.Sources.Sinks): a single graph struct, or a
tuple or vector of them. This allows easily reusing the same data in multiple places in the graph (e.g., values,
labels, hovers). A data function writes only the sinks which are a
[`DataSink`](@extref SomeGraphs SomeGraphs.Sources.DataSink), and a configuration function only those which are a
[`ConfigurationSink`](@extref SomeGraphs SomeGraphs.Sources.ConfigurationSink), so a mixed collection is fine and
either may match nothing at all.

A `fill_` takes the `entries` of the axis to fill (and a matrix one also the `genes` or `row_entries` of its rows). When
these aren't given, it uses the entries the `sinks` are already named after (see
[`get_vector_names_data`](@extref SomeGraphs SomeGraphs.Sources.get_vector_names_data)). Only when the `sinks` have no
names does it use the whole axis. So when the first data source of a graph picks some entries, the rest follow it
without being told which they are. Sinks named after different entries are an error.

A `put_` has a method per leaf it writes and one explicit no-op method for the leaves of its kind it ignores (e.g., a
colors-only configuration ignores an axis). Its walking method takes everything else. A leaf covered by neither is a
`MethodError`; in particular a matrix leaf handed to a vector `put_`, or the other way around.

In general the data and the configuration functions take different sets of parameters (e.g., axis for data and gene
regularization for configuration of log scale). The `title` is a notable exception that gets passed to both sides, as it
is used for computing hover values (data) and as an axis/legend title (configuration).
"""
module DataSources

# The data sources generated from the table of each display bucket (e.g. `fill_n_cells!`) are exported where they are
# defined, by `define_bucket_data_sources`.
export EMPTY_CATEGORY_COLOR
export fill_axes_matrix_data!
export fill_axes_names_data!
export fill_axis_names_data!
export fill_axis_vector_data!
export fill_boolean!
export fill_category!
export fill_column_boolean_annotation!
export fill_column_names_data!
export fill_column_vector_data!
export fill_correlation!
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
export fill_size!
export fill_umap!
export GENE_FRACTION_REGULARIZATION_FOR_GRAPHS
export get_axes_matrix
export get_axis_entries_vector
export get_axis_vector
export get_boolean_vector
export get_category_colors
export get_category_vector
export get_column_vector
export get_correlation_vector
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
export get_size_vector
export get_skeleton_gene_indices
export get_top_marker_gene_indices
export get_umap_vector
export get_vector_query
export MAX_FOLD_FOR_GRAPHS
export put_boolean_configuration!
export put_boolean_values!
export put_category_configuration!
export put_category_values!
export put_coordinates_configuration!
export put_correlation_configuration!
export put_correlation_values!
export put_expression_level_configuration!
export put_fold_configuration!
export put_size_configuration!
export put_size_values!
export put_umap_data!

using DataAxesFormats
using DataFrames
using OrderedCollections
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
The color to give to entities without any value of a category (an empty string), e.g. with no type or no block. A
computed palette takes its colors from the `HSV` palette, which holds only fully saturated colors. The `"grey"` color
isn't one of them.
"""
EMPTY_CATEGORY_COLOR::AbstractString = "grey"

# The keywords which every generated data source takes after its `axis`, as declarations to splice into its signature.
const ENTRIES_KEYWORDS = (
    :(entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing),
    :(via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing),
)

# Define the `get_<property>_vector` and `fill_<property>!` of each property in the `table` of a display bucket. The get
# calls the bucket's `get_function` with its `property`. The fill calls the get, and passes the values to the bucket's
# `put_values_function`. The get takes the `axis`, the `ENTRIES_KEYWORDS` and the bucket's `get_keywords`. The fill also
# takes the `title` and the bucket's `put_keywords`. A table line gives the `property`, its default `axis` and `title`,
# and its `doc`. It may also give its own defaults for any of the bucket's keywords.
function define_bucket_data_sources(;
    get_function::Symbol,
    put_values_function::Symbol,
    get_keywords::Tuple,
    put_keywords::Tuple,
    table::Tuple,
)::Nothing
    for line in table
        get_name = Symbol("get_$(line.property)_vector")
        get_arguments = [
            line_keyword(line, keyword) for
            keyword in (:(axis::AbstractString = $(line.axis)), ENTRIES_KEYWORDS..., get_keywords...)
        ]
        get_names = keyword_names(get_arguments)
        define_property_data_source(;
            line,
            function_name = get_name,
            leading_arguments = (:(daf::DafReader),),
            keyword_arguments = get_arguments,
            result_type = :(AbstractVector),
            doc_suffix = "This is directed at [`$(get_function)`](@ref).",
            body = :(return $get_function(daf; property = $(line.property), $(get_names...))),
        )

        put_arguments = [
            line_keyword(line, keyword) for
            keyword in (:(title::Maybe{AbstractString} = $(line.title)), put_keywords...)
        ]
        put_names = keyword_names(put_arguments)
        fill_arguments =  # NOJET
            vcat(get_arguments, [argument for argument in put_arguments if !(argument.args[1].args[1] in get_names)])
        define_property_data_source(;
            line,
            function_name = Symbol("fill_$(line.property)!"),
            leading_arguments = (:(sinks::VectorDataSinks), :(daf::DafReader)),
            keyword_arguments = fill_arguments,
            result_type = :(Nothing),
            doc_suffix = "This fetches it with [`$(get_name)`](@ref), " *
                         "and shows it with [`$(put_values_function)`](@ref).",
            body = quote
                entries = sinks_entries(sinks, entries)
                value_per_entry = $get_name(daf; $(get_names...))
                return $put_values_function(sinks, daf, value_per_entry; axis, entries, $(put_names...))
            end,
        )
    end
    return nothing
end

# Define and export one generated data source, with its docstring, which ends with the `doc_suffix`.
function define_property_data_source(;
    line::NamedTuple,
    function_name::Symbol,
    leading_arguments::Tuple,
    keyword_arguments::AbstractVector{Expr},
    result_type::Union{Symbol, Expr},
    doc_suffix::AbstractString,
    body::Expr,
)::Nothing
    signature_lines = [
        ["        $(argument)," for argument in leading_arguments[1:(end - 1)]]...,
        "        $(leading_arguments[end]);",
        ["        $(keyword.args[1]) = $(default_text(keyword.args[2]))," for keyword in keyword_arguments]...,
    ]
    docstring = """
                    $(function_name)(
                $(join(signature_lines, "\n"))
                    )::$(result_type)

                $(strip(line.doc)) $(doc_suffix)
                """
    @eval begin
        @doc $docstring function $function_name($(leading_arguments...); $(keyword_arguments...))::$result_type  # NOLINT
            return $body
        end
        export $function_name
    end
    return nothing
end

# The names of keyword arguments.
function keyword_names(keyword_arguments::AbstractVector{Expr})::Vector{Symbol}  # ONLY SEEMS UNTESTED
    return [keyword.args[1].args[1] for keyword in keyword_arguments]
end

# A keyword declaration as a keyword argument, with its default given by the table `line` if it gives one.
function line_keyword(line::NamedTuple, keyword::Expr)::Expr  # ONLY SEEMS UNTESTED
    declaration, default = keyword.args
    name = declaration.args[1]
    if haskey(line, name)
        default = line[name]
    end
    return Expr(:kw, declaration, default)
end

# How a default is written in a docstring.
function default_text(default::Any)::AbstractString  # ONLY SEEMS UNTESTED
    if default isa AbstractString
        return repr(default)
    else
        return string(default)
    end
end

"""
    put_size_configuration!(
        sinks::Sinks;
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Show a size (a count) on a log scale of the `log_base`, regularized by 1, in the `YlOrRd` color scale, named by the
`title`. This is how every count-like quantity is shown: how many cells, how many metacells, how many UMIs, and the
means of these. They all span orders of magnitude, which is what the log scale is for. The regularization gives a zero
count, which is rare, a place on the scale.

The color scale follows the range of the data. It runs from yellow through orange to red, so even the smallest count is
visible. A scale which starts at white would lose it against the background.
"""
function put_size_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector};
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_size_configuration!(sink; title, show_legend, log_base)
    end
    return nothing
end

function put_size_configuration!(
    sizes::SizesConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    sizes.show_legend = show_legend
    put_size_configuration!(sizes.scale; log_base)
    if title !== nothing
        sizes.title = title
    end
    return nothing
end

function put_size_configuration!(
    axis::AxisConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,  # NOLINT
    log_base::LogBase = Log10Base,
)::Nothing
    put_size_configuration!(axis.scale; log_base)
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_size_configuration!(
    colors::ColorsConfiguration;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    colors.palette = "YlOrRd"
    colors.show_legend = show_legend
    put_size_configuration!(colors.scale; log_base)
    if title !== nothing
        colors.title = title
    end
    return nothing
end

function put_size_configuration!(
    scale::ScaleConfiguration;
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = true,  # NOLINT
    log_base::LogBase = Log10Base,
)::Nothing
    scale.log_base = log_base
    scale.log_regularization = 1
    return nothing
end

# The keywords of the size bucket's get and of its `put_size_values!`, other than the ones every generated data source
# takes.
const SIZE_GET_KEYWORDS = (:(empty_value::Real = 1),)
const SIZE_PUT_KEYWORDS = (:(show_legend::Bool = true), :(log_base::LogBase = Log10Base))

"""
    get_size_vector(
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 1,
    )::AbstractVector{<:Real}

Get the size (count) `property` of each of the `entries` of the `daf` `axis`. By default this is the direct property of
the axis entries. However, by specifying `via`, it is possible to access indirect sizes (e.g., the total UMIs of the
block of the metacell of each cell).

An entry whose `via` chain breaks gets the `empty_value`, which is a count like any other, so the log scale has
something to work with.
"""
function get_size_vector(
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 1,
)::AbstractVector{<:Real}
    return get_axis_vector(daf; axis, query_suffix = ": $(escape_value(property))", via, empty_value, entries)
end

"""
    fill_size!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = property,
        empty_value::Real = 1,
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Fill the `sinks` with the size (count) `property` of each of the `entries` of the `daf` `axis` (see
[`get_size_vector`](@ref)), shown by [`put_size_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_size!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = property,
    empty_value::Real = 1,
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_size_vector(daf; axis, property, entries, via, empty_value)
    put_size_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend, log_base)
    return nothing
end

"""
    put_size_values!(
        sinks::VectorDataSinks,
        daf::DafReader,
        value_per_entry::AbstractVector{<:Real};
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Put the size (count) `value_per_entry` of each of the `entries` of the `daf` `axis` into the `sinks`, shown by
[`put_size_configuration!`](@ref), and name the entities after the `entries`. Every size data source ends with this.
"""
function put_size_values!(
    sinks::VectorDataSinks,
    daf::DafReader,
    value_per_entry::AbstractVector{<:Real};
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    put_vector_data!(sinks, value_per_entry; title)
    put_size_configuration!(sinks; title, show_legend, log_base)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

define_bucket_data_sources(;
    get_function = :get_size_vector,
    put_values_function = :put_size_values!,
    get_keywords = SIZE_GET_KEYWORDS,
    put_keywords = SIZE_PUT_KEYWORDS,
    table = (
        (; property = "total_UMIs", axis = "metacell", title = "total UMIs", doc = "The total UMIs of each entry."),
        (;
            property = "n_cells",
            axis = "metacell",
            title = "cells",
            doc = "The number of cells of each entry. Both metacells and blocks have one.",
        ),
        (;
            property = "n_metacells",
            axis = "block",
            title = "metacells",
            doc = "The number of metacells of each entry. Only blocks have one, so the default `axis` is the block.",
        ),
    ),
)

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
    return get_n_cells_vector(daf; axis, entries, via) ./ get_n_metacells_vector(daf; axis, entries, via)  # NOLINT
end

"""
    fill_mean_cells_per_metacell!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean cells per metacell",
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Fill the `sinks` with the mean number of cells per metacell of each of the `entries` of the `daf` `axis`, shown by
[`put_size_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_mean_cells_per_metacell!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean cells per metacell",
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_mean_cells_per_metacell_vector(daf; axis, entries, via)
    put_size_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend, log_base)
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
    return get_total_UMIs_vector(daf; axis, entries, via) ./ get_n_metacells_vector(daf; axis, entries, via)  # NOLINT
end

"""
    fill_mean_total_UMIs_per_metacell!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString = "block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean UMIs per metacell",
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Fill the `sinks` with the mean total UMIs per metacell of each of the `entries` of the `daf` `axis`, shown by
[`put_size_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_mean_total_UMIs_per_metacell!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString = "block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean UMIs per metacell",
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_mean_total_UMIs_per_metacell_vector(daf; axis, entries, via)
    put_size_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend, log_base)
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
    return get_total_UMIs_vector(daf; axis, entries, via) ./ get_n_cells_vector(daf; axis, entries, via)  # NOLINT
end

"""
    fill_mean_total_UMIs_per_cell!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = "mean UMIs per cell",
        show_legend::Bool = true,
        log_base::LogBase = Log10Base,
    )::Nothing

Fill the `sinks` with the mean total UMIs per cell of each of the `entries` of the `daf` `axis`, shown by
[`put_size_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_mean_total_UMIs_per_cell!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = "mean UMIs per cell",
    show_legend::Bool = true,
    log_base::LogBase = Log10Base,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_mean_total_UMIs_per_cell_vector(daf; axis, entries, via)
    put_size_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend, log_base)
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
        sinks::VectorDataSinks,
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
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    type_property::AbstractString = "type",
    empty_value::Real = 0,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    entries = sinks_entries(sinks, entries)
    put_vector_data!(sinks, get_global_flow_order_vector(daf; axis, entries, type_property, empty_value); title)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

# The keywords of the Boolean bucket's get and of its `put_boolean_values!`, other than the ones every generated data
# source takes.
const BOOLEAN_GET_KEYWORDS = (:(empty_value::Bool = false),)
const BOOLEAN_PUT_KEYWORDS = (:(show_legend::Bool = false),)

"""
    get_boolean_vector(
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Bool = false,
    )::AbstractVector{<:AbstractString}

Get the Boolean `property` of each of the `entries` of the `daf` `axis`, as the *strings* `"true"` and `"false"`.
Coloring by a Boolean is coloring by two categories, and a categorical palette is keyed by strings. By specifying `via`,
it is possible to access an indirect property (e.g., whether the anchor gene of each module is a marker). An entry whose
`via` chain breaks gets the `empty_value`.
"""
function get_boolean_vector(
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Bool = false,
)::AbstractVector{<:AbstractString}
    return string.(get_axis_vector(daf; axis, query_suffix = ": $(escape_value(property))", via, empty_value, entries))
end

"""
    put_boolean_configuration!(
        sinks::Sinks;
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = false,
    )::Nothing

Show Boolean *string* values: `true` in black and `false` in light grey, named by the `title`. The values themselves
are put in by [`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!).

Black against light grey has only two values and is named by its title, so `show_legend` is off.
"""
function put_boolean_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector};
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = false,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_boolean_configuration!(sink; title, show_legend)
    end
    return nothing
end

# Booleans are shown only as colors. This method ignores an axis, a scale and sizes.
function put_boolean_configuration!(
    ::Union{AxisConfiguration, ScaleConfiguration, SizesConfiguration};
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = false,  # NOLINT
)::Nothing
    return nothing
end

function put_boolean_configuration!(
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
    fill_boolean!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = property,
        empty_value::Bool = false,
        show_legend::Bool = false,
    )::Nothing

Fill the `sinks` with the Boolean `property` of each of the `entries` of the `daf` `axis` (see
[`get_boolean_vector`](@ref)), shown by [`put_boolean_configuration!`](@ref), and name the entities after the
`entries`.
"""
function fill_boolean!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = property,
    empty_value::Bool = false,
    show_legend::Bool = false,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_boolean_vector(daf; axis, property, entries, via, empty_value)
    put_boolean_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend)
    return nothing
end

"""
    put_boolean_values!(
        sinks::VectorDataSinks,
        daf::DafReader,
        value_per_entry::AbstractVector{<:AbstractString};
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = false,
    )::Nothing

Put the Boolean *string* `value_per_entry` of each of the `entries` of the `daf` `axis` into the `sinks`, shown by
[`put_boolean_configuration!`](@ref), and name the entities after the `entries`. Every Boolean data source ends with
this.
"""
function put_boolean_values!(
    sinks::VectorDataSinks,
    daf::DafReader,
    value_per_entry::AbstractVector{<:AbstractString};
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = false,
)::Nothing
    put_vector_data!(sinks, value_per_entry; title)
    put_boolean_configuration!(sinks; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

define_bucket_data_sources(;
    get_function = :get_boolean_vector,
    put_values_function = :put_boolean_values!,
    get_keywords = BOOLEAN_GET_KEYWORDS,
    put_keywords = BOOLEAN_PUT_KEYWORDS,
    table = (
        (; property = "is_lateral", axis = "gene", title = "lateral", doc = "Whether each gene is lateral."),
        (; property = "is_marker", axis = "gene", title = "marker", doc = "Whether each gene is a marker gene."),
        (;
            property = "is_transcription_factor",
            axis = "gene",
            title = "transcription factor",
            doc = "Whether each gene is a transcription factor (binds to DNA).",
        ),
        (;
            property = "is_regulator",
            axis = "gene",
            title = "regulator",
            doc = "Whether each gene regulates the expression of other genes.",
        ),
        (; property = "is_forbidden", axis = "gene", title = "forbidden", doc = "Whether each gene is forbidden."),
        (; property = "is_skeleton", axis = "gene", title = "skeleton", doc = "Whether each gene is a skeleton gene."),
        (;
            property = "is_correlated_with_skeleton",
            axis = "gene",
            title = "correlated with skeleton",
            doc = "Whether each gene is correlated with some skeleton gene.",
        ),
        (;
            property = "is_mitochondrial",
            axis = "gene",
            title = "mitochondrial",
            doc = "Whether each gene is mitochondrial.",
        ),
        (; property = "is_ribosomal", axis = "gene", title = "ribosomal", doc = "Whether each gene is ribosomal."),
        (; property = "is_excluded", axis = "gene", title = "excluded", doc = "Whether each gene is excluded."),
    ),
)

"""
    fill_gene_expression!(
        sinks::VectorDataSinks,
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
        title::Maybe{AbstractString} = "\$(gene) fraction",
        show_legend::Bool = true,
        log_base::LogBase = Log2Base,
    )::Nothing

Fill the `sinks` with the `gene` expression level (linear fraction) per each of the `entries` of some `daf` `axis`,
shown by [`put_expression_level_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_gene_expression!(
    sinks::VectorDataSinks,
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "$(gene) fraction",
    show_legend::Bool = true,
    log_base::LogBase = Log2Base,
)::Nothing
    entries = sinks_entries(sinks, entries)
    put_vector_data!(sinks, get_gene_expression_vector(daf; gene, axis, entries); title)
    put_expression_level_configuration!(sinks; gene_fraction_regularization, title, show_legend, log_base)
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
        sinks::VectorDataSinks,
        daf::DafReader;
        gene::AbstractString,
        axis::AbstractString = "base_block",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = "correlation",
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the correlation of the `gene` per each of the `entries` of some `daf` `axis`, shown by
[`put_correlation_configuration!`](@ref), and name the entities after the `entries`.
"""
function fill_gene_correlation!(
    sinks::VectorDataSinks,
    daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "correlation",
    show_legend::Bool = true,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_gene_correlation_vector(daf; gene, axis, entries)
    put_correlation_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend)
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
    put_correlation_configuration!(
        sinks::Sinks,
        values::AbstractArray{<:Real};
        title::Maybe{AbstractString} = "correlation",
        show_legend::Bool = true,
    )::Nothing

Show a correlation (or a difference of correlations) on a linear scale, named by the `title`. The `values` are the ones
put in by [`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!) or
[`put_matrix_data!`](@extref SomeGraphs SomeGraphs.Sources.put_matrix_data!), and the colors follow their range.

As colors, zero is white, negative values go towards blue and positive values towards red, each side as far as the
values go. The range always includes zero, even when all the values are positive (or all negative), so white always
means no correlation. Unlike a fold, the range need not be symmetric, since most correlations are positive. If all
the values are zero, the range is -1 to 1.
"""
function put_correlation_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector},
    values::AbstractArray{<:Real};
    title::Maybe{AbstractString} = "correlation",
    show_legend::Bool = true,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_correlation_configuration!(sink, values; title, show_legend)
    end
    return nothing
end

# Correlations are never shown as sizes, and an axis scale of correlations is linear. This method ignores sizes and a
# scale.
function put_correlation_configuration!(
    ::Union{ScaleConfiguration, SizesConfiguration},
    ::AbstractArray{<:Real};
    title::Maybe{AbstractString} = "correlation",  # NOLINT
    show_legend::Bool = true,  # NOLINT
)::Nothing
    return nothing
end

function put_correlation_configuration!(
    axis::AxisConfiguration,
    ::AbstractArray{<:Real};
    title::Maybe{AbstractString} = "correlation",
    show_legend::Bool = true,  # NOLINT
)::Nothing
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_correlation_configuration!(
    colors::ColorsConfiguration,
    values::AbstractArray{<:Real};
    title::Maybe{AbstractString} = "correlation",
    show_legend::Bool = true,
)::Nothing
    minimum_value = min(0, minimum(value for value in values if isfinite(value); init = 0))
    maximum_value = max(0, maximum(value for value in values if isfinite(value); init = 0))
    if minimum_value == maximum_value
        minimum_value = -1
        maximum_value = 1
    end

    color_scale = Tuple{Real, AbstractString}[]
    if minimum_value < 0
        push!(color_scale, (minimum_value, "rgb(5,10,172)"))
    end
    push!(color_scale, (0, "rgb(255,255,255)"))
    if maximum_value > 0
        push!(color_scale, (maximum_value, "rgb(178,10,28)"))
    end

    colors.palette = color_scale
    colors.scale.minimum = minimum_value
    colors.scale.maximum = maximum_value
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

# The keywords of the correlation bucket's get and of its `put_correlation_values!`, other than the ones every
# generated data source takes.
const CORRELATION_GET_KEYWORDS = (:(empty_value::Real = 0),)
const CORRELATION_PUT_KEYWORDS = (:(show_legend::Bool = true),)

"""
    get_correlation_vector(
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::Real = 0,
    )::AbstractVector{<:Real}

Get the correlation `property` of each of the `entries` of the `daf` `axis`. By specifying `via`, it is possible to
access an indirect correlation. An entry whose `via` chain breaks gets the `empty_value`, which by default means no
correlation.
"""
function get_correlation_vector(
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Real = 0,
)::AbstractVector{<:Real}
    return get_axis_vector(daf; axis, query_suffix = ": $(escape_value(property))", via, empty_value, entries)
end

"""
    fill_correlation!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = property,
        empty_value::Real = 0,
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the correlation `property` of each of the `entries` of the `daf` `axis` (see
[`get_correlation_vector`](@ref)), shown by [`put_correlation_configuration!`](@ref), and name the entities after the
`entries`.
"""
function fill_correlation!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = property,
    empty_value::Real = 0,
    show_legend::Bool = true,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_correlation_vector(daf; axis, property, entries, via, empty_value)
    put_correlation_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend)
    return nothing
end

"""
    put_correlation_values!(
        sinks::VectorDataSinks,
        daf::DafReader,
        value_per_entry::AbstractVector{<:Real};
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = true,
    )::Nothing

Put the correlation `value_per_entry` of each of the `entries` of the `daf` `axis` into the `sinks`, shown by
[`put_correlation_configuration!`](@ref), and name the entities after the `entries`. Every correlation data source ends
with this.
"""
function put_correlation_values!(
    sinks::VectorDataSinks,
    daf::DafReader,
    value_per_entry::AbstractVector{<:Real};
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, value_per_entry; title)
    put_correlation_configuration!(sinks, value_per_entry; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

define_bucket_data_sources(;
    get_function = :get_correlation_vector,
    put_values_function = :put_correlation_values!,
    get_keywords = CORRELATION_GET_KEYWORDS,
    put_keywords = CORRELATION_PUT_KEYWORDS,
    table = (
        (;
            property = "correlation_between_cells_and_punctuated_metacells",
            axis = "gene",
            title = "cells correlation",
            doc = "The correlation of each gene between the cells and their punctuated metacells.",
        ),
        (;
            property = "correlation_between_cells_and_projected_punctuated_metacells",
            axis = "gene",
            title = "projected cells correlation",
            doc = "The correlation of each gene between the cells and their projected punctuated metacells.",
        ),
        (;
            property = "correlation_between_cells_and_projected_metacells",
            axis = "gene",
            title = "projection correlation",
            doc = "The correlation of each gene between the query cells and their projected atlas metacells.",
        ),
    ),
)

"""
    fill_gene_correlation_change!(
        sinks::VectorDataSinks,
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
    sinks::VectorDataSinks,
    daf::DafReader,
    base_daf::DafReader;
    gene::AbstractString,
    axis::AbstractString = "base_block",
    base_axis::AbstractString = "base_block",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "correlation change",
    show_legend::Bool = true,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_gene_correlation_change_vector(daf, base_daf; gene, axis, base_axis, entries)
    put_correlation_values!(sinks, daf, value_per_entry; axis, entries, title, show_legend)
    return nothing
end

"""
    fill_genes_expression_matrix!(
        sinks::MatrixDataSinks,
        daf::DafReader;
        genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        axis::AbstractString = "metacell",
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
        title::Maybe{AbstractString} = "fraction",
        show_legend::Bool = true,
        log_base::LogBase = Log2Base,
    )::Nothing

Fill the `sinks` with the `genes` expression level (linear fraction) per each of the `entries` of some `daf` `axis`,
shown by [`put_expression_level_configuration!`](@ref), and name the rows and the columns. The `genes` are the rows
and the `entries` are the columns.
"""
function fill_genes_expression_matrix!(
    sinks::MatrixDataSinks,
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
    log_base::LogBase = Log2Base,
)::Nothing
    genes, entries = sinks_matrix_entries(sinks, genes, entries)
    put_matrix_data!(sinks, get_genes_expression_matrix(daf; genes, axis, entries); title)
    put_expression_level_configuration!(sinks; gene_fraction_regularization, title, show_legend, log_base)
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
    put_expression_level_configuration!(
        sinks::Sinks;
        gene_fraction_regularization::Real = $(GENE_FRACTION_REGULARIZATION_FOR_GRAPHS),
        title::Maybe{AbstractString} = "fraction",
        show_legend::Bool = true,
        log_base::LogBase = Log2Base,
    )::Nothing

Show an expression level (the fraction of the UMIs of a gene, or of a gene module) on a log scale of the `log_base`,
regularized by `gene_fraction_regularization`, named by the `title`. The expression levels themselves are put in by
[`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!) or
[`put_matrix_data!`](@extref SomeGraphs SomeGraphs.Sources.put_matrix_data!).

As colors, the expression levels are shown in the `WtYlOrRdDm` color scale, which starts at zero. Zero expression is
common, and is white, so any expression above it shows a color.

The `show_legend` applies where the expression levels are the colors. An axis is read off its own ticks and ignores it.
"""
function put_expression_level_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector};
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
    log_base::LogBase = Log2Base,
)::Nothing
    @assert gene_fraction_regularization >= 0
    visit_configuration_sinks(sinks) do sink
        return put_expression_level_configuration!(sink; gene_fraction_regularization, title, show_legend, log_base)
    end
    return nothing
end

# Expression levels are never shown as sizes. This method ignores sizes.
function put_expression_level_configuration!(
    ::SizesConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,  # NOLINT
    title::Maybe{AbstractString} = "fraction",  # NOLINT
    show_legend::Bool = true,  # NOLINT
    log_base::LogBase = Log2Base,  # NOLINT
)::Nothing
    return nothing
end

function put_expression_level_configuration!(
    axis::AxisConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,  # NOLINT
    log_base::LogBase = Log2Base,
)::Nothing
    put_expression_level_configuration!(axis.scale; gene_fraction_regularization, log_base)
    if title !== nothing
        axis.title = title
    end
    return nothing
end

function put_expression_level_configuration!(
    colors::ColorsConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fraction",
    show_legend::Bool = true,
    log_base::LogBase = Log2Base,
)::Nothing
    put_expression_level_configuration!(colors.scale; gene_fraction_regularization, log_base)
    colors.palette = "WtYlOrRdDm"
    colors.scale.minimum = 0
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

function put_expression_level_configuration!(
    scale::ScaleConfiguration;
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_GRAPHS,
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = true,  # NOLINT
    log_base::LogBase = Log2Base,
)::Nothing
    scale.log_base = log_base
    scale.log_regularization = gene_fraction_regularization
    return nothing
end

"""
    put_umap_data!(
        sinks::VectorDataSinks,
        coordinate_per_entry::AbstractVector{<:AbstractFloat};
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Put a UMAP `coordinate_per_entry` into the `sinks`, as the values of a role.
"""
function put_umap_data!(
    sinks::Union{AnyContainer, ConfigurationLeaf, Tuple, AbstractVector},
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
# It is not something a heatmap side is arranged by, either.
function put_umap_data!(
    ::Union{VectorEntitiesData, ArrangementData},
    ::AbstractVector{<:AbstractFloat};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    return nothing
end

function put_umap_data!(
    leaf::MatrixDataLeaf,
    ::AbstractVector{<:AbstractFloat};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    return throw(matrix_sink_error(leaf))
end

"""
    fill_umap!(
        sinks::VectorDataSinks,
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
    sinks::VectorDataSinks,
    daf::DafReader;
    coordinate::AbstractString,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = "UMAP $(uppercase(coordinate))",
)::Nothing
    entries = sinks_entries(sinks, entries)
    put_umap_data!(sinks, get_umap_vector(daf; coordinate, axis, entries); title)
    put_coordinates_configuration!(sinks; title)
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
    put_coordinates_configuration!(sinks::Sinks; title::Maybe{AbstractString} = nothing)::Nothing

Show the coordinates of an embedding (e.g. a UMAP) without ticks or a grid, named by the `title`. Only which points are
near which other points means anything, so the values themselves are not worth labelling.
"""
function put_coordinates_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector};
    title::Maybe{AbstractString} = nothing,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_coordinates_configuration!(sink; title)
    end
    return nothing
end

# Coordinates are shown only along an axis. This method ignores a scale, colors and sizes.
function put_coordinates_configuration!(
    ::Union{ScaleConfiguration, ColorsConfiguration, SizesConfiguration};
    title::Maybe{AbstractString} = nothing,  # NOLINT
)::Nothing
    return nothing
end

function put_coordinates_configuration!(axis::AxisConfiguration; title::Maybe{AbstractString} = nothing)::Nothing
    if title !== nothing
        axis.title = title
    end
    axis.show_ticks = false
    axis.show_grid = false
    return nothing
end

"""
    fill_genes_fold_matrix!(
        sinks::MatrixDataSinks,
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
    sinks::MatrixDataSinks,
    daf::DafReader;
    genes::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    axis::AbstractString = "metacell",
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fold from median",
    show_legend::Bool = true,
)::Nothing
    genes, entries = sinks_matrix_entries(sinks, genes, entries)
    put_matrix_data!(sinks, get_genes_fold_matrix(daf; genes, axis, entries); title)  # NOJET
    put_fold_configuration!(sinks; max_fold, title, show_legend)
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
    put_fold_configuration!(
        sinks::Sinks;
        max_fold::Real = MAX_FOLD_FOR_GRAPHS,
        title::Maybe{AbstractString} = "fold from median",
        show_legend::Bool = true,
    )::Nothing

Show a fold (the log base 2 of a ratio, e.g. of a gene's expression to its median) in the range -`max_fold` (blue) to
0 (white) to +`max_fold` (red), named by the `title`, and include it in the legend if `show_legend`. The folds
themselves are put in by [`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!) or
[`put_matrix_data!`](@extref SomeGraphs SomeGraphs.Sources.put_matrix_data!).
"""
function put_fold_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector};
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,
    title::Maybe{AbstractString} = "fold from median",
    show_legend::Bool = true,
)::Nothing
    @assert max_fold > 0
    visit_configuration_sinks(sinks) do sink
        return put_fold_configuration!(sink; max_fold, title, show_legend)
    end
    return nothing
end

# Folds are shown only as colors. This method ignores an axis, a scale and sizes.
function put_fold_configuration!(
    ::Union{AxisConfiguration, ScaleConfiguration, SizesConfiguration};
    max_fold::Real = MAX_FOLD_FOR_GRAPHS,  # NOLINT
    title::Maybe{AbstractString} = "fold from median",  # NOLINT
    show_legend::Bool = true,  # NOLINT
)::Nothing
    return nothing
end

function put_fold_configuration!(
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
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::Nothing

Name the entities of the `sinks` after each of the `entries` of the `daf` `axis`.

This fetches the names with [`get_axis_entries_vector`](@ref), then hands them to
[`put_vector_names_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_names_data!).
"""
function fill_axis_names_data!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::Nothing
    entries = sinks_entries(sinks, entries)
    put_vector_names_data!(sinks, get_axis_entries_vector(daf; axis, entries))
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
    fill_axes_names_data!(
        sinks::MatrixDataSinks,
        daf::DafReader;
        rows_axis::AbstractString,
        columns_axis::AbstractString,
        row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    )::Nothing

Name the rows of the `sinks` after the `row_entries` of the `daf` `rows_axis`, and their columns after the
`column_entries` of its `columns_axis`.

This fetches the names with [`get_axis_entries_vector`](@ref), then hands them to
[`put_matrix_names_data!`](@extref SomeGraphs SomeGraphs.Sources.put_matrix_names_data!).
"""
function fill_axes_names_data!(
    sinks::MatrixDataSinks,
    daf::DafReader;
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
)::Nothing
    row_entries, column_entries = sinks_matrix_entries(sinks, row_entries, column_entries)
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
        sinks::MatrixDataSinks,
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
    sinks::MatrixDataSinks,
    daf::DafReader;
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    query_suffix::AbstractString,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    row_entries, column_entries = sinks_matrix_entries(sinks, row_entries, column_entries)
    put_matrix_data!(
        sinks,
        get_axes_matrix(daf; rows_axis, columns_axis, query_suffix, row_entries, column_entries);
        title,
    )
    fill_axes_names_data!(sinks, daf; rows_axis, columns_axis, row_entries, column_entries)
    return nothing
end

"""
    fill_axis_vector_data!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        query_suffix::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the result of the query `@ axis query_suffix` for each of the `entries` of the `daf` `axis`, and
name the entities after the `entries`.

This fetches the data with [`get_axis_vector`](@ref), then hands it to
[`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!).
"""
function fill_axis_vector_data!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    query_suffix::AbstractString,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::Maybe{StorageScalarBase} = nothing,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    entries = sinks_entries(sinks, entries)
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
        sinks::VectorDataSinks,
        frame::DataFrame;
        column::AbstractString,
        title::Maybe{AbstractString} = nothing,
    )::Nothing

Fill the `sinks` with the `column` of a `frame`, one value per row.

The entities are not named here, unlike [`fill_axis_vector_data!`](@ref); a frame has no one column which identifies
its rows, so say which one does with [`fill_column_names_data!`](@ref).
"""
function fill_column_vector_data!(
    sinks::VectorDataSinks,
    frame::DataFrame;
    column::AbstractString,
    title::Maybe{AbstractString} = nothing,
)::Nothing
    put_vector_data!(sinks, get_column_vector(frame; column); title)
    return nothing
end

"""
    fill_column_names_data!(sinks::VectorDataSinks, frame::DataFrame; column::AbstractString)::Nothing

Name the entities of the `sinks` after the `column` of a `frame`, one name per row.
"""
function fill_column_names_data!(sinks::VectorDataSinks, frame::DataFrame; column::AbstractString)::Nothing
    put_vector_names_data!(sinks, string.(get_column_vector(frame; column)))
    return nothing
end

"""
    fill_column_boolean_annotation!(
        sinks::VectorDataSinks,
        frame::DataFrame;
        column::AbstractString,
        title::Maybe{AbstractString} = column,
        show_legend::Bool = false,
    )::Nothing

Fill the `sinks` with the Boolean `column` of a `frame`, shown as `true` in black and `false` in light grey. The
`DataFrame` twin of [`fill_boolean!`](@ref).
"""
function fill_column_boolean_annotation!(
    sinks::VectorDataSinks,
    frame::DataFrame;
    column::AbstractString,
    title::Maybe{AbstractString} = column,
    show_legend::Bool = false,
)::Nothing
    put_vector_data!(sinks, string.(get_column_vector(frame; column)); title)
    put_boolean_configuration!(sinks; title, show_legend)
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
        sinks::VectorDataSinks,
        frame::DataFrame;
        prefix::AbstractString,
        side_name::AbstractString,
        regulators_count::Integer,
    )::Nothing

Add a hover per row of a gene report `frame` to the `sinks`, saying which regulators the gene is in a module with in
the base blocks of one side. The lines carry their own labels, so nothing is prefixed to them.
"""
function fill_module_regulators_hovers!(
    sinks::VectorDataSinks,
    frame::DataFrame;
    prefix::AbstractString,
    side_name::AbstractString,
    regulators_count::Integer,
)::Nothing
    put_vector_data!(sinks, get_module_regulators_hovers(frame; prefix, side_name, regulators_count))
    return nothing
end

"""
    put_category_configuration!(
        sinks::Sinks,
        color_per_value::CategoricalColors;
        title::Maybe{AbstractString} = nothing,
        show_legend::Bool = true,
    )::Nothing

Show the values of a category (e.g. types or blocks) of the `sinks` using the `color_per_value` palette, named by the
`title`. The values themselves are put in by
[`put_vector_data!`](@extref SomeGraphs SomeGraphs.Sources.put_vector_data!).

The palette is required, and covers every value, including the empty one of entries with no value. Where it comes from
is up to the data source. For example, [`get_category_colors`](@ref) reads the colors given by the repository, or
computes them.

A color says nothing without a legend to read it by, so `show_legend` is on. Turn it off for a category with too many
values to list, or for an annotation of a graph whose colors already name the same values.
"""
function put_category_configuration!(
    sinks::Union{AnyContainer, DataLeaf, Tuple, AbstractVector},
    color_per_value::CategoricalColors;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    visit_configuration_sinks(sinks) do sink
        return put_category_configuration!(sink, color_per_value; title, show_legend)
    end
    return nothing
end

# Categories are shown only as colors. This method ignores an axis, a scale and sizes.
function put_category_configuration!(
    ::Union{AxisConfiguration, ScaleConfiguration, SizesConfiguration},
    ::CategoricalColors;
    title::Maybe{AbstractString} = nothing,  # NOLINT
    show_legend::Bool = true,  # NOLINT
)::Nothing
    return nothing
end

function put_category_configuration!(
    colors::ColorsConfiguration,
    color_per_value::CategoricalColors;
    title::Maybe{AbstractString} = nothing,
    show_legend::Bool = true,
)::Nothing
    colors.palette = color_per_value
    colors.show_legend = show_legend
    if title !== nothing
        colors.title = title
    end
    return nothing
end

"""
    get_category_colors(
        daf::DafReader;
        axis::AbstractString,
        empty_value::AbstractString = "",
        empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
    )::OrderedDict{AbstractString, AbstractString}

Get the palette mapping each entry of the `daf` `axis` of a category (e.g. the types or the blocks) to its color, with
an additional `empty_color` for the `empty_value` of entries with no value. If the `axis` has a `color` vector (as the
types usually do), the palette uses it, in the order of the axis. Otherwise, the colors come from
[`categorical_palette`](@extref SomeGraphs SomeGraphs.Common.categorical_palette), which spreads the entries over the
cyclical `HSV` palette. Either way, the palette covers all the entries of the axis, so each value has the same color in
every graph, whichever entries the graph shows.

A `nothing` `empty_color` leaves the `empty_value` out of the palette, which says you expect every entry to have a
value. An entry which doesn't is then rejected, because its value isn't a key of the palette.
"""
function get_category_colors(
    daf::DafReader;
    axis::AbstractString,
    empty_value::AbstractString = "",
    empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
)::OrderedDict{AbstractString, AbstractString}
    if has_vector(daf, axis, "color")
        color_per_entry = get_vector(daf, axis, "color")
        color_palette = OrderedDict{AbstractString, AbstractString}(zip(names(color_per_entry)[1], color_per_entry))
    else
        color_palette = OrderedDict{AbstractString, AbstractString}(categorical_palette(axis_vector(daf, axis)))  # NOJET
    end
    if empty_color !== nothing
        color_palette[empty_value] = empty_color
    end
    return color_palette
end

# Define and export the `get_<property>_colors` of each line of the category `table`, which gets the palette of the
# line's `colors_axis` (see `get_category_colors`).
function define_category_colors(table::Tuple)::Nothing
    for line in table
        function_name = Symbol("get_$(line.property)_colors")
        docstring = """
                        $(function_name)(
                            daf::DafReader;
                            empty_value::AbstractString = "",
                            empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
                        )::OrderedDict{AbstractString, AbstractString}

                    The palette of the `$(line.colors_axis)` axis, which shows the `$(line.property)`. This is directed
                    at [`get_category_colors`](@ref).
                    """
        @eval begin
            @doc $docstring function $function_name(  # NOLINT
                daf::DafReader;
                empty_value::AbstractString = "",
                empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
            )::OrderedDict{AbstractString, AbstractString}
                return get_category_colors(daf; axis = $(line.colors_axis), empty_value, empty_color)
            end
            export $function_name
        end
    end
    return nothing
end

# The keywords of the category bucket's get and of its `put_category_values!`, other than the ones every generated data
# source takes. Each line of the table gives its own `colors_axis`, since its default (the property) isn't known to the
# generated function.
const CATEGORY_GET_KEYWORDS = (:(empty_value::AbstractString = ""),)
const CATEGORY_PUT_KEYWORDS = (
    :(empty_value::AbstractString = ""),
    :(colors_axis::AbstractString = ""),
    :(empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR),
    :(show_legend::Bool = true),
)

"""
    get_category_vector(
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        empty_value::AbstractString = "",
    )::AbstractVector{<:AbstractString}

Get the category `property` (e.g. the type) of each of the `entries` of the `daf` `axis`. By default this is the direct
property of the axis entries. However, by specifying `via`, it is possible to access indirect values (e.g., the type of
the block of the metacell of each cell). An entry whose `via` chain breaks gets the `empty_value`, which names no value.
"""
function get_category_vector(
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    empty_value::AbstractString = "",
)::AbstractVector{<:AbstractString}
    return get_axis_vector(daf; axis, query_suffix = ": $(escape_value(property))", via, empty_value, entries)
end

"""
    fill_category!(
        sinks::VectorDataSinks,
        daf::DafReader;
        axis::AbstractString,
        property::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
        title::Maybe{AbstractString} = property,
        empty_value::AbstractString = "",
        colors_axis::AbstractString = property,
        empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
        show_legend::Bool = true,
    )::Nothing

Fill the `sinks` with the category `property` of each of the `entries` of the `daf` `axis` (see
[`get_category_vector`](@ref)), shown by [`put_category_configuration!`](@ref), and name the entities after the
`entries`. The palette covers the `colors_axis` (see [`get_category_colors`](@ref)).
"""
function fill_category!(
    sinks::VectorDataSinks,
    daf::DafReader;
    axis::AbstractString,
    property::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    via::Maybe{Union{Tuple{Vararg{AbstractString}}, AbstractVector{<:AbstractString}}} = nothing,
    title::Maybe{AbstractString} = property,
    empty_value::AbstractString = "",
    colors_axis::AbstractString = property,
    empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
    show_legend::Bool = true,
)::Nothing
    entries = sinks_entries(sinks, entries)
    value_per_entry = get_category_vector(daf; axis, property, entries, via, empty_value)
    put_category_values!(
        sinks,
        daf,
        value_per_entry;
        axis,
        entries,
        title,
        empty_value,
        colors_axis,
        empty_color,
        show_legend,
    )
    return nothing
end

"""
    put_category_values!(
        sinks::VectorDataSinks,
        daf::DafReader,
        value_per_entry::AbstractVector{<:AbstractString};
        axis::AbstractString,
        entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
        title::Maybe{AbstractString} = nothing,
        empty_value::AbstractString = "",
        colors_axis::AbstractString,
        empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
        show_legend::Bool = true,
    )::Nothing

Put the category `value_per_entry` of each of the `entries` of the `daf` `axis` into the `sinks`, shown by
[`put_category_configuration!`](@ref) using the palette of the `colors_axis` (see [`get_category_colors`](@ref)), and
name the entities after the `entries`. Every category data source ends with this.
"""
function put_category_values!(
    sinks::VectorDataSinks,
    daf::DafReader,
    value_per_entry::AbstractVector{<:AbstractString};
    axis::AbstractString,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}} = nothing,
    title::Maybe{AbstractString} = nothing,
    empty_value::AbstractString = "",
    colors_axis::AbstractString,
    empty_color::Maybe{AbstractString} = EMPTY_CATEGORY_COLOR,
    show_legend::Bool = true,
)::Nothing
    put_vector_data!(sinks, value_per_entry; title)
    color_palette = get_category_colors(daf; axis = colors_axis, empty_value, empty_color)
    put_category_configuration!(sinks, color_palette; title, show_legend)
    fill_axis_names_data!(sinks, daf; axis, entries)
    return nothing
end

# The properties shown as categories.
const CATEGORY_PROPERTIES = (
    (;
        property = "type",
        axis = "metacell",
        title = "type",
        colors_axis = "type",
        doc = "The type of each entry, shown in the `color` of each type.",
    ),
    (;
        property = "block",
        axis = "metacell",
        title = "block",
        colors_axis = "block",
        show_legend = false,
        doc = """
              The block of each entry. Blocks are shown as a grouping of an axis, as colors, as a hover, or any of
              these. There are too many blocks for a legend to help, so `show_legend` is off.
              """,
    ),
)

define_bucket_data_sources(;
    get_function = :get_category_vector,
    put_values_function = :put_category_values!,
    get_keywords = CATEGORY_GET_KEYWORDS,
    put_keywords = CATEGORY_PUT_KEYWORDS,
    table = CATEGORY_PROPERTIES,
)

define_category_colors(CATEGORY_PROPERTIES)

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
    data = daf[query].array  # NOJET
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

# The `entries` a vector fill uses: the given ones, or else the names its `sinks` already hold. A `nothing` means all of
# them, which is left to the caller to resolve.
function sinks_entries(
    sinks::VectorDataSinks,
    entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}},
)::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}}
    if entries === nothing
        return get_vector_names_data(sinks)
    else
        return entries
    end
end

# The `row_entries` and `column_entries` a matrix fill uses, each chosen the way `sinks_entries` chooses them.
function sinks_matrix_entries(
    sinks::MatrixDataSinks,
    row_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}},
    column_entries::Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}},
)::Tuple{
    Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}},
    Maybe{Union{AbstractVector{<:AbstractString}, AbstractVector{<:Integer}}},
}
    if row_entries === nothing || column_entries === nothing
        name_per_row, name_per_column = get_matrix_names_data(sinks)
        if row_entries === nothing
            row_entries = name_per_row
        end
        if column_entries === nothing
            column_entries = name_per_column
        end
    end
    return (row_entries, column_entries)
end

end  # module
