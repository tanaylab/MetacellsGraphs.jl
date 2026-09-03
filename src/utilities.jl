# Helpers shared by more than one kind of graph. Nothing here is exposed to the user of the package.
module Utilities

export type_colors

using DataAxesFormats
using SomeGraphs
using TanayLabUtilities

# The type of each entry of an axis, and the configuration for coloring by it. Either is `nothing` when the repository
# does not have it: the type is optional, and so are the colors of the type axis. Without the colors there is nothing
# to map a type name to, so a graph which can only show a type as a color shows nothing of a repository which has types
# but no colors for them; a graph which also groups by the type still has the types to group by.
function type_colors(
    daf::DafReader,
    axis::AbstractString;
    show_legend::Bool = false,
)::Tuple{Maybe{AbstractVector{<:AbstractString}}, Maybe{ColorsConfiguration}}
    type_per_entry = get_vector(daf, axis, "type"; default = nothing)
    if type_per_entry === nothing
        return (nothing, nothing)
    elseif !has_axis(daf, "type") || !has_vector(daf, "type", "color")
        return (type_per_entry.array, nothing)
    else
        return (type_per_entry.array, ColorsConfiguration(; palette = get_vector(daf, "type", "color"), show_legend))
    end
end

end  # module
