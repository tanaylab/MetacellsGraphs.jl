using Random
using TOML

using Base.JuliaSyntax: @K_str, children, kind, parseall, source_line, sourcetext, SyntaxNode

PACKAGE_NAME = TOML.parsefile("Project.toml")["name"]

# The module-level assignments of a file which are neither `const` nor typed (`name::Type = value`), as `(line_number,
# name)`. A function which uses such a variable looks it up at run time, and can't infer its type. So a type alias must
# be `const`, and a setting must be typed.
function untyped_module_variables(path::AbstractString)::Vector{Tuple{Int, String}}
    untyped = Tuple{Int, String}[]
    collect_untyped_module_variables!(untyped, parseall(SyntaxNode, read(path, String); filename = path))
    return untyped
end

function collect_untyped_module_variables!(untyped::Vector{Tuple{Int, String}}, node::SyntaxNode)::Nothing
    node_kind = kind(node)
    if node_kind in (K"toplevel", K"module", K"block")
        for child in children(node)
            collect_untyped_module_variables!(untyped, child)
        end
    elseif node_kind == K"doc"
        collect_untyped_module_variables!(untyped, children(node)[end])
    elseif node_kind == K"=" && kind(children(node)[1]) in (K"Identifier", K"curly")
        push!(untyped, (source_line(node), sourcetext(children(node)[1])))
    end
    return nothing
end

# The comprehensions of a file whose element type isn't static, as `(line_number, problem)`. An untyped comprehension
# (`[f(x) for x in xs]`) takes the join of its elements' types, so elements of two types give an abstract element type
# (e.g., `Integer` from `UInt16` and `Int`), and every access to it is untyped. An element type computed at run time
# (`eltype(xs)[f(x) for x in xs]`) makes every element be stored through dynamic dispatch, which boxes it. A static
# element type is a literal type or a type parameter: it has no calls, and its names are capitalized and aren't
# assigned in the file (except as `const`). A run-time element type is made static by a function barrier, a function
# whose parameter's type gives it as a type parameter.
function unstatic_comprehensions(path::AbstractString)::Vector{Tuple{Int, String}}
    unstatic = Tuple{Int, String}[]
    root = parseall(SyntaxNode, read(path, String); filename = path)
    assigned_names = Set{String}()
    collect_assigned_names!(assigned_names, root)
    collect_unstatic_comprehensions!(unstatic, assigned_names, root)
    return unstatic
end

function collect_assigned_names!(assigned_names::Set{String}, node::SyntaxNode)::Nothing
    node_children = children(node)
    if node_children === nothing || kind(node) == K"const"
        return nothing
    end
    if kind(node) == K"=" && kind(node_children[1]) == K"Identifier"
        push!(assigned_names, sourcetext(node_children[1]))
    end
    for child in node_children
        collect_assigned_names!(assigned_names, child)
    end
    return nothing
end

function collect_unstatic_comprehensions!(
    unstatic::Vector{Tuple{Int, String}},
    assigned_names::Set{String},
    node::SyntaxNode,
)::Nothing
    node_children = children(node)
    if node_children === nothing
        return nothing
    end
    text = first(replace(sourcetext(node), r"\s+" => " "), 80)
    if kind(node) == K"comprehension"
        push!(unstatic, (source_line(node), "Comprehension without an element type: $(text)"))
    elseif kind(node) == K"typed_comprehension" && !is_static_type(node_children[1], assigned_names)
        push!(unstatic, (source_line(node), "Comprehension element type is computed at run time: $(text)"))
    end
    for child in node_children
        collect_unstatic_comprehensions!(unstatic, assigned_names, child)
    end
    return nothing
end

function is_static_type(node::SyntaxNode, assigned_names::Set{String})::Bool
    node_kind = kind(node)
    if node_kind == K"Identifier"
        name = sourcetext(node)
        return isuppercase(first(name)) && !(name in assigned_names)
    elseif node_kind == K"Integer"
        return true
    elseif node_kind in (K"curly", K".", K"<:", K">:")
        return all(child -> is_static_type(child, assigned_names), children(node))
    else
        return false
    end
end

# The local variables which closures capture in a `Core.Box`, in the methods of the `modules` defined in the `paths`, as
# `(path, line_number, name)`. A variable is boxed when a closure captures it and it is assigned more than once. Each
# access of a boxed variable is untyped. The line is that of its first assignment in the method's body. A variable which
# is only assigned in closures is reported at the start of its scope.
function boxed_variables(modules::Set{Module}, paths::Set{String})::Vector{Tuple{String, Int, String}}
    boxed = Set{Tuple{String, Int, String}}()
    foreach_package_method(modules, paths) do method, path
        collect_boxed_variables!(boxed, method, path)
        return nothing
    end
    return sort!(collect(boxed))
end

# Call `visit(method, path)` for each method of the `modules` which has a body and is defined in the `paths`.
function foreach_package_method(visit::Function, modules::Set{Module}, paths::Set{String})::Nothing
    path_per_file = Dict{Symbol, String}()
    Base.visit(Core.methodtable) do method
        if method.module in modules && method.source !== nothing
            path = get!(path_per_file, method.file) do
                return isfile(String(method.file)) ? realpath(String(method.file)) : ""
            end
            if path in paths
                visit(method, path)
            end
        end
        return nothing
    end
    return nothing
end

# A boxed variable is created by `slot = Core.Box()` (or `Core.Box(argument)`) at the start of its scope. Each
# assignment to it is a `Core.setfield!(box, :contents, value)`, where the `box` is the slot or an SSA value loading it.
function collect_boxed_variables!(boxed::Set{Tuple{String, Int, String}}, method::Method, path::String)::Nothing
    code = Base.uncompressed_ir(method)
    line_per_slot = Dict{Int, Int}()
    for (position, statement) in enumerate(code.code)
        if statement isa Expr && statement.head == :(=) && is_call_of(statement.args[2], :Box)
            line_per_slot[statement.args[1].id] = statement_line(code, method, position, method.line)
        end
    end
    first_assignment_line_per_slot = Dict{Int, Int}()
    for (position, statement) in enumerate(code.code)
        if is_call_of(statement, :setfield!) && statement.args[3] == QuoteNode(:contents)
            slot = statement.args[2]
            if slot isa Core.SSAValue
                slot = code.code[slot.id]
            end
            if slot isa Core.SlotNumber && haskey(line_per_slot, slot.id)
                line_number = statement_line(code, method, position, line_per_slot[slot.id])
                first_assignment_line_per_slot[slot.id] =
                    min(line_number, get(first_assignment_line_per_slot, slot.id, line_number))
            end
        end
    end
    for (slot, line_number) in line_per_slot
        push!(boxed, (path, get(first_assignment_line_per_slot, slot, line_number), String(code.slotnames[slot])))
    end
    return nothing
end

function is_call_of(statement::Any, name::Symbol)::Bool
    return statement isa Expr && statement.head == :call && statement.args[1] == GlobalRef(Core, name)
end

function statement_line(code::Core.CodeInfo, method::Method, position::Integer, default_line::Integer)::Int
    for line_info in Base.IRShow.buildLineInfoNode(code.debuginfo, method, position)
        if line_info.file == method.file
            return Int(line_info.line)
        end
    end
    return Int(default_line)
end

# The uses of the global random generator in the methods of the `modules` defined in the `paths`, as `(path,
# line_number, problem)`. The global generator's state depends on everything that ran before it, so a result which
# depends on it isn't reproducible. Reported are:
#
#   - A call which doesn't pass a generator to a function which accepts one. A function accepts a generator if some
#     method of it takes an `rng` keyword, or a first positional argument of type `AbstractRNG`. This includes
#     broadcasts (`rand.(...)`) and splatted calls (`rand(arguments...)`, checked when the first argument is given
#     explicitly).
#   - Passing such a function as a value (`map(rand, ...)`), since its calls can't be checked.
#   - Calling `default_rng()`.
#
# A generator is a variable named `rng` or `..._rng`. Not seen: a function which uses the global generator internally
# without accepting one, a generator under another name or in another position, a callee held in a variable, the files
# under subdirectories of `src`, and `@generated` functions.
function global_random_generator_uses(modules::Set{Module}, paths::Set{String})::Vector{Tuple{String, Int, String}}
    uses = Set{Tuple{String, Int, String}}()
    accepts_per_callee = IdDict{Any, Tuple{Bool, Bool}}()
    foreach_package_method(modules, paths) do method, path
        # The methods which apply the keyword arguments' defaults are generated, and call `default_rng()` for them.
        signature = Base.unwrap_unionall(method.sig)
        is_keyword_call_method = signature isa DataType && signature.parameters[1] === typeof(Core.kwcall)
        if !is_keyword_call_method && isempty(Base.kwarg_decl(method))
            collect_global_random_generator_uses!(uses, accepts_per_callee, method, path)
        end
        return nothing
    end
    return sort!(collect(uses))
end

# Keywords forwarded by `kwargs...` are merged into a named tuple whose names aren't known, so such a call isn't
# checked. Nor is the plain call of the same callee which the lowering adds for when the forwarded keywords are empty.
# A keyword function's generated methods pass the defaults to its body (named `#name#...`), so a method which calls
# such a body isn't checked either.
function collect_global_random_generator_uses!(
    uses::Set{Tuple{String, Int, String}},
    accepts_per_callee::IdDict{Any, Tuple{Bool, Bool}},
    method::Method,
    path::String,
)::Nothing
    code = Base.uncompressed_ir(method)
    forwarded_callees = Set{Any}()
    for statement in code.code
        expression = unassigned_expression(statement)
        if !(expression isa Expr && expression.head == :call)
            continue
        end
        callee = resolved_callee(code, expression.args[1])
        if callee isa Function && startswith(String(nameof(callee)), "#")
            return nothing
        end
        if is_call_of(expression, :kwcall) && keyword_names_and_values(code, expression.args[2]) === nothing
            push!(forwarded_callees, resolved_callee(code, expression.args[3]))
        end
    end

    for (position, statement) in enumerate(code.code)
        statement = unassigned_expression(statement)
        if !(statement isa Expr && statement.head == :call)
            continue
        end
        callee_position, keywords, first_argument = call_parts(code, statement)
        callee = resolved_callee(code, statement.args[callee_position])
        line_number = statement_line(code, method, position, method.line)

        if callee === Random.default_rng
            push!(uses, (path, line_number, "Global random generator is used: default_rng"))
        elseif is_checked_callee(callee) && !(callee_position == 1 && callee in forwarded_callees)
            does_accept_keyword, does_accept_positional = cached_accepts_random_generator(accepts_per_callee, callee)
            if (does_accept_keyword || does_accept_positional) &&
               !passes_random_generator(code, keywords, first_argument, does_accept_positional)
                push!(uses, (path, line_number, "Call doesn't pass a random generator: $(nameof(callee))"))
            end
        end

        for argument in statement.args[(callee_position + 1):end]
            passed = resolved_callee(code, argument)
            if is_checked_callee(passed) && any(cached_accepts_random_generator(accepts_per_callee, passed))
                problem = "Function which accepts a random generator is passed as a value: $(nameof(passed))"
                push!(uses, (path, line_number, problem))
            end
        end
    end
    return nothing
end

# The parts of a call: the position of its callee in the call's arguments, its keywords' named tuple (`nothing` if it
# has none), and its first positional argument (`nothing` if it has none, `missing` if it isn't known). A keyword call
# is `Core.kwcall(NamedTuple{names}(Core.tuple(values...)), callee, arguments...)`. A splatted call is
# `Core._apply_iterate(iterate, callee, arguments...)`, whose arguments are collections. A broadcast is
# `Base.broadcasted(callee, arguments...)`.
function call_parts(code::Core.CodeInfo, statement::Expr)::Tuple{Int, Any, Any}
    if is_call_of(statement, :kwcall)
        return (3, statement.args[2], get(statement.args, 4, nothing))
    elseif is_call_of(statement, :_apply_iterate)
        first_collection = get(statement.args, 4, nothing)
        first_collection = first_collection isa Core.SSAValue ? code.code[first_collection.id] : first_collection
        if first_collection isa Tuple && !isempty(first_collection)
            return (3, nothing, first_collection[1])
        elseif is_call_of(first_collection, :tuple) && length(first_collection.args) >= 2
            return (3, nothing, first_collection.args[2])
        else
            return (3, nothing, missing)
        end
    elseif resolved_callee(code, statement.args[1]) === Base.broadcasted
        return (2, nothing, get(statement.args, 3, nothing))
    else
        return (1, nothing, get(statement.args, 2, nothing))
    end
end

# Whether a callee is a function or a data type (a constructor), and isn't the body of a keyword function.
function is_checked_callee(callee::Any)::Bool
    return (callee isa Function || (callee isa Type && Base.unwrap_unionall(callee) isa DataType)) &&
           !startswith(String(nameof(callee)), "#")
end

# Whether a call passes a random generator, as its `rng` keyword or (if the callee accepts it) its first positional
# argument. A call whose keywords' names or first argument aren't known counts as passing one.
function passes_random_generator(
    code::Core.CodeInfo,
    keywords::Any,
    first_argument::Any,
    does_accept_positional::Bool,
)::Bool
    if keywords !== nothing
        names_and_values = keyword_names_and_values(code, keywords)
        if names_and_values === nothing
            return true
        end
        names, values = names_and_values
        rng_position = findfirst(==(:rng), names)
        if rng_position !== nothing
            return is_random_generator(code, values[rng_position])
        end
    end
    if !does_accept_positional
        return false
    end
    return first_argument === missing || (first_argument !== nothing && is_random_generator(code, first_argument))
end

# A statement's expression, without the assignment of its value to a variable if it has one.
function unassigned_expression(statement::Any)::Any
    return statement isa Expr && statement.head == :(=) ? statement.args[2] : statement
end

function cached_accepts_random_generator(
    accepts_per_callee::IdDict{Any, Tuple{Bool, Bool}},
    callee::Union{Function, Type},
)::Tuple{Bool, Bool}
    return get!(accepts_per_callee, callee) do
        return accepts_random_generator(callee)
    end
end

# Whether some method of the `callee` takes an `rng` keyword, and whether some method takes a first positional argument
# of type `AbstractRNG`.
function accepts_random_generator(callee::Union{Function, Type})::Tuple{Bool, Bool}
    does_accept_keyword = false
    does_accept_positional = false
    for callee_method in methods(callee)
        does_accept_keyword |= :rng in Base.kwarg_decl(callee_method)
        signature = Base.unwrap_unionall(callee_method.sig)
        if signature isa DataType && length(signature.parameters) >= 2
            first_parameter = signature.parameters[2]
            if first_parameter isa TypeVar
                first_parameter = first_parameter.ub
            end
            does_accept_positional |= first_parameter === Random.AbstractRNG
        end
    end
    return (does_accept_keyword, does_accept_positional)
end

# The function a call's callee refers to: a global (`f`) or a module's global (`M.f`). Anything else is `nothing`.
function resolved_callee(code::Core.CodeInfo, value::Any)::Any
    if value isa Core.SSAValue
        value = code.code[value.id]
    end
    if value isa GlobalRef
        return isdefined(value.mod, value.name) ? getglobal(value.mod, value.name) : nothing
    end
    if value isa Expr && value.head == :call && value.args[1] == GlobalRef(Base, :getproperty)
        parent = resolved_callee(code, value.args[2])
        name = value.args[3]
        if parent isa Module && name isa QuoteNode && isdefined(parent, name.value)
            return getglobal(parent, name.value)
        end
    end
    return nothing
end

# The names and values of a keyword call's named tuple, or `nothing` if the names aren't a literal tuple.
function keyword_names_and_values(code::Core.CodeInfo, value::Any)::Union{Tuple{Tuple, Vector{Any}}, Nothing}
    construction = value isa Core.SSAValue ? code.code[value.id] : value
    if !(construction isa Expr && construction.head == :call && length(construction.args) == 2)
        return nothing
    end
    named_tuple_type = construction.args[1]
    named_tuple_type = named_tuple_type isa Core.SSAValue ? code.code[named_tuple_type.id] : named_tuple_type
    if !(is_call_of(named_tuple_type, :apply_type) && named_tuple_type.args[2] == GlobalRef(Core, :NamedTuple))
        return nothing
    end
    names = named_tuple_type.args[3]
    names = names isa Core.SSAValue ? code.code[names.id] : names
    names = names isa QuoteNode ? names.value : names
    values = construction.args[2]
    values = values isa Core.SSAValue ? code.code[values.id] : values
    if !(names isa Tuple && is_call_of(values, :tuple))
        return nothing
    end
    return (names, values.args[2:end])
end

# Whether a value is a variable named `rng` or `..._rng`: an argument or a local, or a captured variable of a closure.
function is_random_generator(code::Core.CodeInfo, value::Any)::Bool
    if value isa Core.SSAValue
        value = code.code[value.id]
    end
    if value isa Expr && value.head == :(=)
        value = value.args[1]
    end
    if value isa Core.SlotNumber
        return is_random_generator_name(code.slotnames[value.id])
    end
    if is_call_of(value, :getfield)
        field = value.args[3]
        field = field isa QuoteNode ? field.value : field
        return field isa Symbol && is_random_generator_name(field)
    end
    return false
end

function is_random_generator_name(name::Symbol)::Bool
    return name == :rng || endswith(String(name), "_rng")
end

# Known cases for the self-test. Each boxed variable is marked on the line of its first assignment, which is never the
# first statement of its scope (where the box is created), except for a variable which is only assigned in a closure.
module FlossFixtures

using Random

struct Shown end

function boxed_in_function(values::Vector{Int})::Vector{Int}
    n_values = length(values)
    if n_values > 1
        order = sortperm(values)  # BOXED order
        order = order[2:end]
    else
        order = [1]
    end
    return map(index -> values[order[index]], 1:length(order))
end

function boxed_argument(values::Vector{Int})::Vector{Int}
    n_values = length(values)
    if n_values > 1
        values = sort(values)  # BOXED values
    end
    return map(index -> values[index], 1:n_values)
end

function boxed_in_closure(values::Vector{Int})::Vector{Int}
    return map(1:2) do offset
        n_values = length(values)
        if offset > 1
            shifted = values .+ offset  # BOXED shifted
        else
            shifted = values
        end
        return sum(map(index -> shifted[index], 1:n_values))
    end
end

function assigned_only_in_closure(values::Vector{Int})::Int
    n_values = length(values)  # BOXED total
    local total
    foreach(values) do value
        total = value
        return nothing
    end
    return total + n_values
end

function Base.show(io::IO, ::Shown)::Nothing
    n_marks = 1
    if n_marks > 0
        text = "shown"  # BOXED text
    else
        text = ""
    end
    foreach(_ -> print(io, text), 1:n_marks)
    return nothing
end

function not_boxed(values::Vector{Int})::Vector{Int}
    sorted = sortperm(values)
    order = sorted[2:end]
    return map(index -> values[order[index]], 1:length(order))
end

function random_scale(value::Int; rng::AbstractRNG = Random.default_rng())::Float64
    return rand(rng) * value
end

function passes_rng(values::Vector{Int}; rng::AbstractRNG = Random.default_rng())::Vector{Float64}
    shuffled = shuffle(rng, values)
    return map(shuffled) do value
        return random_scale(value; rng)
    end
end

function drops_keyword_rng(value::Int; rng::AbstractRNG = Random.default_rng())::Float64
    return random_scale(value) + rand(rng)  # RANDOM random_scale
end

function drops_positional_rng(values::Vector{Int})::Vector{Int}
    return shuffle(values)  # RANDOM shuffle
end

function drops_rng_in_closure(values::Vector{Int}; rng::AbstractRNG = Random.default_rng())::Vector{Float64}
    return map(values) do value
        return rand() * random_scale(value; rng)  # RANDOM rand
    end
end

function passes_local_rng(value::Int)::Float64
    local_rng = Xoshiro(value)
    return random_scale(value; rng = local_rng)
end

function forwards_keywords(value::Int; keywords...)::Float64
    return random_scale(value; keywords...)
end

function drops_rng_in_splat(counts::Vector{Int})::Array{Float64}
    return rand(Float64, counts...)  # RANDOM rand
end

function passes_rng_in_splat(counts::Vector{Int}; rng::AbstractRNG = Random.default_rng())::Array{Float64}
    return rand(rng, Float64, counts...)
end

function drops_rng_in_broadcast(values::Vector{Int})::Vector{Int}
    return rand.(values)  # RANDOM rand
end

function passes_random_function(values::Vector{Vector{Int}})::Vector{Vector{Int}}
    return map(shuffle, values)  # RANDOM shuffle
end

function uses_global_rng(value::Int)::Float64
    global_rng = Random.default_rng()  # RANDOM default_rng
    return random_scale(value; rng = global_rng)
end

end  # module

# Verify the rules find exactly the known cases. The boxed variables are found through Julia internals, which a new
# Julia version may change. Without this, such a change would make the rule silently find nothing.
function self_test()::Nothing
    path = realpath(@__FILE__)
    expected_boxed = Tuple{Int, String}[]
    for (line_number, line_text) in enumerate(readlines(path))
        boxed_match = match(r"# BOXED (\w+)$", line_text)
        if boxed_match !== nothing
            push!(expected_boxed, (line_number, boxed_match.captures[1]))
        end
    end
    found_boxed = [
        (line_number, name) for
        (_, line_number, name) in boxed_variables(Set{Module}([FlossFixtures]), Set{String}([path]))
    ]
    if found_boxed != expected_boxed || length(expected_boxed) != 5
        error("floss self-test failed: boxed variables: $(found_boxed) expected: $(expected_boxed)")
    end

    expected_random = Tuple{Int, String}[]
    for (line_number, line_text) in enumerate(readlines(path))
        random_match = match(r"# RANDOM (\w+)$", line_text)
        if random_match !== nothing
            push!(expected_random, (line_number, random_match.captures[1]))
        end
    end
    found_random = Tuple{Int, String}[
        (line_number, String(last(split(problem)))) for
        (_, line_number, problem) in global_random_generator_uses(Set{Module}([FlossFixtures]), Set{String}([path]))
    ]
    if found_random != expected_random || length(expected_random) != 7
        error("floss self-test failed: global random generator uses: $(found_random) expected: $(expected_random)")
    end

    found_untyped = mktemp() do module_path, io
        write(io, "module M\nUNTYPED = 1\nconst ALIAS = Int\nTYPED::Int = 1\nend\n")
        close(io)
        return untyped_module_variables(module_path)
    end
    if found_untyped != [(2, "UNTYPED")]
        error("floss self-test failed: module-level variables: $(found_untyped) expected: [(2, \"UNTYPED\")]")
    end

    found_comprehensions = mktemp() do module_path, io
        write(
            io,
            """
            module M
            indices(dict, names) = [get(dict, name, 0) for name in names]
            typed_indices(dict, names) = UInt32[get(dict, name, 0) for name in names]
            total(values) = sum(value for value in values)
            pairs(values) = [(left, right) for left in values for right in values]
            copied(values) = eltype(values)[value for value in values]
            nested(values) = Vector{eltype(values)}[[value] for value in values]
            parametric(values::Vector{T}) where {T} = T[value for value in values]
            named(values) = (value_type = eltype(values); value_type[value for value in values])
            function assigned(values)
                Element = eltype(values)
                return Element[value for value in values]
            end
            bounded(sets) = AbstractSet{<:AbstractString}[set for set in sets]
            dotted(values) = Base.RefValue{Int}[Ref(value) for value in values]
            end
            """,
        )
        close(io)
        return Int[line_number for (line_number, _) in unstatic_comprehensions(module_path)]
    end
    if found_comprehensions != [2, 5, 6, 7, 9, 12]
        error("floss self-test failed: unstatic comprehensions: $(found_comprehensions) expected: [2, 5, 6, 7, 9, 12]")
    end

    return nothing
end

# The modules of the package: its own and all its submodules.
function package_modules(root_module::Module)::Set{Module}
    modules = Set{Module}([root_module])
    pending = [root_module]
    while !isempty(pending)
        parent_module = pop!(pending)
        for name in names(parent_module; all = true)
            if isdefined(parent_module, name)
                value = getfield(parent_module, name)
                if value isa Module && value !== parent_module && parentmodule(value) === parent_module
                    push!(modules, value)
                    push!(pending, value)
                end
            end
        end
    end
    return modules
end

function floss(root_module::Module)::Nothing
    relative_path_per_path = Dict{String, String}()
    for file_name in readdir("src")
        if endswith(file_name, ".jl")
            relative_path_per_path[realpath("src/$(file_name)")] = "src/$(file_name)"
        end
    end

    problems_per_path = Dict(path => Tuple{Int, String}[] for path in keys(relative_path_per_path))
    for (path, problems) in problems_per_path
        for (line_number, name) in untyped_module_variables(path)
            push!(problems, (line_number, "Module-level variable is neither const nor typed: $(name)"))
        end
        for (line_number, problem) in unstatic_comprehensions(path)
            push!(problems, (line_number, problem))
        end
    end
    modules = package_modules(root_module)
    paths = Set{String}(keys(relative_path_per_path))
    for (path, line_number, name) in boxed_variables(modules, paths)
        push!(problems_per_path[path], (line_number, "Captured variable is boxed: $(name)"))
    end
    for (path, line_number, problem) in global_random_generator_uses(modules, paths)
        push!(problems_per_path[path], (line_number, problem))
    end

    errors = 0
    skipped = 0
    unused = 0
    for path in sort!(collect(keys(problems_per_path)); by = path -> relative_path_per_path[path])
        relative_path = relative_path_per_path[path]
        no_floss_lines = Set{Int}()
        for (line_number, line_text) in enumerate(readlines(path))
            if contains(line_text, "NOFLOSS")
                push!(no_floss_lines, line_number)
            end
        end
        unused_no_floss_lines = copy(no_floss_lines)
        for (line_number, problem) in sort!(problems_per_path[path])
            if line_number in no_floss_lines
                delete!(unused_no_floss_lines, line_number)
                skipped += 1
            else
                println("$(relative_path):$(line_number): $(problem)")
                errors += 1
            end
        end
        for line_number in sort!(collect(unused_no_floss_lines))
            println("$(relative_path):$(line_number): unused NOFLOSS directive")
            unused += 1
        end
    end

    counts = String[]
    if errors > 0
        push!(counts, "$(errors) errors")
    end
    if skipped > 0
        push!(counts, "$(skipped) skipped")
    end
    if unused > 0
        push!(counts, "$(unused) unused")
    end
    if isempty(counts)
        println("Floss: clean!")
    else
        println("Floss: $(join(counts, ", "))")
    end

    if errors + unused > 0
        exit(1)
    end
    return nothing
end

self_test()
@eval using $(Symbol(PACKAGE_NAME))
floss(getfield(Main, Symbol(PACKAGE_NAME)))
