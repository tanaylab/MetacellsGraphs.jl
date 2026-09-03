push!(LOAD_PATH, ".")

using Aqua
using MetacellsGraphs
Aqua.test_ambiguities([MetacellsGraphs])
Aqua.test_all(MetacellsGraphs; ambiguities = false, unbound_args = false, deps_compat = false)
