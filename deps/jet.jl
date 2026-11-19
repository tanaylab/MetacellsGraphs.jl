using JET
using TOML

# The package in the current directory comes first, so JET checks it and not another copy in the default environment.
pushfirst!(LOAD_PATH, ".")

PACKAGE_NAME = TOML.parsefile("Project.toml")["name"]
@eval using $(Symbol(PACKAGE_NAME))

println(report_package(getfield(Main, Symbol(PACKAGE_NAME))))
