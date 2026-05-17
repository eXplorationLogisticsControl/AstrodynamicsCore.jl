import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))
if !any(d -> d.name == "AstrodynamicsCore", values(Pkg.dependencies()))
    Pkg.develop(Pkg.PackageSpec(path=PKG_ROOT))
end

using Documenter
using AstrodynamicsCore

makedocs(
    modules = [AstrodynamicsCore],
    sitename = "AstrodynamicsCore.jl",
    authors = "Yuri <yuri.shimane@gmail.com>",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "API" => [
            "Orbital Elements" => "api/elements.md",
            "Anomaly Conversions" => "api/anomalies.md",
            "Propagation" => "api/propagation.md",
            "Mean Elements (J2)" => "api/mean_elements.md",
            "Planets & Transfers" => "api/planet_transfer.md",
            "Utilities" => "api/utilities.md",
        ],
    ],
    checkdocs = :exports,
    warnonly = [:cross_references, :missing_docs],
)

if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    deploydocs(
        repo = "github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl.git",
        devbranch = "main",
        push = true,
    )
end
