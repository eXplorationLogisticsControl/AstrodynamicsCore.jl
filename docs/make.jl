import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

const SRC = normpath(joinpath(@__DIR__, "..", "src"))
include(joinpath(SRC, "AstrodynamicsCore.jl"))

using Documenter

makedocs(
    modules = [AstrodynamicsCore],
    sitename = "AstrodynamicsCore.jl",
    authors = "Yuri <yuri.shimane@gmail.com>",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = nothing,
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => [
            "Orbital Elements" => "tutorials/elements.md",
            "Mean Elements (J2)" => "tutorials/mean_elements.md",
            "Kepler Problem" => "tutorials/kepler.md",
            "Lambert Problem" => "tutorials/lambert.md",
        ],
        "API" => [
            "Orbital Elements API" => "api/elements.md",
            "Anomaly Conversions" => "api/anomalies.md",
            "Propagation API" => "api/propagation.md",
            "Mean Elements (J2) API" => "api/mean_elements.md",
            "Planets & Transfers" => "api/planet_transfer.md",
            "Utilities" => "api/utilities.md",
        ],
    ],
    checkdocs = :exports,
    warnonly = [:cross_references, :missing_docs],
)
