using Pkg
using MoleculeFlow
using Documenter

DocMeta.setdocmeta!(MoleculeFlow, :DocTestSetup, :(using MoleculeFlow); recursive = true)

makedocs(;
    modules = [MoleculeFlow],
    authors = "Renée Gil",
    repo = "https://github.com/MoleculeHub/MoleculeFlow.jl/blob/{commit}{path}#{line}",
    sitename = "MoleculeFlow.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://moleculehub.github.io/MoleculeFlow.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Practical Examples" => "examples.md",
        "API Reference" => [
            "Basic I/O" => "api/io.md",
            "Molecular Operations" => "api/operations.md",
            "Drawing & Visualization" => "api/drawing.md",
            "Molecular Descriptors" => "api/descriptors.md",
            "Fingerprints" => "api/fingerprints.md",
            "Substructure Search" => "api/substructure.md",
            "Atom Operations" => "api/atoms.md",
            "Bond Operations" => "api/bonds.md",
            "Similarity" => "api/similarity.md",
            "Standardization" => "api/standardization.md",
            "Conformers" => "api/conformers.md",
            "Fragmentation" => "api/fragmentation.md",
            "Graph Operations" => "api/graph.md",
            "Progress Tracking" => "api/progress.md",
        ],
    ],
    checkdocs = :none,  # Skip docstring checks
    warnonly = [:docs_block, :missing_docs],  # Only warn about missing docs
)

deploydocs(;
    repo = "github.com/MoleculeHub/MoleculeFlow.jl", devbranch = "main", push_preview = true
)
