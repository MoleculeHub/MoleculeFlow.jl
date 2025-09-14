using Test
using MoleculeFlow
using Graphs

@testset "MoleculeFlow.jl Tests" begin
    include("test_molecules.jl")
    include("test_descriptors.jl")
    include("test_fingerprints.jl")
    include("test_similarity.jl")
    include("test_atoms.jl")
    include("test_bonds.jl")
    include("test_substructure.jl")
    include("test_graph.jl")
    include("test_drawing.jl")
    include("test_standardization.jl")
    include("test_conformers.jl")
    include("test_error_handling.jl")
    include("test_aqua.jl")
end
