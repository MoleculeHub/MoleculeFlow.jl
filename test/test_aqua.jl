using Test
using MoleculeFlow
using Aqua

@testset "Aqua.jl Quality Checks" begin
    Aqua.test_all(MoleculeFlow)
end