using Test
using MoleculeFlow

@testset "Fingerprints" begin
    mol = mol_from_smiles("CCO")

    morgan_fp = morgan_fingerprint(mol)
    @test isa(morgan_fp, Vector{Bool})
    @test length(morgan_fp) == 2048

    rdk_fp = rdk_fingerprint(mol)
    @test isa(rdk_fp, Vector{Bool})
    @test length(rdk_fp) == 2048

    maccs_fp = maccs_fingerprint(mol)
    @test isa(maccs_fp, Vector{Bool})
    @test length(maccs_fp) == 167

    morgan_fp_256 = morgan_fingerprint(mol; nbits = 256)
    @test length(morgan_fp_256) == 256
end
