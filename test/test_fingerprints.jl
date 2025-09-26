using Test
using MoleculeFlow

@testset "Fingerprints" begin
    mol = mol_from_smiles("CCO")
    @test mol.valid

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

    # Test with invalid molecule
    invalid_mol = mol_from_smiles("invalid_smiles")
    @test !invalid_mol.valid
    @test morgan_fingerprint(invalid_mol) === missing
    @test rdk_fingerprint(invalid_mol) === missing
    @test maccs_fingerprint(invalid_mol) === missing
end
