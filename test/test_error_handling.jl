using Test
using MoleculeFlow

@testset "Error Handling" begin
    bad_mol = mol_from_smiles("invalid")
    @test ismissing(molecular_weight(bad_mol))
    @test ismissing(morgan_fingerprint(bad_mol))
    @test ismissing(get_atoms(bad_mol))

    good_mol = mol_from_smiles("CCO")
    @test ismissing(tanimoto_similarity(good_mol, bad_mol))
    @test ismissing(tanimoto_similarity(bad_mol, good_mol))
end
