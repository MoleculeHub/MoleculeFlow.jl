using Test
using MoleculeFlow

@testset "Similarity Calculations" begin
    mol1 = mol_from_smiles("CCO")  # ethanol
    mol2 = mol_from_smiles("CCC")  # propane
    mol3 = mol_from_smiles("CCO")  # ethanol again

    @test mol1.valid
    @test mol2.valid
    @test mol3.valid

    # Test Tanimoto similarity
    sim_different = tanimoto_similarity(mol1, mol2)
    sim_identical = tanimoto_similarity(mol1, mol3)

    @test isa(sim_different, Float64)
    @test isa(sim_identical, Float64)
    @test sim_identical > sim_different
    @test sim_identical ≈ 1.0 atol = 0.01

    # Test other similarity metrics
    @test isa(dice_similarity(mol1, mol2), Float64)
    @test isa(cosine_similarity(mol1, mol2), Float64)

    # Test bulk similarity
    target_mols = [mol1, mol2, mol3]
    similarities = bulk_similarity(mol1, target_mols)
    @test length(similarities) == 3
    @test similarities[1] ≈ 1.0 atol = 0.01  # identical to self
    @test similarities[3] ≈ 1.0 atol = 0.01  # identical to mol3

    # Test with invalid molecule
    invalid_mol = mol_from_smiles("invalid_smiles")
    @test !invalid_mol.valid
    @test tanimoto_similarity(mol1, invalid_mol) === missing
    @test dice_similarity(mol1, invalid_mol) === missing
    @test cosine_similarity(mol1, invalid_mol) === missing
    @test all(ismissing, bulk_similarity(invalid_mol, [mol1, mol2]))
end
