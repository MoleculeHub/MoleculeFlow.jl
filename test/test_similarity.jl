using Test
using MoleculeFlow

@testset "Similarity Calculations" begin
    mol1 = mol_from_smiles("CCO")  # ethanol
    mol2 = mol_from_smiles("CCC")  # propane
    mol3 = mol_from_smiles("CCO")  # ethanol again

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
end
