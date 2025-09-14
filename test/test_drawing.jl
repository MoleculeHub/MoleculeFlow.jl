using Test
using MoleculeFlow

@testset "Drawing Functions" begin
    mol = mol_from_smiles("CCO")

    img = mol_to_image(mol)
    @test !isnothing(img)

    mols = mol_from_smiles(["CCO", "CCC", "CC(C)O"])
    grid_img = mols_to_grid_image(mols)
    @test !isnothing(grid_img)
end
