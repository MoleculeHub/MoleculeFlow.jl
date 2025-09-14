using Test
using MoleculeFlow

@testset "Molecular Descriptors" begin
    mol = mol_from_smiles("CC(C)O")

    @test molecular_weight(mol) ≈ 60.096 atol = 0.01
    @test heavy_atom_count(mol) == 4
    @test num_hbd(mol) == 1  # OH group
    @test num_hba(mol) == 1  # O atom
    @test num_rotatable_bonds(mol) == 0  # Isopropanol has no rotatable bonds by RDKit definition

    benzene = mol_from_smiles("c1ccccc1")
    @test num_rings(benzene) == 1
    @test num_aromatic_rings(benzene) == 1

    mols = mol_from_smiles(["CCO", "CCC", "CC(C)O"])
    weights = molecular_weight(mols)
    @test length(weights) == 3
    @test all(w -> isa(w, Float64), weights)
end
