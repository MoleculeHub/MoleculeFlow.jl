using Test
using MoleculeFlow

@testset "Basic Molecule Operations" begin
    # Test molecule creation from SMILES
    mol = mol_from_smiles("CCCO")  # ethanol
    @test mol.valid == true
    @test mol.source == "CCCO"

    # Test invalid SMILES
    bad_mol = mol_from_smiles("invalid_smiles")
    @test bad_mol.valid == false

    # Test SMILES conversion
    smiles = mol_to_smiles(mol)
    @test isa(smiles, String)
    @test length(smiles) > 0

    # Test kekule_smiles
    benzene = mol_from_smiles("c1ccccc1")
    @test mol_to_smiles(benzene; kekule_smiles = true) == "C1=CC=CC=C1"
    @test mol_to_smiles(benzene; kekule_smiles = false) == "c1ccccc1"

    # Test all_bonds_explicit
    @test mol_to_smiles(benzene; kekule_smiles = true, all_bonds_explicit = true) ==
        "C1=C-C=C-C=C-1"

    # Test InChI conversion  
    inchi = mol_to_inchi(mol)
    @test isa(inchi, String)
    @test startswith(inchi, "InChI=")

    # Test dynamic properties
    mol.logP = 0.25
    @test mol.logP == 0.25

    mol.custom_prop = "test_value"
    @test mol.custom_prop == "test_value"
end

@testset "Vectorized Operations" begin
    smiles_list = ["CCO", "CCC", "invalid", "CC(C)O"]
    mols = mol_from_smiles(smiles_list)

    @test length(mols) == 4
    @test mols[1].valid == true
    @test mols[2].valid == true
    @test mols[3].valid == false  # invalid SMILES
    @test mols[4].valid == true

    # Test conversion back to SMILES
    smiles_back = mol_to_smiles(mols)
    @test length(smiles_back) == 4
    @test !ismissing(smiles_back[1])
    @test !ismissing(smiles_back[2])
    @test ismissing(smiles_back[3])  # invalid molecule
    @test !ismissing(smiles_back[4])
end
