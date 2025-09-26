using Test
using MoleculeFlow

@testset "Atom Operations" begin
    mol = mol_from_smiles("CCO")
    @test mol.valid

    # Test getting atoms
    atoms = get_atoms(mol)
    @test length(atoms) == 3  # C, C, O
    @test all(atom -> isa(atom, Atom), atoms)

    # Test individual atom access
    atom1 = get_atom(mol, 1)
    @test isa(atom1, Atom)

    # Test atom descriptors
    @test get_symbol(atom1) == "C"
    @test get_atomic_number(atom1) == 6
    @test isa(get_degree(atom1), Int)
    @test isa(get_formal_charge(atom1), Int)

    # Test oxygen atom (should be at index 3)
    oxygen = get_atom(mol, 3)
    @test get_symbol(oxygen) == "O"
    @test get_atomic_number(oxygen) == 8

    # Test with invalid molecule
    invalid_mol = mol_from_smiles("invalid_smiles")
    @test !invalid_mol.valid
    @test get_atoms(invalid_mol) === missing
    @test get_atom(invalid_mol, 1) === missing
end
