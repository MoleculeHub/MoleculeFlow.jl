using Test
using MoleculeFlow

@testset "Bond Operations" begin
    mol = mol_from_smiles("CCO")

    # Test getting bonds from atom
    bonds = get_bonds_from_atom(mol, 1)  # First carbon
    @test length(bonds) >= 1
    @test all(bond -> isa(bond, Bond), bonds)

    if length(bonds) > 0
        bond = bonds[1]
        @test isa(get_begin_atom_idx(bond), Int)
        @test isa(get_end_atom_idx(bond), Int)
        @test isa(get_bond_type(bond), String)
        @test isa(is_aromatic(bond), Bool)
    end
end
