using Test
using MoleculeFlow

@testset "All Descriptors" begin
    @testset "calc_all_descriptors" begin
        mol = mol_from_smiles("CCO")
        descriptors = calc_all_descriptors(mol)

        @test isa(descriptors, Dict{Symbol, Any})
        @test length(descriptors) > 200  

        @test haskey(descriptors, :MolWt)
        @test haskey(descriptors, :TPSA)
        @test haskey(descriptors, :NumHDonors)
        @test haskey(descriptors, :NumHAcceptors)
        @test haskey(descriptors, :MolLogP)

        @test abs(descriptors[:MolWt] - 46.069) < 0.01
        @test abs(descriptors[:TPSA] - 20.23) < 0.01
        @test descriptors[:NumHDonors] == 1
        @test descriptors[:NumHAcceptors] == 1

        @test isa(descriptors[:MolWt], Float64)
        @test isa(descriptors[:NumHDonors], Real)  # Could be Int or Float64
        @test isa(descriptors[:NumHAcceptors], Real)
    end

    @testset "calc_all_descriptors with invalid molecule" begin
        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid_smiles")
        descriptors = calc_all_descriptors(bad_mol)
        @test descriptors === missing
    end

    @testset "calc_all_descriptors vectorized" begin
        smiles_list = ["CCO", "CCC", "c1ccccc1"]
        mols = mol_from_smiles(smiles_list)
        all_descriptors = calc_all_descriptors(mols)

        @test length(all_descriptors) == 3
        @test isa(all_descriptors[1], Dict{Symbol, Any})
        @test isa(all_descriptors[2], Dict{Symbol, Any})
        @test isa(all_descriptors[3], Dict{Symbol, Any})

        # Check that different molecules have different values
        @test all_descriptors[1][:MolWt] != all_descriptors[2][:MolWt]
        @test all_descriptors[1][:MolWt] != all_descriptors[3][:MolWt]

        # Test with mixed valid/invalid molecules
        mixed_smiles = ["CCO", "invalid", "CCC"]
        mixed_mols = mol_from_smiles(mixed_smiles)
        mixed_descriptors = calc_all_descriptors(mixed_mols)

        @test length(mixed_descriptors) == 3
        @test isa(mixed_descriptors[1], Dict{Symbol, Any})
        @test mixed_descriptors[2] === missing
        @test isa(mixed_descriptors[3], Dict{Symbol, Any})
    end
end