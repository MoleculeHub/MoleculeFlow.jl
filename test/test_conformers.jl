using Test
using MoleculeFlow

@testset "Conformer Generation" begin
    @testset "Single 3D Conformer Generation" begin
        mol = mol_from_smiles("CCO")  # ethanol
        @test mol.valid == true

        conformers = generate_3d_conformers(mol, 1; optimize = true, force_field = :mmff)
        @test !isempty(conformers)
        @test length(conformers) == 1
        @test conformers[1] isa ConformerMolecule
        @test conformers[1].conformer_result.optimized == true
        @test conformers[1].conformer_result.energy < Inf
        @test haskey(conformers[1].molecule.props, :coordinates_3d)

        mol2 = mol_from_smiles("CCO")
        conformers_unopt = generate_3d_conformers(
            mol2, 1; optimize = false, force_field = :mmff
        )
        @test !isempty(conformers_unopt)
        @test conformers_unopt[1].conformer_result.optimized == false
        @test haskey(conformers_unopt[1].molecule.props, :coordinates_3d)
    end

    @testset "Multiple 3D Conformer Generation" begin
        mol = mol_from_smiles("CCCCCCCC")
        @test mol.valid == true

        conformers = generate_3d_conformers(
            mol, 5; optimize = true, force_field = :mmff, random_seed = 42
        )
        @test !isempty(conformers)
        @test length(conformers) <= 5  # May be fewer due to pruning

        if length(conformers) > 1
            for i in eachindex(conformers)
                @test conformers[i].conformer_result.energy >=
                    conformers[i - 1].conformer_result.energy
            end
        end

        for conf in conformers
            @test conf.conformer_result.energy < Inf
            @test conf.conformer_result.optimized == true
            @test haskey(conf.molecule.props, :coordinates_3d)
            coords = conf.molecule.props[:coordinates_3d]
            @test size(coords, 2) == 3  # Should have x, y, z coordinates
        end
    end

    @testset "Force Field Options" begin
        mol = mol_from_smiles("c1ccccc1")
        @test mol.valid == true

        conformers_mmff = generate_3d_conformers(
            mol, 1; optimize = true, force_field = :mmff
        )
        @test !isempty(conformers_mmff)
        @test conformers_mmff[1].conformer_result.energy < Inf

        mol2 = mol_from_smiles("c1ccccc1")
        conformers_uff = generate_3d_conformers(
            mol2, 1; optimize = true, force_field = :uff
        )
        @test !isempty(conformers_uff)
        @test conformers_uff[1].conformer_result.energy < Inf
    end

    @testset "Error Handling" begin
        bad_mol = mol_from_smiles("invalid_smiles")
        @test bad_mol.valid == false

        conformers_3d = generate_3d_conformers(bad_mol, 1)
        @test isempty(conformers_3d)
    end

    @testset "2D Conformer Generation" begin
        mol = mol_from_smiles("CCO")
        @test mol.valid == true

        conformers_2d = generate_2d_conformers(mol)
        @test !isempty(conformers_2d)
        @test length(conformers_2d) == 1
        @test conformers_2d[1] isa ConformerMolecule
        @test haskey(conformers_2d[1].molecule.props, :coordinates_2d)
        @test haskey(conformers_2d[1].molecule.props, :conformer_type)
        @test conformers_2d[1].molecule.props[:conformer_type] == "2d"

        coords_2d = conformers_2d[1].molecule.props[:coordinates_2d]
        @test coords_2d isa Matrix{Float64}
        @test size(coords_2d, 2) == 2
        @test size(coords_2d, 1) > 0   # Should have atoms

        benzene = mol_from_smiles("c1ccccc1")
        conformers_benzene = generate_2d_conformers(benzene)
        @test !isempty(conformers_benzene)
        coords_benzene = conformers_benzene[1].molecule.props[:coordinates_2d]
        @test size(coords_benzene, 1) == 6  # 6 carbon atoms
        @test size(coords_benzene, 2) == 2  # x, y coordinates

        bad_mol = mol_from_smiles("invalid_smiles")
        @test bad_mol.valid == false
        conformers_2d_bad = generate_2d_conformers(bad_mol)
        @test isempty(conformers_2d_bad)
    end

    @testset "Display Functions" begin
        mol = mol_from_smiles("CC")
        conformers = generate_3d_conformers(mol, 1)

        if !isempty(conformers)
            conf_result_str = string(conformers[1].conformer_result)
            @test occursin("ConformerResult", conf_result_str)
            @test occursin("energy=", conf_result_str)
            @test occursin("kcal/mol", conf_result_str)

            conf_mol_str = string(conformers[1])
            @test occursin("ConformerMolecule", conf_mol_str)
            @test occursin("ConformerResult", conf_mol_str)
        end
    end
end
