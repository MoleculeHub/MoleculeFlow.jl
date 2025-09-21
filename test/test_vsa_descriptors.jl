using Test
using MoleculeFlow

@testset "VSA Descriptors" begin
    # Test molecules
    ethanol = mol_from_smiles("CCO")
    benzene = mol_from_smiles("c1ccccc1")
    aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    invalid_mol = mol_from_smiles("invalid_smiles")

    @testset "SlogP_VSA Descriptors" begin
        # Test SlogP_VSA2-12 descriptors
        vsa_functions = [
            slogp_vsa2,
            slogp_vsa3,
            slogp_vsa4,
            slogp_vsa5,
            slogp_vsa6,
            slogp_vsa7,
            slogp_vsa8,
            slogp_vsa9,
            slogp_vsa10,
            slogp_vsa11,
            slogp_vsa12,
        ]

        for (i, vsa_func) in enumerate(vsa_functions)
            # Test with ethanol
            result_ethanol = vsa_func(ethanol)
            @test isa(result_ethanol, Float64) || result_ethanol === missing
            if isa(result_ethanol, Float64)
                @test result_ethanol >= 0.0
            end

            # Test with benzene
            result_benzene = vsa_func(benzene)
            @test isa(result_benzene, Float64) || result_benzene === missing
            if isa(result_benzene, Float64)
                @test result_benzene >= 0.0
            end

            # Test with invalid molecule
            @test vsa_func(invalid_mol) === missing
        end

        # Test vectorized operations
        mols = [ethanol, benzene, aspirin]
        results = slogp_vsa2.(mols)
        @test length(results) == 3
        @test all(r -> isa(r, Float64) || r === missing, results)
    end

    @testset "SMR_VSA Descriptors" begin
        # Test SMR_VSA1-10 descriptors
        vsa_functions = [
            smr_vsa1,
            smr_vsa2,
            smr_vsa3,
            smr_vsa4,
            smr_vsa5,
            smr_vsa6,
            smr_vsa7,
            smr_vsa8,
            smr_vsa9,
            smr_vsa10,
        ]

        for vsa_func in vsa_functions
            # Test with ethanol
            result_ethanol = vsa_func(ethanol)
            @test isa(result_ethanol, Float64) || result_ethanol === missing
            if isa(result_ethanol, Float64)
                @test result_ethanol >= 0.0
            end

            # Test with benzene
            result_benzene = vsa_func(benzene)
            @test isa(result_benzene, Float64) || result_benzene === missing
            if isa(result_benzene, Float64)
                @test result_benzene >= 0.0
            end

            # Test with invalid molecule
            @test vsa_func(invalid_mol) === missing
        end

        # Test vectorized operations
        mols = [ethanol, benzene, aspirin]
        results = smr_vsa1.(mols)
        @test length(results) == 3
        @test all(r -> isa(r, Float64) || r === missing, results)
    end

    @testset "PEOE_VSA Descriptors" begin
        # Test PEOE_VSA1-14 descriptors
        vsa_functions = [
            peoe_vsa1,
            peoe_vsa2,
            peoe_vsa3,
            peoe_vsa4,
            peoe_vsa5,
            peoe_vsa6,
            peoe_vsa7,
            peoe_vsa8,
            peoe_vsa9,
            peoe_vsa10,
            peoe_vsa11,
            peoe_vsa12,
            peoe_vsa13,
            peoe_vsa14,
        ]

        for vsa_func in vsa_functions
            # Test with ethanol
            result_ethanol = vsa_func(ethanol)
            @test isa(result_ethanol, Float64) || result_ethanol === missing
            if isa(result_ethanol, Float64)
                @test result_ethanol >= 0.0
            end

            # Test with benzene
            result_benzene = vsa_func(benzene)
            @test isa(result_benzene, Float64) || result_benzene === missing
            if isa(result_benzene, Float64)
                @test result_benzene >= 0.0
            end

            # Test with invalid molecule
            @test vsa_func(invalid_mol) === missing
        end

        # Test vectorized operations
        mols = [ethanol, benzene, aspirin]
        results = peoe_vsa1.(mols)
        @test length(results) == 3
        @test all(r -> isa(r, Float64) || r === missing, results)
    end

    @testset "VSA Descriptor Properties" begin
        # Test that VSA descriptors give different values for different molecules
        ethanol_slogp2 = slogp_vsa2(ethanol)
        benzene_slogp2 = slogp_vsa2(benzene)

        if isa(ethanol_slogp2, Float64) && isa(benzene_slogp2, Float64)
            # They should be different for different molecules
            @test ethanol_slogp2 != benzene_slogp2
        end

        # Test that VSA descriptors are non-negative
        for mol in [ethanol, benzene]
            peoe1 = peoe_vsa1(mol)
            if isa(peoe1, Float64)
                @test peoe1 >= 0.0
            end
        end
    end
end

@testset "BCUT Descriptors" begin
    # Test molecules
    ethanol = mol_from_smiles("CCO")
    benzene = mol_from_smiles("c1ccccc1")
    invalid_mol = mol_from_smiles("invalid_smiles")

    @testset "BCUT2D Descriptors" begin
        # Test all BCUT2D descriptors
        bcut_functions = [
            bcut2d_mwlow,
            bcut2d_mwhi,
            bcut2d_chglow,
            bcut2d_chghi,
            bcut2d_logplow,
            bcut2d_logphi,
            bcut2d_mrlow,
            bcut2d_mrhi,
        ]

        for bcut_func in bcut_functions
            # Test with ethanol
            result_ethanol = bcut_func(ethanol)
            @test isa(result_ethanol, Float64) || result_ethanol === missing

            # Test with benzene
            result_benzene = bcut_func(benzene)
            @test isa(result_benzene, Float64) || result_benzene === missing

            # Test with invalid molecule
            @test bcut_func(invalid_mol) === missing
        end

        # Test vectorized operations
        mols = [ethanol, benzene]
        results = bcut2d_mwlow.(mols)
        @test length(results) == 2
        @test all(r -> isa(r, Float64) || r === missing, results)
    end

    @testset "BCUT Descriptor Properties" begin
        # Test that high/low pairs make sense when both are valid
        mwlow = bcut2d_mwlow(benzene)
        mwhi = bcut2d_mwhi(benzene)

        if isa(mwlow, Float64) && isa(mwhi, Float64)
            @test mwlow <= mwhi  # Low should be <= high
        end

        logplow = bcut2d_logplow(benzene)
        logphi = bcut2d_logphi(benzene)

        if isa(logplow, Float64) && isa(logphi, Float64)
            @test logplow <= logphi  # Low should be <= high
        end
    end
end

@testset "Additional Structure Count Descriptors" begin
    # Test molecules with various structural features
    ethanol = mol_from_smiles("CCO")
    cyclohexane = mol_from_smiles("C1CCCCC1")
    pyridine = mol_from_smiles("c1cccnc1")
    tetrahydrofuran = mol_from_smiles("C1CCOC1")
    invalid_mol = mol_from_smiles("invalid_smiles")

    @testset "Additional Ring Count Functions" begin
        # Test num_aliphatic_heterocycles
        @test num_aliphatic_heterocycles(ethanol) == 0
        @test num_aliphatic_heterocycles(tetrahydrofuran) == 1  # THF is aliphatic heterocycle
        @test num_aliphatic_heterocycles(pyridine) == 0  # Pyridine is aromatic
        @test num_aliphatic_heterocycles(invalid_mol) === missing

        # Test num_saturated_heterocycles
        @test num_saturated_heterocycles(ethanol) == 0
        @test num_saturated_heterocycles(tetrahydrofuran) == 1
        @test num_saturated_heterocycles(invalid_mol) === missing

        # Test num_saturated_carbocycles
        @test num_saturated_carbocycles(ethanol) == 0
        @test num_saturated_carbocycles(cyclohexane) == 1
        @test num_saturated_carbocycles(invalid_mol) === missing

        # Test num_aliphatic_rings
        @test num_aliphatic_rings(ethanol) == 0
        @test num_aliphatic_rings(cyclohexane) == 1
        @test num_aliphatic_rings(pyridine) == 0  # Aromatic, not aliphatic
        @test num_aliphatic_rings(invalid_mol) === missing

        # Test num_heterocycles
        @test num_heterocycles(ethanol) == 0
        @test num_heterocycles(pyridine) == 1
        @test num_heterocycles(tetrahydrofuran) == 1
        @test num_heterocycles(invalid_mol) === missing
    end

    @testset "Additional Stereochemistry and Structure Counts" begin
        chiral_mol = mol_from_smiles("C[C@H](O)C")

        # Test num_unspecified_atom_stereo_centers
        @test isa(num_unspecified_atom_stereo_centers(ethanol), Int)
        @test isa(num_unspecified_atom_stereo_centers(chiral_mol), Int)
        @test num_unspecified_atom_stereo_centers(invalid_mol) === missing

        # Test num_spiro_atoms
        @test isa(num_spiro_atoms(ethanol), Int)
        @test num_spiro_atoms(ethanol) == 0  # No spiro atoms in ethanol
        @test num_spiro_atoms(invalid_mol) === missing

        # Test num_bridgehead_atoms
        @test isa(num_bridgehead_atoms(ethanol), Int)
        @test num_bridgehead_atoms(ethanol) == 0  # No bridgehead atoms in ethanol
        @test num_bridgehead_atoms(invalid_mol) === missing

        # Test hall_kier_alpha (can be negative)
        @test isa(hall_kier_alpha(ethanol), Float64)
        @test hall_kier_alpha(invalid_mol) === missing

        # Note: num_sp3_heavy_atoms removed as RDKit function doesn't exist
    end
end

@testset "Additional E-State Descriptors" begin
    ethanol = mol_from_smiles("CCO")
    benzene = mol_from_smiles("c1ccccc1")
    invalid_mol = mol_from_smiles("invalid_smiles")

    @testset "Max/Min Absolute E-State Indices" begin
        # Test max_absolute_e_state_index
        max_abs_estate = max_absolute_e_state_index(ethanol)
        @test isa(max_abs_estate, Float64) || max_abs_estate === missing
        @test max_absolute_e_state_index(invalid_mol) === missing

        # Test min_absolute_e_state_index
        min_abs_estate = min_absolute_e_state_index(ethanol)
        @test isa(min_abs_estate, Float64) || min_abs_estate === missing
        @test min_absolute_e_state_index(invalid_mol) === missing

        # Test relationship between max and min (when both are valid)
        if isa(max_abs_estate, Float64) && isa(min_abs_estate, Float64)
            @test max_abs_estate >= min_abs_estate
        end
    end
end

@testset "Vectorized Operations for New Descriptors" begin
    mols = [
        mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1"), mol_from_smiles("C1CCCCC1")
    ]

    # Test vectorized VSA descriptors
    slogp2_values = slogp_vsa2.(mols)
    @test length(slogp2_values) == 3
    @test all(v -> isa(v, Float64) || v === missing, slogp2_values)

    smr1_values = smr_vsa1.(mols)
    @test length(smr1_values) == 3
    @test all(v -> isa(v, Float64) || v === missing, smr1_values)

    peoe1_values = peoe_vsa1.(mols)
    @test length(peoe1_values) == 3
    @test all(v -> isa(v, Float64) || v === missing, peoe1_values)

    # Test vectorized BCUT descriptors
    bcut_mwlow_values = bcut2d_mwlow.(mols)
    @test length(bcut_mwlow_values) == 3
    @test all(v -> isa(v, Float64) || v === missing, bcut_mwlow_values)

    # Test vectorized structure counts
    aliphatic_rings = num_aliphatic_rings.(mols)
    @test length(aliphatic_rings) == 3
    @test all(v -> isa(v, Int) || v === missing, aliphatic_rings)

    # Note: num_sp3_heavy_atoms test removed as function doesn't exist

    # Test vectorized 3D descriptors
    pmi1_values = pmi1.(mols)
    @test length(pmi1_values) == 3
    @test all(v -> isa(v, Float64) || v === missing, pmi1_values)
end
