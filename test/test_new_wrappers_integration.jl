using Test
using MoleculeFlow

@testset "New Wrappers Integration Tests" begin
    # Test molecules covering different structural features
    test_molecules = [
        ("ethanol", "CCO"),
        ("benzene", "c1ccccc1"),
        ("cyclohexane", "C1CCCCC1"),
        ("pyridine", "c1cccnc1"),
        ("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("chiral_molecule", "C[C@H](O)C"),
        ("tetrahydrofuran", "C1CCOC1"),
        ("acetamide", "CC(=O)N"),
    ]

    @testset "Comprehensive Descriptor Coverage" begin
        for (name, smiles) in test_molecules
            mol = mol_from_smiles(smiles)
            @test mol.valid

            # Test all new VSA descriptors
            vsa_descriptors = [
                slogp_vsa2, slogp_vsa3, slogp_vsa4, slogp_vsa5, slogp_vsa6,
                slogp_vsa7, slogp_vsa8, slogp_vsa9, slogp_vsa10, slogp_vsa11, slogp_vsa12,
                smr_vsa1, smr_vsa2, smr_vsa3, smr_vsa4, smr_vsa5,
                smr_vsa6, smr_vsa7, smr_vsa8, smr_vsa9, smr_vsa10,
                peoe_vsa1, peoe_vsa2, peoe_vsa3, peoe_vsa4, peoe_vsa5, peoe_vsa6, peoe_vsa7,
                peoe_vsa8, peoe_vsa9, peoe_vsa10, peoe_vsa11, peoe_vsa12, peoe_vsa13, peoe_vsa14,
            ]

            for desc_func in vsa_descriptors
                result = desc_func(mol)
                @test isa(result, Float64) || result === missing
                if isa(result, Float64)
                    @test result >= 0.0  # VSA descriptors should be non-negative
                end
            end

            # Test all new BCUT descriptors
            bcut_descriptors = [
                bcut2d_mwlow, bcut2d_mwhi, bcut2d_chglow, bcut2d_chghi,
                bcut2d_logplow, bcut2d_logphi, bcut2d_mrlow, bcut2d_mrhi,
            ]

            for desc_func in bcut_descriptors
                result = desc_func(mol)
                @test isa(result, Float64) || result === missing
            end

            # Test 3D descriptors (may be missing without coordinates)
            threeD_descriptors = [pmi1, pmi2, pmi3]
            for desc_func in threeD_descriptors
                result = desc_func(mol)
                @test isa(result, Float64) || result === missing
            end

            # Test additional structure count descriptors
            structure_descriptors = [
                num_aliphatic_heterocycles, num_saturated_heterocycles, num_saturated_carbocycles,
                num_unspecified_atom_stereo_centers, num_spiro_atoms, num_bridgehead_atoms,
                num_aliphatic_rings, num_heterocycles,
            ]

            for desc_func in structure_descriptors
                result = desc_func(mol)
                @test isa(result, Int) || result === missing
                if isa(result, Int)
                    @test result >= 0  # Counts should be non-negative
                end
            end

            # Test additional descriptors
            @test isa(hall_kier_alpha(mol), Float64) || hall_kier_alpha(mol) === missing
            @test isa(max_absolute_e_state_index(mol), Float64) || max_absolute_e_state_index(mol) === missing
            @test isa(min_absolute_e_state_index(mol), Float64) || min_absolute_e_state_index(mol) === missing
        end
    end

    @testset "Molecular Operations Integration" begin
        mol = mol_from_smiles("c1ccccc1CCO")  # Benzyl alcohol

        # Test ring information
        ring_info = get_ring_info(mol)
        @test ring_info !== missing

        # Test canonical ranking
        ranks = canonical_atom_ranks(mol)
        @test isa(ranks, Vector{Int}) || ranks === missing

        # Test atom environment
        if isa(ranks, Vector{Int}) && length(ranks) > 3
            env = find_atom_environment(mol, 1, 0)
            @test isa(env, Vector{Int}) || env === missing
        end

        # Test fragment operations
        if heavy_atom_count(mol) >= 3
            frag_smarts = mol_fragment_to_smarts(mol, [0, 1, 2])
            @test isa(frag_smarts, String) || frag_smarts === missing

            frag_smiles = mol_fragment_to_smiles(mol, [0, 1])
            @test isa(frag_smiles, String) || frag_smiles === missing
        end

        # Test molecular editing operations
        mol_copy = mol_from_smiles("CCO")

        # Test in-place operations
        sanitize_mol!(mol_copy)
        @test mol_copy.valid

        compute_2d_coords!(mol_copy)
        @test mol_copy.valid

        fast_find_rings!(mol_copy)
        @test mol_copy.valid

        remove_stereochemistry!(mol_copy)
        @test mol_copy.valid

        # Test atom renumbering
        if heavy_atom_count(mol_copy) == 3
            new_mol = renumber_atoms(mol_copy, [2, 1, 0])
            @test new_mol.valid
        end
    end

    @testset "Error Handling Consistency" begin
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test !invalid_mol.valid

        # Test that all new descriptors handle invalid molecules consistently
        @test slogp_vsa2(invalid_mol) === missing
        @test smr_vsa1(invalid_mol) === missing
        @test peoe_vsa1(invalid_mol) === missing
        @test bcut2d_mwlow(invalid_mol) === missing
        @test pmi1(invalid_mol) === missing
        @test num_aliphatic_heterocycles(invalid_mol) === missing
        @test hall_kier_alpha(invalid_mol) === missing

        # Test that all new operations handle invalid molecules consistently
        @test get_ring_info(invalid_mol) === missing
        @test canonical_atom_ranks(invalid_mol) == Int[]
        @test find_atom_environment(invalid_mol, 1, 0) === missing
        @test mol_fragment_to_smarts(invalid_mol, [0, 1]) === missing
        @test mol_fragment_to_smiles(invalid_mol, [0, 1]) === missing

        # In-place operations should return invalid molecule unchanged
        @test !sanitize_mol!(invalid_mol).valid
        @test !compute_2d_coords!(invalid_mol).valid
        @test fast_find_rings!(invalid_mol) == false
        @test !remove_stereochemistry!(invalid_mol).valid
        @test !renumber_atoms(invalid_mol, [0, 1, 2]).valid
    end

    @testset "Vectorized Operations" begin
        mols = [mol_from_smiles(smiles) for (_, smiles) in test_molecules]

        # Test vectorized VSA descriptors
        slogp2_values = slogp_vsa2.(mols)
        @test length(slogp2_values) == length(mols)
        @test all(v -> isa(v, Float64) || v === missing, slogp2_values)

        peoe1_values = peoe_vsa1.(mols)
        @test length(peoe1_values) == length(mols)
        @test all(v -> isa(v, Float64) || v === missing, peoe1_values)

        # Test vectorized BCUT descriptors
        bcut_values = bcut2d_mwlow.(mols)
        @test length(bcut_values) == length(mols)
        @test all(v -> isa(v, Float64) || v === missing, bcut_values)

        # Test vectorized structure counts
        aliphatic_counts = num_aliphatic_rings.(mols)
        @test length(aliphatic_counts) == length(mols)
        @test all(v -> isa(v, Int) || v === missing, aliphatic_counts)

        # Test vectorized 3D descriptors
        pmi1_values = pmi1.(mols)
        @test length(pmi1_values) == length(mols)
        @test all(v -> isa(v, Float64) || v === missing, pmi1_values)
    end

    @testset "Descriptor Value Sanity Checks" begin
        # Test that descriptors give sensible values for known molecules
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")
        cyclohexane = mol_from_smiles("C1CCCCC1")

        # Test ring counting consistency
        @test num_aliphatic_rings(ethanol) == 0
        @test num_aliphatic_rings(cyclohexane) == 1
        @test num_aliphatic_rings(benzene) == 0  # Aromatic, not aliphatic

        @test num_aromatic_rings(benzene) == 1
        @test num_aromatic_rings(cyclohexane) == 0
        @test num_aromatic_rings(ethanol) == 0

        # Test heterocycle counting
        pyridine = mol_from_smiles("c1cccnc1")
        @test num_heterocycles(pyridine) == 1
        @test num_heterocycles(benzene) == 0
        @test num_aromatic_heterocycles(pyridine) == 1

        # Test structure features
        @test num_spiro_atoms(ethanol) == 0
        @test num_bridgehead_atoms(ethanol) == 0

        # Test that BCUT low/high relationships make sense (when both are valid)
        mwlow = bcut2d_mwlow(benzene)
        mwhi = bcut2d_mwhi(benzene)
        if isa(mwlow, Float64) && isa(mwhi, Float64)
            @test mwlow <= mwhi
        end
    end

    @testset "Performance and Memory" begin
        # Test that descriptor calculations don't cause memory issues
        mol = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin

        # Calculate many descriptors in sequence
        descriptors = [
            slogp_vsa2, slogp_vsa5, slogp_vsa10,
            smr_vsa3, smr_vsa7,
            peoe_vsa2, peoe_vsa8, peoe_vsa14,
            bcut2d_mwlow, bcut2d_logphi,
            num_aliphatic_rings, num_heterocycles,
            hall_kier_alpha,
        ]

        for desc in descriptors
            result = desc(mol)
            @test isa(result, Number) || result === missing
        end

        # Test operations don't cause issues
        get_ring_info(mol)
        canonical_atom_ranks(mol)
        sanitize_mol!(mol)
        compute_2d_coords!(mol)

        @test mol.valid  # Should still be valid after all operations
    end
end

@testset "Backward Compatibility" begin
    # Ensure new wrappers don't break existing functionality
    mol = mol_from_smiles("CCO")

    # Test that original descriptors still work
    @test isa(molecular_weight(mol), Float64)
    @test isa(heavy_atom_count(mol), Int)
    @test isa(num_hbd(mol), Int)
    @test isa(num_hba(mol), Int)
    @test isa(logp(mol), Float64)
    @test isa(tpsa(mol), Float64)

    # Test that original operations still work
    @test isa(add_hs(mol), Molecule)
    @test isa(remove_hs(mol), Molecule)
    @test isa(mol_to_smiles(mol), String)
    @test isa(mol_to_inchi_key(mol), String)

    # Test fingerprints still work
    @test isa(morgan_fingerprint(mol), Vector{Bool})
    @test isa(rdk_fingerprint(mol), Vector{Bool})
    @test isa(maccs_fingerprint(mol), Vector{Bool})
end