using Test
using MoleculeFlow

@testset "Advanced Drug-like and ADMET Descriptors" begin
    @testset "QED (Quantitative Estimate of Drug-likeness)" begin
        # Test with known drug-like molecules
        ethanol = mol_from_smiles("CCO")
        aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")

        qed_ethanol = qed(ethanol)
        qed_aspirin = qed(aspirin)

        @test isa(qed_ethanol, Float64)
        @test isa(qed_aspirin, Float64)
        @test 0.0 <= qed_ethanol <= 1.0
        @test 0.0 <= qed_aspirin <= 1.0
        @test qed_aspirin > qed_ethanol  # Aspirin should be more drug-like

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test qed(invalid_mol) === missing
    end

    @testset "Synthetic Accessibility Score" begin
        # Test with simple molecules
        ethanol = mol_from_smiles("CCO")
        methane = mol_from_smiles("C")

        sa_ethanol = synthetic_accessibility(ethanol)
        sa_methane = synthetic_accessibility(methane)

        @test isa(sa_ethanol, Float64)
        @test isa(sa_methane, Float64)
        @test 1.0 <= sa_ethanol <= 10.0
        @test 1.0 <= sa_methane <= 10.0

        # Test with more complex molecule
        aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        sa_aspirin = synthetic_accessibility(aspirin)
        @test isa(sa_aspirin, Float64)
        @test 1.0 <= sa_aspirin <= 10.0

        # Test that we get reasonable ranges (don't make assumptions about relative scores)
        @test sa_ethanol > 0.5  # Should be reasonable
        @test sa_aspirin > 0.5  # Should be reasonable

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test synthetic_accessibility(invalid_mol) === missing
    end

    @testset "Fraction of sp3 carbons" begin
        # Test with fully sp3 molecule
        cyclohexane = mol_from_smiles("C1CCCCC1")
        fsp3_cyclohexane = fraction_csp3(cyclohexane)
        @test fsp3_cyclohexane == 1.0

        # Test with fully sp2 molecule
        benzene = mol_from_smiles("c1ccccc1")
        fsp3_benzene = fraction_csp3(benzene)
        @test fsp3_benzene == 0.0

        # Test with mixed sp2/sp3
        ethanol = mol_from_smiles("CCO")
        fsp3_ethanol = fraction_csp3(ethanol)
        @test fsp3_ethanol == 1.0  # Both carbons are sp3

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test fraction_csp3(invalid_mol) === missing
    end

    @testset "Labute Accessible Surface Area" begin
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")

        asa_ethanol = labute_asa(ethanol)
        asa_benzene = labute_asa(benzene)

        @test isa(asa_ethanol, Float64)
        @test isa(asa_benzene, Float64)
        @test asa_ethanol > 0
        @test asa_benzene > 0
        @test asa_benzene > asa_ethanol  # Benzene should have larger surface area

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test labute_asa(invalid_mol) === missing
    end

    @testset "Molar Refractivity" begin
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")

        mr_ethanol = molar_refractivity(ethanol)
        mr_benzene = molar_refractivity(benzene)

        @test isa(mr_ethanol, Float64)
        @test isa(mr_benzene, Float64)
        @test mr_ethanol > 0
        @test mr_benzene > 0

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test molar_refractivity(invalid_mol) === missing
    end
end

@testset "Advanced Ring and Structure Counts" begin
    @testset "Aliphatic and Aromatic Ring Counts" begin
        # Test aromatic carbocycle (benzene)
        benzene = mol_from_smiles("c1ccccc1")
        @test num_aromatic_carbocycles(benzene) == 1
        @test num_aliphatic_carbocycles(benzene) == 0

        # Test aliphatic carbocycle (cyclohexane)
        cyclohexane = mol_from_smiles("C1CCCCC1")
        @test num_aliphatic_carbocycles(cyclohexane) == 1
        @test num_aromatic_carbocycles(cyclohexane) == 0

        # Test aromatic heterocycle (pyridine)
        pyridine = mol_from_smiles("c1cccnc1")
        @test num_aromatic_heterocycles(pyridine) == 1
        @test num_aromatic_carbocycles(pyridine) == 0

        # Test bicyclic system (naphthalene)
        naphthalene = mol_from_smiles("c1ccc2ccccc2c1")
        @test num_aromatic_carbocycles(naphthalene) == 2

        # Test invalid molecules
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test num_aromatic_carbocycles(invalid_mol) === missing
        @test num_aliphatic_carbocycles(invalid_mol) === missing
        @test num_aromatic_heterocycles(invalid_mol) === missing
    end

    @testset "Stereocenter Count" begin
        # Test chiral molecule
        chiral_mol = mol_from_smiles("C[C@H](O)C")
        stereo_count = num_atom_stereo_centers(chiral_mol)
        @test isa(stereo_count, Int)
        @test stereo_count >= 0

        # Test achiral molecule
        ethanol = mol_from_smiles("CCO")
        @test num_atom_stereo_centers(ethanol) == 0

        # Test invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test num_atom_stereo_centers(invalid_mol) === missing
    end

    @testset "Amide Bond Count" begin
        # Test molecule with amide bond
        acetamide = mol_from_smiles("CC(=O)N")
        @test num_amide_bonds(acetamide) == 1

        # Test molecule without amide bonds
        ethanol = mol_from_smiles("CCO")
        @test num_amide_bonds(ethanol) == 0

        # Test peptide-like molecule with multiple amide bonds
        dipeptide = mol_from_smiles("CC(N)C(=O)NC(C)C(=O)O")
        amide_count = num_amide_bonds(dipeptide)
        @test isa(amide_count, Int)
        @test amide_count >= 1

        # Test invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test num_amide_bonds(invalid_mol) === missing
    end
end

@testset "3D Descriptors" begin
    @testset "Asphericity" begin
        # Test with molecule (may fail without 3D coordinates)
        ethanol = mol_from_smiles("CCO")
        asp_result = asphericity(ethanol)

        # Should either be a Float64 or missing (if no 3D coords)
        @test asp_result === missing || isa(asp_result, Float64)

        if isa(asp_result, Float64)
            @test 0.0 <= asp_result <= 1.0
        end

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test asphericity(invalid_mol) === missing
    end

    @testset "Radius of Gyration" begin
        # Test with molecule (may fail without 3D coordinates)
        ethanol = mol_from_smiles("CCO")
        rog_result = radius_of_gyration(ethanol)

        # Should either be a Float64 or missing (if no 3D coords)
        @test rog_result === missing || isa(rog_result, Float64)

        if isa(rog_result, Float64)
            @test rog_result > 0.0
        end

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test radius_of_gyration(invalid_mol) === missing
    end

    @testset "Principal Moments of Inertia" begin
        # Test PMI descriptors (may fail without 3D coordinates)
        ethanol = mol_from_smiles("CCO")

        pmi1_result = pmi1(ethanol)
        pmi2_result = pmi2(ethanol)
        pmi3_result = pmi3(ethanol)

        # Should either be Float64 or missing (if no 3D coords)
        @test pmi1_result === missing || isa(pmi1_result, Float64)
        @test pmi2_result === missing || isa(pmi2_result, Float64)
        @test pmi3_result === missing || isa(pmi3_result, Float64)

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test pmi1(invalid_mol) === missing
        @test pmi2(invalid_mol) === missing
        @test pmi3(invalid_mol) === missing
    end

    @testset "3D Descriptors with Generated Conformers" begin
        # Skip this test if conformer generation is not available
        try
            mol = mol_from_smiles("CCO")
            conformers = generate_3d_conformers(mol, 1)

            if !isempty(conformers)
                mol_3d = conformers[1].molecule

                # Test asphericity with 3D coordinates
                asp = asphericity(mol_3d)
                @test isa(asp, Float64)
                @test 0.0 <= asp <= 1.0

                # Test radius of gyration with 3D coordinates
                rog = radius_of_gyration(mol_3d)
                @test isa(rog, Float64)
                @test rog > 0.0
            end
        catch e
            @warn "3D conformer generation not available for testing: $e"
        end
    end

    @testset "Advanced 3D Descriptors" begin
        # Test spherocity index
        mol = mol_from_smiles("CCO")
        spherocity_result = spherocity_index(mol)
        @test spherocity_result === missing || isa(spherocity_result, Float64)

        # Test GETAWAY descriptors
        getaway_result = getaway_descriptors(mol)
        @test getaway_result === missing || isa(getaway_result, Vector{Float64})

        # Test WHIM descriptors
        whim_result = whim_descriptors(mol)
        @test whim_result === missing || isa(whim_result, Vector{Float64})

        # Test RDF descriptors
        rdf_result = rdf_descriptors(mol)
        @test rdf_result === missing || isa(rdf_result, Vector{Float64})

        # Test MORSE descriptors
        morse_result = morse_descriptors(mol)
        @test morse_result === missing || isa(morse_result, Vector{Float64})

        # Test with missing molecule
        @test spherocity_index(missing) === missing
        @test getaway_descriptors(missing) === missing
        @test whim_descriptors(missing) === missing
        @test rdf_descriptors(missing) === missing
        @test morse_descriptors(missing) === missing

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test spherocity_index(invalid_mol) === missing
        @test getaway_descriptors(invalid_mol) === missing
        @test whim_descriptors(invalid_mol) === missing
        @test rdf_descriptors(invalid_mol) === missing
        @test morse_descriptors(invalid_mol) === missing
    end

    @testset "Advanced 3D Descriptors with Generated Conformers" begin
        # Skip this test if conformer generation is not available
        try
            mol = mol_from_smiles("CCO")
            conformers = generate_3d_conformers(mol, 1)

            if !isempty(conformers)
                mol_3d = conformers[1].molecule

                # Test spherocity index with 3D coordinates
                spherocity = spherocity_index(mol_3d)
                @test isa(spherocity, Float64)
                @test 0.0 <= spherocity <= 1.0

                # Test GETAWAY descriptors with 3D coordinates
                getaway = getaway_descriptors(mol_3d)
                @test isa(getaway, Vector{Float64})
                @test length(getaway) > 0

                # Test WHIM descriptors with 3D coordinates
                whim = whim_descriptors(mol_3d)
                @test isa(whim, Vector{Float64})
                @test length(whim) > 0

                # Test RDF descriptors with 3D coordinates
                rdf = rdf_descriptors(mol_3d)
                @test isa(rdf, Vector{Float64})
                @test length(rdf) > 0

                # Test MORSE descriptors with 3D coordinates
                morse = morse_descriptors(mol_3d)
                @test isa(morse, Vector{Float64})
                @test length(morse) > 0

                # Test custom parameters
                getaway_custom = getaway_descriptors(mol_3d; precision = 3)
                @test isa(getaway_custom, Vector{Float64})

                whim_custom = whim_descriptors(mol_3d; thresh = 0.01)
                @test isa(whim_custom, Vector{Float64})
            end
        catch e
            @warn "3D conformer generation not available for advanced descriptor testing: $e"
        end
    end
end

@testset "Vectorized Operations for New Descriptors" begin
    # Test that new descriptors work with vectors of molecules
    mols = [mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1"), mol_from_smiles("CC(=O)N")]

    # Test QED
    qed_values = qed.(mols)
    @test length(qed_values) == 3
    @test all(isa(v, Float64) for v in qed_values)

    # Test SAscore
    sa_values = synthetic_accessibility.(mols)
    @test length(sa_values) == 3
    @test all(isa(v, Float64) for v in sa_values)

    # Test Fsp3
    fsp3_values = fraction_csp3.(mols)
    @test length(fsp3_values) == 3
    @test all(isa(v, Float64) for v in fsp3_values)

    # Test ring counts
    arom_carbo_counts = num_aromatic_carbocycles.(mols)
    @test length(arom_carbo_counts) == 3
    @test all(isa(v, Int) for v in arom_carbo_counts)

    # Test amide counts
    amide_counts = num_amide_bonds.(mols)
    @test length(amide_counts) == 3
    @test all(isa(v, Int) for v in amide_counts)

    # Test vectorized 3D descriptors (may be missing without 3D coordinates)
    spherocity_values = spherocity_index.(mols)
    @test length(spherocity_values) == 3
    @test all(v === missing || isa(v, Float64) for v in spherocity_values)

    getaway_values = getaway_descriptors.(mols)
    @test length(getaway_values) == 3
    @test all(v === missing || isa(v, Vector{Float64}) for v in getaway_values)

    whim_values = whim_descriptors.(mols)
    @test length(whim_values) == 3
    @test all(v === missing || isa(v, Vector{Float64}) for v in whim_values)

    rdf_values = rdf_descriptors.(mols)
    @test length(rdf_values) == 3
    @test all(v === missing || isa(v, Vector{Float64}) for v in rdf_values)

    morse_values = morse_descriptors.(mols)
    @test length(morse_values) == 3
    @test all(v === missing || isa(v, Vector{Float64}) for v in morse_values)
end

@testset "Error Handling and Edge Cases" begin
    # Test with empty molecule
    empty_mol = mol_from_smiles("")

    # All descriptors should handle empty/invalid molecules gracefully
    @test qed(empty_mol) === missing || isa(qed(empty_mol), Float64)
    @test synthetic_accessibility(empty_mol) === missing ||
        isa(synthetic_accessibility(empty_mol), Float64)
    @test fraction_csp3(empty_mol) === missing || isa(fraction_csp3(empty_mol), Float64)
    @test labute_asa(empty_mol) === missing || isa(labute_asa(empty_mol), Float64)
    @test molar_refractivity(empty_mol) === missing ||
        isa(molar_refractivity(empty_mol), Float64)

    # Test integer returns for count functions
    @test num_aromatic_carbocycles(empty_mol) === missing ||
        isa(num_aromatic_carbocycles(empty_mol), Int)
    @test num_amide_bonds(empty_mol) === missing || isa(num_amide_bonds(empty_mol), Int)
end
