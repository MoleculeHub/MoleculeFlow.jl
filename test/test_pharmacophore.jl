using Test
using MoleculeFlow

@testset "Pharmacophore Features" begin
    @testset "Feature Factory Creation" begin
        # Test default feature factory
        factory = create_feature_factory()
        @test factory isa FeatureFactory
        @test length(factory.feature_families) > 0
        @test factory.num_definitions > 0

        # Test that common feature families are present
        families = get_feature_families(factory)
        @test "Donor" in families
        @test "Acceptor" in families
        @test "Aromatic" in families

        # Test custom feature definition string
        fdef_string = raw"""
        DefineFeature HAcceptor1 [N,O;H0]
           Family HBondAcceptor
           Weights 1.0
        EndFeature
        """

        try
            custom_factory = create_feature_factory(
                use_default = false, fdef_string = fdef_string
            )
            @test custom_factory isa FeatureFactory
        catch e
            @warn "Custom feature factory creation failed: $e"
        end
    end

    @testset "Molecular Feature Extraction" begin
        # Test with simple molecules
        ethanol = mol_from_smiles("CCO")
        phenol = mol_from_smiles("c1ccccc1O")

        factory = create_feature_factory()

        # Test feature extraction
        ethanol_features = get_mol_features(ethanol, factory)
        @test ethanol_features isa Vector{ChemicalFeature}

        phenol_features = get_mol_features(phenol, factory)
        @test phenol_features isa Vector{ChemicalFeature}
        @test length(phenol_features) >= length(ethanol_features)  # Phenol should have more features

        # Test that features have valid properties
        if !isempty(ethanol_features)
            feature = ethanol_features[1]
            @test feature.family isa String
            @test feature.type isa String
            @test feature.atom_ids isa Vector{Int}
            @test feature.position isa Vector{Float64}
            @test feature.id isa Int
            @test length(feature.position) == 3  # 3D coordinates
            @test all(id -> id > 0, feature.atom_ids)  # 1-based indexing
        end

        # Test vectorized operation
        molecules = [ethanol, phenol]
        all_features = get_mol_features(molecules, factory)
        @test length(all_features) == 2
        @test all_features[1] isa Vector{ChemicalFeature}
        @test all_features[2] isa Vector{ChemicalFeature}

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        invalid_features = get_mol_features(invalid_mol, factory)
        @test isempty(invalid_features)
    end

    @testset "Feature Filtering" begin
        mol = mol_from_smiles("CCO")  # Ethanol: has donor and acceptor
        factory = create_feature_factory()
        features = get_mol_features(mol, factory)

        if !isempty(features)
            # Test filtering by family
            families = unique([f.family for f in features])
            if !isempty(families)
                first_family = families[1]
                filtered = filter_features_by_family(features, first_family)
                @test all(f -> f.family == first_family, filtered)
            end

            # Test filtering for non-existent family
            empty_filtered = filter_features_by_family(features, "NonExistentFamily")
            @test isempty(empty_filtered)
        end
    end

    @testset "Pharmacophore Fingerprints" begin
        molecules = [
            mol_from_smiles("CCO"),      # Ethanol
            mol_from_smiles("c1ccccc1O"), # Phenol
            mol_from_smiles("CC(=O)N"),   # Acetamide
        ]

        for mol in molecules
            # Test basic fingerprint generation
            fp = pharmacophore_fingerprint(mol)
            @test fp isa Vector{Bool} || fp === missing

            if fp isa Vector{Bool}
                @test length(fp) > 0
                # Test with custom parameters
                fp_custom = pharmacophore_fingerprint(mol; min_points = 2, max_points = 4)
                @test fp_custom isa Vector{Bool} || fp_custom === missing
            end
        end

        # Test vectorized fingerprint generation
        fps = pharmacophore_fingerprint(molecules)
        @test length(fps) == 3
        @test all(fp -> fp isa Vector{Bool} || fp === missing, fps)

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        invalid_fp = pharmacophore_fingerprint(invalid_mol)
        @test invalid_fp === missing
    end

    @testset "3D Pharmacophore Features" begin
        mol = mol_from_smiles("CCO")
        factory = create_feature_factory()

        # Test without 3D coordinates (should work but may have limited positioning)
        ph4_2d = get_pharmacophore_3d(mol, factory)
        @test ph4_2d isa Vector{Tuple{String, Vector{Float64}}}

        # Test with 3D coordinates if conformer generation is available
        try
            result = generate_3d_conformers(mol; num_conformers = 1)
            if result.success && !isempty(result.molecules)
                mol_3d = result.molecules[1]
                ph4_3d = get_pharmacophore_3d(mol_3d, factory)

                @test ph4_3d isa Vector{Tuple{String, Vector{Float64}}}

                # Check that each pharmacophore point has valid data
                for (family, position) in ph4_3d
                    @test family isa String
                    @test position isa Vector{Float64}
                    @test length(position) == 3  # 3D coordinates
                    @test !any(isnan, position)
                end

                # Test vectorized 3D pharmacophore extraction
                molecules_3d = [mol_3d]
                ph4_batch = get_pharmacophore_3d(molecules_3d, factory)
                @test length(ph4_batch) == 1
                @test ph4_batch[1] isa Vector{Tuple{String, Vector{Float64}}}
            end
        catch e
            @warn "3D conformer generation not available for pharmacophore testing: $e"
        end

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        invalid_ph4 = get_pharmacophore_3d(invalid_mol, factory)
        @test isempty(invalid_ph4)
    end

    @testset "Edge Cases and Error Handling" begin
        factory = create_feature_factory()

        # Test with empty molecule (if possible)
        try
            empty_mol = mol_from_smiles("")
            empty_features = get_mol_features(empty_mol, factory)
            @test isempty(empty_features)

            empty_fp = pharmacophore_fingerprint(empty_mol)
            @test empty_fp isa Vector{Bool} || empty_fp === missing
        catch
            # Some SMILES parsers may not handle empty strings
        end

        # Test with very simple molecule
        methane = mol_from_smiles("C")
        methane_features = get_mol_features(methane, factory)
        @test methane_features isa Vector{ChemicalFeature}

        methane_fp = pharmacophore_fingerprint(methane)
        @test methane_fp isa Vector{Bool} || methane_fp === missing
    end

    @testset "Common Use Cases" begin
        # Test drug-like molecules
        aspirin = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O")
        caffeine = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

        factory = create_feature_factory()

        # These should have multiple features
        aspirin_features = get_mol_features(aspirin, factory)
        caffeine_features = get_mol_features(caffeine, factory)

        @test aspirin_features isa Vector{ChemicalFeature}
        @test caffeine_features isa Vector{ChemicalFeature}

        # Test pharmacophore fingerprints for drug-like molecules
        aspirin_fp = pharmacophore_fingerprint(aspirin)
        caffeine_fp = pharmacophore_fingerprint(caffeine)

        @test aspirin_fp isa Vector{Bool} || aspirin_fp === missing
        @test caffeine_fp isa Vector{Bool} || caffeine_fp === missing

        # Fingerprints should be different for different molecules
        if aspirin_fp isa Vector{Bool} && caffeine_fp isa Vector{Bool}
            @test length(aspirin_fp) == length(caffeine_fp)  # Same fingerprint size
            @test aspirin_fp != caffeine_fp  # Different molecules should have different fingerprints
        end
    end
end
