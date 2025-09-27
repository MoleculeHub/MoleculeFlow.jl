using Test
using MoleculeFlow

@testset "Molecular Standardization" begin
    @testset "Salt Stripping" begin
        # Test salt stripping with ethanol + chloride
        salt_mol = mol_from_smiles("CCO.Cl")
        @test salt_mol.valid == true

        stripped = strip_salts(salt_mol)
        @test stripped.valid == true

        # Stripped molecule should be just ethanol
        stripped_smiles = mol_to_smiles(stripped)
        @test stripped_smiles == "CCO" || stripped_smiles == "OCC"

        # Test with sodium acetate
        sodium_acetate = mol_from_smiles("CC(=O)[O-].[Na+]")
        if sodium_acetate.valid
            stripped_acetate = strip_salts(sodium_acetate)
            @test stripped_acetate.valid == true
            # Should keep the acetate part
            acetate_smiles = mol_to_smiles(stripped_acetate)
            @test occursin("CC", acetate_smiles) && occursin("O", acetate_smiles)
        end

        # Test with molecule that has no salts
        pure_mol = mol_from_smiles("CCO")
        stripped_pure = strip_salts(pure_mol)
        @test stripped_pure.valid == true
        @test mol_to_smiles(stripped_pure) == mol_to_smiles(pure_mol)
    end

    @testset "Tautomer Enumeration" begin
        # Test with acetylacetone (classic keto-enol tautomerism)
        acetylacetone = mol_from_smiles("CC(=O)CC(=O)C")
        @test acetylacetone.valid == true

        tautomers = enumerate_tautomers(acetylacetone; max_tautomers = 5)
        @test length(tautomers) >= 1
        @test all(t -> t.valid, tautomers)

        # All tautomers should be valid molecules
        for taut in tautomers
            smiles = mol_to_smiles(taut)
            @test !ismissing(smiles)
            @test length(smiles) > 0
        end

        # Test with a molecule that has no reasonable tautomers
        methane = mol_from_smiles("C")
        @test methane.valid
        tautomers_methane = enumerate_tautomers(methane)
        @test length(tautomers_methane) >= 1  # Should at least return the original

        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid")
        tautomers_bad = enumerate_tautomers(bad_mol)
        @test length(tautomers_bad) == 1
        @test !tautomers_bad[1].valid
    end

    @testset "Canonical Tautomer" begin
        # Test enol -> keto conversion
        enol_form = mol_from_smiles("CC(O)=CC(=O)C")
        if enol_form.valid
            canonical = canonical_tautomer(enol_form)
            @test canonical.valid == true

            canonical_smiles = mol_to_smiles(canonical)
            @test !ismissing(canonical_smiles)
        end

        # Test with simple molecule (should return same)
        ethanol = mol_from_smiles("CCO")
        @test ethanol.valid
        canonical_ethanol = canonical_tautomer(ethanol)
        @test canonical_ethanol.valid == true
        @test mol_to_smiles(canonical_ethanol) == mol_to_smiles(ethanol)

        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid")
        canonical_bad = canonical_tautomer(bad_mol)
        @test !canonical_bad.valid
    end

    @testset "Charge Neutralization" begin
        # Test with acetate anion
        acetate = mol_from_smiles("CC(=O)[O-]")
        if acetate.valid
            neutralized = neutralize_charges(acetate)
            @test neutralized.valid == true

            # Should become acetic acid
            neutral_smiles = mol_to_smiles(neutralized)
            @test occursin("CC", neutral_smiles) && occursin("O", neutral_smiles)
        end

        # Test with neutral molecule (should remain unchanged)
        ethanol = mol_from_smiles("CCO")
        neutral_ethanol = neutralize_charges(ethanol)
        @test neutral_ethanol.valid == true
        @test mol_to_smiles(neutral_ethanol) == mol_to_smiles(ethanol)

        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid")
        neutral_bad = neutralize_charges(bad_mol)
        @test !neutral_bad.valid
    end

    @testset "Molecule Normalization" begin
        # Test with quinone structure
        quinone = mol_from_smiles("O=C1C=CC(=O)C=C1")
        if quinone.valid
            normalized = normalize_molecule(quinone)
            @test normalized.valid == true

            normalized_smiles = mol_to_smiles(normalized)
            @test !ismissing(normalized_smiles)
        end

        # Test with normal molecule
        ethanol = mol_from_smiles("CCO")
        normalized_ethanol = normalize_molecule(ethanol)
        @test normalized_ethanol.valid == true

        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid")
        normalized_bad = normalize_molecule(bad_mol)
        @test !normalized_bad.valid
    end

    @testset "Comprehensive Standardization" begin
        # Test with complex molecule: enol form with sodium salt
        complex_mol = mol_from_smiles("CC(O)=CC(=O)C.[Na+]")
        if complex_mol.valid
            standardized = standardize_molecule(complex_mol)
            @test standardized.valid == true

            standardized_smiles = mol_to_smiles(standardized)
            @test !ismissing(standardized_smiles)
            # Should not contain sodium
            @test !occursin("Na", standardized_smiles)
        end

        # Test with various standardization options
        test_mol = mol_from_smiles("CCO.Cl")
        if test_mol.valid
            # With salt stripping
            std_with_salt = standardize_molecule(test_mol; strip_salts_flag = true)
            @test std_with_salt.valid == true

            # Without salt stripping
            std_no_salt = standardize_molecule(test_mol; strip_salts_flag = false)
            @test std_no_salt.valid == true

            # The salt-stripped version should be different
            @test mol_to_smiles(std_with_salt) != mol_to_smiles(std_no_salt)
        end

        # Test with stereochemistry removal
        chiral_mol = mol_from_smiles("C[C@H](O)C")  # Chiral center
        if chiral_mol.valid
            # Keep stereochemistry
            std_keep_stereo = standardize_molecule(
                chiral_mol; remove_stereochemistry = false
            )
            @test std_keep_stereo.valid == true

            # Remove stereochemistry
            std_no_stereo = standardize_molecule(chiral_mol; remove_stereochemistry = true)
            @test std_no_stereo.valid == true

            # Check if stereochemistry markers are removed
            stereo_smiles = mol_to_smiles(std_keep_stereo)
            no_stereo_smiles = mol_to_smiles(std_no_stereo)
            # The no-stereo version should not contain @ symbols
            if no_stereo_smiles !== missing
                @test !occursin("@", no_stereo_smiles)
            end
        end

        # Test with invalid molecule
        bad_mol = mol_from_smiles("invalid")
        standardized_bad = standardize_molecule(bad_mol)
        @test !standardized_bad.valid
    end

    @testset "Standardization Error Handling" begin
        # Test all functions with invalid molecules
        bad_mol = mol_from_smiles("invalid")

        @test !strip_salts(bad_mol).valid
        @test length(enumerate_tautomers(bad_mol)) == 1
        @test !enumerate_tautomers(bad_mol)[1].valid
        @test !canonical_tautomer(bad_mol).valid
        @test !neutralize_charges(bad_mol).valid
        @test !normalize_molecule(bad_mol).valid
        @test !standardize_molecule(bad_mol).valid
    end
end
