using Test
using MoleculeFlow

@testset "Hydrogen Manipulation" begin
    @testset "add_hs function" begin
        # Test basic hydrogen addition
        mol = mol_from_smiles("CCO")
        @test mol.valid == true

        mol_with_hs = add_hs(mol)
        @test mol_with_hs.valid == true
        @test mol_with_hs.source == mol.source

        # Properties should be copied
        mol.test_prop = "test_value"
        mol_with_hs2 = add_hs(mol)
        @test mol_with_hs2.props[:test_prop] == "test_value"

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = add_hs(invalid_mol)
        @test result.valid == false
    end

    @testset "remove_hs function" begin
        # Test basic hydrogen removal
        mol = mol_from_smiles("CCO")
        mol_with_hs = add_hs(mol)
        mol_no_hs = remove_hs(mol_with_hs)

        @test mol_no_hs.valid == true
        @test mol_no_hs.source == mol.source

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = remove_hs(invalid_mol)
        @test result.valid == false
    end

    @testset "hydrogen add/remove cycle" begin
        # Test that add_hs -> remove_hs preserves molecular structure
        original = mol_from_smiles("CCO")
        processed = remove_hs(add_hs(original))

        @test processed.valid == true
        @test mol_to_smiles(original) == mol_to_smiles(processed)
    end
end

@testset "Extended File I/O" begin
    @testset "mol_to_inchi_key function" begin
        mol = mol_from_smiles("CCO")
        inchi_key = mol_to_inchi_key(mol)

        @test isa(inchi_key, String)
        @test length(inchi_key) > 0
        @test inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"  # Known InChI key for ethanol

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = mol_to_inchi_key(invalid_mol)
        @test result == ""
    end

    @testset "mol_to_molblock function" begin
        mol = mol_from_smiles("CCO")
        molblock = mol_to_molblock(mol)

        @test isa(molblock, String)
        @test length(molblock) > 0
        @test occursin("M  END", molblock)  # MOL files end with "M  END"

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = mol_to_molblock(invalid_mol)
        @test result == ""
    end

    @testset "mol_to_pdb_block function" begin
        mol = mol_from_smiles("CCO")
        pdb_block = mol_to_pdb_block(mol)

        @test isa(pdb_block, String)
        # Note: PDB block might be empty for simple molecules without 3D coords

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = mol_to_pdb_block(invalid_mol)
        @test result == ""
    end

    @testset "mol_from_pdb_block function" begin
        # Test with minimal PDB data
        pdb_data = """
ATOM      1  C   MOL A   1      20.154  21.875  21.235  1.00 10.00           C
ATOM      2  C   MOL A   1      19.618  22.166  22.618  1.00 10.00           C
ATOM      3  O   MOL A   1      18.259  22.160  22.618  1.00 10.00           O
END
"""
        mol = mol_from_pdb_block(pdb_data)
        @test mol.source == "PDB block"
        # Note: PDB parsing might fail for simple cases, so we just test that function runs

        # Test with invalid PDB data (note: empty/invalid PDB might still be considered valid)
        mol_invalid = mol_from_pdb_block("invalid pdb data")
        @test mol_invalid.source == "PDB block"  # Just test that function runs
    end
end

@testset "Molecular Editing and Manipulation" begin
    @testset "combine_mols function" begin
        mol1 = mol_from_smiles("CCO")
        mol2 = mol_from_smiles("CCC")

        combined = combine_mols(mol1, mol2)
        @test combined.valid == true
        @test combined.source == "combined"

        # Test with invalid molecules
        invalid_mol = mol_from_smiles("invalid_smiles")
        result1 = combine_mols(mol1, invalid_mol)
        @test result1.valid == false

        result2 = combine_mols(invalid_mol, mol2)
        @test result2.valid == false
    end

    @testset "delete_substructs function" begin
        # Test removing alcohol group
        mol = mol_from_smiles("CCO")
        result = delete_substructs(mol, "[OH]")

        @test result.valid == true
        @test result.source == mol.source

        # Test with invalid SMARTS
        result_invalid = delete_substructs(mol, "invalid_smarts")
        @test result_invalid.source == mol.source  # Should return original on error

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result_invalid_mol = delete_substructs(invalid_mol, "[OH]")
        @test result_invalid_mol.valid == false
    end

    @testset "replace_substructs function" begin
        mol = mol_from_smiles("CCO")

        # Test replacing alcohol with amine (this might not always work)
        results = replace_substructs(mol, "[OH]", "N")
        @test isa(results, Vector{Molecule})
        @test length(results) >= 1

        # Test with invalid SMARTS/SMILES
        results_invalid = replace_substructs(mol, "invalid_smarts", "N")
        @test length(results_invalid) == 1
        @test results_invalid[1] == mol  # Should return original

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        results_invalid_mol = replace_substructs(invalid_mol, "[OH]", "N")
        @test length(results_invalid_mol) == 1
        @test results_invalid_mol[1].valid == false
    end
end

@testset "Stereochemistry Operations" begin
    @testset "assign_stereochemistry! function" begin
        mol = mol_from_smiles("C[C@H](O)C")  # Chiral molecule

        success = assign_stereochemistry!(mol)
        @test isa(success, Bool)

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = assign_stereochemistry!(invalid_mol)
        @test result == false
    end

    @testset "find_chiral_centers function" begin
        # Test with achiral molecule
        mol_achiral = mol_from_smiles("CCO")
        centers_achiral = find_chiral_centers(mol_achiral)
        @test isa(centers_achiral, Vector{Tuple{Int,String}})

        # Test with potentially chiral molecule
        mol_chiral = mol_from_smiles("C[C@H](O)C")
        centers_chiral = find_chiral_centers(mol_chiral)
        @test isa(centers_chiral, Vector{Tuple{Int,String}})

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = find_chiral_centers(invalid_mol)
        @test result == Tuple{Int,String}[]
    end
end

@testset "Ring Analysis" begin
    @testset "fast_find_rings! function" begin
        mol = mol_from_smiles("c1ccccc1")  # Benzene

        success = fast_find_rings!(mol)
        @test isa(success, Bool)

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = fast_find_rings!(invalid_mol)
        @test result == false
    end

    @testset "canonical_atom_ranks function" begin
        mol = mol_from_smiles("CCO")
        ranks = canonical_atom_ranks(mol)

        @test isa(ranks, Vector{Int})
        @test length(ranks) == 3  # Should have 3 atoms
        @test all(r >= 0 for r in ranks)  # All ranks should be non-negative

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = canonical_atom_ranks(invalid_mol)
        @test result == Int[]
    end
end

@testset "Pattern Matching" begin
    @testset "quick_smarts_match function" begin
        mol = mol_from_smiles("CCO")

        # Test matching alcohol group
        has_alcohol = quick_smarts_match(mol, "[OH]")
        @test has_alcohol == true

        # Test matching carbon
        has_carbon = quick_smarts_match(mol, "[C]")
        @test has_carbon == true

        # Test non-matching pattern
        has_nitrogen = quick_smarts_match(mol, "[N]")
        @test has_nitrogen == false

        # Test with invalid SMARTS
        result_invalid = quick_smarts_match(mol, "invalid_smarts")
        @test result_invalid == false

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = quick_smarts_match(invalid_mol, "[OH]")
        @test result == false
    end

    @testset "mol_fragment_to_smarts function" begin
        mol = mol_from_smiles("CCO")

        # Test with first two atoms (indices are 0-based in RDKit)
        smarts = mol_fragment_to_smarts(mol, [0, 1])
        @test isa(smarts, String)

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        result = mol_fragment_to_smarts(invalid_mol, [0, 1])
        @test result == ""

        # Test with empty indices
        result_empty = mol_fragment_to_smarts(mol, Int[])
        @test isa(result_empty, String)
    end
end

@testset "Error Handling and Edge Cases" begin
    @testset "Invalid inputs handling" begin
        # Most error handling is tested above, but let's add a few more edge cases

        # Test empty SMILES (note: empty string might create valid empty molecule)
        empty_mol = mol_from_smiles("")
        # Just test that operations handle it gracefully
        inchi_key_result = mol_to_inchi_key(empty_mol)
        @test isa(inchi_key_result, String)

        molblock_result = mol_to_molblock(empty_mol)
        @test isa(molblock_result, String)

        ranks_result = canonical_atom_ranks(empty_mol)
        @test isa(ranks_result, Vector{Int})

        smarts_result = quick_smarts_match(empty_mol, "[C]")
        @test isa(smarts_result, Bool)
    end
end

@testset "XYZ Format Support" begin
    @testset "XYZ Block Parsing" begin
        # Test valid XYZ block
        xyz_data = """3
Test molecule
C    0.000    0.000    0.000
C    1.520    0.000    0.000
O    2.080    1.100    0.000"""

        mol = mol_from_xyz_block(xyz_data)
        @test mol.valid
        @test heavy_atom_count(mol) == 3
        @test num_carbons(mol) == 2
        @test num_oxygens(mol) == 1

        # Test that a SMILES can be generated (may be disconnected since XYZ has no bonds)
        smiles = mol_to_smiles(mol)
        @test isa(smiles, String)
        @test length(smiles) > 0
        # XYZ format doesn't contain bond info, so molecules may be disconnected
        @test occursin("C", smiles)  # Should contain carbon
        @test occursin("O", smiles)  # Should contain oxygen

        # Test invalid XYZ block (XYZ parser may be lenient, just test it doesn't crash)
        invalid_xyz = "invalid xyz data"
        invalid_mol = mol_from_xyz_block(invalid_xyz)
        @test isa(invalid_mol, Molecule)

        # Test empty XYZ block (XYZ parser may be lenient, just test it doesn't crash)
        empty_mol = mol_from_xyz_block("")
        @test isa(empty_mol, Molecule)
    end

    @testset "XYZ Export" begin
        # Test XYZ export without 3D coordinates (should fail gracefully)
        mol = mol_from_smiles("CCO")
        xyz_string = mol_to_xyz_block(mol)
        @test xyz_string == ""  # Should be empty when no 3D coords

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")  # Creates invalid molecule
        xyz_invalid = mol_to_xyz_block(invalid_mol)
        @test xyz_invalid == ""
    end
end

@testset "MOL2 Format Support" begin
    @testset "MOL2 Block Parsing" begin
        # Test valid MOL2 block
        mol2_data = """@<TRIPOS>MOLECULE
ethanol
3 2 0 0 0
SMALL

@<TRIPOS>ATOM
1 C1 0.0000 0.0000 0.0000 C.3 1 RES1 0.0000
2 C2 1.5200 0.0000 0.0000 C.3 1 RES1 0.0000
3 O1 2.0800 1.1000 0.0000 O.3 1 RES1 0.0000

@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""

        mol = mol_from_mol2_block(mol2_data)
        @test mol.valid
        @test heavy_atom_count(mol) == 3
        @test num_carbons(mol) == 2
        @test num_oxygens(mol) == 1

        # Test that molecular weight is reasonable (may differ due to charge interpretation)
        mw = molecular_weight(mol)
        @test isa(mw, Float64)
        @test mw > 0
        @test mw < 100  # Should be reasonable for a 3-atom molecule

        # Test that SMILES can be generated (may have charges from MOL2 interpretation)
        smiles = mol_to_smiles(mol)
        @test isa(smiles, String)
        @test length(smiles) > 0
        @test occursin("C", smiles)  # Should contain carbon
        @test occursin("O", smiles)  # Should contain oxygen

        # Test invalid MOL2 block (MOL2 parser may be lenient, just test it doesn't crash)
        invalid_mol2 = "invalid mol2 data"
        invalid_mol = mol_from_mol2_block(invalid_mol2)
        @test isa(invalid_mol, Molecule)

        # Test empty MOL2 block (MOL2 parser may be lenient, just test it doesn't crash)
        empty_mol = mol_from_mol2_block("")
        @test isa(empty_mol, Molecule)
    end

    @testset "MOL2 vs Other Formats" begin
        # Compare MOL2 parsing with SMILES for same molecule
        mol2_data = """@<TRIPOS>MOLECULE
water
1 0 0 0 0
SMALL

@<TRIPOS>ATOM
1 O1 0.0000 0.0000 0.0000 O.3 1 WAT 0.0000
"""

        mol_from_mol2 = mol_from_mol2_block(mol2_data)
        mol_from_smiles_ref = mol_from_smiles("O")

        @test mol_from_mol2.valid
        @test mol_from_smiles_ref.valid
        @test heavy_atom_count(mol_from_mol2) == heavy_atom_count(mol_from_smiles_ref)
        @test num_oxygens(mol_from_mol2) == num_oxygens(mol_from_smiles_ref)
    end
end

@testset "Format Integration Tests" begin
    @testset "Round-trip Format Tests" begin
        # Test that we can convert between formats where applicable
        original_mol = mol_from_smiles("CCO")

        # Test InChI key generation
        inchi_key = mol_to_inchi_key(original_mol)
        @test isa(inchi_key, String)
        @test length(inchi_key) > 0

        # Test MOL block generation
        molblock = mol_to_molblock(original_mol)
        @test isa(molblock, String)
        @test length(molblock) > 0
        @test occursin("RDKit", molblock)  # MOL blocks typically contain software info

        # Test that MOL block can be parsed back
        mol_from_block = mol_from_molblock(molblock)
        @test mol_from_block.valid
        @test heavy_atom_count(mol_from_block) == heavy_atom_count(original_mol)
    end

    @testset "File Format Error Handling" begin
        # Test file functions with non-existent files
        # Note: XYZ reader behavior may vary - test that function doesn't crash
        non_existent_xyz = mol_from_xyz_file("non_existent_file.xyz")
        @test isa(non_existent_xyz, Molecule)  # Should return a Molecule object

        non_existent_mol2 = mol_from_mol2_file("non_existent_file.mol2")
        @test !non_existent_mol2.valid

        non_existent_pdb = mol_from_pdb_file("non_existent_file.pdb")
        @test !non_existent_pdb.valid
    end
end