using Test
using MoleculeFlow

@testset "Molecular Fragmentation" begin

    @testset "BRICS Decomposition" begin
        mol = mol_from_smiles("CCCOCCc1ccccc1")
        fragments = brics_decompose(mol)

        @test isa(fragments, Vector{String})
        @test length(fragments) > 0

        fragments_limited = brics_decompose(mol, min_fragment_size=2, max_fragment_size=6)
        @test isa(fragments_limited, Vector{String})

        bad_mol = mol_from_smiles("invalid_smiles")
        @test brics_decompose(bad_mol) === missing
    end

    @testset "RECAP Decomposition" begin
        mol = mol_from_smiles("CC(=O)NCCc1ccccc1")
        fragments = recap_decompose(mol)

        @test isa(fragments, Vector{String})

        fragments_limited = recap_decompose(mol, min_fragment_size=3)
        @test isa(fragments_limited, Vector{String})

        bad_mol = mol_from_smiles("invalid_smiles")
        @test recap_decompose(bad_mol) === missing
    end

    @testset "Murcko Scaffolds" begin
        mol = mol_from_smiles("CCCOCCc1ccccc1")
        scaffold = get_murcko_scaffold(mol)

        @test isa(scaffold, String)
        @test scaffold == "c1ccccc1"  # Should extract the benzene ring

        generic_scaffold = get_generic_scaffold(mol)
        @test isa(generic_scaffold, String)
        @test generic_scaffold == "C1CCCCC1" 

        acyclic_mol = mol_from_smiles("CCCCCC")
        acyclic_scaffold = get_murcko_scaffold(acyclic_mol)
        @test acyclic_scaffold === missing || acyclic_scaffold == ""

        bad_mol = mol_from_smiles("invalid_smiles")
        @test get_murcko_scaffold(bad_mol) === missing
        @test get_generic_scaffold(bad_mol) === missing
    end

    @testset "Fragment Counting and Splitting" begin
        mol = mol_from_smiles("CCO.CCC")

        count = get_fragment_count(mol)
        @test count == 2

        fragments = split_fragments(mol)
        @test isa(fragments, Vector)
        @test length(fragments) == 2
        @test all(f -> f.valid, fragments)

        # Check that we can convert fragments back to SMILES
        frag_smiles = [mol_to_smiles(f) for f in fragments]
        @test all(s -> !ismissing(s), frag_smiles)
        @test Set(frag_smiles) == Set(["CCO", "CCC"])

        largest = get_largest_fragment(mol)
        @test largest !== missing
        @test largest.valid
        # CCO and CCC have same number of heavy atoms, so either could be "largest"
        largest_smiles = mol_to_smiles(largest)
        @test largest_smiles in ["CCO", "CCC"]

        connected_mol = mol_from_smiles("CCCc1ccccc1")
        @test get_fragment_count(connected_mol) == 1
        @test length(split_fragments(connected_mol)) == 1
        @test mol_to_smiles(get_largest_fragment(connected_mol)) == mol_to_smiles(connected_mol)

        bad_mol = mol_from_smiles("invalid_smiles")
        @test get_fragment_count(bad_mol) === missing
        @test split_fragments(bad_mol) === missing
        @test get_largest_fragment(bad_mol) === missing
    end

    @testset "Bond-based Fragmentation" begin
        mol = mol_from_smiles("CCCC")

        # Fragment at bond index 1 (second bond, 0-based indexing)
        fragments = fragment_by_bonds(mol, [1])
        @test isa(fragments, Vector{String})
        @test length(fragments) >= 2  # Should produce at least 2 fragments

        bad_mol = mol_from_smiles("invalid_smiles")
        @test fragment_by_bonds(bad_mol, [0]) === missing
    end

    @testset "Vectorized Operations" begin
        # Test vectorized BRICS decomposition
        smiles_list = ["CCCOCCc1ccccc1", "CC(=O)NCCc1ccccc1", "CCCC"]
        mols = mol_from_smiles(smiles_list)

        all_brics = brics_decompose(mols)
        @test length(all_brics) == 3
        @test all(frags -> isa(frags, Vector{String}), all_brics)

        # Test vectorized scaffold extraction
        all_scaffolds = get_murcko_scaffold(mols)
        @test length(all_scaffolds) == 3
        @test all_scaffolds[1] == "c1ccccc1"  # First molecule has benzene scaffold
        @test all_scaffolds[2] == "c1ccccc1"  # Second molecule has benzene scaffold
        @test all_scaffolds[3] === missing || all_scaffolds[3] == ""  # Third molecule has no rings

        # Test with mixed valid/invalid molecules
        mixed_smiles = ["CCO.CCC", "invalid", "CCCC"]
        mixed_mols = mol_from_smiles(mixed_smiles)
        mixed_counts = get_fragment_count(mixed_mols)

        @test length(mixed_counts) == 3
        @test mixed_counts[1] == 2
        @test mixed_counts[2] === missing
        @test mixed_counts[3] == 1
    end
end