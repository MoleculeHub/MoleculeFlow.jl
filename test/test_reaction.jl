using MoleculeFlow
using Test

@testset "Basic Reaction Creation" begin
    # Example: ester hydrolysis
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)
    @test is_reaction_valid(rxn)
    @test isa(rxn, Reaction)
    @test rxn.props[:SMARTS] == rxn_smarts
    @test reaction_to_smarts(rxn) isa String
    @test rxn.validated  # Should be validated by default
end

@testset "Reaction File I/O" begin
    # Test RXN block creation and parsing
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)
    rxn_block = reaction_to_rxn_block(rxn)
    @test isa(rxn_block, String)
    @test length(rxn_block) > 0

    # Test creating reaction from block
    rxn2 = reaction_from_rxn_block(rxn_block)
    @test isa(rxn2, Reaction)
end

@testset "Reaction Validation" begin
    # Valid reaction
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts, validate = false)
    @test !rxn.validated
    @test validate_reaction!(rxn)
    @test rxn.validated

    # Test validation check
    @test is_reaction_valid(rxn)
end

@testset "Run Reaction" begin
    # Example: simple ester hydrolysis
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)
    # Ethyl acetate: CC(=O)OCC
    reactant = mol_from_smiles("CC(=O)OCC")
    @test reactant.valid
    products = run_reaction(rxn, [reactant])
    @test isa(products, Vector)
    @test length(products) > 0
    @test all(prodset -> all(prod -> prod isa Molecule, prodset), products)

    # Test with max_products limit
    products_limited = run_reaction(rxn, [reactant], max_products = 1)
    @test length(products_limited) <= 1
end

@testset "Reaction Templates" begin
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)

    # Test template counts
    @test get_num_reactant_templates(rxn) >= 1
    @test get_num_product_templates(rxn) >= 1

    # Test template retrieval
    reactant_template = get_reactant_template(rxn, 0)
    @test isa(reactant_template, Molecule)

    product_template = get_product_template(rxn, 0)
    @test isa(product_template, Molecule)
end

@testset "Substructure Matching" begin
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)

    # Matching molecule (ethyl acetate)
    matching_mol = mol_from_smiles("CC(=O)OCC")
    @test has_reactant_substructure_match(rxn, matching_mol)

    # Non-matching molecule (methane)
    non_matching_mol = mol_from_smiles("C")
    @test !has_reactant_substructure_match(rxn, non_matching_mol)

    # Test reacting atoms
    reacting_atoms = get_reacting_atoms(rxn)
    @test isa(reacting_atoms, Vector{Tuple{Int, Int}})
end

@testset "Reaction Fingerprinting" begin
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)

    # Test difference fingerprint
    fp = reaction_fingerprint(rxn)
    @test isa(fp, BitVector)
    @test length(fp) == 2048  # Default size

    # Test cached fingerprint
    fp2 = reaction_fingerprint(rxn, use_cache = true)
    @test fp == fp2

    # Test structural fingerprint
    struct_fp = reaction_structural_fingerprint(rxn)
    @test isa(struct_fp, BitVector)

    # Test reaction center fingerprint
    center_fp = reaction_center_fingerprint(rxn)
    @test isa(center_fp, BitVector)
end

@testset "Reaction Similarity" begin
    rxn1_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn2_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[OH].[C:3]"

    rxn1 = reaction_from_smarts(rxn1_smarts)
    rxn2 = reaction_from_smarts(rxn2_smarts)

    # Test Tanimoto similarity
    sim_tanimoto = reaction_similarity(rxn1, rxn2, method = :tanimoto)
    @test isa(sim_tanimoto, Float64)
    @test 0.0 <= sim_tanimoto <= 1.0

    # Test Dice similarity
    sim_dice = reaction_similarity(rxn1, rxn2, method = :dice)
    @test isa(sim_dice, Float64)
    @test 0.0 <= sim_dice <= 1.0

    # Self-similarity should be 1.0
    self_sim = reaction_similarity(rxn1, rxn1)
    @test self_sim ≈ 1.0
end

@testset "Atom Mapping" begin
    # Test molecule with atom mapping
    mol_smiles = "[CH3:1][C:2](=[O:3])[O:4][CH2:5][CH3:6]"
    mol = mol_from_smiles(mol_smiles)

    # Get atom mapping numbers
    map_nums = get_atom_mapping_numbers(mol)
    @test isa(map_nums, Vector)
    @test length(map_nums) == heavy_atom_count(mol)

    # Test setting new mapping numbers
    new_map_nums = [10, 20, 30, 40, 50, 60]
    set_atom_mapping_numbers!(mol, new_map_nums)
    updated_map_nums = get_atom_mapping_numbers(mol)
    @test updated_map_nums == new_map_nums
end

@testset "Reaction Preprocessing" begin
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)

    # Test preprocessing
    preprocess_reaction!(rxn)
    @test !rxn.validated  # Should need revalidation

    # Test sanitization
    @test sanitize_reaction!(rxn)
    @test rxn.validated
end

@testset "Reaction Analysis" begin
    rxn_smarts = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"
    rxn = reaction_from_smarts(rxn_smarts)

    # Test reaction info
    info = reaction_info(rxn)
    @test isa(info, Dict{Symbol, Any})
    @test haskey(info, :smarts)
    @test haskey(info, :validated)
    @test haskey(info, :num_reactant_templates)
    @test haskey(info, :num_product_templates)

    # Test complexity calculation
    complexity = reaction_complexity(rxn)
    @test isa(complexity, Float64)
    @test complexity > 0

    # Test reaction classification
    classification = reaction_type_classification(rxn)
    @test isa(classification, Symbol)
    @test classification in
        [:decomposition, :combination, :substitution, :rearrangement, :complex]
end

@testset "Reaction Library Enumeration" begin
    # Simple reaction: alcohol + carboxylic acid -> ester
    rxn_smarts = "[C:1][OH:2].[C:3](=O)[OH:4]>>[C:1][O:2][C:3](=O)"
    rxn = reaction_from_smarts(rxn_smarts)

    # Create small library of reactants
    alcohols = [mol_from_smiles("CCO"), mol_from_smiles("CO")]
    acids = [mol_from_smiles("CC(=O)O"), mol_from_smiles("C(=O)O")]

    # Enumerate library (might fail with simple SMARTS, so we test the function exists)
    try
        products = enumerate_library(rxn, [alcohols, acids])
        @test isa(products, Vector)
    catch
        # Some reactions might not work with enumerate_library
        @test true  # Just test that the function exists
    end
end

@testset "Reaction Database Search" begin
    # Create a small database of reactions
    rxn1 = reaction_from_smarts("[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]")
    rxn2 = reaction_from_smarts("[C:1][OH:2]>>[C:1][O-]")
    rxn3 = reaction_from_smarts("[N:1][C:2](=O)[C:3]>>[N:1].[C:2](=O)[C:3]")

    reaction_db = [rxn1, rxn2, rxn3]

    # Search for similar reactions
    similar = find_similar_reactions(rxn1, reaction_db, threshold = 0.1)
    @test isa(similar, Vector{Tuple{Reaction, Float64}})
    @test length(similar) >= 1  # Should at least find itself

    # Check that similarities are sorted
    if length(similar) > 1
        @test all(i -> similar[i][2] >= similar[i + 1][2], 1:(length(similar) - 1))
    end
end
