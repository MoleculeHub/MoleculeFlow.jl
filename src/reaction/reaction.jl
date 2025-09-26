#######################################################
# Chemical Reactions 
#######################################################

"""
    Reaction

Represents a chemical reaction with comprehensive functionality

# Fields

  - `_rdkit_rxn::Py`: The underlying RDKit reaction object
  - `props::Dict{Symbol, Any}`: Dictionary for storing additional reaction properties
  - `validated::Bool`: Whether the reaction has been validated
  - `fingerprint::Union{Nothing, BitVector}`: Cached reaction fingerprint
"""
@kwdef mutable struct Reaction
    _rdkit_rxn::Py
    props::Dict{Symbol, Any} = Dict()
    validated::Bool = false
    fingerprint::Union{Nothing, BitVector} = nothing
end

function Base.show(io::IO, rxn::Reaction)
    smarts = get(rxn.props, :SMARTS, "unknown")
    status = rxn.validated ? "validated" : "unvalidated"
    print(io, "Reaction($smarts, $status)")
end

"""
    reaction_from_smarts(smarts::String; validate::Bool=true) -> Reaction

Create a Reaction object from a SMARTS string.

# Arguments

  - `smarts::String`: A valid reaction SMARTS string with atom mapping (e.g., `"[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]"`)
  - `validate::Bool=true`: Whether to validate and sanitize the reaction after creation

# Returns

  - `Reaction`: A validated reaction object ready for use

# Examples

```julia
# Basic ester hydrolysis reaction
rxn = reaction_from_smarts(
    "\\[C:1\\](=O)\\[O:2\\]\\[C:3\\]>>\\[C:1\\](=O)\\[O-\\].\\[C:3\\]\\[O+\\]"
)

# Create without validation (faster, but potentially unsafe)
rxn = reaction_from_smarts(smarts; validate = false)
```

# Notes

  - Atom mapping (numbers after colons) is highly recommended for predictable results
  - Invalid SMARTS strings will throw an error
  - Validation includes sanitization and structural checks
"""
function reaction_from_smarts(smarts::String; validate::Bool = true)
    rxn = _reaction_from_smarts(smarts)
    if pynot(rxn)
        error("Invalid reaction SMARTS: $smarts")
    end
    reaction = Reaction(; _rdkit_rxn = rxn, props = Dict(:SMARTS => smarts))
    if validate
        validate_reaction!(reaction)
    end
    return reaction
end

"""
    reaction_from_rxn_file(filename::String; validate::Bool=true) -> Reaction

Create a Reaction object from an RXN file.
"""
function reaction_from_rxn_file(filename::String; validate::Bool = true)
    rxn = _reaction_from_rxn_file(filename)
    if pynot(rxn)
        error("Invalid RXN file: $filename")
    end
    reaction = Reaction(; _rdkit_rxn = rxn, props = Dict(:source => filename))
    if validate
        validate_reaction!(reaction)
    end
    return reaction
end

"""
    reaction_from_rxn_block(rxnblock::String; validate::Bool=true) -> Reaction

Create a Reaction object from an RXN block string.
"""
function reaction_from_rxn_block(rxnblock::String; validate::Bool = true)
    rxn = _reaction_from_rxn_block(rxnblock)
    if pynot(rxn)
        error("Invalid RXN block")
    end
    reaction = Reaction(; _rdkit_rxn = rxn, props = Dict(:source => "rxn_block"))
    if validate
        validate_reaction!(reaction)
    end
    return reaction
end

"""
    run_reaction(rxn::Reaction, reactants::Vector{Molecule}; max_products::Int=1000) -> Vector{Vector{Molecule}}

Apply a reaction to a set of reactants and generate products.

# Arguments

  - `rxn::Reaction`: The reaction to apply
  - `reactants::Vector{Molecule}`: Vector of reactant molecules
  - `max_products::Int=1000`: Maximum number of product sets to generate

# Returns

  - `Vector{Vector{Molecule}}`: Vector of product sets, where each set represents one possible reaction outcome

# Examples

```julia
# Apply ester hydrolysis to ethyl acetate
rxn = reaction_from_smarts(
    "\\[C:1\\](=O)\\[O:2\\]\\[C:3\\]>>\\[C:1\\](=O)\\[O-\\].\\[C:3\\]\\[O+\\]"
)
reactant = mol_from_smiles("CC(=O)OCC")
products = run_reaction(rxn, [reactant])

# Process all product sets
for (i, product_set) in enumerate(products)
    println("Product set \$i:")
    for product in product_set
        println("  \$(mol_to_smiles(product))")
    end
end
```

# Notes

  - Returns empty vector if no products can be formed
  - Each product set represents one possible stereochemical or regioisomeric outcome
  - Unvalidated reactions will show a warning but still execute
  - Use `has_reactant_substructure_match()` to pre-check compatibility
"""
function run_reaction(rxn::Reaction, reactants::Vector{Molecule}; max_products::Int = 1000)
    if !rxn.validated
        @warn "Running unvalidated reaction. Consider calling validate_reaction! first."
    end

    py_reactants = [r._rdkit_mol for r in reactants]
    product_sets = _reaction_run_reactants_inline_properties(
        rxn._rdkit_rxn, py_reactants, max_products
    )
    return [
        [
            Molecule(; _rdkit_mol = prod, valid = !pynot(prod), source = "reaction") for
            prod in set
        ] for set in product_sets
    ]
end

"""
    reaction_to_smarts(rxn::Reaction) -> String

Export a Reaction object to a SMARTS string.
"""
function reaction_to_smarts(rxn::Reaction)
    return string(_reaction_to_smarts(rxn._rdkit_rxn))
end

"""
    reaction_to_rxn_block(rxn::Reaction) -> String

Export a Reaction object to an RXN block string.
"""
function reaction_to_rxn_block(rxn::Reaction)
    return string(_reaction_to_rxn_block(rxn._rdkit_rxn))
end

#######################################################
# Reaction Validation
#######################################################

"""
    validate_reaction!(rxn::Reaction) -> Bool

Validate and sanitize a reaction. Returns true if successful.
"""
function validate_reaction!(rxn::Reaction)
    try
        _reaction_validate(rxn._rdkit_rxn)
        rxn.validated = true
        return true
    catch e
        rxn.validated = false
        error("Reaction validation failed: $e")
    end
end

"""
    is_reaction_valid(rxn::Reaction) -> Bool

Check if a reaction is valid without modifying it.
"""
function is_reaction_valid(rxn::Reaction)
    try
        # Create a copy and try to validate it
        test_rxn = _reaction_from_smarts(reaction_to_smarts(rxn))
        _reaction_validate(test_rxn)
        return true
    catch
        return false
    end
end

#######################################################
# Reaction Templates and Matching
#######################################################

"""
    get_num_reactant_templates(rxn::Reaction) -> Int

Get the number of reactant templates in the reaction.
"""
function get_num_reactant_templates(rxn::Reaction)
    return pyconvert(Int, _reaction_get_num_reactant_templates(rxn._rdkit_rxn))
end

"""
    get_num_product_templates(rxn::Reaction) -> Int

Get the number of product templates in the reaction.
"""
function get_num_product_templates(rxn::Reaction)
    return pyconvert(Int, _reaction_get_num_product_templates(rxn._rdkit_rxn))
end

"""
    get_reactant_template(rxn::Reaction, idx::Int) -> Molecule

Get a reactant template by index (0-based).
"""
function get_reactant_template(rxn::Reaction, idx::Int)
    template = _reaction_get_reactant_template(rxn._rdkit_rxn, idx)
    return Molecule(;
        _rdkit_mol = template, valid = !pynot(template), source = "reaction_template"
    )
end

"""
    get_product_template(rxn::Reaction, idx::Int) -> Molecule

Get a product template by index (0-based).
"""
function get_product_template(rxn::Reaction, idx::Int)
    template = _reaction_get_product_template(rxn._rdkit_rxn, idx)
    return Molecule(;
        _rdkit_mol = template, valid = !pynot(template), source = "reaction_template"
    )
end

"""
    has_reactant_substructure_match(rxn::Reaction, mol::Molecule) -> Bool

Check if a molecule matches any reactant template in the reaction.
"""
function has_reactant_substructure_match(rxn::Reaction, mol::Molecule)
    return pyconvert(
        Bool, _reaction_has_reactant_substructure_match(rxn._rdkit_rxn, mol._rdkit_mol)
    )
end

"""
    get_reacting_atoms(rxn::Reaction) -> Vector{Tuple{Int, Int}}

Get the indices of atoms that participate in the reaction for each reactant template.
Returns a vector of tuples where each tuple contains (template_idx, atom_idx).
"""
function get_reacting_atoms(rxn::Reaction)
    reacting_atoms = _reaction_get_reacting_atoms(rxn._rdkit_rxn)
    result = Tuple{Int, Int}[]
    for template_atoms in reacting_atoms
        for atom_idx in template_atoms
            push!(result, (length(result), pyconvert(Int, atom_idx)))
        end
    end
    return result
end

#######################################################
# Reaction Fingerprinting
#######################################################

"""
    reaction_fingerprint(rxn::Reaction; fp_size::Int=2048, use_cache::Bool=true) -> BitVector

Generate a difference fingerprint for the reaction, capturing the chemical transformation.

# Arguments

  - `rxn::Reaction`: The reaction to fingerprint
  - `fp_size::Int=2048`: Size of the fingerprint bit vector
  - `use_cache::Bool=true`: Whether to cache the fingerprint for reuse

# Returns

  - `BitVector`: Fingerprint representing the chemical transformation

# Examples

```julia
rxn1 = reaction_from_smarts(
    "\\[C:1\\](=O)\\[O:2\\]\\[C:3\\]>>\\[C:1\\](=O)\\[O-\\].\\[C:3\\]\\[O+\\]"
)
rxn2 = reaction_from_smarts("\\[C:1\\]\\[OH:2\\]>>\\[C:1\\]\\[O-\\]")

fp1 = reaction_fingerprint(rxn1)
fp2 = reaction_fingerprint(rxn2)

# Calculate Tanimoto similarity
similarity = sum(fp1 .& fp2) / sum(fp1 .| fp2)
```

# Notes

  - Difference fingerprints focus on the atoms and bonds that change during the reaction
  - Cached fingerprints are automatically retrieved on subsequent calls
  - Use for reaction similarity analysis and database searching
  - Different from structural fingerprints which capture overall molecular features
"""
function reaction_fingerprint(rxn::Reaction; fp_size::Int = 2048, use_cache::Bool = true)
    if use_cache && rxn.fingerprint !== nothing && length(rxn.fingerprint) == fp_size
        return rxn.fingerprint
    end

    fp = _reaction_fingerprint(rxn._rdkit_rxn, fp_size)

    # Handle different fingerprint types
    if pyhasattr(fp, "GetBit")
        # ExplicitBitVect - can access bits directly
        fp_bits = BitVector([pyconvert(Bool, fp.GetBit(i)) for i in 0:(fp_size - 1)])
    else
        # UIntSparseIntVect - convert to dense representation
        fp_bits = BitVector(zeros(Bool, fp_size))
        nz_elements = fp.GetNonzeroElements()
        py_dict = pyconvert(Dict, nz_elements)
        for (idx, val) in py_dict
            bit_idx = idx + 1  # Convert to 1-based indexing
            if bit_idx <= fp_size && val != 0
                fp_bits[bit_idx] = true
            end
        end
    end

    if use_cache
        rxn.fingerprint = fp_bits
    end

    return fp_bits
end

"""
    reaction_structural_fingerprint(rxn::Reaction; fp_size::Int=2048) -> BitVector

Generate a structural fingerprint for the reaction.
"""
function reaction_structural_fingerprint(rxn::Reaction; fp_size::Int = 2048)
    fp = _reaction_structural_fingerprint(rxn._rdkit_rxn, fp_size)
    return BitVector([pyconvert(Bool, fp.GetBit(i)) for i in 0:(fp_size - 1)])
end

"""
    reaction_center_fingerprint(rxn::Reaction; fp_size::Int=2048) -> BitVector

Generate a reaction center fingerprint focusing on the reaction core.
"""
function reaction_center_fingerprint(rxn::Reaction; fp_size::Int = 2048)
    fp = _compute_reaction_center_fingerprint(rxn._rdkit_rxn, fp_size)

    # Handle different fingerprint types
    if pyhasattr(fp, "GetBit")
        # ExplicitBitVect - can access bits directly
        return BitVector([pyconvert(Bool, fp.GetBit(i)) for i in 0:(fp_size - 1)])
    else
        # UIntSparseIntVect - convert to dense representation
        fp_bits = BitVector(zeros(Bool, fp_size))
        nz_elements = fp.GetNonzeroElements()
        py_dict = pyconvert(Dict, nz_elements)
        for (idx, val) in py_dict
            bit_idx = idx + 1  # Convert to 1-based indexing
            if bit_idx <= fp_size && val != 0
                fp_bits[bit_idx] = true
            end
        end
        return fp_bits
    end
end

"""
    reaction_similarity(rxn1::Reaction, rxn2::Reaction; method::Symbol=:tanimoto) -> Float64

Calculate similarity between two reactions using their fingerprints.

# Arguments

  - `rxn1::Reaction`: First reaction
  - `rxn2::Reaction`: Second reaction
  - `method::Symbol=:tanimoto`: Similarity method (`:tanimoto` or `:dice`)

# Returns

  - `Float64`: Similarity score between 0.0 (no similarity) and 1.0 (identical)

# Examples

```julia
rxn1 = reaction_from_smarts(
    "\\[C:1\\](=O)\\[O:2\\]\\[C:3\\]>>\\[C:1\\](=O)\\[O-\\].\\[C:3\\]\\[O+\\]"
)
rxn2 = reaction_from_smarts(
    "\\[C:1\\](=O)\\[O:2\\]\\[C:3\\]>>\\[C:1\\](=O)\\[OH\\].\\[C:3\\]"
)

# Tanimoto similarity
tanimoto_sim = reaction_similarity(rxn1, rxn2; method = :tanimoto)

# Dice similarity
dice_sim = reaction_similarity(rxn1, rxn2; method = :dice)

println("Tanimoto: \$tanimoto_sim, Dice: \$dice_sim")
```

# Notes

  - Tanimoto similarity: intersection / union (most common)
  - Dice similarity: 2 * intersection / (|A| + |B|) (gives higher scores)
  - Uses cached fingerprints when available for better performance
  - Self-similarity should equal 1.0 for identical reactions
"""
function reaction_similarity(rxn1::Reaction, rxn2::Reaction; method::Symbol = :tanimoto)
    fp1 = reaction_fingerprint(rxn1)
    fp2 = reaction_fingerprint(rxn2)

    if method == :tanimoto
        intersection = sum(fp1 .& fp2)
        union = sum(fp1 .| fp2)
        return union > 0 ? intersection / union : 0.0
    elseif method == :dice
        intersection = sum(fp1 .& fp2)
        total = sum(fp1) + sum(fp2)
        return total > 0 ? 2 * intersection / total : 0.0
    else
        error("Unknown similarity method: $method")
    end
end

#######################################################
# Reaction Library and Enumeration
#######################################################

"""
    enumerate_library(rxn::Reaction, reactant_lists::Vector{Vector{Molecule}}) -> Vector{Vector{Molecule}}

Enumerate a library of products from lists of reactants.
"""
function enumerate_library(rxn::Reaction, reactant_lists::Vector{Vector{Molecule}})
    py_reactant_lists = [[mol._rdkit_mol for mol in mols] for mols in reactant_lists]
    product_sets = _reaction_enumerate_library_from_reaction(
        rxn._rdkit_rxn, py_reactant_lists
    )
    return [
        [Molecule(; _rdkit_mol = prod, valid = true, source = "library") for prod in set]
        for set in product_sets
    ]
end

#######################################################
# Reaction Analysis
#######################################################

"""
    reaction_info(rxn::Reaction) -> Dict{Symbol, Any}

Get comprehensive information about a reaction.
"""
function reaction_info(rxn::Reaction)
    info = Dict{Symbol, Any}(
        :smarts => reaction_to_smarts(rxn),
        :validated => rxn.validated,
        :num_reactant_templates => get_num_reactant_templates(rxn),
        :num_product_templates => get_num_product_templates(rxn),
        :has_fingerprint => rxn.fingerprint !== nothing,
        :properties => copy(rxn.props),
    )
    return info
end

"""
    is_balanced(rxn::Reaction) -> Bool

Check if a reaction is atom-balanced (same atoms on both sides).
Note: This is a simplified implementation.
"""
function is_balanced(rxn::Reaction)
    try
        # Get reactant and product templates
        num_reactants = get_num_reactant_templates(rxn)
        num_products = get_num_product_templates(rxn)

        # Simple check: same number of heavy atoms on both sides
        reactant_atoms = 0
        product_atoms = 0

        # Count reactant atoms
        for i in 0:(num_reactants - 1)
            template = get_reactant_template(rxn, i)
            reactant_atoms += heavy_atom_count(template)
        end

        # Count product atoms
        for i in 0:(num_products - 1)
            template = get_product_template(rxn, i)
            product_atoms += heavy_atom_count(template)
        end

        return reactant_atoms == product_atoms
    catch
        return false
    end
end

#######################################################
# Atom Mapping and Stereochemistry
#######################################################

"""
    get_atom_mapping_numbers(mol::Molecule) -> Vector{Int}

Get atom mapping numbers for all atoms in a molecule.
"""
function get_atom_mapping_numbers(mol::Molecule)
    py_result = _get_atom_mapping_numbers(mol._rdkit_mol)
    return [pyconvert(Int, x) for x in py_result]
end

"""
    set_atom_mapping_numbers!(mol::Molecule, map_nums::Vector{Int})

Set atom mapping numbers for all atoms in a molecule.
"""
function set_atom_mapping_numbers!(mol::Molecule, map_nums::Vector{Int})
    _set_atom_mapping_numbers(mol._rdkit_mol, map_nums)
end

"""
    remove_unmapped_reactant_templates!(rxn::Reaction; mode::Int=1)

Remove reactant templates that don't have mapped atoms.
"""
function remove_unmapped_reactant_templates!(rxn::Reaction; mode::Int = 1)
    _reaction_remove_unmapped_reactant_templates(rxn._rdkit_rxn, mode)
    rxn.validated = false  # Need to revalidate after modification
end

"""
    remove_unmapped_product_templates!(rxn::Reaction; mode::Int=1)

Remove product templates that don't have mapped atoms.
"""
function remove_unmapped_product_templates!(rxn::Reaction; mode::Int = 1)
    _reaction_remove_unmapped_product_templates(rxn._rdkit_rxn, mode)
    rxn.validated = false  # Need to revalidate after modification
end

"""
    preprocess_reaction!(rxn::Reaction)

Preprocess a reaction to standardize it.
"""
function preprocess_reaction!(rxn::Reaction)
    _reaction_preprocess(rxn._rdkit_rxn)
    rxn.validated = false  # Need to revalidate after modification
end

"""
    compute_atom_mapping!(rxn::Reaction)

Compute and apply atom mapping to the reaction.
"""
function compute_atom_mapping!(rxn::Reaction)
    _reaction_compute_atom_mapping(rxn._rdkit_rxn)
    rxn.validated = false  # Need to revalidate after modification
end

"""
    is_template_molecule_agent(mol::Molecule) -> Bool

Check if a molecule is an agent in a reaction template.
"""
function is_template_molecule_agent(mol::Molecule)
    return pyconvert(Bool, _reaction_is_template_molecule_agent(mol._rdkit_mol))
end

"""
    sanitize_reaction!(rxn::Reaction; ops::Int=15) -> Bool

Sanitize a reaction with specific operations. Returns true if successful.
"""
function sanitize_reaction!(rxn::Reaction; ops::Int = 15)
    try
        _reaction_sanitize_reaction(rxn._rdkit_rxn, ops)
        rxn.validated = true
        return true
    catch e
        rxn.validated = false
        error("Reaction sanitization failed: $e")
    end
end

#######################################################
# Advanced Reaction Analysis
#######################################################

"""
    reaction_complexity(rxn::Reaction) -> Float64

Calculate a complexity score for the reaction based on various factors.
"""
function reaction_complexity(rxn::Reaction)
    num_reactants = get_num_reactant_templates(rxn)
    num_products = get_num_product_templates(rxn)

    # Base complexity from number of components
    complexity = (num_reactants + num_products) / 2.0

    # Add complexity from molecular size
    total_atoms = 0
    for i in 0:(num_reactants - 1)
        template = get_reactant_template(rxn, i)
        total_atoms += heavy_atom_count(template)
    end
    for i in 0:(num_products - 1)
        template = get_product_template(rxn, i)
        total_atoms += heavy_atom_count(template)
    end

    complexity += total_atoms / 20.0  # Normalize by typical molecule size

    return complexity
end

"""
    reaction_type_classification(rxn::Reaction) -> Symbol

Classify the reaction type based on structural changes.
"""
function reaction_type_classification(rxn::Reaction)
    num_reactants = get_num_reactant_templates(rxn)
    num_products = get_num_product_templates(rxn)

    if num_reactants == 1 && num_products == 2
        return :decomposition
    elseif num_reactants == 2 && num_products == 1
        return :combination
    elseif num_reactants == 2 && num_products == 2
        return :substitution
    elseif num_reactants == 1 && num_products == 1
        return :rearrangement
    else
        return :complex
    end
end

"""
    find_similar_reactions(target_rxn::Reaction, reaction_db::Vector{Reaction};
    					  threshold::Float64=0.7, method::Symbol=:tanimoto) -> Vector{Tuple{Reaction, Float64}}

Find reactions similar to the target reaction from a database.
"""
function find_similar_reactions(
    target_rxn::Reaction,
    reaction_db::Vector{Reaction};
    threshold::Float64 = 0.7,
    method::Symbol = :tanimoto,
)
    similar_reactions = Tuple{Reaction, Float64}[]

    for rxn in reaction_db
        similarity = reaction_similarity(target_rxn, rxn; method = method)
        if similarity >= threshold
            push!(similar_reactions, (rxn, similarity))
        end
    end

    # Sort by similarity (descending)
    sort!(similar_reactions; by = x->x[2], rev = true)

    return similar_reactions
end
