#######################################################
# Molecular Editing and Manipulation
#######################################################

"""
    combine_mols(mol1::Molecule, mol2::Molecule) -> Molecule

Combine two molecules into a single molecule.

# Arguments
- `mol1::Molecule`: First molecule
- `mol2::Molecule`: Second molecule

# Returns
- `Molecule`: Combined molecule

# Example
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
combined = combine_mols(mol1, mol2)
```
"""
function combine_mols(mol1::Molecule, mol2::Molecule)
    if !mol1.valid || !mol2.valid
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = "combined")
    end

    try
        combined_mol = _combine_mols(mol1._rdkit_mol, mol2._rdkit_mol)
        return Molecule(_rdkit_mol = combined_mol, valid = true, source = "combined")
    catch e
        @warn "Error combining molecules: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = "combined")
    end
end

"""
    delete_substructs(mol::Molecule, query::String) -> Molecule

Delete substructures matching a SMARTS pattern from a molecule.

# Arguments
- `mol::Molecule`: Input molecule
- `query::String`: SMARTS pattern to match and delete

# Returns
- `Molecule`: Molecule with matching substructures removed

# Example
```julia
mol = mol_from_smiles("CCO")
# Remove alcohol groups
mol_modified = delete_substructs(mol, "[OH]")
```
"""
function delete_substructs(mol::Molecule, query::String)
    if !mol.valid
        return mol
    end

    try
        query_mol = _mol_from_smarts(query)
        if pyconvert(Bool, query_mol === pybuiltins.None)
            @warn "Invalid SMARTS pattern: $query"
            return mol
        end

        modified_mol = _delete_substructs(mol._rdkit_mol, query_mol)
        return Molecule(_rdkit_mol = modified_mol, valid = true, source = mol.source)
    catch e
        @warn "Error deleting substructures: $e"
        return mol
    end
end

"""
    replace_substructs(mol::Molecule, query::String, replacement::String) -> Vector{Molecule}

Replace substructures matching a SMARTS pattern with a replacement structure.

# Arguments
- `mol::Molecule`: Input molecule
- `query::String`: SMARTS pattern to match
- `replacement::String`: SMILES string for replacement

# Returns
- `Vector{Molecule}`: Vector of molecules with replacements made

# Example
```julia
mol = mol_from_smiles("CCO")
# Replace alcohol with amine
replacements = replace_substructs(mol, "[OH]", "N")
```
"""
# Ring analysis functions
"""
    get_ring_info(mol::Molecule) -> Union{Py, Missing}

Get ring information for a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{Py, Missing}`: RingInfo object or missing if molecule is invalid

# Example
```julia
mol = mol_from_smiles("c1ccccc1")
ring_info = get_ring_info(mol)
```
"""
function get_ring_info(mol::Molecule)
    !mol.valid && return missing
    try
        ring_info = _get_ring_info(mol._rdkit_mol)
        return Dict(
            :num_rings => pyconvert(Int, ring_info.NumRings()),
            :atom_rings => [pyconvert(Vector{Int}, ring) .+ 1 for ring in ring_info.AtomRings()],
            :bond_rings => [pyconvert(Vector{Int}, ring) .+ 1 for ring in ring_info.BondRings()],
        )
    catch e
        @warn "Error getting ring info: $e"
        return missing
    end
end

"""
    canonical_atom_ranks(mol::Molecule) -> Union{Vector{Int}, Missing}

Get canonical atom ranks for a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{Vector{Int}, Missing}`: Vector of atom ranks or missing if molecule is invalid

# Example
```julia
mol = mol_from_smiles("CCO")
ranks = canonical_atom_ranks(mol)
```
"""
function canonical_atom_ranks(mol::Molecule)
    !mol.valid && return Int[]
    try
        py_ranks = _canonical_rank_atoms(mol._rdkit_mol)
        return pyconvert(Vector{Int}, py_ranks)
    catch e
        @warn "Error getting canonical atom ranks: $e"
        return Int[]
    end
end

"""
    find_atom_environment(mol::Molecule, radius::Int, atom_idx::Int) -> Union{Vector{Int}, Missing}

Find atom environment of specified radius around an atom.

# Arguments
- `mol::Molecule`: Input molecule
- `radius::Int`: Radius to search
- `atom_idx::Int`: Index of central atom (0-based)

# Returns
- `Union{Vector{Int}, Missing}`: Vector of bond indices in environment or missing if invalid

# Example
```julia
mol = mol_from_smiles("CCCCC")
env = find_atom_environment(mol, 2, 1)
```
"""
function find_atom_environment(mol::Molecule, radius::Int, atom_idx::Int)
    !mol.valid && return missing
    try
        py_env = _find_atom_environment_of_radius_n(mol._rdkit_mol, radius, atom_idx)
        return pyconvert(Vector{Int}, py_env)
    catch e
        @warn "Error finding atom environment: $e"
        return missing
    end
end

"""
    mol_fragment_to_smarts(mol::Molecule, atom_indices::Vector{Int}) -> Union{String, Missing}

Convert a molecular fragment to SMARTS pattern.

# Arguments
- `mol::Molecule`: Input molecule
- `atom_indices::Vector{Int}`: Indices of atoms to include in fragment

# Returns
- `Union{String, Missing}`: SMARTS pattern or missing if invalid

# Example
```julia
mol = mol_from_smiles("CCCO")
smarts = mol_fragment_to_smarts(mol, [0, 1])  # First two carbons
```
"""
function mol_fragment_to_smarts(mol::Molecule, atom_indices::Vector{Int})
    !mol.valid && return missing
    try
        py_smarts = _mol_fragment_to_smarts(mol._rdkit_mol, atom_indices)
        return pyconvert(String, py_smarts)
    catch e
        @warn "Error converting fragment to SMARTS: $e"
        return missing
    end
end

"""
    mol_fragment_to_smiles(mol::Molecule, atom_indices::Vector{Int}) -> Union{String, Missing}

Convert a molecular fragment to SMILES string.

# Arguments
- `mol::Molecule`: Input molecule
- `atom_indices::Vector{Int}`: Indices of atoms to include in fragment

# Returns
- `Union{String, Missing}`: SMILES string or missing if invalid

# Example
```julia
mol = mol_from_smiles("CCCO")
smiles = mol_fragment_to_smiles(mol, [0, 1])  # First two carbons
```
"""
function mol_fragment_to_smiles(mol::Molecule, atom_indices::Vector{Int})
    !mol.valid && return missing
    try
        py_smiles = _mol_fragment_to_smiles(mol._rdkit_mol, atom_indices)
        return pyconvert(String, py_smiles)
    catch e
        @warn "Error converting fragment to SMILES: $e"
        return missing
    end
end

"""
    renumber_atoms(mol::Molecule, new_order::Vector{Int}) -> Molecule

Renumber atoms in a molecule according to a new ordering.

# Arguments
- `mol::Molecule`: Input molecule
- `new_order::Vector{Int}`: New atom ordering (0-based indices)

# Returns
- `Molecule`: Molecule with renumbered atoms

# Example
```julia
mol = mol_from_smiles("CCO")
# Reverse atom order
new_mol = renumber_atoms(mol, [2, 1, 0])
```
"""
function renumber_atoms(mol::Molecule, new_order::Vector{Int})
    !mol.valid && return mol
    try
        new_mol = _renumber_atoms(mol._rdkit_mol, new_order)
        return Molecule(_rdkit_mol = new_mol, valid = true, source = mol.source)
    catch e
        @warn "Error renumbering atoms: $e"
        return mol
    end
end

"""
    remove_stereochemistry!(mol::Molecule) -> Molecule

Remove stereochemistry information from a molecule (in-place operation).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Modified molecule

# Example
```julia
mol = mol_from_smiles("C[C@H](O)C")
remove_stereochemistry!(mol)
```
"""
function remove_stereochemistry!(mol::Molecule)
    !mol.valid && return mol
    try
        _remove_stereochemistry(mol._rdkit_mol)
        return mol
    catch e
        @warn "Error removing stereochemistry: $e"
        return mol
    end
end

"""
    sanitize_mol!(mol::Molecule) -> Molecule

Sanitize a molecule (in-place operation).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Sanitized molecule

# Example
```julia
mol = mol_from_smiles("CCO")
sanitize_mol!(mol)
```
"""
function sanitize_mol!(mol::Molecule)
    !mol.valid && return mol
    try
        _sanitize_mol(mol._rdkit_mol)
        return mol
    catch e
        @warn "Error sanitizing molecule: $e"
        return mol
    end
end

"""
    fast_find_rings!(mol::Molecule) -> Molecule

Find rings in a molecule using the fast algorithm (in-place operation).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Molecule with ring information computed

# Example
```julia
mol = mol_from_smiles("c1ccccc1")
fast_find_rings!(mol)
```
"""
function fast_find_rings!(mol::Molecule)
    !mol.valid && return false
    try
        _fast_find_rings(mol._rdkit_mol)
        return true
    catch e
        @warn "Error finding rings: $e"
        return false
    end
end

"""
    compute_2d_coords!(mol::Molecule) -> Molecule

Compute 2D coordinates for a molecule (in-place operation).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Molecule with 2D coordinates

# Example
```julia
mol = mol_from_smiles("CCO")
compute_2d_coords!(mol)
```
"""
function compute_2d_coords!(mol::Molecule)
    !mol.valid && return mol
    try
        _compute_2d_coords(mol._rdkit_mol)
        return mol
    catch e
        @warn "Error computing 2D coordinates: $e"
        return mol
    end
end

function replace_substructs(mol::Molecule, query::String, replacement::String)
    if !mol.valid
        return [mol]
    end

    try
        query_mol = _mol_from_smarts(query)
        replacement_mol = _mol_from_smiles(replacement)

        if pyconvert(Bool, query_mol === pybuiltins.None) || pyconvert(Bool, replacement_mol === pybuiltins.None)
            @warn "Invalid SMARTS/SMILES pattern"
            return [mol]
        end

        result_mols = _replace_substructs(mol._rdkit_mol, query_mol, replacement_mol)

        # Convert Python tuple/list to Julia vector
        molecules = Molecule[]
        for rdkit_mol in result_mols
            push!(molecules, Molecule(_rdkit_mol = rdkit_mol, valid = true, source = mol.source))
        end

        return molecules
    catch e
        @warn "Error replacing substructures: $e"
        return [mol]
    end
end

#######################################################
# Stereochemistry Operations
#######################################################

"""
    assign_stereochemistry!(mol::Molecule; clean_it::Bool=true, force::Bool=false) -> Bool

Assign stereochemistry to a molecule in place.

# Arguments
- `mol::Molecule`: Input molecule (modified in place)
- `clean_it::Bool`: Whether to clean up the molecule
- `force::Bool`: Whether to force assignment

# Returns
- `Bool`: Success status

# Example
```julia
mol = mol_from_smiles("C[C@H](O)C")
success = assign_stereochemistry!(mol)
```
"""
function assign_stereochemistry!(mol::Molecule; clean_it::Bool=true, force::Bool=false)
    if !mol.valid
        return false
    end

    try
        _assign_stereochemistry(mol._rdkit_mol; cleanIt=clean_it, force=force)
        return true
    catch e
        @warn "Error assigning stereochemistry: $e"
        return false
    end
end

"""
    find_chiral_centers(mol::Molecule; include_unassigned::Bool=false) -> Vector{Tuple{Int,String}}

Find chiral centers in a molecule.

# Arguments
- `mol::Molecule`: Input molecule
- `include_unassigned::Bool`: Whether to include unassigned chiral centers

# Returns
- `Vector{Tuple{Int,String}}`: Vector of (atom_index, chirality) tuples

# Example
```julia
mol = mol_from_smiles("C[C@H](O)C")
chiral_centers = find_chiral_centers(mol)
```
"""
function find_chiral_centers(mol::Molecule; include_unassigned::Bool=false)
    if !mol.valid
        return Tuple{Int,String}[]
    end

    try
        centers = _find_mol_chiral_centers(mol._rdkit_mol; includeUnassigned=include_unassigned)
        result = Tuple{Int,String}[]

        for center in centers
            atom_idx = pyconvert(Int, center[0])
            chirality = pyconvert(String, center[1])
            push!(result, (atom_idx, chirality))
        end

        return result
    catch e
        @warn "Error finding chiral centers: $e"
        return Tuple{Int,String}[]
    end
end

#######################################################
# Ring Analysis
#######################################################


#######################################################
# Pattern Matching
#######################################################

"""
    quick_smarts_match(mol::Molecule, smarts::String) -> Bool

Quick check if a molecule matches a SMARTS pattern.

# Arguments
- `mol::Molecule`: Input molecule
- `smarts::String`: SMARTS pattern

# Returns
- `Bool`: Whether the molecule matches the pattern

# Example
```julia
mol = mol_from_smiles("CCO")
has_alcohol = quick_smarts_match(mol, "[OH]")
```
"""
function quick_smarts_match(mol::Molecule, smarts::String)
    if !mol.valid
        return false
    end

    try
        # Use the standard substructure matching instead of QuickSmartsMatch
        query_mol = _mol_from_smarts(smarts)
        if pyconvert(Bool, query_mol === pybuiltins.None)
            @warn "Invalid SMARTS pattern: $smarts"
            return false
        end

        return pyconvert(Bool, mol._rdkit_mol.HasSubstructMatch(query_mol))
    catch e
        @warn "Error in SMARTS matching: $e"
        return false
    end
end





#######################################################
# Additional Missing Operations
#######################################################



"""
    split_mol_by_pdb_chain_id(mol::Molecule) -> Union{Vector{Molecule}, Missing}

Split a PDB molecule by chain ID.

# Returns
- `Union{Vector{Molecule}, Missing}`: Vector of molecules (one per chain) or missing if operation fails
"""
function split_mol_by_pdb_chain_id(mol::Molecule)
    !mol.valid && return missing
    try
        mol_dict = _split_mol_by_pdb_chain_id(mol._rdkit_mol)
        molecules = Molecule[]
        for (chain_id, rdkit_mol) in mol_dict
            push!(molecules, Molecule(_rdkit_mol = rdkit_mol, valid = true, source = "PDB_chain_$chain_id"))
        end
        return molecules
    catch e
        @warn "Error splitting molecule by chain ID: $e"
        return missing
    end
end

"""
    split_mol_by_pdb_residues(mol::Molecule) -> Union{Vector{Molecule}, Missing}

Split a PDB molecule by residues.

# Returns
- `Union{Vector{Molecule}, Missing}`: Vector of molecules (one per residue) or missing if operation fails
"""
function split_mol_by_pdb_residues(mol::Molecule)
    !mol.valid && return missing
    try
        mol_list = _split_mol_by_pdb_residues(mol._rdkit_mol)
        molecules = Molecule[]
        for (i, rdkit_mol) in enumerate(mol_list)
            push!(molecules, Molecule(_rdkit_mol = rdkit_mol, valid = true, source = "PDB_residue_$i"))
        end
        return molecules
    catch e
        @warn "Error splitting molecule by residues: $e"
        return missing
    end
end

"""
    assign_stereochemistry_from_3d!(mol::Molecule; conf_id::Int = -1) -> Bool

Assign stereochemistry from 3D coordinates.

# Arguments
- `mol::Molecule`: Input molecule (modified in place)
- `conf_id::Int`: Conformer ID to use (-1 for default)

# Returns
- `Bool`: true if successful, false otherwise
"""
function assign_stereochemistry_from_3d!(mol::Molecule; conf_id::Int = -1)
    !mol.valid && return false
    try
        _assign_stereochemistry_from_3d(mol._rdkit_mol; confId = conf_id)
        return true
    catch e
        @warn "Error assigning stereochemistry from 3D: $e"
        return false
    end
end

"""
    detect_bond_stereochemistry(mol::Molecule, bond_idx::Int) -> Union{String, Missing}

Detect the stereochemistry of a specific bond.

# Arguments
- `mol::Molecule`: Input molecule
- `bond_idx::Int`: Bond index

# Returns
- `Union{String, Missing}`: Stereochemistry designation or missing if detection fails
"""
function detect_bond_stereochemistry(mol::Molecule, bond_idx::Int)
    !mol.valid && return missing
    try
        return pyconvert(String, _detect_bond_stereochemistry(mol._rdkit_mol, bond_idx))
    catch e
        @warn "Error detecting bond stereochemistry: $e"
        return missing
    end
end

"""
    find_potential_stereo(mol::Molecule) -> Union{Vector, Missing}

Find potential stereo centers in a molecule.

# Returns
- `Union{Vector, Missing}`: Vector of potential stereo information or missing if detection fails
"""
function find_potential_stereo(mol::Molecule)
    !mol.valid && return missing
    try
        return pyconvert(Vector, _find_potential_stereo(mol._rdkit_mol))
    catch e
        @warn "Error finding potential stereo: $e"
        return missing
    end
end

"""
    canonicalize_enhanced_stereo!(mol::Molecule) -> Bool

Canonicalize enhanced stereochemistry information.

# Returns
- `Bool`: true if successful, false otherwise
"""
function canonicalize_enhanced_stereo!(mol::Molecule)
    !mol.valid && return false
    try
        _canonicalize_enhanced_stereo(mol._rdkit_mol)
        return true
    catch e
        @warn "Error canonicalizing enhanced stereo: $e"
        return false
    end
end

"""
    set_bond_stereo_from_directions!(mol::Molecule) -> Bool

Set bond stereochemistry from directional information.

# Returns
- `Bool`: true if successful, false otherwise
"""
function set_bond_stereo_from_directions!(mol::Molecule)
    !mol.valid && return false
    try
        _set_bond_stereo_from_directions(mol._rdkit_mol)
        return true
    catch e
        @warn "Error setting bond stereo from directions: $e"
        return false
    end
end

"""
    wedge_mol_bonds!(mol::Molecule; wedge_bonds::Bool = true) -> Bool

Add wedge information to molecule bonds for visualization.

# Arguments
- `mol::Molecule`: Input molecule (modified in place)
- `wedge_bonds::Bool`: Whether to add wedge bonds

# Returns
- `Bool`: true if successful, false otherwise
"""
function wedge_mol_bonds!(mol::Molecule; wedge_bonds::Bool = true)
    !mol.valid && return false
    try
        _wedge_mol_bonds(mol._rdkit_mol; wedgeBonds = wedge_bonds)
        return true
    catch e
        @warn "Error wedging mol bonds: $e"
        return false
    end
end

"""
    find_ring_families(mol::Molecule) -> Union{Vector, Missing}

Find ring families in the molecule.

# Returns
- `Union{Vector, Missing}`: Vector of ring families or missing if operation fails
"""
function find_ring_families(mol::Molecule)
    !mol.valid && return missing
    try
        return pyconvert(Vector, _find_ring_families(mol._rdkit_mol))
    catch e
        @warn "Error finding ring families: $e"
        return missing
    end
end


"""
    canonical_rank_atoms_in_fragment(mol::Molecule, atoms_to_use::Vector{Int}) -> Union{Vector{Int}, Missing}

Get canonical ranking of atoms in a molecular fragment.

# Arguments
- `mol::Molecule`: Input molecule
- `atoms_to_use::Vector{Int}`: Atom indices to include in ranking

# Returns
- `Union{Vector{Int}, Missing}`: Canonical ranks or missing if operation fails
"""
function canonical_rank_atoms_in_fragment(mol::Molecule, atoms_to_use::Vector{Int})
    !mol.valid && return missing
    try
        ranks = _canonical_rank_atoms_in_fragment(mol._rdkit_mol, atoms_to_use)
        return pyconvert(Vector{Int}, ranks)
    catch e
        @warn "Error ranking atoms in fragment: $e"
        return missing
    end
end

"""
    find_atom_environment_of_radius_n(mol::Molecule, radius::Int, atom_idx::Int) -> Union{Vector{Int}, Missing}

Find atoms within a specified radius of a given atom.

# Arguments
- `mol::Molecule`: Input molecule
- `radius::Int`: Search radius
- `atom_idx::Int`: Central atom index

# Returns
- `Union{Vector{Int}, Missing}`: Atom indices within radius or missing if operation fails
"""
function find_atom_environment_of_radius_n(mol::Molecule, radius::Int, atom_idx::Int)
    !mol.valid && return missing
    try
        env = _find_atom_environment_of_radius_n(mol._rdkit_mol, radius, atom_idx)
        return pyconvert(Vector{Int}, env)
    catch e
        @warn "Error finding atom environment: $e"
        return missing
    end
end



"""
    get_most_substituted_core_match(mol::Molecule, core::Molecule) -> Union{Vector{Int}, Missing}

Get the most substituted core match in a molecule.

# Arguments
- `mol::Molecule`: Input molecule
- `core::Molecule`: Core pattern to match

# Returns
- `Union{Vector{Int}, Missing}`: Atom indices of the match or missing if no match found
"""
function get_most_substituted_core_match(mol::Molecule, core::Molecule)
    !mol.valid && return missing
    !core.valid && return missing
    try
        match = _get_most_substituted_core_match(mol._rdkit_mol, core._rdkit_mol)
        return pyconvert(Vector{Int}, match)
    catch e
        @warn "Error finding most substituted core match: $e"
        return missing
    end
end



