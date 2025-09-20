#######################################################
# Additional Molecular Operations
#######################################################

#######################################################
# File I/O Extensions
#######################################################

"""
    mol_from_pdb_block(pdb_block::String) -> Molecule

Create a molecule from a PDB block string.

# Arguments
- `pdb_block::String`: PDB format string

# Returns
- `Molecule`: Parsed molecule or invalid molecule if parsing fails

# Example
```julia
pdb_data = "ATOM      1  C   MOL A   1      20.154  21.875  21.235  1.00 10.00           C"
mol = mol_from_pdb_block(pdb_data)
```
"""
function mol_from_pdb_block(pdb_block::String)
    try
        rdkit_mol = _mol_from_pdb_block(pdb_block)
        if pyconvert(Bool, rdkit_mol !== pybuiltins.None)
            return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = "PDB block")
        else
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = "PDB block")
        end
    catch e
        @warn "Error parsing PDB block: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = "PDB block")
    end
end

"""
    mol_from_pdb_file(filename::String) -> Molecule

Create a molecule from a PDB file.

# Arguments
- `filename::String`: Path to PDB file

# Returns
- `Molecule`: Parsed molecule or invalid molecule if parsing fails

# Example
```julia
mol = mol_from_pdb_file("protein.pdb")
```
"""
function mol_from_pdb_file(filename::String)
    try
        rdkit_mol = _mol_from_pdb_file(filename)
        if pyconvert(Bool, rdkit_mol !== pybuiltins.None)
            return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = filename)
        else
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = filename)
        end
    catch e
        @warn "Error reading PDB file: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = filename)
    end
end

"""
    mol_to_pdb_block(mol::Molecule; confId::Int=-1) -> String

Convert a molecule to PDB block format.

# Arguments
- `mol::Molecule`: Input molecule
- `confId::Int`: Conformer ID to use (-1 for default)

# Returns
- `String`: PDB format string or empty string if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
pdb_block = mol_to_pdb_block(mol)
```
"""
function mol_to_pdb_block(mol::Molecule; confId::Int=-1)
    if !mol.valid
        return ""
    end

    try
        return pyconvert(String, _mol_to_pdb_block(mol._rdkit_mol))
    catch e
        @warn "Error converting to PDB block: $e"
        return ""
    end
end

"""
    mol_to_inchi_key(mol::Molecule) -> String

Generate InChI key for a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `String`: InChI key or empty string if generation fails

# Example
```julia
mol = mol_from_smiles("CCO")
inchi_key = mol_to_inchi_key(mol)
```
"""
function mol_to_inchi_key(mol::Molecule)
    if !mol.valid
        return ""
    end

    try
        return pyconvert(String, _mol_to_inchi_key(mol._rdkit_mol))
    catch e
        @warn "Error generating InChI key: $e"
        return ""
    end
end

"""
    mol_to_molblock(mol::Molecule) -> String

Convert a molecule to MOL block format.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `String`: MOL block string or empty string if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
molblock = mol_to_molblock(mol)
```
"""
function mol_to_molblock(mol::Molecule)
    if !mol.valid
        return ""
    end

    try
        return pyconvert(String, _mol_to_molblock(mol._rdkit_mol))
    catch e
        @warn "Error converting to MOL block: $e"
        return ""
    end
end

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

"""
    fast_find_rings!(mol::Molecule) -> Bool

Perform fast ring finding on a molecule in place.

# Arguments
- `mol::Molecule`: Input molecule (modified in place)

# Returns
- `Bool`: Success status

# Example
```julia
mol = mol_from_smiles("c1ccccc1")
success = fast_find_rings!(mol)
```
"""
function fast_find_rings!(mol::Molecule)
    if !mol.valid
        return false
    end

    try
        _fast_find_rings(mol._rdkit_mol)
        return true
    catch e
        @warn "Error finding rings: $e"
        return false
    end
end

"""
    canonical_atom_ranks(mol::Molecule) -> Vector{Int}

Get canonical ranking of atoms in a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Vector{Int}`: Canonical ranks for each atom

# Example
```julia
mol = mol_from_smiles("CCO")
ranks = canonical_atom_ranks(mol)
```
"""
function canonical_atom_ranks(mol::Molecule)
    if !mol.valid
        return Int[]
    end

    try
        ranks = _canonical_rank_atoms(mol._rdkit_mol)
        return [pyconvert(Int, rank) for rank in ranks]
    catch e
        @warn "Error getting canonical ranks: $e"
        return Int[]
    end
end

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

"""
    mol_fragment_to_smarts(mol::Molecule, atom_indices::Vector{Int}) -> String

Convert a molecular fragment to SMARTS representation.

# Arguments
- `mol::Molecule`: Input molecule
- `atom_indices::Vector{Int}`: Atom indices to include in fragment (0-based)

# Returns
- `String`: SMARTS representation of the fragment

# Example
```julia
mol = mol_from_smiles("CCO")
smarts = mol_fragment_to_smarts(mol, [0, 1])  # First two atoms
```
"""
function mol_fragment_to_smarts(mol::Molecule, atom_indices::Vector{Int})
    if !mol.valid
        return ""
    end

    try
        return pyconvert(String, _mol_fragment_to_smarts(mol._rdkit_mol, atom_indices))
    catch e
        @warn "Error converting fragment to SMARTS: $e"
        return ""
    end
end

# XYZ format support
"""
    mol_from_xyz_file(filename::String) -> Molecule

Read a molecule from an XYZ file.

XYZ is a simple file format for storing molecular coordinates. The format typically
contains the number of atoms, an optional comment line, and atom coordinates.

# Arguments
- `filename`: Path to the XYZ file

# Returns
- `Molecule`: The parsed molecule object
- Returns invalid molecule if parsing fails

# Examples
```julia
mol = mol_from_xyz_file("molecule.xyz")
if mol.valid
    println("Successfully loaded molecule with ", heavy_atom_count(mol), " heavy atoms")
end
```

# Notes
- XYZ files contain only 3D coordinates, no bond information
- RDKit may infer bonds based on interatomic distances
- Not all XYZ files will produce chemically meaningful molecules
"""
function mol_from_xyz_file(filename::String)
    try
        rdkit_mol = _mol_from_xyz_file(filename)
        if pyconvert(Bool, rdkit_mol === pybuiltins.None)
            @warn "Failed to parse XYZ file: $filename"
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = filename, props = Dict{Symbol, Any}())
        end
        return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = filename, props = Dict{Symbol, Any}())
    catch e
        @warn "Error reading XYZ file $filename: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = filename, props = Dict{Symbol, Any}())
    end
end

"""
    mol_from_xyz_block(xyz_block::String) -> Molecule

Parse a molecule from an XYZ format string.

# Arguments
- `xyz_block`: XYZ format string containing atom coordinates

# Returns
- `Molecule`: The parsed molecule object
- Returns invalid molecule if parsing fails

# Examples
```julia
xyz_data = \"\"\"3
Ethanol molecule
C    0.000    0.000    0.000
C    1.520    0.000    0.000
O    2.080    1.100    0.000\"\"\"
mol = mol_from_xyz_block(xyz_data)
```

# Notes
- First line: number of atoms
- Second line: comment (molecule name/description)
- Following lines: element symbol and x, y, z coordinates
"""
function mol_from_xyz_block(xyz_block::String)
    try
        rdkit_mol = _mol_from_xyz_block(xyz_block)
        if pyconvert(Bool, rdkit_mol === pybuiltins.None)
            @warn "Failed to parse XYZ block"
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = "xyz_block", props = Dict{Symbol, Any}())
        end
        return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = "xyz_block", props = Dict{Symbol, Any}())
    catch e
        @warn "Error parsing XYZ block: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = "xyz_block", props = Dict{Symbol, Any}())
    end
end

"""
    mol_to_xyz_block(mol::Molecule) -> String

Convert a molecule to XYZ format string.

Exports the 3D coordinates of the molecule in XYZ format. The molecule must have
3D coordinates assigned for this to work properly.

# Arguments
- `mol`: Molecule object with 3D coordinates

# Returns
- `String`: XYZ format representation
- Empty string if conversion fails or molecule has no 3D coordinates

# Examples
```julia
mol = mol_from_smiles("CCO")
# Generate 3D coordinates first
conformers = generate_3d_conformers(mol, 1)
if !isempty(conformers)
    mol_3d = conformers[1].molecule
    xyz_string = mol_to_xyz_block(mol_3d)
    println(xyz_string)
end
```

# Notes
- Requires 3D coordinates to be present in the molecule
- Only exports heavy atoms (no hydrogens unless explicit)
- Bond information is lost in XYZ format
"""
function mol_to_xyz_block(mol::Molecule)
    !mol.valid && return ""
    try
        return pyconvert(String, _mol_to_xyz_block(mol._rdkit_mol))
    catch e
        @warn "Error converting molecule to XYZ: $e"
        return ""
    end
end

# MOL2 format support
"""
    mol_from_mol2_file(filename::String) -> Molecule

Read a molecule from a MOL2 file.

MOL2 (Sybyl format) is a comprehensive molecular file format that includes
atoms, bonds, and additional molecular information like partial charges.

# Arguments
- `filename`: Path to the MOL2 file

# Returns
- `Molecule`: The parsed molecule object
- Returns invalid molecule if parsing fails

# Examples
```julia
mol = mol_from_mol2_file("ligand.mol2")
if mol.valid
    println("Loaded molecule: ", mol_to_smiles(mol))
    println("Molecular weight: ", molecular_weight(mol))
end
```

# Notes
- MOL2 format preserves bond orders and formal charges
- May contain multiple molecules (only first one is loaded)
- Partial charges are preserved if present in the file
"""
function mol_from_mol2_file(filename::String)
    try
        rdkit_mol = _mol_from_mol2_file(filename)
        if pyconvert(Bool, rdkit_mol === pybuiltins.None)
            @warn "Failed to parse MOL2 file: $filename"
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = filename, props = Dict{Symbol, Any}())
        end
        return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = filename, props = Dict{Symbol, Any}())
    catch e
        @warn "Error reading MOL2 file $filename: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = filename, props = Dict{Symbol, Any}())
    end
end

"""
    mol_from_mol2_block(mol2_block::String) -> Molecule

Parse a molecule from a MOL2 format string.

# Arguments
- `mol2_block`: MOL2 format string

# Returns
- `Molecule`: The parsed molecule object
- Returns invalid molecule if parsing fails

# Examples
```julia
mol2_data = \"\"\"@<TRIPOS>MOLECULE
ethanol
3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 C1 0.0000 0.0000 0.0000 C.3 1 RES1 0.0000
2 C2 1.5200 0.0000 0.0000 C.3 1 RES1 0.0000
3 O1 2.0800 1.1000 0.0000 O.3 1 RES1 0.0000

@<TRIPOS>BOND
1 1 2 1
2 2 3 1
\"\"\"
mol = mol_from_mol2_block(mol2_data)
```

# Notes
- MOL2 format uses TRIPOS record types
- Supports various atom and bond types
- May include substructure and partial charge information
"""
function mol_from_mol2_block(mol2_block::String)
    try
        rdkit_mol = _mol_from_mol2_block(mol2_block)
        if pyconvert(Bool, rdkit_mol === pybuiltins.None)
            @warn "Failed to parse MOL2 block"
            return Molecule(_rdkit_mol = rdkit_mol, valid = false, source = "mol2_block", props = Dict{Symbol, Any}())
        end
        return Molecule(_rdkit_mol = rdkit_mol, valid = true, source = "mol2_block", props = Dict{Symbol, Any}())
    catch e
        @warn "Error parsing MOL2 block: $e"
        return Molecule(_rdkit_mol = pybuiltins.None, valid = false, source = "mol2_block", props = Dict{Symbol, Any}())
    end
end