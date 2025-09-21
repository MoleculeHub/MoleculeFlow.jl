#######################################################
# Writing Mols
#######################################################

"""
    mol_to_smiles(mol::Union{Molecule,Missing}) -> Union{String,Missing}

Convert a Molecule object to a SMILES string.

# Arguments

  - `mol::Union{Molecule,Missing}`: A Molecule object or missing

# Returns

  - `Union{String,Missing}`: SMILES string, or missing if the molecule is invalid or missing

# Examples

```julia
mol = mol_from_smiles("CCO")
smiles = mol_to_smiles(mol)  # "CCO"

# Handles missing values gracefully
smiles = mol_to_smiles(missing)  # missing
```
"""
function mol_to_smiles(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    return pyconvert(String, _mol_to_smiles(mol._rdkit_mol))
end

function mol_to_smiles(mol_list::Vector{Union{Molecule, Missing}})
    results = Vector{Union{String, Missing}}(undef, length(mol_list))

    @inbounds for i in eachindex(mol_list)
        smiles = mol_to_smiles(mol_list[i])
        results[i] = smiles
    end

    return pyconvert(Vector{Union{String, Missing}}, results)
end

function mol_to_smiles(mol_list::Vector{Molecule})
    results = Vector{Union{String, Missing}}(undef, length(mol_list))

    @inbounds for i in eachindex(mol_list)
        smiles = mol_to_smiles(mol_list[i])
        results[i] = smiles
    end

    return pyconvert(Vector{Union{String, Missing}}, results)
end

"""
    mol_to_inchi(mol::Union{Molecule,Missing}) -> Union{String,Missing}

Convert a Molecule object to an InChI (International Chemical Identifier) string.

# Arguments

  - `mol::Union{Molecule,Missing}`: A Molecule object or missing

# Returns

  - `Union{String,Missing}`: InChI string, or missing if the molecule is invalid or missing

# Examples

```julia
mol = mol_from_smiles("CCO")
inchi = mol_to_inchi(mol)  # "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
```

# Notes

  - InChI provides a unique identifier for chemical structures
  - More verbose than SMILES but provides standardized representation
"""
function mol_to_inchi(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    return pyconvert(String, _mol_to_inchi(mol._rdkit_mol))
end

function mol_to_inchi(mol_list::Vector{Union{Molecule, Missing}})
    results = Vector{Union{String, Missing}}(undef, length(mol_list))

    @inbounds for i in eachindex(mol_list)
        inchi = mol_to_inchi(mol_list[i])
        results[i] = inchi
    end

    return pyconvert(Vector{Union{String, Missing}}, results)
end

"""
    mol_to_inchi_key(mol::Molecule) -> Union{String, Missing}

Convert a molecule to InChI Key format.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{String, Missing}`: InChI Key or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
key = mol_to_inchi_key(mol)  # "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
```
"""
function mol_to_inchi_key(mol::Molecule)
    !mol.valid && return ""
    try
        return pyconvert(String, _mol_to_inchi_key(mol._rdkit_mol))
    catch e
        @warn "Error converting to InChI key: $e"
        return ""
    end
end

"""
    mol_to_inchi_and_aux_info(mol::Molecule) -> Union{Tuple{String, String}, Missing}

Convert a molecule to InChI string with auxiliary information.

# Returns
- `Union{Tuple{String, String}, Missing}`: (InChI, AuxInfo) or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
inchi, aux = mol_to_inchi_and_aux_info(mol)
```
"""
function mol_to_inchi_and_aux_info(mol::Molecule)
    !mol.valid && return missing
    try
        result = _mol_to_inchi_and_aux_info(mol._rdkit_mol)
        inchi = pyconvert(String, result[0])
        aux_info = pyconvert(String, result[1])
        return (inchi, aux_info)
    catch e
        @warn "Error converting to InChI with aux info: $e"
        return missing
    end
end

"""
    mol_to_molblock(mol::Molecule) -> Union{String, Missing}

Convert a molecule to MOL block format.

# Returns
- `Union{String, Missing}`: MOL block or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
molblock = mol_to_molblock(mol)
```
"""
function mol_to_molblock(mol::Molecule)
    !mol.valid && return ""
    try
        return pyconvert(String, _mol_to_molblock(mol._rdkit_mol))
    catch e
        @warn "Error converting to MOL block: $e"
        return ""
    end
end

"""
    mol_to_v3k_molblock(mol::Molecule) -> Union{String, Missing}

Convert a molecule to V3000 MOL block format.

# Returns
- `Union{String, Missing}`: V3000 MOL block or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
v3k_block = mol_to_v3k_molblock(mol)
```
"""
function mol_to_v3k_molblock(mol::Molecule)
    !mol.valid && return missing
    try
        return pyconvert(String, _mol_to_v3k_molblock(mol._rdkit_mol))
    catch e
        @warn "Error converting to V3K MOL block: $e"
        return missing
    end
end

"""
    mol_to_pdb_block(mol::Molecule; conf_id::Int=-1) -> Union{String, Missing}

Convert a molecule to PDB block format.

# Arguments
- `mol::Molecule`: Input molecule
- `conf_id::Int`: Conformer ID to use (-1 for default)

# Returns
- `Union{String, Missing}`: PDB block or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
pdb_block = mol_to_pdb_block(mol)
```
"""
function mol_to_pdb_block(mol::Molecule; conf_id::Int=-1)
    !mol.valid && return ""
    try
        return pyconvert(String, _mol_to_pdb_block(mol._rdkit_mol))
    catch e
        @warn "Error converting to PDB block: $e"
        return ""
    end
end

"""
    mol_to_xyz_block(mol::Molecule) -> Union{String, Missing}

Convert a molecule to XYZ block format.

# Returns
- `Union{String, Missing}`: XYZ block or missing if conversion fails

# Example
```julia
mol = mol_from_smiles("CCO")
xyz_block = mol_to_xyz_block(mol)
```
"""
function mol_to_xyz_block(mol::Molecule)
    !mol.valid && return ""
    try
        return pyconvert(String, _mol_to_xyz_block(mol._rdkit_mol))
    catch e
        @warn "Error converting to XYZ block: $e"
        return ""
    end
end

"""
    mol_to_pdb_file(mol::Molecule, filename::String; conf_id::Int=-1) -> Bool

Write a molecule to a PDB file.

# Arguments
- `mol::Molecule`: Input molecule
- `filename::String`: Output file path
- `conf_id::Int`: Conformer ID to use (-1 for default)

# Returns
- `Bool`: true if successful, false otherwise

# Example
```julia
mol = mol_from_smiles("CCO")
success = mol_to_pdb_file(mol, "ethanol.pdb")
```
"""
function mol_to_pdb_file(mol::Molecule, filename::String; conf_id::Int=-1)
    !mol.valid && return false
    try
        _mol_to_pdb_file(mol._rdkit_mol, filename)
        return true
    catch e
        @warn "Error writing PDB file: $e"
        return false
    end
end

"""
    mol_to_xyz_file(mol::Molecule, filename::String) -> Bool

Write a molecule to an XYZ file.

# Arguments
- `mol::Molecule`: Input molecule
- `filename::String`: Output file path

# Returns
- `Bool`: true if successful, false otherwise

# Example
```julia
mol = mol_from_smiles("CCO")
success = mol_to_xyz_file(mol, "ethanol.xyz")
```
"""
function mol_to_xyz_file(mol::Molecule, filename::String)
    !mol.valid && return false
    try
        _mol_to_xyz_file(mol._rdkit_mol, filename)
        return true
    catch e
        @warn "Error writing XYZ file: $e"
        return false
    end
end
