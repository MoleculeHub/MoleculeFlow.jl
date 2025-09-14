#######################################################
# Writing Mols
#######################################################
_mol_to_smiles(mol::Py) = @pyconst(pyimport("rdkit.Chem").MolToSmiles)(mol)
_mol_to_inchi(mol::Py) = @pyconst(pyimport("rdkit.Chem").MolToInchi)(mol)
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
