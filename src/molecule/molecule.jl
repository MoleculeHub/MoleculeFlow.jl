#######################################################
# Molecule structs and properties
#######################################################
@kwdef struct Molecule
    _rdkit_mol::Py
    valid::Bool
    source::String
    props::Dict{Symbol, Any} = Dict()
end

function Base.show(io::IO, mol::Molecule)
    if mol.valid
        print(io, "Molecule($(mol.source))")
    else
        print(io, "BadMolecule($(mol.source))")
    end
end

function Base.setproperty!(mol::Molecule, name::Symbol, value)
    if !(name in (:_rdkit_mol, :valid, :source, :props))
        mol.props[name] = value
    else
        throw(ArgumentError("Cannot set this property directly."))
    end
end

function Base.getproperty(mol::Molecule, name::Symbol)
    if name in fieldnames(Molecule)
        return getfield(mol, name)
    else
        return getfield(mol, :props)[name]
    end
end

#######################################################
# Hydrogen manipulation
#######################################################

"""
    add_hs(mol::Molecule) -> Molecule

Add explicit hydrogens to a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: New molecule with explicit hydrogens added

# Example
```julia
mol = mol_from_smiles("CCO")
mol_with_hs = add_hs(mol)
```
"""
function add_hs(mol::Molecule)
    if !mol.valid
        return mol
    end

    try
        new_rdkit_mol = _add_hs(mol._rdkit_mol)
        return Molecule(
            _rdkit_mol = new_rdkit_mol,
            valid = true,
            source = mol.source,
            props = copy(mol.props)
        )
    catch e
        @warn "Error adding hydrogens: $e"
        return Molecule(
            _rdkit_mol = mol._rdkit_mol,
            valid = false,
            source = mol.source,
            props = copy(mol.props)
        )
    end
end

"""
    remove_hs(mol::Molecule) -> Molecule

Remove explicit hydrogens from a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: New molecule with explicit hydrogens removed

# Example
```julia
mol_with_hs = add_hs(mol_from_smiles("CCO"))
mol_no_hs = remove_hs(mol_with_hs)
```
"""
function remove_hs(mol::Molecule)
    if !mol.valid
        return mol
    end

    try
        new_rdkit_mol = _remove_hs(mol._rdkit_mol)
        return Molecule(
            _rdkit_mol = new_rdkit_mol,
            valid = true,
            source = mol.source,
            props = copy(mol.props)
        )
    catch e
        @warn "Error removing hydrogens: $e"
        return Molecule(
            _rdkit_mol = mol._rdkit_mol,
            valid = false,
            source = mol.source,
            props = copy(mol.props)
        )
    end
end
