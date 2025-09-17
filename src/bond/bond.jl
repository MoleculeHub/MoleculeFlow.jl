#######################################################
# Bond structs and properties
#######################################################
"""
    Bond

A structure representing a chemical bond between two atoms in a molecule.

# Fields
- `_rdkit_bond::Py`: The underlying RDKit bond object
- `props::Dict{Symbol, Any}`: Dictionary for storing additional bond properties

# Examples
```julia
mol = mol_from_smiles("CCO")
atoms = get_atoms(mol)
bonds = get_bonds_from_atom(mol, 1)  # Bonds from first atom
bond = bonds[1]

# Access bond properties
bond_type = get_bond_type(bond)
begin_idx = get_begin_atom_idx(bond)
end_idx = get_end_atom_idx(bond)
```
"""
@kwdef struct Bond
    _rdkit_bond::Py
    props::Dict{Symbol, Any} = Dict()
end

function Base.show(io::IO, bond::Bond)
    msg = join(["$k => $(v)" for (k, v) in bond.props], ", ")
    print(io, "Bond($msg)")
end

function Base.setproperty!(bond::Bond, name::Symbol, value)
    if !(name in (:_rdkit_bond, :props))
        bond.props[name] = value
    else
        throw(ArgumentError("Cannot set this property directly."))
    end
end

function Base.getproperty(bond::Bond, name::Symbol)
    if name in fieldnames(Bond)
        return getfield(bond, name)
    else
        return getfield(bond, :props)[name]
    end
end

"""
    get_bond_type(bond::Bond) -> String

Get the bond type as a string.

# Arguments
- `bond`: A Bond object

# Returns
- `String`: The bond type ("SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", etc.)

# Example
```julia
mol = mol_from_smiles("C=C")
bonds = get_bonds(mol)
bond_type = get_bond_type(bonds[1])  # Returns "DOUBLE"
```
"""
function get_bond_type(bond::Bond)
    return pyconvert(String, bond._rdkit_bond.GetBondType().__str__())
end

"""
    get_begin_atom_idx(bond::Bond) -> Int

Get the index of the first atom in the bond (1-based indexing).

# Arguments
- `bond`: A Bond object

# Returns
- `Int`: The 1-based index of the beginning atom

# Example
```julia
mol = mol_from_smiles("CCO")
bonds = get_bonds(mol)
begin_idx = get_begin_atom_idx(bonds[1])  # Returns 1
```
"""
function get_begin_atom_idx(bond::Bond)
    return pyconvert(Int, bond._rdkit_bond.GetBeginAtomIdx()) + 1  # Convert to 1-based indexing
end

"""
    get_end_atom_idx(bond::Bond) -> Int

Get the index of the second atom in the bond (1-based indexing).

# Arguments
- `bond`: A Bond object

# Returns
- `Int`: The 1-based index of the ending atom

# Example
```julia
mol = mol_from_smiles("CCO")
bonds = get_bonds(mol)
end_idx = get_end_atom_idx(bonds[1])  # Returns 2
```
"""
function get_end_atom_idx(bond::Bond)
    return pyconvert(Int, bond._rdkit_bond.GetEndAtomIdx()) + 1  # Convert to 1-based indexing
end

"""
    is_aromatic(bond::Bond) -> Bool

Check if a bond is aromatic.

# Arguments
- `bond`: A Bond object

# Returns
- `Bool`: true if the bond is aromatic, false otherwise

# Example
```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
bonds = get_bonds(mol)
is_aromatic(bonds[1])  # Returns true for aromatic C-C bond
```
"""
function is_aromatic(bond::Bond)
    return pyconvert(Bool, bond._rdkit_bond.GetIsAromatic())
end

"""
    is_conjugated(bond::Bond) -> Bool

Check if a bond is conjugated (part of a conjugated system).

# Arguments
- `bond`: A Bond object

# Returns
- `Bool`: true if the bond is conjugated, false otherwise

# Example
```julia
mol = mol_from_smiles("C=C-C=C")  # Conjugated diene
bonds = get_bonds(mol)
is_conjugated(bonds[2])  # Returns true for the central C-C bond
```
"""
function is_conjugated(bond::Bond)
    return pyconvert(Bool, bond._rdkit_bond.GetIsConjugated())
end

"""
    is_in_ring(bond::Bond) -> Bool

Check if a bond is part of a ring.

# Arguments
- `bond`: A Bond object

# Returns
- `Bool`: true if the bond is in a ring, false otherwise

# Example
```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
bonds = get_bonds(mol)
is_in_ring(bonds[1])  # Returns true
```
"""
function is_in_ring(bond::Bond)
    return pyconvert(Bool, bond._rdkit_bond.IsInRing())
end
