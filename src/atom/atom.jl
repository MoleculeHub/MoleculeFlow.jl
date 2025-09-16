"""
    Atom

A structure representing an individual atom within a molecule.

# Fields
- `_rdkit_atom::Py`: The underlying RDKit atom object
- `props::Dict{Symbol, Any}`: Dictionary for storing additional atom properties

# Examples
```julia
mol = mol_from_smiles("CCO")
atoms = get_atoms(mol)
atom = atoms[1]  # First atom (carbon)

# Access atom properties
symbol = get_symbol(atom)
atomic_num = get_atomic_number(atom)
```
"""
@kwdef struct Atom
    _rdkit_atom::Py
    props::Dict{Symbol, Any} = Dict()
end

function Base.show(io::IO, atom::Atom)
    msg = join(["$k => $(v)" for (k, v) in atom.props], ", ")
    print(io, "Atom($msg)")
end

function Base.setproperty!(atom::Atom, name::Symbol, value)
    if !(name in (:_rdkit_atom, :props))
        atom.props[name] = value
    else
        throw(ArgumentError("Cannot set this property directly."))
    end
end

function Base.getproperty(atom::Atom, name::Symbol)
    if name in fieldnames(Atom)
        return getfield(atom, name)
    else
        return getfield(atom, :props)[name]
    end
end
