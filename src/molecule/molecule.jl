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
