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
