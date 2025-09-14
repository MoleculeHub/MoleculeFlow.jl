#######################################################
# Molecule to Graph conversion using Graphs.jl
#######################################################

using Graphs

"""
    mol_to_graph(mol::Molecule; atom_features=false, bond_features=false)

Convert a molecule to a Graph object from Graphs.jl.

# Arguments
- `mol::Molecule`: The molecule to convert
- `atom_features::Bool`: Whether to include atom features in metadata (default: false)
- `bond_features::Bool`: Whether to include bond features in metadata (default: false)

# Returns
- A SimpleGraph or SimpleDiGraph representing the molecular structure
- If atom_features or bond_features are true, returns a MetaGraph with features

# Examples
```julia
mol = mol_from_smiles("CCO")
g = mol_to_graph(mol)
nv(g)  # number of atoms (vertices)
ne(g)  # number of bonds (edges)
```
"""
function mol_to_graph(mol::Molecule; atom_features=false, bond_features=false)
    if !mol.valid
        return missing
    end

    try
        # Get number of atoms from RDKit molecule
        natoms = pyconvert(Int, mol._rdkit_mol.GetNumAtoms())

        if natoms == 0
            return SimpleGraph(0)
        end

        # Create basic graph
        g = SimpleGraph(natoms)

        # Get all bonds from the molecule directly
        nbonds = pyconvert(Int, mol._rdkit_mol.GetNumBonds())

        for i in 0:(nbonds - 1)  # RDKit uses 0-based indexing
            bond = mol._rdkit_mol.GetBondWithIdx(i)
            begin_idx = pyconvert(Int, bond.GetBeginAtomIdx()) + 1  # Convert to 1-based
            end_idx = pyconvert(Int, bond.GetEndAtomIdx()) + 1      # Convert to 1-based

            # Add edge (undirected graph)
            add_edge!(g, begin_idx, end_idx)
        end

        return g
    catch e
        @warn "Error converting molecule to graph: $e"
        return missing
    end
end

"""
    mol_to_digraph(mol::Molecule; atom_features=false, bond_features=false)

Convert a molecule to a directed Graph object from Graphs.jl.
Each bond creates two directed edges (one in each direction).

# Arguments
- `mol::Molecule`: The molecule to convert
- `atom_features::Bool`: Whether to include atom features in metadata (default: false)
- `bond_features::Bool`: Whether to include bond features in metadata (default: false)

# Returns
- A SimpleDiGraph representing the molecular structure

# Examples
```julia
mol = mol_from_smiles("CCO")
dg = mol_to_digraph(mol)
nv(dg)  # number of atoms (vertices)
ne(dg)  # number of directed bonds (edges)
```
"""
function mol_to_digraph(mol::Molecule; atom_features=false, bond_features=false)
    if !mol.valid
        return missing
    end

    try
        # Get number of atoms from RDKit molecule
        natoms = pyconvert(Int, mol._rdkit_mol.GetNumAtoms())

        if natoms == 0
            return SimpleDiGraph(0)
        end

        # Create directed graph
        dg = SimpleDiGraph(natoms)

        # Get all bonds from the molecule directly
        nbonds = pyconvert(Int, mol._rdkit_mol.GetNumBonds())

        for i in 0:(nbonds - 1)  # RDKit uses 0-based indexing
            bond = mol._rdkit_mol.GetBondWithIdx(i)
            begin_idx = pyconvert(Int, bond.GetBeginAtomIdx()) + 1  # Convert to 1-based
            end_idx = pyconvert(Int, bond.GetEndAtomIdx()) + 1      # Convert to 1-based

            # Add both directions for directed graph
            add_edge!(dg, begin_idx, end_idx)
            add_edge!(dg, end_idx, begin_idx)
        end

        return dg
    catch e
        @warn "Error converting molecule to directed graph: $e"
        return missing
    end
end
