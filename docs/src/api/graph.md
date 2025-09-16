```@meta
CurrentModule = MoleculeFlow
```

# Graph Operations

Functions for converting molecules to graph representations and analyzing molecular graphs.

## Graph Conversion

```@docs
mol_to_graph
mol_to_digraph
```

## Graph Analysis

These functions are re-exported from [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) for convenience when working with molecular graphs:

- `nv(g)` - Number of vertices (atoms) in the graph
- `ne(g)` - Number of edges (bonds) in the graph
- `density(g)` - Edge density of the graph
- `is_connected(g)` - Check if the graph is connected
- `diameter(g)` - Diameter of the graph
- `radius(g)` - Radius of the graph

See the [Graphs.jl documentation](https://juliagraphs.org/Graphs.jl/stable/) for complete function documentation.

## Graph Matrices

These functions are also re-exported from Graphs.jl:

- `adjacency_matrix(g)` - Adjacency matrix representation
- `laplacian_matrix(g)` - Laplacian matrix representation

## Examples

```julia
using MoleculeFlow

mol = mol_from_smiles("c1ccccc1")  # Benzene
g = mol_to_graph(mol)

# Basic graph properties
println("Number of atoms: ", nv(g))     # 6
println("Number of bonds: ", ne(g))     # 6
println("Is connected: ", is_connected(g))  # true
println("Density: ", density(g))        # 0.33

# Graph matrices
adj_mat = adjacency_matrix(g)
lap_mat = laplacian_matrix(g)
```