```@meta
CurrentModule = MoleculeFlow
```

# MoleculeFlow.jl

MoleculeFlow.jl is a comprehensive Julia library for cheminformatics, providing tools for molecular and reaction analysis. It enables common molecular manipulation tasks, such as reading/writing, descriptor calculation, fingerprint generation, similarity analysis, chemical reaction processing, and others.

!!! info
    This library is based on the [rdkit](https://github.com/rdkit/rdkit) and [PythonCall](https://github.com/JuliaPy/PythonCall.jl). Author sees little sense in trying to implement things from scratch, as rdkit has been, is, and very likely will continue being the benchmark cheminformatics library in the future.
    The point of the library is two-fold: 1) Provide a pleasant and user-friendly julian interface that is very easy to extend and/or build upon using a reliable source library. 2) Provide cheminformatics functionality that is not straightforward to implement in in rdkit.
    
    
## Installation

```julia
using Pkg
Pkg.add("MoleculeFlow")
```

## Next Steps

- Check the [Getting Started](getting-started.md) for a basic overview
- See [Examples](examples.md) for some practical use cases
- Browse the [API Reference](api/io.md) for complete function documentation
