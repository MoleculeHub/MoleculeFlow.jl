# For Developers

## Architecture Overview

MoleculeFlow.jl is built as a Julia wrapper around 
[RDKit](https://github.com/rdkit/rdkit) cheminformatics library. 
Rather than reimplementing complex
cheminformatics algorithms from scratch, we leverage RDKit's mature and
battle-tested implementations while providing a clean, Julia-native
interface, which also allows to build additional functionality on top.

!!! info 
    Knowledge of Python RDKit is assumed.
    This background knowledge is essential because MoleculeFlow's functions are direct mappings to RDKit functionality.

## Development Pattern

### 1. Low-Level Caching Layer

All RDKit functionality is accessed through caching functions in `src/rdkit.jl`. This file contains functions that:
- Import and cache Python RDKit modules and functions using `@pyconst`
- Provide direct access to RDKit functions
- Handle the Python-Julia interface via PythonCall.jl

!!! info
    The `@pyconst` macro from PythonCall.jl is used to cache Python function calls for better performance. This ensures that the Python import and function lookup only happens once, and subsequent calls reuse the cached function object.

Example from `rdkit.jl`:
```julia
# RDKit function caching with @pyconst
_mol_from_smiles(smiles::String) = @pyconst(pyimport("rdkit.Chem").MolFromSmiles)(smiles)
_mol_to_smiles(mol::Py) = @pyconst(pyimport("rdkit.Chem").MolToSmiles)(mol)

function _mol_wt(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MolWt)(mol)
end

function _heavy_atom_count(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").HeavyAtomCount)(mol)
end
```

### 2. Higher-Level Julia Wrappers

The main API functions are implemented as Julia wrappers around the cached RDKit functions. 

Example wrapper pattern:
```julia
function molecular_weight(mol::Molecule)
    !mol.valid && return missing

    try
        result = _mol_wt(mol._rdkit_mol)
        return pyconvert(Float64, result)
    catch
        return missing
    end
end

function get_atom(mol::Molecule, idx::Int)
    !mol.valid && throw(ArgumentError("Invalid molecule"))

    # Convert indices (Julia 1-based → Python 0-based)
    python_idx = idx - 1

    rdkit_atom = mol._rdkit_mol.GetAtomWithIdx(python_idx)
    return Atom(rdkit_atom, mol)
end
```

## Critical Index Conversion

!!! danger
    Index Conversion is Critical: Julia uses 1-based indexing while Python uses 0-based indexing. 

### Examples of Index Handling
```julia

function get_atom(mol::Molecule, idx::Int)
    rdkit_mol = mol._rdkit_mol
    rdkit_idx = idx - 1  
    return get_chem().GetAtomWithIdx(rdkit_mol, rdkit_idx)
end
```

## Error Handling Patterns

### Molecule Validation
Always check molecule validity:
```julia
function some_function(mol::Molecule)
    !mol.valid && throw(ArgumentError("Invalid molecule"))
    # ... proceed with function
end
```

### Graceful Degradation
For operations that may fail, return appropriate fallback values:
```julia
function calculate_descriptor(mol::Molecule)
    !mol.valid && return missing

    try
        result = rdkit_calculation(mol._rdkit_mol)
        return convert_result(result)
    catch
        return missing  # Or appropriate fallback
    end
end
```

## Contributing New Functions

When adding new functionality:

1. **Add RDKit module caching** to `src/rdkit.jl` if needed
2. **Implement the wrapper** in the appropriate module file
3. **Export the function** in `src/MoleculeFlow.jl`
4. **Add comprehensive tests** following existing patterns
5. **Document with examples** in the appropriate API docs file
6. **Handle index conversion** carefully
7. **Follow existing error handling patterns**

## Best Practices

1. **Follow existing naming conventions** (snake_case for functions)
2. **Use keyword arguments** for optional parameters
3. **Provide sensible defaults** where possible
4. **Handle missing/invalid inputs gracefully**
5. **Add type annotations** where helpful
