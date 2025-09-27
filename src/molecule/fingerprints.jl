#######################################################
# Molecular fingerprints
#######################################################

# Morgan fingerprints (ECFP)
"""
    morgan_fingerprint(mol::Molecule; radius::Int=2, nbits::Int=2048) -> Union{Vector{Bool},Missing}

Generate Morgan (ECFP - Extended Connectivity Fingerprint) for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `radius::Int=2`: Radius for the fingerprint (ECFP4 uses radius=2, ECFP6 uses radius=3)
  - `nbits::Int=2048`: Length of the fingerprint bit vector

# Returns

  - `Union{Vector{Bool},Missing}`: Binary fingerprint vector, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = morgan_fingerprint(mol)  # ECFP4 with 2048 bits
fp6 = morgan_fingerprint(mol; radius = 3)  # ECFP6
```

# Notes

  - Morgan fingerprints are circular fingerprints based on atom environments
  - ECFP4 (radius=2) is most commonly used for similarity searches
  - Higher radius captures larger molecular fragments    # Use the new MorganGenerator API to avoid deprecation warning
"""
function morgan_fingerprint(mol::Molecule; radius::Int = 2, nbits::Int = 2048)
    !mol.valid && return missing
    # Use the new MorganGenerator API to avoid deprecation warning
    morgan_gen = _get_morgan_generator(; radius = radius, fpSize = nbits)
    fp = morgan_gen.GetFingerprint(mol._rdkit_mol)
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(nbits - 1)])
end

function morgan_fingerprint(
    mols::Vector{Union{Molecule, Missing}}; radius::Int = 2, nbits::Int = 2048
)
    # Bulk optimization: process all valid molecules at once
    valid_indices = findall(mol -> mol !== missing && mol.valid, mols)
    if isempty(valid_indices)
        return [missing for _ in mols]
    end

    # Extract valid RDKit molecules
    valid_rdkit_mols = [mols[i]._rdkit_mol for i in valid_indices]

    # Create generator once
    morgan_gen = _get_morgan_generator(; radius = radius, fpSize = nbits)

    # Bulk fingerprint generation
    bulk_fps = _get_bulk_morgan_fingerprints(valid_rdkit_mols, morgan_gen)

    # Convert to Julia format
    result = Vector{Union{Vector{Bool}, Missing}}(undef, length(mols))
    for i in eachindex(mols)
        if mols[i] !== missing && mols[i].valid
            # Find position in valid_indices
            valid_pos = findfirst(==(i), valid_indices)
            if valid_pos !== nothing
                fp = bulk_fps[valid_pos - 1]  # Python 0-indexed
                result[i] = pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(j)) for j in 0:(nbits - 1)])
            else
                result[i] = missing
            end
        else
            result[i] = missing
        end
    end

    return result
end

function morgan_fingerprint(mols::Vector{Molecule}; radius::Int = 2, nbits::Int = 2048)
    # Bulk optimization: process all molecules at once
    valid_indices = findall(mol -> mol.valid, mols)
    if isempty(valid_indices)
        return [missing for _ in mols]
    end

    # Extract valid RDKit molecules
    valid_rdkit_mols = [mols[i]._rdkit_mol for i in valid_indices]

    # Create generator once
    morgan_gen = _get_morgan_generator(; radius = radius, fpSize = nbits)

    # Bulk fingerprint generation
    bulk_fps = _get_bulk_morgan_fingerprints(valid_rdkit_mols, morgan_gen)

    # Convert to Julia format
    result = Vector{Union{Vector{Bool}, Missing}}(undef, length(mols))
    for i in eachindex(mols)
        if mols[i].valid
            # Find position in valid_indices
            valid_pos = findfirst(==(i), valid_indices)
            if valid_pos !== nothing
                fp = bulk_fps[valid_pos - 1]  # Python 0-indexed
                result[i] = pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(j)) for j in 0:(nbits - 1)])
            else
                result[i] = missing
            end
        else
            result[i] = missing
        end
    end

    return result
end

# RDK fingerprints
"""
    rdk_fingerprint(mol::Molecule; nbits::Int=2048, min_path::Int=1, max_path::Int=7) -> Union{Vector{Bool}, Missing}

Generate an RDKit fingerprint for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `nbits::Int=2048`: Number of bits in the fingerprint (will be folded/extended if different from default)
  - `min_path::Int=1`: Minimum path length (currently not used in newer RDKit)
  - `max_path::Int=7`: Maximum path length (currently not used in newer RDKit)

# Returns

  - `Vector{Bool}`: Binary fingerprint vector
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = rdk_fingerprint(mol)
length(fp)  # 2048
```

# Notes

  - RDKit fingerprints encode structural features as bit patterns
  - Based on linear and branched subgraphs of molecules
  - Different from Morgan fingerprints in their structural encoding    # GetRDKFingerprint in newer RDKit doesn't support custom parameters like fpSize
"""
function rdk_fingerprint(
    mol::Molecule; nbits::Int = 2048, min_path::Int = 1, max_path::Int = 7
)
    !mol.valid && return missing
    # GetRDKFingerprint in newer RDKit doesn't support custom parameters like fpSize
    # Using the default fingerprint size
    fp = _get_rdk_fingerprint(mol._rdkit_mol)
    fp_size = pyconvert(Int, fp.GetNumBits())

    # If user wants a different size, we'll need to fold or extend the fingerprint
    bits = pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(fp_size - 1)])

    if nbits != fp_size
        if nbits < fp_size
            # Fold the fingerprint by XORing bits
            folded = fill(false, nbits)
            for i in 1:fp_size
                folded[((i - 1) % nbits) + 1] = folded[((i - 1) % nbits) + 1] ⊻ bits[i]
            end
            return folded
        else
            # Extend by padding with zeros
            return vcat(bits, fill(false, nbits - fp_size))
        end
    end

    return bits
end

function rdk_fingerprint(
    mols::Vector{Union{Molecule, Missing}};
    nbits::Int = 2048,
    min_path::Int = 1,
    max_path::Int = 7,
)
    return [
        rdk_fingerprint(mol; nbits = nbits, min_path = min_path, max_path = max_path) for
        mol in mols
    ]
end

function rdk_fingerprint(
    mols::Vector{Molecule}; nbits::Int = 2048, min_path::Int = 1, max_path::Int = 7
)
    return [
        rdk_fingerprint(mol; nbits = nbits, min_path = min_path, max_path = max_path) for
        mol in mols
    ]
end

# MACCS keys
"""
    maccs_fingerprint(mol::Molecule) -> Union{Vector{Bool},Missing}

Generate MACCS (Molecular ACCess System) keys fingerprint for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Vector{Bool},Missing}`: 167-bit MACCS fingerprint, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = maccs_fingerprint(mol)  # 167-bit vector
```

# Notes

  - MACCS keys are a fixed set of 167 predefined structural patterns
  - Each bit represents presence/absence of a specific substructure
  - Widely used for similarity searching and clustering
  - More interpretable than Morgan fingerprints
"""
function maccs_fingerprint(mol::Molecule)
    !mol.valid && return missing
    fp = _maccs_keys(mol._rdkit_mol)
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:166])  # MACCS has 167 bits
end

function maccs_fingerprint(mols::Vector{Union{Molecule, Missing}})
    return [maccs_fingerprint(mol) for mol in mols]
end

function maccs_fingerprint(mols::Vector{Molecule})
    return [maccs_fingerprint(mol) for mol in mols]
end

"""
    atom_pair_fingerprint(mol::Molecule; nbits::Int=2048) -> Union{Vector{Bool}, Missing}

Generate an atom pair fingerprint for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `nbits::Int=2048`: Number of bits in the fingerprint

# Returns

  - `Vector{Bool}`: Binary fingerprint vector
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = atom_pair_fingerprint(mol)
length(fp)  # 2048
```

# Notes

  - Atom pair fingerprints encode pairs of atoms and the distance between them
  - Based on the concept that similar molecules have similar atom pair patterns
  - Useful for scaffold hopping and diverse similarity searching
"""
function atom_pair_fingerprint(mol::Molecule; nbits::Int = 2048)
    !mol.valid && return missing
    fp = _get_hashed_atom_pair_fingerprint_as_bit_vect(mol._rdkit_mol; nBits = nbits)
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(nbits - 1)])
end

function atom_pair_fingerprint(mols::Vector{Union{Molecule, Missing}}; nbits::Int = 2048)
    return [atom_pair_fingerprint(mol; nbits = nbits) for mol in mols]
end

"""
    topological_torsion_fingerprint(mol::Molecule; nbits::Int=2048) -> Union{Vector{Bool}, Missing}

Generate a topological torsion fingerprint for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `nbits::Int=2048`: Number of bits in the fingerprint

# Returns

  - `Vector{Bool}`: Binary fingerprint vector
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = topological_torsion_fingerprint(mol)
length(fp)  # 2048
```

# Notes

  - Topological torsion fingerprints encode four-atom paths in molecules
  - Based on the torsion angles and atom types in these paths
  - Useful for 3D pharmacophore-like similarity searching
"""
function topological_torsion_fingerprint(mol::Molecule; nbits::Int = 2048)
    !mol.valid && return missing
    fp = _get_hashed_topological_torsion_fingerprint_as_bit_vect(
        mol._rdkit_mol; nBits = nbits
    )
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(nbits - 1)])
end

function topological_torsion_fingerprint(
    mols::Vector{Union{Molecule, Missing}}; nbits::Int = 2048
)
    return [topological_torsion_fingerprint(mol; nbits = nbits) for mol in mols]
end

"""
    fcfp_fingerprint(mol::Molecule; radius::Int=2, nbits::Int=2048) -> Union{Vector{Bool}, Missing}

Generate a Functional-Class FingerPrint (FCFP) for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `radius::Int=2`: Radius of the circular fingerprint
  - `nbits::Int=2048`: Number of bits in the fingerprint

# Returns

  - `Vector{Bool}`: Binary fingerprint vector
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = fcfp_fingerprint(mol; radius = 3)
length(fp)  # 2048
```

# Notes

  - FCFP fingerprints are similar to Morgan fingerprints but use functional class atom invariants
  - Groups atoms by pharmacophoric features rather than exact atom types
  - Better for finding molecules with similar biological activity    # Use the new MorganGenerator API with feature atom invariants generator to avoid deprecation warning
"""
function fcfp_fingerprint(mol::Molecule; radius::Int = 2, nbits::Int = 2048)
    !mol.valid && return missing
    # Use the new MorganGenerator API with feature atom invariants generator to avoid deprecation warning
    feature_inv_gen = _get_morgan_feature_atom_inv_gen()
    morgan_gen = _get_morgan_generator(;
        radius = radius, fpSize = nbits, atomInvariantsGenerator = feature_inv_gen
    )
    fp = morgan_gen.GetFingerprint(mol._rdkit_mol)
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(nbits - 1)])
end

function fcfp_fingerprint(
    mols::Vector{Union{Molecule, Missing}}; radius::Int = 2, nbits::Int = 2048
)
    return [fcfp_fingerprint(mol; radius = radius, nbits = nbits) for mol in mols]
end

"""
    pattern_fingerprint(mol::Molecule; nbits::Int=2048) -> Union{Vector{Bool}, Missing}

Generate a pattern fingerprint for a molecule.

# Arguments

  - `mol::Molecule`: Input molecule
  - `nbits::Int=2048`: Number of bits in the fingerprint

# Returns

  - `Vector{Bool}`: Binary fingerprint vector
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
fp = pattern_fingerprint(mol)
length(fp)  # 2048
```

# Notes

  - Pattern fingerprints encode molecular substructures as bit patterns
  - Based on predefined structural patterns and motifs
  - Useful for substructure-based similarity searching
"""
function pattern_fingerprint(mol::Molecule; nbits::Int = 2048)
    !mol.valid && return missing
    fp = _pattern_fingerprint(mol._rdkit_mol; fpSize = nbits)
    return pyconvert(Vector{Bool}, [pyconvert(Bool, fp.GetBit(i)) for i in 0:(nbits - 1)])
end

function pattern_fingerprint(mols::Vector{Union{Molecule, Missing}}; nbits::Int = 2048)
    return [pattern_fingerprint(mol; nbits = nbits) for mol in mols]
end
