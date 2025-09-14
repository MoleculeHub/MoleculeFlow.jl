#######################################################
# Molecular similarity calculations
#######################################################

"""
    tanimoto_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    tanimoto_similarity(mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...)

Calculate the Tanimoto similarity coefficient between two molecules or fingerprints.

The Tanimoto coefficient is defined as the size of the intersection divided by the size of the union of two sets.
For fingerprints, this translates to the number of bits set in both fingerprints divided by the number of bits set in either fingerprint.

# Arguments
- `fp1`, `fp2`: Boolean vectors representing molecular fingerprints
- `mol1`, `mol2`: Molecule objects to compare
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`, `:atom_pair`, `:topological_torsion`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Float64`: Similarity score between 0.0 (no similarity) and 1.0 (identical)
- `missing`: If either molecule is invalid

# Examples
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
similarity = tanimoto_similarity(mol1, mol2)  # Using Morgan fingerprints by default
similarity_rdk = tanimoto_similarity(mol1, mol2, fingerprint_type=:rdk)
```
"""
function tanimoto_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must have the same length"))
    end

    intersection = sum(fp1 .& fp2)
    union = sum(fp1 .| fp2)

    return union == 0 ? 0.0 : intersection / union
end

function tanimoto_similarity(
    mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...
)
    if !mol1.valid || !mol2.valid
        return missing
    end

    if fingerprint_type == :morgan
        fp1 = morgan_fingerprint(mol1; kwargs...)
        fp2 = morgan_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :rdk
        fp1 = rdk_fingerprint(mol1; kwargs...)
        fp2 = rdk_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :maccs
        fp1 = maccs_fingerprint(mol1)
        fp2 = maccs_fingerprint(mol2)
    elseif fingerprint_type == :atom_pair
        fp1 = atom_pair_fingerprint(mol1; kwargs...)
        fp2 = atom_pair_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :topological_torsion
        fp1 = topological_torsion_fingerprint(mol1; kwargs...)
        fp2 = topological_torsion_fingerprint(mol2; kwargs...)
    else
        throw(ArgumentError("Unsupported fingerprint type: $fingerprint_type"))
    end

    return tanimoto_similarity(fp1, fp2)
end

"""
    dice_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    dice_similarity(mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...)

Calculate the Dice similarity coefficient (also known as Sørensen-Dice coefficient) between two molecules or fingerprints.

The Dice coefficient is defined as twice the size of the intersection divided by the sum of the sizes of the two sets.
For fingerprints, this translates to 2 × (number of bits set in both) / (total bits set in both fingerprints).

# Arguments
- `fp1`, `fp2`: Boolean vectors representing molecular fingerprints
- `mol1`, `mol2`: Molecule objects to compare
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`, `:atom_pair`, `:topological_torsion`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Float64`: Similarity score between 0.0 (no similarity) and 1.0 (identical)
- `missing`: If either molecule is invalid

# Examples
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
similarity = dice_similarity(mol1, mol2)
```
"""
function dice_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must have the same length"))
    end

    intersection = sum(fp1 .& fp2)
    sum_bits = sum(fp1) + sum(fp2)

    return sum_bits == 0 ? 0.0 : 2 * intersection / sum_bits
end

function dice_similarity(
    mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...
)
    if !mol1.valid || !mol2.valid
        return missing
    end

    if fingerprint_type == :morgan
        fp1 = morgan_fingerprint(mol1; kwargs...)
        fp2 = morgan_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :rdk
        fp1 = rdk_fingerprint(mol1; kwargs...)
        fp2 = rdk_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :maccs
        fp1 = maccs_fingerprint(mol1)
        fp2 = maccs_fingerprint(mol2)
    else
        throw(ArgumentError("Unsupported fingerprint type: $fingerprint_type"))
    end

    return dice_similarity(fp1, fp2)
end

"""
    cosine_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    cosine_similarity(mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...)

Calculate the cosine similarity between two molecules or fingerprints.

The cosine similarity measures the cosine of the angle between two vectors in a multi-dimensional space.
For binary fingerprints, this is computed as the dot product divided by the product of the magnitudes.

# Arguments
- `fp1`, `fp2`: Boolean vectors representing molecular fingerprints
- `mol1`, `mol2`: Molecule objects to compare
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Float64`: Similarity score between 0.0 (orthogonal) and 1.0 (identical direction)
- `missing`: If either molecule is invalid

# Examples
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
similarity = cosine_similarity(mol1, mol2)
```
"""
function cosine_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must have the same length"))
    end

    dot_product = sum(fp1 .& fp2)
    norm1 = sqrt(sum(fp1))
    norm2 = sqrt(sum(fp2))

    return (norm1 == 0 || norm2 == 0) ? 0.0 : dot_product / (norm1 * norm2)
end

function cosine_similarity(
    mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...
)
    if !mol1.valid || !mol2.valid
        return missing
    end

    if fingerprint_type == :morgan
        fp1 = morgan_fingerprint(mol1; kwargs...)
        fp2 = morgan_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :rdk
        fp1 = rdk_fingerprint(mol1; kwargs...)
        fp2 = rdk_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :maccs
        fp1 = maccs_fingerprint(mol1)
        fp2 = maccs_fingerprint(mol2)
    else
        throw(ArgumentError("Unsupported fingerprint type: $fingerprint_type"))
    end

    return cosine_similarity(fp1, fp2)
end

"""
    sokal_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    sokal_similarity(mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...)

Calculate the Sokal similarity coefficient between two molecules or fingerprints.

The Sokal similarity is defined as intersection / (2 × union - intersection).
This metric gives less weight to shared features compared to Tanimoto similarity.

# Arguments
- `fp1`, `fp2`: Boolean vectors representing molecular fingerprints
- `mol1`, `mol2`: Molecule objects to compare
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Float64`: Similarity score between 0.0 (no similarity) and 1.0 (identical)
- `missing`: If either molecule is invalid

# Examples
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
similarity = sokal_similarity(mol1, mol2)
```
"""
function sokal_similarity(fp1::Vector{Bool}, fp2::Vector{Bool})
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must have the same length"))
    end

    intersection = sum(fp1 .& fp2)
    union = sum(fp1 .| fp2)

    return union == 0 ? 0.0 : intersection / (2 * union - intersection)
end

function sokal_similarity(
    mol1::Molecule, mol2::Molecule; fingerprint_type=:morgan, kwargs...
)
    if !mol1.valid || !mol2.valid
        return missing
    end

    if fingerprint_type == :morgan
        fp1 = morgan_fingerprint(mol1; kwargs...)
        fp2 = morgan_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :rdk
        fp1 = rdk_fingerprint(mol1; kwargs...)
        fp2 = rdk_fingerprint(mol2; kwargs...)
    elseif fingerprint_type == :maccs
        fp1 = maccs_fingerprint(mol1)
        fp2 = maccs_fingerprint(mol2)
    else
        throw(ArgumentError("Unsupported fingerprint type: $fingerprint_type"))
    end

    return sokal_similarity(fp1, fp2)
end

"""
    bulk_similarity(query_mol::Molecule, target_mols::Vector{Molecule};
                   similarity_function=tanimoto_similarity, fingerprint_type=:morgan, kwargs...)

Calculate similarity between one query molecule and a vector of target molecules.

This function efficiently computes pairwise similarities between a single query molecule
and multiple target molecules, useful for similarity searches and ranking.

# Arguments
- `query_mol`: The reference molecule to compare against
- `target_mols`: Vector of molecules to compare with the query
- `similarity_function`: Function to use for similarity calculation (default: `tanimoto_similarity`)
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Vector{Union{Float64,Missing}}`: Vector of similarity scores, `missing` for invalid molecules

# Examples
```julia
query = mol_from_smiles("CCO")
targets = [mol_from_smiles("CCC"), mol_from_smiles("CC"), mol_from_smiles("CCCO")]
similarities = bulk_similarity(query, targets)
most_similar_idx = argmax(similarities)
```
"""
function bulk_similarity(
    query_mol::Molecule,
    target_mols::Vector{Union{Molecule, Missing}};
    similarity_function=tanimoto_similarity,
    fingerprint_type=:morgan,
    kwargs...,
)
    return [
        similarity_function(query_mol, mol; fingerprint_type=fingerprint_type, kwargs...)
        for mol in target_mols
    ]
end

function bulk_similarity(
    query_mol::Molecule,
    target_mols::Vector{Molecule};
    similarity_function=tanimoto_similarity,
    fingerprint_type=:morgan,
    kwargs...,
)
    return [
        similarity_function(query_mol, mol; fingerprint_type=fingerprint_type, kwargs...)
        for mol in target_mols
    ]
end

"""
    similarity_matrix(mols::Vector{Molecule};
                     similarity_function=tanimoto_similarity, fingerprint_type=:morgan, kwargs...)

Compute a symmetric similarity matrix for a vector of molecules.

This function creates an n×n symmetric matrix where element (i,j) represents the
similarity between molecules i and j. The diagonal elements are always 1.0 (self-similarity).

# Arguments
- `mols`: Vector of molecules to compare
- `similarity_function`: Function to use for similarity calculation (default: `tanimoto_similarity`)
- `fingerprint_type`: Type of fingerprint to use (`:morgan`, `:rdk`, `:maccs`)
- `kwargs...`: Additional arguments passed to fingerprint generation functions

# Returns
- `Matrix{Union{Float64,Missing}}`: Symmetric similarity matrix, `missing` for invalid molecules

# Examples
```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("CCC"), mol_from_smiles("CCCO")]
sim_matrix = similarity_matrix(mols)
# Access similarity between molecules 1 and 3
sim_1_3 = sim_matrix[1, 3]
```
"""
function similarity_matrix(
    mols::Vector{Union{Molecule, Missing}};
    similarity_function=tanimoto_similarity,
    fingerprint_type=:morgan,
    kwargs...,
)
    n = length(mols)
    sim_matrix = Matrix{Union{Float64, Missing}}(undef, n, n)

    for i in 1:n
        for j in 1:n
            if i == j
                sim_matrix[i, j] = 1.0
            elseif i < j
                sim = similarity_function(
                    mols[i], mols[j]; fingerprint_type=fingerprint_type, kwargs...
                )
                sim_matrix[i, j] = sim
                sim_matrix[j, i] = sim  # Symmetric matrix
            end
        end
    end

    return sim_matrix
end

function similarity_matrix(
    mols::Vector{Molecule};
    similarity_function=tanimoto_similarity,
    fingerprint_type=:morgan,
    kwargs...,
)
    n = length(mols)
    sim_matrix = Matrix{Float64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            if i == j
                sim_matrix[i, j] = 1.0
            elseif i < j
                sim = similarity_function(
                    mols[i], mols[j]; fingerprint_type=fingerprint_type, kwargs...
                )
                sim_matrix[i, j] = sim
                sim_matrix[j, i] = sim  # Symmetric matrix
            end
        end
    end

    return sim_matrix
end
