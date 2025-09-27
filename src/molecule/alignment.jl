#######################################################
# Molecular Alignment Functions
#######################################################
using LinearAlgebra

"""
    align_mol(probe_mol::Molecule, ref_mol::Molecule;
             probe_conf_id::Int=-1, ref_conf_id::Int=-1,
             atom_map=nothing, weights=nothing, reflect::Bool=false,
             max_iterations::Int=50) -> Float64

Optimally align a probe molecule to a reference molecule to minimize RMSD.

# Arguments

  - `probe_mol::Molecule`: The molecule to be aligned (will be modified in place)
  - `ref_mol::Molecule`: The reference molecule to align to (remains unchanged)
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `atom_map=nothing`: Optional atom mapping as vector of tuples `[(probe_idx, ref_idx), ...]` using 1-based indexing
  - `weights=nothing`: Optional vector of weights for atoms during alignment (same order as atoms)
  - `reflect::Bool=false`: Whether to allow reflection (improper rotations) during alignment
  - `max_iterations::Int=50`: Maximum iterations for iterative alignment algorithms

# Returns

  - `Float64`: The root-mean-square deviation (RMSD) after alignment, or `Inf` if alignment fails

# Throws

  - `ArgumentError`: If either molecule is invalid

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - The probe molecule is modified in place - its coordinates are transformed
  - Atom mapping uses 1-based Julia indexing (not 0-based Python indexing)
  - Returns `Inf` if no valid alignment can be found

# Examples

```julia
# Basic alignment
mol1_base = mol_from_smiles("CCO")
mol2_base = mol_from_smiles("CCO")
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule
rmsd = align_mol(mol1, mol2)

# Alignment with atom mapping (align first 3 atoms)
atom_map = [(1, 1), (2, 2), (3, 3)]
rmsd = align_mol(mol1, mol2; atom_map = atom_map)

# Alignment allowing reflection
rmsd = align_mol(mol1, mol2; reflect = true)
```
"""
function align_mol(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    atom_map = nothing,
    weights = nothing,
    reflect::Bool = false,
    max_iterations::Int = 50,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        # Use RDKit's AlignMol function
        if atom_map === nothing && weights === nothing
            # Simple alignment without atom mapping or weights
            rmsd = _align_mol(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                reflect = reflect,
                maxIters = max_iterations,
            )
        else
            # Alignment with atom mapping and/or weights
            atom_map_py =
                atom_map === nothing ? pylist([]) :
                pylist([(i - 1, j - 1) for (i, j) in atom_map])
            weights_py = weights === nothing ? pylist([]) : pylist(weights)

            rmsd = _align_mol_with_map(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol,
                atom_map_py,
                weights_py;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                reflect = reflect,
                maxIters = max_iterations,
            )
        end

        return pyconvert(Float64, rmsd)
    catch e
        @warn "Error in align_mol: $e"
        return Inf
    end
end

"""
    calc_rms(probe_mol::Molecule, ref_mol::Molecule;
            probe_conf_id::Int=-1, ref_conf_id::Int=-1,
            weights=nothing, transform_probe::Bool=false) -> Float64

Calculate the root-mean-square distance between two molecules without performing alignment.

This function computes the RMS distance between corresponding atoms in two molecules. Unlike
`align_mol`, this function does not perform any alignment - it calculates the RMS distance
between the molecules in their current orientations. Optionally, it can consider molecular
symmetry and transform the probe molecule to find the best RMS.

# Arguments

  - `probe_mol::Molecule`: The probe molecule for comparison
  - `ref_mol::Molecule`: The reference molecule for comparison
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `weights=nothing`: Optional vector of weights for atoms during calculation
  - `transform_probe::Bool=false`: If `true`, optimally transform the probe to minimize RMS (similar to alignment)

# Returns

  - `Float64`: The root-mean-square distance between molecules in Å, or `Inf` if calculation fails

# Throws

  - `ArgumentError`: If either molecule is invalid

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - If `transform_probe=false`, calculates RMS in current orientations (no alignment)
  - If `transform_probe=true`, finds optimal transformation to minimize RMS
  - Can handle molecular symmetry automatically
  - For explicit atom mapping, use the `align_mol` function instead

# Examples

```julia
# Basic RMS calculation (no alignment)
mol1_base = mol_from_smiles("CCO")
mol2_base = mol_from_smiles("CCO")
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule
rms = calc_rms(mol1, mol2)

# RMS with optimal transformation
rms_aligned = calc_rms(mol1, mol2; transform_probe = true)

# RMS with weights
weights = [1.0, 1.0, 2.0]  # Give more weight to the third atom
rms_weighted = calc_rms(mol1, mol2; weights = weights)
```
"""
function calc_rms(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    weights = nothing,
    transform_probe::Bool = false,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        if weights === nothing
            rms = _calc_rms(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                transform = transform_probe,
            )
        else
            weights_py = pylist(weights)
            rms = _calc_rms_with_weights(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol,
                weights_py;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                transform = transform_probe,
            )
        end

        return pyconvert(Float64, rms)
    catch e
        @warn "Error in calc_rms: $e"
        return Inf
    end
end

"""
    get_best_rms(probe_mol::Molecule, ref_mol::Molecule;
                probe_conf_id::Int=-1, ref_conf_id::Int=-1,
                weights=nothing, max_matches::Int=1000000) -> Float64

Calculate the best possible RMS between two molecules considering molecular symmetry.

This function finds the optimal RMS distance by considering all possible symmetry-equivalent
alignments between the two molecules. It explores different atom mappings that arise from
molecular symmetry to find the minimum possible RMS distance. This is particularly useful
for symmetric molecules like benzene where multiple equivalent alignments are possible.

# Arguments

  - `probe_mol::Molecule`: The probe molecule for comparison
  - `ref_mol::Molecule`: The reference molecule for comparison
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `weights=nothing`: Optional vector of weights for atoms during calculation
  - `max_matches::Int=1000000`: Maximum number of symmetry-equivalent matches to consider

# Returns

  - `Float64`: The minimum RMS distance considering all symmetry-equivalent alignments in Å, or `Inf` if calculation fails

# Throws

  - `ArgumentError`: If either molecule is invalid

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - This function is computationally more expensive than `calc_rms` as it explores multiple alignments
  - Particularly useful for symmetric molecules (e.g., benzene, cyclohexane)
  - The probe molecule is transformed to find the best possible alignment
  - May take longer for highly symmetric molecules
  - For explicit atom mapping, use the `align_mol` function instead

# Examples

```julia
# Best RMS for symmetric molecules (benzene)
mol1_base = mol_from_smiles("c1ccccc1")  # benzene
mol2_base = mol_from_smiles("c1ccccc1")  # benzene
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule

# Regular RMS (may not be optimal due to symmetry)
regular_rms = calc_rms(mol1, mol2)

# Best RMS considering symmetry
best_rms = get_best_rms(mol1, mol2)
# best_rms <= regular_rms due to symmetry considerations

# Limit the number of symmetry matches to consider
best_rms_limited = get_best_rms(mol1, mol2; max_matches = 100)
```
"""
function get_best_rms(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    weights = nothing,
    max_matches::Int = 1000000,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        if weights === nothing
            rms = _get_best_rms(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                maxMatches = max_matches,
            )
        else
            weights_py = pylist(weights)
            rms = _get_best_rms_with_weights(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol,
                weights_py;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                maxMatches = max_matches,
            )
        end

        return pyconvert(Float64, rms)
    catch e
        @warn "Error in get_best_rms: $e"
        return Inf
    end
end

"""
    get_alignment_transform(probe_mol::Molecule, ref_mol::Molecule;
                           probe_conf_id::Int=-1, ref_conf_id::Int=-1,
                           weights=nothing, reflect::Bool=false) -> Matrix{Float64}

Compute the transformation matrix required to align a probe molecule to a reference molecule.

# Arguments

  - `probe_mol::Molecule`: The molecule to be aligned
  - `ref_mol::Molecule`: The reference molecule
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default)
  - `weights=nothing`: Optional weights for atoms during alignment
  - `reflect::Bool=false`: Whether to allow reflection during alignment

# Returns

  - `Matrix{Float64}`: 4x4 transformation matrix for the alignment

# Notes

  - For explicit atom mapping, use the `align_mol` function instead
  - The transformation matrix can be applied using `apply_transform`

# Example

```julia
mol1_base = mol_from_smiles("CCO")
mol2_base = mol_from_smiles("CCO")
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule
transform = get_alignment_transform(mol1, mol2)
```
"""
function get_alignment_transform(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    weights = nothing,
    reflect::Bool = false,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        if weights === nothing
            transform_py = _get_alignment_transform(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                reflect = reflect,
            )
        else
            weights_py = pylist(weights)
            transform_py = _get_alignment_transform_with_weights(
                probe_mol._rdkit_mol,
                ref_mol._rdkit_mol,
                weights_py;
                prbCid = probe_conf_id,
                refCid = ref_conf_id,
                reflect = reflect,
            )
        end

        # RDKit returns a tuple (RMSD, transform_matrix)
        # Extract just the transformation matrix (second element)
        transform_matrix = transform_py[1]  # Python 0-based indexing, so [1] is second element
        transform_array = pyconvert(Array, transform_matrix)
        return reshape(transform_array, 4, 4)
    catch e
        @warn "Error in get_alignment_transform: $e"
        return Matrix{Float64}(I, 4, 4)  # Return identity matrix on error
    end
end

"""
    random_transform(mol::Molecule; conf_id::Int=-1,
                    seed::Union{Int,Nothing}=nothing) -> Molecule

Perform a random rigid body transformation (rotation + translation) on a molecule.

# Arguments

  - `mol::Molecule`: The molecule to transform
  - `conf_id::Int=-1`: Conformer ID to transform (-1 for default)
  - `seed::Union{Int,Nothing}=nothing`: Random seed for reproducibility

# Returns

  - `Molecule`: The transformed molecule (original molecule is also modified)

# Example

```julia
mol_base = mol_from_smiles("CCO")
conformers = generate_3d_conformers(mol_base, 1)
mol = conformers[1].molecule
transformed_mol = random_transform(mol; seed = 42)
```
"""
function random_transform(
    mol::Molecule; conf_id::Int = -1, seed::Union{Int, Nothing} = nothing
)
    if !mol.valid
        throw(ArgumentError("Molecule must be valid"))
    end

    try
        # Pass seed directly to _random_transform (RDKit's RandomTransform accepts seed parameter)
        seed_value = seed === nothing ? -1 : seed
        _random_transform(mol._rdkit_mol; confId = conf_id, seed = seed_value)

        return mol
    catch e
        @warn "Error in random_transform: $e"
        return mol
    end
end

"""
    apply_transform(mol::Molecule, transform::Matrix{Float64}; conf_id::Int=-1) -> Molecule

Apply a transformation matrix to a molecule's coordinates.

# Arguments

  - `mol::Molecule`: The molecule to transform
  - `transform::Matrix{Float64}`: 4x4 transformation matrix
  - `conf_id::Int=-1`: Conformer ID to transform (-1 for default)

# Returns

  - `Molecule`: The transformed molecule (original molecule is also modified)

# Example

```julia
mol_base = mol_from_smiles("CCO")
conformers = generate_3d_conformers(mol_base, 1)
mol = conformers[1].molecule
# Get identity transformation as example
transform = Matrix{Float64}(I, 4, 4)
transformed_mol = apply_transform(mol, transform)
```
"""
function apply_transform(mol::Molecule, transform::Matrix{Float64}; conf_id::Int = -1)
    if !mol.valid
        throw(ArgumentError("Molecule must be valid"))
    end

    if size(transform) != (4, 4)
        throw(ArgumentError("Transform matrix must be 4x4"))
    end

    try
        # Get the specific conformer
        conf = mol._rdkit_mol.GetConformer(conf_id)

        # Convert Julia matrix to numpy array
        numpy = pyimport("numpy")
        transform_matrix = numpy.array(transform; dtype = numpy.float64)

        # Apply the transformation directly
        rdMolTransforms = pyimport("rdkit.Chem.rdMolTransforms")
        rdMolTransforms.TransformConformer(conf, transform_matrix)

        return mol
    catch e
        @warn "Error in apply_transform: $e"
        return mol
    end
end

"""
    O3AResult

Structure to hold results from O3A (Open3DAlign) alignment operations.

# Fields

  - `score::Float64`: O3A alignment score (higher values indicate better alignment)
  - `rmsd::Float64`: Root Mean Square Deviation after alignment in Å
  - `transform::Matrix{Float64}`: 4x4 transformation matrix applied to the probe molecule
  - `matched_atoms::Vector{Tuple{Int,Int}}`: Pairs of matched atom indices using 1-based indexing (probe_idx, ref_idx)

# Notes

  - Failed alignments return `score = -1.0` and `rmsd = Inf`
  - The transformation matrix is in homogeneous coordinates format
  - Atom indices in `matched_atoms` use Julia's 1-based indexing convention
"""
struct O3AResult
    score::Float64
    rmsd::Float64
    transform::Matrix{Float64}
    matched_atoms::Vector{Tuple{Int, Int}}
end

function Base.show(io::IO, result::O3AResult)
    print(
        io,
        "O3AResult(score=$(round(result.score, digits=3)), rmsd=$(round(result.rmsd, digits=3)), $(length(result.matched_atoms)) matched atoms)",
    )
end

"""
    get_o3a(probe_mol::Molecule, ref_mol::Molecule;
           probe_conf_id::Int=-1, ref_conf_id::Int=-1,
           reflect::Bool=false, accuracy::Float64=0.0001,
           attempt_generic_features::Bool=true,
           prune_conformers::Bool=true) -> O3AResult

Perform Open3DAlign (O3A) alignment using MMFF atom types for molecular overlay.
This method aligns molecules based on their 3D pharmacophoric features using MMFF
molecular properties to define feature points.

# Arguments

  - `probe_mol::Molecule`: The molecule to be aligned (will be modified in place)
  - `ref_mol::Molecule`: The reference molecule to align to (remains unchanged)
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `reflect::Bool=false`: Whether to allow reflection during alignment
  - `accuracy::Float64=0.0001`: Accuracy threshold for feature matching (currently not used)
  - `attempt_generic_features::Bool=true`: Whether to use generic pharmacophoric features (currently not used)
  - `prune_conformers::Bool=true`: Whether to prune conformers during alignment (currently not used)

# Returns

  - `O3AResult`: Structure containing:

      + `score::Float64`: O3A alignment score (higher is better)
      + `rmsd::Float64`: RMSD after alignment in Å
      + `transform::Matrix{Float64}`: 4x4 transformation matrix applied to probe
      + `matched_atoms::Vector{Tuple{Int,Int}}`: Pairs of matched atom indices (probe, reference)

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - The probe molecule is modified in place during alignment
  - MMFF molecular properties are automatically computed for feature generation
  - Returns failed result (score=-1.0, rmsd=Inf) if alignment fails

# Example

```julia
mol1_base = mol_from_smiles("c1ccc(cc1)CCN")  # phenethylamine
mol2_base = mol_from_smiles("c1ccc(cc1)CCNC") # N-methylphenethylamine
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule
result = get_o3a(mol1, mol2)
```
"""
function get_o3a(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    reflect::Bool = false,
    accuracy::Float64 = 0.0001,
    attempt_generic_features::Bool = true,
    prune_conformers::Bool = true,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        # Get O3A alignment object
        o3a = _get_o3a(
            probe_mol._rdkit_mol,
            ref_mol._rdkit_mol;
            prbCid = probe_conf_id,
            refCid = ref_conf_id,
            reflect = reflect,
            accuracy = accuracy,
            attemptGenericFeatures = attempt_generic_features,
            pruneConfs = prune_conformers,
        )

        # Perform alignment and get RMSD
        rmsd = pyconvert(Float64, o3a.Align())

        # Get score
        score = pyconvert(Float64, o3a.Score())

        # Get transformation matrix
        transform_py = o3a.Trans()
        # Trans() returns a tuple (RMSD, transform_matrix), get the matrix (second element)
        transform_matrix = transform_py[1]  # Python 0-based indexing, so [1] is second element
        transform_array = pyconvert(Array, transform_matrix)
        transform = transform_array  # Already a 4x4 matrix, no reshape needed

        # Get matched atoms
        matches = o3a.Matches()
        matched_atoms = Tuple{Int, Int}[]
        for match in matches
            probe_idx = pyconvert(Int, match[0]) + 1  # Convert to 1-based indexing
            ref_idx = pyconvert(Int, match[1]) + 1
            push!(matched_atoms, (probe_idx, ref_idx))
        end

        return O3AResult(score, rmsd, transform, matched_atoms)
    catch e
        @warn "Error in get_o3a: $e"
        return O3AResult(-1.0, Inf, Matrix{Float64}(I, 4, 4), Tuple{Int, Int}[])
    end
end

"""
    get_crippen_o3a(probe_mol::Molecule, ref_mol::Molecule;
                    probe_conf_id::Int=-1, ref_conf_id::Int=-1,
                    reflect::Bool=false, accuracy::Float64=0.0001,
                    attempt_generic_features::Bool=true,
                    prune_conformers::Bool=true) -> O3AResult

Perform Open3DAlign (O3A) alignment using Crippen atom contributions for molecular overlay.
This method uses Crippen LogP and molar refractivity contributions to define feature points
for pharmacophore-based alignment.

# Arguments

  - `probe_mol::Molecule`: The molecule to be aligned (will be modified in place)
  - `ref_mol::Molecule`: The reference molecule to align to (remains unchanged)
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `reflect::Bool=false`: Whether to allow reflection during alignment
  - `accuracy::Float64=0.0001`: Accuracy threshold for feature matching (currently not used)
  - `attempt_generic_features::Bool=true`: Whether to use generic pharmacophoric features (currently not used)
  - `prune_conformers::Bool=true`: Whether to prune conformers during alignment (currently not used)

# Returns

  - `O3AResult`: Structure containing:

      + `score::Float64`: O3A alignment score (higher is better)
      + `rmsd::Float64`: RMSD after alignment in Å
      + `transform::Matrix{Float64}`: 4x4 transformation matrix applied to probe
      + `matched_atoms::Vector{Tuple{Int,Int}}`: Pairs of matched atom indices (probe, reference)

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - The probe molecule is modified in place during alignment
  - Crippen contributions (LogP and molar refractivity) are automatically computed for feature generation
  - Returns failed result (score=-1.0, rmsd=Inf) if alignment fails
  - Generally more robust than MMFF-based alignment for diverse molecule types

# Example

```julia
mol1_base = mol_from_smiles("CCc1ccccc1")  # ethylbenzene
mol2_base = mol_from_smiles("CCc1ccc(O)cc1")  # 4-ethylphenol
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule
result = get_crippen_o3a(mol1, mol2)
```
"""
function get_crippen_o3a(
    probe_mol::Molecule,
    ref_mol::Molecule;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    reflect::Bool = false,
    accuracy::Float64 = 0.0001,
    attempt_generic_features::Bool = true,
    prune_conformers::Bool = true,
)
    if !probe_mol.valid || !ref_mol.valid
        throw(ArgumentError("Both molecules must be valid"))
    end

    try
        # Get Crippen O3A alignment object
        o3a = _get_crippen_o3a(
            probe_mol._rdkit_mol,
            ref_mol._rdkit_mol;
            prbCid = probe_conf_id,
            refCid = ref_conf_id,
            reflect = reflect,
            accuracy = accuracy,
            attemptGenericFeatures = attempt_generic_features,
            pruneConfs = prune_conformers,
        )

        # Perform alignment and get RMSD
        rmsd = pyconvert(Float64, o3a.Align())

        # Get score
        score = pyconvert(Float64, o3a.Score())

        # Get transformation matrix
        transform_py = o3a.Trans()
        # Trans() returns a tuple (RMSD, transform_matrix), get the matrix (second element)
        transform_matrix = transform_py[1]  # Python 0-based indexing, so [1] is second element
        transform_array = pyconvert(Array, transform_matrix)
        transform = transform_array  # Already a 4x4 matrix, no reshape needed

        # Get matched atoms
        matches = o3a.Matches()
        matched_atoms = Tuple{Int, Int}[]
        for match in matches
            probe_idx = pyconvert(Int, match[0]) + 1  # Convert to 1-based indexing
            ref_idx = pyconvert(Int, match[1]) + 1
            push!(matched_atoms, (probe_idx, ref_idx))
        end

        return O3AResult(score, rmsd, transform, matched_atoms)
    catch e
        @warn "Error in get_crippen_o3a: $e"
        return O3AResult(-1.0, Inf, Matrix{Float64}(I, 4, 4), Tuple{Int, Int}[])
    end
end

"""
    o3a_align!(probe_mol::Molecule, ref_mol::Molecule, alignment_method::Symbol=:mmff;
              probe_conf_id::Int=-1, ref_conf_id::Int=-1,
              reflect::Bool=false, accuracy::Float64=0.0001) -> O3AResult

Convenience function to perform O3A (Open3DAlign) alignment and modify the probe molecule in place.
This function provides a simple interface to choose between MMFF-based or Crippen-based alignment.

# Arguments

  - `probe_mol::Molecule`: The molecule to be aligned (will be modified in place)
  - `ref_mol::Molecule`: The reference molecule to align to (remains unchanged)
  - `alignment_method::Symbol=:mmff`: Alignment method (:mmff or :crippen)

      + `:mmff`: Uses MMFF molecular properties for feature generation (calls `get_o3a`)
      + `:crippen`: Uses Crippen LogP/molar refractivity contributions (calls `get_crippen_o3a`)
  - `probe_conf_id::Int=-1`: Conformer ID of the probe molecule (-1 for default/first conformer)
  - `ref_conf_id::Int=-1`: Conformer ID of the reference molecule (-1 for default/first conformer)
  - `reflect::Bool=false`: Whether to allow reflection during alignment
  - `accuracy::Float64=0.0001`: Accuracy threshold for feature matching (currently not used)

# Returns

  - `O3AResult`: Structure containing:

      + `score::Float64`: O3A alignment score (higher is better)
      + `rmsd::Float64`: RMSD after alignment in Å
      + `transform::Matrix{Float64}`: 4x4 transformation matrix applied to probe
      + `matched_atoms::Vector{Tuple{Int,Int}}`: Pairs of matched atom indices (probe, reference)

# Notes

  - Both molecules must have 3D coordinates (conformers)
  - The probe molecule is modified in place during alignment
  - Returns failed result (score=-1.0, rmsd=Inf) if alignment fails
  - For most cases, `:crippen` method is more robust than `:mmff`

# Example

```julia
mol1_base = mol_from_smiles("c1ccccc1CCN")
mol2_base = mol_from_smiles("c1ccccc1CCNC")
conformers1 = generate_3d_conformers(mol1_base, 1)
conformers2 = generate_3d_conformers(mol2_base, 1)
mol1 = conformers1[1].molecule
mol2 = conformers2[1].molecule

# MMFF-based alignment
result1 = o3a_align!(mol1, mol2, :mmff)
println("MMFF result: \$(result1)")

# Crippen-based alignment (generally more robust)
result2 = o3a_align!(mol1, mol2, :crippen)
println("Crippen result: \$(result2)")
```
"""
function o3a_align!(
    probe_mol::Molecule,
    ref_mol::Molecule,
    alignment_method::Symbol = :mmff;
    probe_conf_id::Int = -1,
    ref_conf_id::Int = -1,
    reflect::Bool = false,
    accuracy::Float64 = 0.0001,
)
    if alignment_method == :mmff
        return get_o3a(
            probe_mol,
            ref_mol;
            probe_conf_id = probe_conf_id,
            ref_conf_id = ref_conf_id,
            reflect = reflect,
            accuracy = accuracy,
        )
    elseif alignment_method == :crippen
        return get_crippen_o3a(
            probe_mol,
            ref_mol;
            probe_conf_id = probe_conf_id,
            ref_conf_id = ref_conf_id,
            reflect = reflect,
            accuracy = accuracy,
        )
    else
        throw(ArgumentError("alignment_method must be :mmff or :crippen"))
    end
end
