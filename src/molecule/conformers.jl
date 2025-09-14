#######################################################
# 3D Conformer Generation and Optimization
#######################################################

function _rdkit_conformer()
    return @pyconst(pyimport("rdkit.Chem.rdDistGeom"))
end

function _rdkit_forcefield()
    return @pyconst(pyimport("rdkit.Chem.rdForceFieldHelpers"))
end

function _rdkit_rdmolops()
    return @pyconst(pyimport("rdkit.Chem.rdMolOps"))
end

function _rdkit_allchem()
    return @pyconst(pyimport("rdkit.Chem.AllChem"))
end

"""
    ConformerResult

Structure to hold conformer information including ID, energy, and optimization status.

# Fields
- `id::Int`: Conformer ID (0-based)
- `energy::Float64`: Energy in kcal/mol
- `optimized::Bool`: Whether the conformer has been optimized
- `converged::Bool`: Whether optimization converged (if optimized)
"""
struct ConformerResult
    id::Int
    energy::Float64
    optimized::Bool
    converged::Bool
end

function Base.show(io::IO, conf::ConformerResult)
    status = conf.optimized ? (conf.converged ? "optimized" : "opt_failed") : "unoptimized"
    print(
        io,
        "ConformerResult(id=$(conf.id), energy=$(round(conf.energy, digits=2)) kcal/mol, $(status))",
    )
end

"""
    ConformerMolecule

Structure to hold both conformer metadata and the actual molecule with 3D coordinates.

# Fields
- `molecule::Molecule`: The molecule with 3D coordinates stored as properties
- `conformer_result::ConformerResult`: Metadata about the conformer (ID, energy, etc.)
"""
struct ConformerMolecule
    molecule::Molecule
    conformer_result::ConformerResult
end

function Base.show(io::IO, conf_mol::ConformerMolecule)
    print(io, "ConformerMolecule($(conf_mol.conformer_result))")
end

"""
    _extract_3d_coordinates(mol_rdkit, conformer_id::Int)

Extract 3D coordinates from an RDKit molecule conformer as a Matrix{Float64}.
Returns a matrix of size (num_atoms, 3) where each row contains [x, y, z] coordinates.
"""
function _extract_3d_coordinates(mol_rdkit, conformer_id::Int)
    try
        conf = mol_rdkit.GetConformer(conformer_id)
        num_atoms = pyconvert(Int, mol_rdkit.GetNumAtoms())
        coords = Matrix{Float64}(undef, num_atoms, 3)

        for i in 0:(num_atoms - 1)
            pos = conf.GetAtomPosition(i)
            coords[i + 1, 1] = pyconvert(Float64, pos.x)
            coords[i + 1, 2] = pyconvert(Float64, pos.y)
            coords[i + 1, 3] = pyconvert(Float64, pos.z)
        end

        return coords
    catch e
        @warn "Error extracting 3D coordinates: $e"
        return Matrix{Float64}(undef, 0, 3)
    end
end

"""
    generate_3d_conformers(mol::Molecule, num_conformers::Int=1; 
                          optimize::Bool=true, force_field::Symbol=:mmff,
                          max_attempts::Int=1000, random_seed::Union{Int,Nothing}=nothing,
                          prune_rms_thresh::Float64=1.0, max_iterations::Int=200)

Generate and optionally optimize 3D conformers for a molecule, returning both the conformer metadata 
and actual molecules with 3D coordinates stored as properties.

# Arguments
- `mol::Molecule`: The molecule to generate conformers for
- `num_conformers::Int=1`: Number of conformers to generate
- `optimize::Bool=true`: Whether to optimize conformers after generation
- `force_field::Symbol=:mmff`: Force field for optimization (:mmff or :uff)
- `max_attempts::Int=1000`: Maximum attempts for conformer generation
- `random_seed::Union{Int,Nothing}=nothing`: Random seed for reproducibility
- `prune_rms_thresh::Float64=1.0`: RMS threshold for pruning similar conformers
- `max_iterations::Int=200`: Maximum optimization iterations

# Returns
- `Vector{ConformerMolecule}`: Conformer molecules sorted by energy (lowest first)

# Example
```julia
mol = mol_from_smiles("CCO")
conformer_mols = generate_3d_conformers(mol, 10)
println("Best conformer: \$(conformer_mols[1].conformer_result)")
coords = conformer_mols[1].molecule.props[:coordinates_3d]
```
"""
function generate_3d_conformers(
    mol::Molecule,
    num_conformers::Int=1;
    optimize::Bool=true,
    force_field::Symbol=:mmff,
    max_attempts::Int=1000,
    random_seed::Union{Int, Nothing}=nothing,
    prune_rms_thresh::Float64=1.0,
    max_iterations::Int=200,
)
    if !mol.valid
        return ConformerMolecule[]
    end

    try
        # Work with a copy of the molecule to avoid modifying the original during conformer generation
        mol_copy = _rdkit_allchem().AddHs(mol._rdkit_mol)

        # Set up parameters for conformer generation
        params = _rdkit_conformer().ETKDGv3()
        if random_seed !== nothing
            params.randomSeed = random_seed
        end
        params.pruneRmsThresh = prune_rms_thresh
        # Use numThreads instead of maxAttempts for recent RDKit versions
        try
            params.maxAttempts = max_attempts
        catch
            # For newer RDKit versions, use a different approach
            params.numThreads = 1  # Use single thread for consistent results
        end

        # Generate conformers
        if num_conformers == 1
            conf_id = _rdkit_conformer().EmbedMolecule(mol_copy, params)
            conf_ids = pyconvert(Int, conf_id) == -1 ? Int[] : [pyconvert(Int, conf_id)]
        else
            try
                conf_ids = _rdkit_conformer().EmbedMultipleConfs(
                    mol_copy; numConfs=num_conformers, params=params
                )
                conf_ids = [pyconvert(Int, cid) for cid in conf_ids]
            catch e
                # Try alternative method for newer RDKit versions
                conf_ids = _rdkit_conformer().EmbedMultipleConfs(
                    mol_copy, num_conformers, params
                )
                conf_ids = [pyconvert(Int, cid) for cid in conf_ids]
            end
        end

        if isempty(conf_ids)
            return ConformerMolecule[]
        end

        # Create results vector
        results = Vector{ConformerMolecule}()

        for (i, conf_id) in enumerate(conf_ids)
            # Extract 3D coordinates
            coords_3d = _extract_3d_coordinates(mol_copy, conf_id)

            # Create a new molecule with the same properties but add 3D coordinates
            mol_props = copy(mol.props)
            mol_props[:coordinates_3d] = coords_3d
            mol_props[:conformer_id] = conf_id

            # Create a copy of the original molecule and add the conformer
            new_rdkit_mol = @pyconst(pyimport("rdkit.Chem").Mol)(mol._rdkit_mol)
            new_rdkit_mol.RemoveAllConformers()

            # Get conformer from molecule without hydrogens to match original structure
            mol_no_hs = _rdkit_allchem().RemoveHs(mol_copy)
            conf_no_hs = mol_no_hs.GetConformer(conf_id)
            new_rdkit_mol.AddConformer(conf_no_hs; assignId=true)

            # Create new Molecule object
            new_mol = Molecule(;
                _rdkit_mol=new_rdkit_mol, valid=true, source=mol.source, props=mol_props
            )

            # Calculate energy and create ConformerResult
            if optimize
                success, energy = _optimize_conformer_with_mol!(
                    mol_copy, conf_id, force_field, max_iterations
                )
                conf_result = ConformerResult(conf_id, energy, true, Bool(success))
            else
                energy = _calculate_conformer_energy_with_mol(
                    mol_copy, conf_id, force_field
                )
                energy = energy === nothing ? Inf : energy
                conf_result = ConformerResult(conf_id, energy, false, false)
            end

            push!(results, ConformerMolecule(new_mol, conf_result))
        end

        # Sort by energy (lowest first)
        sort!(results; by=x -> x.conformer_result.energy)

        return results
    catch e
        @warn "Error in conformer generation: $e"
        return ConformerMolecule[]
    end
end

"""
    _optimize_conformer_with_mol!(mol_rdkit, conformer_id::Int, force_field::Symbol, max_iterations::Int)

Internal function to optimize a single conformer using the provided RDKit molecule.
"""
function _optimize_conformer_with_mol!(
    mol_rdkit, conformer_id::Int, force_field::Symbol, max_iterations::Int
)
    try
        if force_field == :mmff
            ff = _rdkit_forcefield().MMFFGetMoleculeForceField(
                mol_rdkit,
                _rdkit_forcefield().MMFFGetMoleculeProperties(mol_rdkit);
                confId=conformer_id,
            )
        elseif force_field == :uff
            ff = _rdkit_forcefield().UFFGetMoleculeForceField(
                mol_rdkit; confId=conformer_id
            )
        else
            error("Unsupported force field: $force_field")
        end

        if ff === nothing
            return false, Inf
        end

        # Optimize
        converged = ff.Minimize(; maxIts=max_iterations)
        final_energy = ff.CalcEnergy()

        return pyconvert(Int, converged) == 0, pyconvert(Float64, final_energy)
    catch
        return false, Inf
    end
end

"""
    _optimize_conformer!(mol::Molecule, conformer_id::Int, force_field::Symbol, max_iterations::Int)

Internal function to optimize a single conformer.
"""
function _optimize_conformer!(
    mol::Molecule, conformer_id::Int, force_field::Symbol, max_iterations::Int
)
    try
        mol_with_hs = _rdkit_allchem().AddHs(mol._rdkit_mol)

        if force_field == :mmff
            ff = _rdkit_forcefield().MMFFGetMoleculeForceField(
                mol_with_hs,
                _rdkit_forcefield().MMFFGetMoleculeProperties(mol_with_hs);
                confId=conformer_id,
            )
        elseif force_field == :uff
            ff = _rdkit_forcefield().UFFGetMoleculeForceField(
                mol_with_hs; confId=conformer_id
            )
        else
            error("Unsupported force field: $force_field")
        end

        if ff === nothing
            return false, Inf
        end

        # Optimize
        converged = ff.Minimize(; maxIts=max_iterations)
        final_energy = ff.CalcEnergy()

        # Copy the optimized conformer back
        conf = mol_with_hs.GetConformer(conformer_id)
        mol._rdkit_mol.RemoveConformer(conformer_id)
        mol._rdkit_mol.AddConformer(conf; assignId=true)

        return pyconvert(Int, converged) == 0, pyconvert(Float64, final_energy)
    catch
        return false, Inf
    end
end

"""
    _calculate_conformer_energy_with_mol(mol_rdkit, conformer_id::Int, force_field::Symbol)

Internal function to calculate conformer energy without optimization using the provided RDKit molecule.
"""
function _calculate_conformer_energy_with_mol(
    mol_rdkit, conformer_id::Int, force_field::Symbol
)
    try
        if force_field == :mmff
            ff = _rdkit_forcefield().MMFFGetMoleculeForceField(
                mol_rdkit,
                _rdkit_forcefield().MMFFGetMoleculeProperties(mol_rdkit);
                confId=conformer_id,
            )
        elseif force_field == :uff
            ff = _rdkit_forcefield().UFFGetMoleculeForceField(
                mol_rdkit; confId=conformer_id
            )
        else
            return nothing
        end

        if ff === nothing
            return nothing
        end

        return pyconvert(Float64, ff.CalcEnergy())
    catch
        return nothing
    end
end

"""
    _calculate_conformer_energy(mol::Molecule, conformer_id::Int, force_field::Symbol)

Internal function to calculate conformer energy without optimization.
"""
function _calculate_conformer_energy(mol::Molecule, conformer_id::Int, force_field::Symbol)
    try
        mol_with_hs = _rdkit_allchem().AddHs(mol._rdkit_mol)

        if force_field == :mmff
            ff = _rdkit_forcefield().MMFFGetMoleculeForceField(
                mol_with_hs,
                _rdkit_forcefield().MMFFGetMoleculeProperties(mol_with_hs);
                confId=conformer_id,
            )
        elseif force_field == :uff
            ff = _rdkit_forcefield().UFFGetMoleculeForceField(
                mol_with_hs; confId=conformer_id
            )
        else
            return nothing
        end

        if ff === nothing
            return nothing
        end

        return pyconvert(Float64, ff.CalcEnergy())
    catch
        return nothing
    end
end

"""
    generate_2d_conformers(mol::Molecule) -> Vector{ConformerMolecule}

Generate 2D coordinates for a molecule using RDKit's 2D coordinate generation.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Vector{ConformerMolecule}`: Vector containing a single ConformerMolecule with 2D coordinates
- `ConformerMolecule[]`: Empty vector if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO")
conformers_2d = generate_2d_conformers(mol)
coords = conformers_2d[1].molecule.props[:coordinates_2d]  # 2D coordinates matrix
```

# Notes
- Generates standardized 2D layout for visualization
- Coordinates are stored in the molecule's properties as `:coordinates_2d`
- Always returns exactly one conformer (or empty vector if invalid)
- Used for 2D molecular visualization and plotting
"""
function generate_2d_conformers(mol::Molecule)
    if !mol.valid
        return ConformerMolecule[]
    end

    mol_copy = @pyconst(pyimport("rdkit.Chem").Mol)(mol._rdkit_mol)
    rdDepictor = @pyconst(pyimport("rdkit.Chem.rdDepictor"))
    rdDepictor.Compute2DCoords(mol_copy)

    conf = mol_copy.GetConformer()
    num_atoms = pyconvert(Int, mol_copy.GetNumAtoms())
    coords_2d = Matrix{Float64}(undef, num_atoms, 2)

    for i in 0:(num_atoms - 1)
        pos = conf.GetAtomPosition(i)
        coords_2d[i + 1, 1] = pyconvert(Float64, pos.x)
        coords_2d[i + 1, 2] = pyconvert(Float64, pos.y)
    end

    mol_props = copy(mol.props)
    mol_props[:coordinates_2d] = coords_2d
    mol_props[:conformer_id] = 0
    mol_props[:conformer_type] = "2d"

    new_mol = Molecule(;
        _rdkit_mol=mol_copy, valid=true, source=mol.source, props=mol_props
    )
    conf_result = ConformerResult(0, 0.0, false, false)

    return [ConformerMolecule(new_mol, conf_result)]
end
