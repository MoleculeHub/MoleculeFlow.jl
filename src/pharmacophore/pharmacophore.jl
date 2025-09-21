#######################################################
# Pharmacophore Features and Analysis
#######################################################

"""
    FeatureFactory

Struct to hold a chemical feature factory for pharmacophore analysis.

# Fields
- `_rdkit_factory::Py`: The underlying RDKit feature factory object
- `feature_families::Vector{String}`: List of available feature families
- `num_definitions::Int`: Number of feature definitions

# Example
```julia
factory = create_feature_factory()
println("Available families: ", factory.feature_families)
```
"""
@kwdef struct FeatureFactory
    _rdkit_factory::Py
    feature_families::Vector{String} = String[]
    num_definitions::Int = 0
end

function Base.show(io::IO, factory::FeatureFactory)
    nfamilies = length(factory.feature_families)
    ndefs = factory.num_definitions

    print(io, "FeatureFactory(")
    print(io, "$nfamilies families, $ndefs definitions")

    if nfamilies > 0
        print(io, ": ")
        if nfamilies <= 6
            # Show all families if 6 or fewer
            print(io, join(factory.feature_families, ", "))
        else
            # Show first 5 and indicate there are more
            print(io, join(factory.feature_families[1:5], ", "))
            print(io, ", ... ($(nfamilies-5) more)")
        end
    end

    print(io, ")")
end

"""
    ChemicalFeature

Represents a chemical feature identified in a molecule.

# Fields
- `family::String`: Feature family (e.g., "Donor", "Acceptor", "Aromatic")
- `type::String`: Specific feature type
- `atom_ids::Vector{Int}`: Atom indices involved in the feature (1-based)
- `position::Vector{Float64}`: 3D coordinates of the feature center
- `id::Int`: Feature identifier

# Example
```julia
features = get_mol_features(mol, factory)
for feature in features
    println("Feature: ", feature.family, " at position ", feature.position)
end
```
"""
@kwdef struct ChemicalFeature
    family::String
    type::String
    atom_ids::Vector{Int}
    position::Vector{Float64}
    id::Int
end

"""
    create_feature_factory(;use_default::Bool = true, fdef_file::Union{String, Nothing} = nothing, fdef_string::Union{String, Nothing} = nothing) -> FeatureFactory

Create a chemical feature factory for pharmacophore analysis.

# Arguments
- `use_default::Bool = true`: Whether to use RDKit's default feature definitions
- `fdef_file::Union{String, Nothing} = nothing`: Path to custom feature definition file
- `fdef_string::Union{String, Nothing} = nothing`: Custom feature definition string

# Returns
- `FeatureFactory`: A feature factory for identifying molecular features

# Examples
```julia
# Use default RDKit feature definitions
factory = create_feature_factory()

# Use custom feature definition file
factory = create_feature_factory(use_default=false, fdef_file="custom.fdef")

# Use custom feature definition string
fdef = \"\"\"
DefineFeature HAcceptor1 [N,O;H0]
   Family HBondAcceptor
   Weights 1.0
EndFeature
\"\"\"
factory = create_feature_factory(use_default=false, fdef_string=fdef)
```

# Notes
- Default factory includes: Donor, Acceptor, NegIonizable, PosIonizable, ZnBinder, Aromatic, Hydrophobe, LumpedHydrophobe
- Custom definitions use SMARTS patterns
"""
function create_feature_factory(;use_default::Bool = true, fdef_file::Union{String, Nothing} = nothing, fdef_string::Union{String, Nothing} = nothing)
    try
        if use_default
            # Use RDKit's default BaseFeatures.fdef
            data_dir = _get_rdconfig_data_dir()
            default_fdef = joinpath(pyconvert(String, data_dir), "BaseFeatures.fdef")
            rdkit_factory = _build_feature_factory(default_fdef)
        elseif fdef_file !== nothing
            rdkit_factory = _build_feature_factory(fdef_file)
        elseif fdef_string !== nothing
            rdkit_factory = _build_feature_factory_from_string(fdef_string)
        else
            throw(ArgumentError("Must specify either use_default=true, fdef_file, or fdef_string"))
        end

        # Extract feature information
        families = pyconvert(Vector{String}, _get_feature_families(rdkit_factory))
        num_defs = pyconvert(Int, _get_num_feature_defs(rdkit_factory))

        return FeatureFactory(
            _rdkit_factory = rdkit_factory,
            feature_families = families,
            num_definitions = num_defs
        )
    catch e
        @warn "Error creating feature factory: $e"
        return FeatureFactory(
            _rdkit_factory = pybuiltins.None,
            feature_families = String[],
            num_definitions = 0
        )
    end
end

"""
    get_mol_features(mol::Molecule, factory::FeatureFactory; conf_id::Int = -1) -> Vector{ChemicalFeature}

Extract chemical features from a molecule using the specified feature factory.

# Arguments
- `mol::Molecule`: Input molecule (must have 2D or 3D conformer with coordinates)
- `factory::FeatureFactory`: Feature factory for feature identification
- `conf_id::Int = -1`: Conformer ID to use (-1 for default)

# Returns
- `Vector{ChemicalFeature}`: List of identified chemical features

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
factory = create_feature_factory()

# Generate coordinates (required for feature extraction)
conformers_2d = generate_2d_conformers(mol)
if !isempty(conformers_2d)
    mol_2d = conformers_2d[1].molecule
    features = get_mol_features(mol_2d, factory)

    println("Found ", length(features), " features:")
    for feature in features
        println("  ", feature.family, " at atoms ", feature.atom_ids)
    end
end
```

# Notes
- **Coordinate Requirement**: Molecule must have 2D or 3D conformer coordinates. Use `generate_2d_conformers()` or `generate_3d_conformers()` first.
- Returns empty vector for molecules without conformers or invalid molecules
- 3D coordinates provide more accurate feature positioning than 2D
"""
function get_mol_features(mol::Molecule, factory::FeatureFactory; conf_id::Int = -1)
    !mol.valid && return ChemicalFeature[]

    try
        rdkit_features = _get_features_for_mol(factory._rdkit_factory, mol._rdkit_mol; conf_id = conf_id)
        features = ChemicalFeature[]

        for (i, rdkit_feature) in enumerate(rdkit_features)
            family = pyconvert(String, rdkit_feature.GetFamily())
            feature_type = pyconvert(String, rdkit_feature.GetType())
            atom_ids = pyconvert(Vector{Int}, rdkit_feature.GetAtomIds()) .+ 1  # Convert to 1-based

            # Extract position coordinates from Point3D object
            pos_obj = rdkit_feature.GetPos()
            pos_coords = try
                [pyconvert(Float64, pos_obj.x), pyconvert(Float64, pos_obj.y), pyconvert(Float64, pos_obj.z)]
            catch
                # If direct access fails, use default position
                [0.0, 0.0, 0.0]
            end

            feature_id = pyconvert(Int, rdkit_feature.GetId())

            push!(features, ChemicalFeature(
                family = family,
                type = feature_type,
                atom_ids = atom_ids,
                position = pos_coords,
                id = feature_id
            ))
        end

        return features
    catch e
        @warn "Error extracting molecular features: $e"
        return ChemicalFeature[]
    end
end

"""
    pharmacophore_fingerprint(mol::Molecule; min_points::Int = 2, max_points::Int = 3, factory::Union{FeatureFactory, Nothing} = nothing) -> Union{Vector{Bool}, Missing}

Generate a pharmacophore fingerprint for a molecule.

# Arguments
- `mol::Molecule`: Input molecule (no conformer coordinates required)
- `min_points::Int = 2`: Minimum number of pharmacophore points
- `max_points::Int = 3`: Maximum number of pharmacophore points
- `factory::Union{FeatureFactory, Nothing} = nothing`: Feature factory (default will be created if not provided)

# Returns
- `Union{Vector{Bool}, Missing}`: Binary pharmacophore fingerprint, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("c1ccccc1O")  # Phenol
fp = pharmacophore_fingerprint(mol)

# Custom parameters
fp = pharmacophore_fingerprint(mol; min_points=2, max_points=4)
```

# Notes
- **No Conformer Required**: Works directly with molecular structure using 2D topological distances
- Combines chemical features with distance information
- Fingerprint length depends on feature combinations and distance bins
- More efficient than 3D approaches as no coordinate generation is needed
"""
function pharmacophore_fingerprint(mol::Molecule; min_points::Int = 2, max_points::Int = 3, factory::Union{FeatureFactory, Nothing} = nothing)
    !mol.valid && return missing

    try
        # Create factory if not provided
        if factory === nothing
            factory = create_feature_factory()
        end

        # Create signature factory
        sig_factory = _create_sig_factory(factory._rdkit_factory; min_point_count = min_points, max_point_count = max_points)

        # Generate fingerprint
        fp = _get_pharmacophore_fingerprint(mol._rdkit_mol, factory._rdkit_factory, sig_factory)

        # Convert IntSparseIntVect to Julia vector
        fp_size = pyconvert(Int, fp.GetLength())
        result = zeros(Bool, fp_size)

        # Get non-zero indices and set them to true
        non_zero_indices = fp.GetNonzeroElements()
        for (idx, val) in non_zero_indices.items()
            if pyconvert(Int, val) > 0
                result[pyconvert(Int, idx) + 1] = true  # Convert to 1-based indexing
            end
        end

        return result
    catch e
        @warn "Error generating pharmacophore fingerprint: $e"
        return missing
    end
end

"""
    get_feature_families(factory::FeatureFactory) -> Vector{String}

Get the list of feature families defined in the feature factory.

# Arguments
- `factory::FeatureFactory`: The feature factory

# Returns
- `Vector{String}`: List of feature family names

# Example
```julia
factory = create_feature_factory()
families = get_feature_families(factory)
println("Available families: ", families)
```
"""
function get_feature_families(factory::FeatureFactory)
    return factory.feature_families
end

"""
    filter_features_by_family(features::Vector{ChemicalFeature}, family::String) -> Vector{ChemicalFeature}

Filter chemical features by their family type.

# Arguments
- `features::Vector{ChemicalFeature}`: List of chemical features (from `get_mol_features`)
- `family::String`: Feature family to filter by (e.g., "Donor", "Acceptor")

# Returns
- `Vector{ChemicalFeature}`: Filtered list of features

# Example
```julia
mol = mol_from_smiles("CCO")
factory = create_feature_factory()

# Generate coordinates for feature extraction
conformers_2d = generate_2d_conformers(mol)
if !isempty(conformers_2d)
    mol_2d = conformers_2d[1].molecule
    features = get_mol_features(mol_2d, factory)

    # Get only hydrogen bond donors
    donors = filter_features_by_family(features, "Donor")
    println("Found ", length(donors), " hydrogen bond donors")
end
```
"""
function filter_features_by_family(features::Vector{ChemicalFeature}, family::String)
    return filter(f -> f.family == family, features)
end

"""
    get_pharmacophore_3d(mol::Molecule, factory::FeatureFactory; conf_id::Int = -1) -> Vector{Tuple{String, Vector{Float64}}}

Extract 3D pharmacophore points from a molecule.

# Arguments
- `mol::Molecule`: Input molecule (must have 3D conformer with coordinates)
- `factory::FeatureFactory`: Feature factory for feature identification
- `conf_id::Int = -1`: Conformer ID to use (-1 for default)

# Returns
- `Vector{Tuple{String, Vector{Float64}}}`: List of (feature_family, 3d_position) tuples

# Example
```julia
mol = mol_from_smiles("CCO")
conformers_3d = generate_3d_conformers(mol, 1)
if !isempty(conformers_3d)
    mol_3d = conformers_3d[1].molecule
    factory = create_feature_factory()
    ph4_points = get_pharmacophore_3d(mol_3d, factory)

    for (family, pos) in ph4_points
        println("  \$(family) at position [\$(pos[1]), \$(pos[2]), \$(pos[3])]")
    end
end
```

# Notes
- **3D Conformer Required**: Molecule must have 3D coordinates. Use `generate_3d_conformers()` first.
- Returns empty vector for molecules without 3D conformers
- Provides spatial positioning of pharmacophore features for 3D analysis
"""
function get_pharmacophore_3d(mol::Molecule, factory::FeatureFactory; conf_id::Int = -1)
    !mol.valid && return Tuple{String, Vector{Float64}}[]

    try
        feature_list = _explicit_pharmacophore_from_mol(mol._rdkit_mol, factory._rdkit_factory; conf_id = conf_id)
        return feature_list
    catch e
        @warn "Error extracting 3D pharmacophore: $e"
        return Tuple{String, Vector{Float64}}[]
    end
end

# Vectorized operations
function get_mol_features(mols::Vector{Molecule}, factory::FeatureFactory; conf_id::Int = -1)
    return [get_mol_features(mol, factory; conf_id = conf_id) for mol in mols]
end

function pharmacophore_fingerprint(mols::Vector{Molecule}; min_points::Int = 2, max_points::Int = 3, factory::Union{FeatureFactory, Nothing} = nothing)
    return [pharmacophore_fingerprint(mol; min_points = min_points, max_points = max_points, factory = factory) for mol in mols]
end

function get_pharmacophore_3d(mols::Vector{Molecule}, factory::FeatureFactory; conf_id::Int = -1)
    return [get_pharmacophore_3d(mol, factory; conf_id = conf_id) for mol in mols]
end