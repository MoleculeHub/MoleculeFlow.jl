```@meta
CurrentModule = MoleculeFlow
```

# Pharmacophore Features

Functions for pharmacophore analysis, chemical feature identification, and pharmacophore fingerprint generation.

## Overview

Pharmacophores represent the spatial arrangement of chemical features that are essential for molecular recognition. MoleculeFlow provides comprehensive pharmacophore functionality including:

- Chemical feature identification (donors, acceptors, aromatic rings, etc.)
- 3D pharmacophore point extraction
- Pharmacophore fingerprint generation
- Feature-based molecular analysis

## Types

### FeatureFactory

```@docs
FeatureFactory
```

### ChemicalFeature

```@docs
ChemicalFeature
```

## Factory Creation

### Default Feature Factory

```@docs
create_feature_factory
```

## Feature Extraction

### Molecular Features

```@docs
get_mol_features
```

### Feature Families

```@docs
get_feature_families
```

### Feature Filtering

```@docs
filter_features_by_family
```

## Pharmacophore Analysis

### 3D Pharmacophore Points

```@docs
get_pharmacophore_3d
```

### Pharmacophore Fingerprints

```@docs
pharmacophore_fingerprint
```

## Usage Examples

### Basic Feature Identification

```julia
using MoleculeFlow

# Create a feature factory with default definitions
factory = create_feature_factory()

# Analyze a molecule (features require coordinates)
mol = mol_from_smiles("c1ccccc1O")  # Phenol
conformers_2d = generate_2d_conformers(mol)  # Generate 2D coordinates

if !isempty(conformers_2d)
    mol_2d = conformers_2d[1].molecule
    features = get_mol_features(mol_2d, factory)

    println("Found $(length(features)) features:")
    for feature in features
        println("  $(feature.family) at atoms $(feature.atom_ids)")
    end
end
```

### Pharmacophore Fingerprints

```julia
# Generate pharmacophore fingerprints for multiple molecules
molecules = [
    mol_from_smiles("CCO"),           # Ethanol
    mol_from_smiles("c1ccccc1O"),     # Phenol
    mol_from_smiles("CC(=O)N")        # Acetamide
]

# Note: Pharmacophore fingerprints work with molecules as-is
# (they generate 2D topological features)
fingerprints = pharmacophore_fingerprint(molecules)

for (i, fp) in enumerate(fingerprints)
    if fp !== missing
        println("Molecule $(i): $(count(fp)) bits set in fingerprint of length $(length(fp))")
    end
end
```

### 3D Pharmacophore Analysis

```julia
# Analyze 3D pharmacophore points (requires 3D coordinates)
mol = mol_from_smiles("c1ccccc1O")  # Phenol

# Generate 3D conformer
conformers_3d = generate_3d_conformers(mol, 1)
if !isempty(conformers_3d)
    mol_3d = conformers_3d[1].molecule
    factory = create_feature_factory()

    # Extract 3D pharmacophore points
    ph4_points = get_pharmacophore_3d(mol_3d, factory)

    println("3D Pharmacophore points:")
    for (family, position) in ph4_points
        println("$(family) at [$(position[1]), $(position[2]), $(position[3])]")           
    end
end
```

### Feature Family Analysis

```julia
# Analyze specific feature families
mol = mol_from_smiles("CC(=O)Nc1ccccc1O")  # N-acetyl-p-hydroxybenzamide

# Generate 3D conformer for accurate feature positioning
conformers = generate_3d_conformers(mol, 1)
if !isempty(conformers)
    mol_3d = conformers[1].molecule
    factory = create_feature_factory()
    features = get_mol_features(mol_3d, factory)

    # Get all hydrogen bond donors
    donors = filter_features_by_family(features, "Donor")
    println("Hydrogen bond donors: $(length(donors))")

    # Get all hydrogen bond acceptors
    acceptors = filter_features_by_family(features, "Acceptor")
    println("Hydrogen bond acceptors: $(length(acceptors))")

    # Get all aromatic features
    aromatics = filter_features_by_family(features, "Aromatic")
    println("Aromatic rings: $(length(aromatics))")
end
```

### Custom Feature Definitions

```julia
# Define custom pharmacophore features
custom_fdef = raw"""
DefineFeature Carboxyl [CX3](=O)[OX2H1]
   Family AcidicGroup
   Weights 1.0,1.0,1.0
EndFeature

DefineFeature Amine [NX3;H2,H1;!$(NC=O)]
   Family BasicGroup
   Weights 1.0
EndFeature
"""

# Create factory with custom definitions
custom_factory = create_feature_factory(
    use_default=false,
    fdef_string=custom_fdef
)

# Available families in custom factory
families = get_feature_families(custom_factory)
println("Custom feature families: $(families)")
```

### Drug-like Molecule Analysis

```julia
# Analyze pharmacophore features of drug molecules
aspirin = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O")
ibuprofen = mol_from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")

factory = create_feature_factory()

# Compare pharmacophore fingerprints
aspirin_fp = pharmacophore_fingerprint(aspirin)
ibuprofen_fp = pharmacophore_fingerprint(ibuprofen)

if aspirin_fp isa Vector{Bool} && ibuprofen_fp isa Vector{Bool}
    # Calculate Tanimoto similarity
    similarity = tanimoto_similarity(aspirin_fp, ibuprofen_fp)
    println("Pharmacophore similarity: $(similarity)")
end

# Analyze individual features (need coordinates for feature extraction)
aspirin_2d = generate_2d_conformers(aspirin)
ibuprofen_2d = generate_2d_conformers(ibuprofen)

if !isempty(aspirin_2d) && !isempty(ibuprofen_2d)
    aspirin_features = get_mol_features(aspirin_2d[1].molecule, factory)
    ibuprofen_features = get_mol_features(ibuprofen_2d[1].molecule, factory)

    println("Aspirin features:")
    for feature in aspirin_features
        println("  $(feature.family)")
    end

    println("Ibuprofen features:")
    for feature in ibuprofen_features
        println("  $(feature.family)")
    end
end
```

## Available Feature Families

The default feature factory includes the following pharmacophore feature families:

- **Donor**: Hydrogen bond donors
- **Acceptor**: Hydrogen bond acceptors
- **NegIonizable**: Negatively ionizable groups
- **PosIonizable**: Positively ionizable groups
- **ZnBinder**: Zinc binding groups
- **Aromatic**: Aromatic ring systems
- **Hydrophobe**: Hydrophobic regions
- **LumpedHydrophobe**: Larger hydrophobic regions

## Notes

- **Coordinate Requirements**: Feature extraction (`get_mol_features`) requires molecules to have 2D or 3D coordinates. Use `generate_2d_conformers()` or `generate_3d_conformers()` before feature extraction.
- **Pharmacophore Fingerprints**: Work directly with molecules (no coordinates required) and use 2D topological distances by default
- **3D Pharmacophore Analysis**: Requires molecules with 3D coordinates for accurate spatial positioning
- **Feature Identification**: Uses SMARTS patterns defined in the feature factory
- **Custom Definitions**: Custom feature definitions can be created using the RDKit feature definition language
- **Vectorized Operations**: All functions support vectorized operations for processing multiple molecules efficiently
