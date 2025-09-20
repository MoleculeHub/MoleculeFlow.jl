```@meta
CurrentModule = MoleculeFlow
```

# Molecular Descriptors

Functions for calculating molecular properties and descriptors.

## Molecular Weight

```@docs
molecular_weight
```

## Exact Molecular Weight

```@docs
exact_molecular_weight
```

## Heavy Atom Count

```@docs
heavy_atom_count
```

## Number of Heteroatoms

```@docs
num_heteroatoms
```

## LogP

```@docs
logp
```

## TPSA

```@docs
tpsa
```

## Number of Hydrogen Bond Donors

```@docs
num_hbd
```

## Number of Hydrogen Bond Acceptors

```@docs
num_hba
```

## Number of Rotatable Bonds

```@docs
num_rotatable_bonds
```

## Number of Rings

```@docs
num_rings
```

## Number of Aromatic Rings

```@docs
num_aromatic_rings
```

## Number of Saturated Rings

```@docs
num_saturated_rings
```

## Bertz CT

```@docs
bertz_ct
```

## Balaban J

```@docs
balaban_j
```

## Chi0v

```@docs
chi0v
```

## Kappa1

```@docs
kappa1
```

## SlogP VSA

```@docs
slogp_vsa
```

## Calculate All Descriptors

```@docs
calc_all_descriptors
```

## Advanced Drug-like and ADMET Descriptors

These descriptors are particularly important for drug discovery and ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) prediction.

### Drug-likeness Scores

```@docs
qed
synthetic_accessibility
```

### Molecular Complexity and 3D Character

```@docs
fraction_csp3
labute_asa
molar_refractivity
```

## Advanced Ring and Structure Counts

More detailed structural analysis beyond basic ring counts.

```@docs
num_aliphatic_carbocycles
num_aromatic_carbocycles
num_aromatic_heterocycles
num_atom_stereo_centers
num_amide_bonds
```

## 3D Descriptors

These descriptors require 3D coordinates to be present in the molecule.

```@docs
asphericity
radius_of_gyration
```

## Examples

### Drug Discovery Workflow

```julia
using MoleculeFlow

# Analyze drug-like properties
compounds = [
    mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O"),  # Aspirin
    mol_from_smiles("CCO"),                         # Ethanol
    mol_from_smiles("c1ccc2c(c1)nnnc2N")           # Drug-like compound
]

for (i, mol) in enumerate(compounds)
    println("Compound $i:")
    println("  QED score: ", qed(mol))
    println("  SAscore: ", synthetic_accessibility(mol))
    println("  Fsp3: ", fraction_csp3(mol))
    println("  Aromatic rings: ", num_aromatic_carbocycles(mol))
    println("  Stereocenters: ", num_atom_stereo_centers(mol))
    println()
end
```

### 3D Shape Analysis

```julia
# Generate 3D conformer and analyze shape
mol = mol_from_smiles("CCO")
conformers = generate_3d_conformers(mol, 1)

if !isempty(conformers)
    mol_3d = conformers[1].molecule

    println("3D Shape Descriptors:")
    println("  Asphericity: ", asphericity(mol_3d))
    println("  Radius of gyration: ", radius_of_gyration(mol_3d))
end
```

### Comprehensive Analysis

```julia
# Calculate multiple advanced descriptors at once
mol = mol_from_smiles("CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F")

println("Comprehensive Descriptor Analysis:")
println("  Drug-likeness (QED): ", qed(mol))
println("  Synthetic accessibility: ", synthetic_accessibility(mol))
println("  3D character (Fsp3): ", fraction_csp3(mol))
println("  Surface area (LabuteASA): ", labute_asa(mol))
println("  Molar refractivity: ", molar_refractivity(mol))
println("  Aromatic carbocycles: ", num_aromatic_carbocycles(mol))
println("  Aromatic heterocycles: ", num_aromatic_heterocycles(mol))
println("  Amide bonds: ", num_amide_bonds(mol))
```

## Additional Molecular Connectivity Descriptors

These descriptors provide detailed information about molecular connectivity and shape.

### Chi Connectivity Indices

```@docs
chi0n
chi1n
chi2n
chi3n
chi4n
chi1v
chi2v
chi3v
chi4v
```

### Kappa Shape Descriptors

```@docs
kappa2
kappa3
```

### E-State Descriptors

```@docs
max_e_state_index
min_e_state_index
```

### Information Content

```@docs
ipc
```

## Atom Counts

Simple but useful atom counting functions for chemical analysis.

```@docs
num_carbons
num_nitrogens
num_oxygens
num_sulfurs
num_halogens
```

## Get Address

```@docs
get_address
```

## Examples

### Molecular Connectivity Analysis

```julia
using MoleculeFlow

# Analyze molecular connectivity for drug-like molecules
compounds = [
    mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O"),  # Aspirin
    mol_from_smiles("CCO"),                         # Ethanol
    mol_from_smiles("c1ccc2c(c1)nnnc2N")           # Triazole compound
]

for (i, mol) in enumerate(compounds)
    println("Compound $i:")
    println("Chi0n: ", chi0n(mol))
    println("Chi1v: ", chi1v(mol))
    println("Kappa2: ", kappa2(mol))
    println("Max E-state: ", max_e_state_index(mol))
    println("Carbons: ", num_carbons(mol))
    println("Nitrogens: ", num_nitrogens(mol))
    println("Halogens: ", num_halogens(mol))
    println()
end
```

### Working with Multiple Molecules

Most descriptor functions in MoleculeFlow can accept vectors of molecules directly for efficient batch processing:

```julia
# Create a vector of molecules
molecules = [mol_from_smiles(smi) for smi in ["CCO", "c1ccccc1", "CC(=O)N", "CCCl"]]

# Descriptor functions accept vectors directly
mw_values = molecular_weight(molecules)
logp_values = logp(molecules)
chi_values = chi1v(molecules)
carbon_counts = num_carbons(molecules)

println("Molecular weights: ", mw_values)
println("LogP values: ", logp_values)
println("Chi1v values: ", chi_values)
println("Carbon counts: ", carbon_counts)
```
