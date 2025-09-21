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

## Quantitative Estimate of Drug-likeness (QED)

```@docs
qed
```

## Synthetic Accessibility

```@docs
synthetic_accessibility
```

## Molecular Complexity and 3D Character

```@docs
fraction_csp3
labute_asa
molar_refractivity
```

## Number of Aliphatic Carbocycles

```@docs
num_aliphatic_carbocycles
```

## Number of Aromatic Carbocycles

```@docs
num_aromatic_carbocycles
```

## Number of Aromatic Heterocycles

```@docs
num_aromatic_heterocycles
```

## Number of Stereo Centres

```@docs
num_amide_bonds
```

## Number of Amide Bonds

```@docs
num_amide_bonds
```

## Additional Ring Type Counts

```@docs
num_aliphatic_heterocycles
```

```@docs
num_saturated_heterocycles
```

```@docs
num_saturated_carbocycles
```

```@docs
num_aliphatic_rings
```

```@docs
num_heterocycles
```

## Additional Stereochemistry and Structure Counts

```@docs
num_unspecified_atom_stereo_centers
```

```@docs
num_spiro_atoms
```

```@docs
num_bridgehead_atoms
```

## Molecular Complexity

```@docs
hall_kier_alpha
```

## 3D Descriptors

These descriptors require 3D coordinates to be present in the molecule.

```@docs
asphericity
radius_of_gyration
eccentricity
inertial_shape_factor
```

## Principal Moments of Inertia

```@docs
pmi1
```

```@docs
pmi2
```

```@docs
pmi3
```

## Chi Connectivity Indices

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

## Kappa Shape Descriptors

```@docs
kappa2
kappa3
```

## E-State Descriptors

```@docs
max_e_state_index
min_e_state_index
```

```@docs
max_absolute_e_state_index
```

```@docs
min_absolute_e_state_index
```

## Information Content

```@docs
ipc
```

## SlogP_VSA Descriptors

Surface area contributions based on SlogP values.

```@docs
slogp_vsa2
```

```@docs
slogp_vsa3
```

```@docs
slogp_vsa4
```

```@docs
slogp_vsa5
```

```@docs
slogp_vsa6
```

```@docs
slogp_vsa7
```

```@docs
slogp_vsa8
```

```@docs
slogp_vsa9
```

```@docs
slogp_vsa10
```

```@docs
slogp_vsa11
```

```@docs
slogp_vsa12
```

## SMR_VSA Descriptors

Surface area contributions based on Molar Refractivity values.

```@docs
smr_vsa1
```

```@docs
smr_vsa2
```

```@docs
smr_vsa3
```

```@docs
smr_vsa4
```

```@docs
smr_vsa5
```

```@docs
smr_vsa6
```

```@docs
smr_vsa7
```

```@docs
smr_vsa8
```

```@docs
smr_vsa9
```

```@docs
smr_vsa10
```

## PEOE_VSA Descriptors

Surface area contributions based on Partial Equalization of Orbital Electronegativities (PEOE) charges.

```@docs
peoe_vsa1
```

```@docs
peoe_vsa2
```

```@docs
peoe_vsa3
```

```@docs
peoe_vsa4
```

```@docs
peoe_vsa5
```

```@docs
peoe_vsa6
```

```@docs
peoe_vsa7
```

```@docs
peoe_vsa8
```

```@docs
peoe_vsa9
```

```@docs
peoe_vsa10
```

```@docs
peoe_vsa11
```

```@docs
peoe_vsa12
```

```@docs
peoe_vsa13
```

```@docs
peoe_vsa14
```

## BCUT2D Molecular Weight

```@docs
bcut2d_mwlow
```

```@docs
bcut2d_mwhi
```

## BCUT2D Partial Charge

```@docs
bcut2d_chglow
```

```@docs
bcut2d_chghi
```

## BCUT2D LogP

```@docs
bcut2d_logplow
```

```@docs
bcut2d_logphi
```

## BCUT2D Molar Refractivity

```@docs
bcut2d_mrlow
```

```@docs
bcut2d_mrhi
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

## Working with Multiple Molecules

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

