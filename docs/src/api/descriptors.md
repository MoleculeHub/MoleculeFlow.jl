```@meta
CurrentModule = MoleculeFlow
```

# Molecular Descriptors

Functions for calculating molecular properties and descriptors.

## Basic Properties

```@docs
molecular_weight
exact_molecular_weight
heavy_atom_count
num_heteroatoms
```

## Lipinski Descriptors

```@docs
logp
tpsa
num_hbd
num_hba
num_rotatable_bonds
```

## Ring Information

```@docs
num_rings
num_aromatic_rings
num_saturated_rings
```

## Advanced Descriptors

```@docs
bertz_ct
balaban_j
chi0v
kappa1
slogp_vsa
```

## All Descriptors at Once

```@docs
calc_all_descriptors
```

## Utility Functions

```@docs
get_address
```