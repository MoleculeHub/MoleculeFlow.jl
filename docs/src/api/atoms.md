```@meta
CurrentModule = MoleculeFlow
```

# Atom Operations

Functions for working with individual atoms in molecules.

## Atom Type

```@docs
Atom
```

## Atom Access

```@docs
get_atoms
get_atom
```

## Basic Properties

```@docs
get_atomic_number
get_symbol
get_mass
get_isotope
```

## Bonding Information

```@docs
get_degree
get_valence
get_formal_charge
get_hybridization
```

## Hydrogen Counts

```@docs
get_num_explicit_hs
get_num_implicit_hs
get_total_num_hs
```

## Structural Properties

```@docs
is_aromatic
is_in_ring
is_in_ring_size
get_chiral_tag
```

## Connectivity

```@docs
get_neighbors
get_bonds_from_atom
```

## Extended Properties

```@docs
get_num_radical_electrons
is_hetero
is_hydrogen_donor
is_hydrogen_acceptor
is_chiral_center
get_cip_code
get_ring_size
```

## Contribution Properties

```@docs
get_crippen_log_p_contribution
get_crippen_molar_refractivity_contribution
get_tpsa_contribution
get_labute_asa_contribution
```

## Comprehensive Atom Data

```@docs
get_all_atom_properties
```

## Partial Charges

```@docs
compute_gasteiger_charges!
get_gasteiger_charge
```