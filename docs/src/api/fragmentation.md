```@meta
CurrentModule = MoleculeFlow
```

# Molecular Fragmentation

Functions for fragmenting molecules using various computational approaches.

## BRICS Fragmentation

BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures) uses retrosynthetic rules to identify breakable bonds in molecules.

## BRICS Decompose

```@docs
brics_decompose
```

## RECAP Fragmentation

RECAP (Retrosynthetic Combinatorial Analysis Procedure) uses different fragmentation rules focused on bonds commonly formed in synthetic chemistry.

## RECAP Decompose

```@docs
recap_decompose
```

## Scaffold Extraction

Scaffold extraction methods for identifying core molecular frameworks.

## Get Murcko Scaffold

```@docs
get_murcko_scaffold
```

## Get Generic Scaffold

```@docs
get_generic_scaffold
```

## Custom Fragmentation

Tools for custom fragmentation strategies and fragment manipulation.

## Fragment by Bonds

```@docs
fragment_by_bonds
```

## Get Fragment Count

```@docs
get_fragment_count
```

## Split Fragments

```@docs
split_fragments
```

## Get Largest Fragment

```@docs
get_largest_fragment
```

## Usage Examples

### Basic Fragmentation

```julia
using MoleculeFlow

# Create a molecule
mol = mol_from_smiles("CCCOCCc1ccccc1")

# BRICS fragmentation
brics_frags = brics_decompose(mol)
println("BRICS fragments: ", brics_frags)

# RECAP fragmentation
recap_frags = recap_decompose(mol)
println("RECAP fragments: ", recap_frags)

# Extract Murcko scaffold
scaffold = get_murcko_scaffold(mol)
println("Murcko scaffold: ", scaffold)
```

### Working with Disconnected Molecules

```julia
# Molecule with disconnected fragments
mol = mol_from_smiles("CCO.CCC.c1ccccc1")

# Count fragments
count = get_fragment_count(mol)
println("Number of fragments: ", count)

# Split into individual molecules
fragments = split_fragments(mol)
for (i, frag) in enumerate(fragments)
    println("Fragment $i: ", mol_to_smiles(frag))
end

# Get the largest fragment
largest = get_largest_fragment(mol)
println("Largest fragment: ", mol_to_smiles(largest))
```

### Advanced Fragmentation

```julia
# Custom bond-based fragmentation
mol = mol_from_smiles("CCCCCCCC")

# Fragment by breaking specific bonds (0-based indexing)
fragments = fragment_by_bonds(mol, [2, 5])
println("Custom fragments: ", fragments)

# Size-constrained BRICS fragmentation
fragments = brics_decompose(mol, min_fragment_size=3, max_fragment_size=8)
println("Size-constrained fragments: ", fragments)
```

### Vectorized Operations

```julia
# Process multiple molecules at once
smiles_list = ["CCCOCCc1ccccc1", "CC(=O)NCCc1ccccc1", "CCCCCC"]
mols = mol_from_smiles(smiles_list)

# Get all BRICS fragments
all_brics = brics_decompose(mols)
for (i, frags) in enumerate(all_brics)
    println("Molecule $i BRICS: ", frags)
end

# Get all scaffolds
all_scaffolds = get_murcko_scaffold(mols)
for (i, scaffold) in enumerate(all_scaffolds)
    println("Molecule $i scaffold: ", scaffold)
end
```
