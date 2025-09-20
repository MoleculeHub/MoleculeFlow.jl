```@meta
CurrentModule = MoleculeFlow
```

# Molecular Operations

Advanced molecular manipulation, editing, and analysis functions.

## Hydrogen Manipulation

Functions for adding and removing explicit hydrogens from molecules.

### Adding Hydrogens

```@docs
add_hs
```

### Removing Hydrogens

```@docs
remove_hs
```

## Molecular Editing

Functions for modifying molecular structures programmatically.

### Combining Molecules

```@docs
combine_mols
```

### Substructure Modifications

```@docs
delete_substructs
replace_substructs
```

## Stereochemistry Operations

Functions for analyzing and manipulating molecular stereochemistry.

### Stereochemistry Assignment

```@docs
assign_stereochemistry!
```

### Chiral Center Analysis

```@docs
find_chiral_centers
```

## Ring Analysis

Functions for analyzing ring systems and molecular topology.

### Ring Detection

```@docs
fast_find_rings!
```

### Atom Ranking

```@docs
canonical_atom_ranks
```

## Pattern Matching

Advanced pattern matching and substructure analysis.

### SMARTS Matching

```@docs
quick_smarts_match
```

### Fragment Analysis

```@docs
mol_fragment_to_smarts
```

## Examples

### Working with Hydrogens

```julia
using MoleculeFlow

# Start with a simple molecule
mol = mol_from_smiles("CCO")
println("Original: ", mol_to_smiles(mol))

# Add explicit hydrogens
mol_with_hs = add_hs(mol)
println("With H's: ", mol_to_smiles(mol_with_hs))

# Remove explicit hydrogens
mol_clean = remove_hs(mol_with_hs)
println("Clean: ", mol_to_smiles(mol_clean))
```

### Molecular Editing

```julia
# Combine two molecules
ethanol = mol_from_smiles("CCO")
propane = mol_from_smiles("CCC")
combined = combine_mols(ethanol, propane)

# Remove alcohol groups
benzyl_alcohol = mol_from_smiles("c1ccc(CO)cc1")
benzene_derivative = delete_substructs(benzyl_alcohol, "[OH]")
```

### Pattern Matching

```julia
# Check for functional groups
mol = mol_from_smiles("CCO")

has_alcohol = quick_smarts_match(mol, "[OH]")  # true
has_amine = quick_smarts_match(mol, "[NH2]")   # false
has_carbon = quick_smarts_match(mol, "[C]")    # true
```

### Stereochemistry Analysis

```julia
# Work with chiral molecules
chiral_mol = mol_from_smiles("C[C@H](O)C")

# Assign stereochemistry
assign_stereochemistry!(chiral_mol)

# Find chiral centers
centers = find_chiral_centers(chiral_mol)
for (atom_idx, chirality) in centers
    println("Atom $atom_idx: $chirality")
end
```

### Advanced I/O

```julia
# Generate InChI keys for database storage
mol = mol_from_smiles("CCO")
inchi_key = mol_to_inchi_key(mol)
println("InChI Key: ", inchi_key)

# Export as MOL block
molblock = mol_to_molblock(mol)
println("MOL Block:")
println(molblock)
```

### Working with 3D Coordinates

```julia
# Load molecule from XYZ file with 3D coordinates
mol_3d = mol_from_xyz_file("molecule.xyz")

# Or create XYZ format from existing molecule with 3D coordinates
mol = mol_from_smiles("CCO")
conformers = generate_3d_conformers(mol, 1)
if !isempty(conformers)
    mol_with_coords = conformers[1].molecule
    xyz_string = mol_to_xyz_block(mol_with_coords)
    println("XYZ Format:")
    println(xyz_string)
end
```

### MOL2 Format for Pharmacophore Modeling

```julia
# Load ligand from MOL2 file (preserves charges and atom types)
ligand = mol_from_mol2_file("drug_compound.mol2")

if ligand.valid
    println("Loaded ligand: ", mol_to_smiles(ligand))
    println("Heavy atoms: ", heavy_atom_count(ligand))

    # MOL2 format preserves partial charges for QSAR studies
    atoms = get_atoms(ligand)
    for atom in atoms[1:3]  # Show first 3 atoms
        println("Atom ", get_symbol(atom), " charge: ", get_formal_charge(atom))
    end
end
```

## Notes

- Many functions modify molecules in-place (indicated by `!` suffix) for performance
- Functions return `false` or empty results when operations fail
- Invalid molecules are handled gracefully without throwing exceptions
- All new functions integrate seamlessly with existing MoleculeFlow.jl workflows