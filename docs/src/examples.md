# Practical Examples

Here are comprehensive practical examples of using MoleculeFlow.jl for molecular analysis, manipulation, and visualization.

## Table of Contents

1. [Basic Molecule Management](#basic-molecule-management)
2. [Working with Atoms](#working-with-atoms)
3. [Working with Bonds](#working-with-bonds)
4. [Molecular Properties and Descriptors](#molecular-properties-and-descriptors)
5. [Substructure Analysis](#substructure-analysis)
6. [Similarity and Fingerprints](#similarity-and-fingerprints)
7. [Molecular Visualization](#molecular-visualization)
8. [Dataset Analysis](#dataset-analysis)
9. [Molecular Standardization](#molecular-standardization)
10. [3D Conformer Generation](#3d-conformer-generation)

## Basic Molecule Management

### Creating and Validating Molecules

```julia
using MoleculeFlow

# Create molecules from SMILES
ethanol = mol_from_smiles("CCO")
benzene = mol_from_smiles("c1ccccc1")
invalid_mol = mol_from_smiles("invalid_smiles")

# Check validity
println("Ethanol valid: $(ethanol.valid)")
println("Benzene valid: $(benzene.valid)")
println("Invalid mol valid: $(invalid_mol.valid)")

# Access molecule source (original SMILES)
println("Ethanol source: $(ethanol.source)")

# Batch creation from SMILES list
smiles_list = ["CCO", "CCC", "c1ccccc1", "CC(=O)O"]
molecules = mol_from_smiles(smiles_list)
valid_molecules = filter(mol -> mol.valid, molecules)
```

### Converting Molecules Back to SMILES

```julia
# Convert single molecule
ethanol_smiles = mol_to_smiles(ethanol)
println("Ethanol SMILES: $ethanol_smiles")

# Convert multiple molecules
output_smiles = mol_to_smiles(molecules)
println("Output SMILES: $output_smiles")

# Generate InChI representation
ethanol_inchi = mol_to_inchi(ethanol)
println("Ethanol InChI: $ethanol_inchi")
```

### Reading from SDF Files

```julia
# Read all molecules from an SDF file
# Note: This example assumes you have an SDF file
try
    molecules = read_sdf("compounds.sdf")
    println("Read $(length(molecules)) molecules from SDF")

    # Filter valid molecules
    valid_molecules = filter(mol -> mol.valid, molecules)
    println("$(length(valid_molecules)) valid molecules")

    # Access SDF properties
    for mol in valid_molecules[1:min(3, length(valid_molecules))]
        println("Molecule: $(mol_to_smiles(mol))")
        println("  SDF Index: $(mol.props[:sdf_index])")
        println("  Source File: $(mol.props[:source_file])")

        # Display other SDF properties
        for (key, value) in mol.props
            if key ∉ [:sdf_index, :source_file]
                println("  $key: $value")
            end
        end
        println()
    end

catch e
    println("SDF file not found or error reading: $e")
end

# Read only first 100 molecules for large files
try
    molecules = read_sdf("large_database.sdf", max_mols=100)
    println("Read first 100 molecules from large SDF")
catch e
    println("Large SDF file not found: $e")
end

# Memory-efficient reading for huge files
try
    next_mol = read_sdf_lazy("huge_database.sdf")
    mol_count = 0

    while mol_count < 10  # Process first 10 molecules
        mol = next_mol()
        if mol === nothing
            break  # End of file
        end

        mol_count += 1
        if mol.valid
            smiles = mol_to_smiles(mol)
            println("Molecule $mol_count: $smiles")
        end
    end
catch e
    println("Huge SDF file not found: $e")
end
```

### Reading MOL Blocks

```julia
# Example MOL block (ethane)
molblock = """

  RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""

mol = mol_from_molblock(molblock)
if mol.valid
    println("Successfully parsed MOL block")
    println("SMILES: $(mol_to_smiles(mol))")
else
    println("Failed to parse MOL block")
end
```

## Working with Atoms

### Basic Atom Access and Properties

```julia
# Get all atoms from a molecule
caffeine = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
atoms = get_atoms(caffeine)
println("Caffeine has $(length(atoms)) atoms")

# Get specific atom by index (1-based)
first_atom = get_atom(caffeine, 1)
println("First atom symbol: $(get_symbol(first_atom))")
println("First atom atomic number: $(get_atomic_number(first_atom))")

# Iterate over all atoms
println("\nAll atoms in caffeine:")
for (i, atom) in enumerate(atoms)
    symbol = get_symbol(atom)
    atomic_num = get_atomic_number(atom)
    degree = get_degree(atom)
    println("Atom $i: $symbol (Z=$atomic_num, degree=$degree)")
end
```

### Advanced Atom Properties

```julia
# Get detailed atom properties
for (i, atom) in enumerate(atoms[1:5])  # First 5 atoms
    symbol = get_symbol(atom)
    degree = get_degree(atom)
    valence = get_valence(atom)
    formal_charge = get_formal_charge(atom)
    hybridization = get_hybridization(atom)
    is_aromatic_atom = is_aromatic(atom)
    is_in_ring_atom = is_in_ring(atom)

    println("Atom $i ($symbol):")
    println("Degree: $degree, Valence: $valence")
    println("Formal charge: $formal_charge")
    println("Hybridization: $hybridization")
    println("Aromatic: $is_aromatic_atom, In ring: $is_in_ring_atom")
    println()
end
```

### Atom Neighbors and Connectivity

```julia
# Get neighbors of specific atoms
println("Analyzing atom connectivity:")
for i in 1:min(5, length(atoms))
    neighbors = get_neighbors(caffeine, i)
    neighbor_symbols = [get_symbol(get_atom(caffeine, j)) for j in neighbors]

    atom_symbol = get_symbol(get_atom(caffeine, i))
    println("Atom $i ($atom_symbol) connected to: $(join(neighbor_symbols, ", "))")
end

# Get bonds involving specific atoms
println("\nBonds from first atom:")
bonds_from_atom1 = get_bonds_from_atom(caffeine, 1)
for (i, bond) in enumerate(bonds_from_atom1)
    bond_type = get_bond_type(bond)
    begin_atom = get_begin_atom_idx(bond)
    end_atom = get_end_atom_idx(bond)
    println("Bond $i: Atom $begin_atom - Atom $end_atom ($bond_type)")
end
```

### Gasteiger Partial Charges

```julia
charged_mol = mol_from_smiles("[NH3+]CCO")  # Protonated ethylamine
compute_gasteiger_charges!(charged_mol)

atoms = get_atoms(charged_mol)
println("Gasteiger partial charges:")
for (i, atom) in enumerate(atoms)
    symbol = get_symbol(atom)
    charge = get_gasteiger_charge(atom)
    println("Atom $i ($symbol): $charge")
end
```

## Working with Bonds

### Basic Bond Analysis

```julia
molecule = mol_from_smiles("CC=C(C)C#N")

println("Bond analysis for CC=C(C)C#N:")
atoms = get_atoms(molecule)
for i in eachindex(atoms)
    bonds = get_bonds_from_atom(molecule, i)
    if !isempty(bonds)
        atom_symbol = get_symbol(get_atom(molecule, i))
        println("\nAtom $i ($atom_symbol) bonds:")

        for bond in bonds
            bond_type = get_bond_type(bond)
            begin_idx = get_begin_atom_idx(bond)
            end_idx = get_end_atom_idx(bond)
            is_aromatic_bond = is_aromatic(bond)
            is_in_ring_bond = is_in_ring(bond)

            other_atom_idx = begin_idx == i ? end_idx : begin_idx
            other_symbol = get_symbol(get_atom(molecule, other_atom_idx))

            println(" to Atom $other_atom_idx ($other_symbol): $bond_type" *
                   (is_aromatic_bond ? " (aromatic)" : "") *
                   (is_in_ring_bond ? " (in ring)" : ""))
        end
    end
end
```

### Bond Type Classification

```julia
complex_mol = mol_from_smiles("c1ccc2c(c1)ccc(=O)c2=O")  # Naphthoquinone

single_bonds = 0
double_bonds = 0
aromatic_bonds = 0
ring_bonds = 0

atoms = get_atoms(complex_mol)
processed_bonds = Set()  # To avoid counting bonds twice

for i in eachindex(atoms)
    bonds = get_bonds_from_atom(complex_mol, i)
    for bond in bonds
        begin_idx = get_begin_atom_idx(bond)
        end_idx = get_end_atom_idx(bond)
        bond_key = (min(begin_idx, end_idx), max(begin_idx, end_idx))

        if bond_key ∉ processed_bonds
            push!(processed_bonds, bond_key)

            bond_type = get_bond_type(bond)
            if bond_type == "SINGLE"
                single_bonds += 1
            elseif bond_type == "DOUBLE"
                double_bonds += 1
            elseif bond_type == "AROMATIC"
                aromatic_bonds += 1
            end

            if is_in_ring(bond)
                ring_bonds += 1
            end
        end
    end
end

println("Bond analysis for naphthoquinone:")
println("Single bonds: $single_bonds")
println("Double bonds: $double_bonds")
println("Aromatic bonds: $aromatic_bonds")
println("Ring bonds: $ring_bonds")
```

## Molecular Properties and Descriptors

### Comprehensive Property Analysis

```julia
mol = mol_from_smiles("CCO")

println("Molecular Properties Analysis:")
mw = molecular_weight(mol)
logp_val = logp(mol)
tpsa_val = tpsa(mol)
hbd = num_hbd(mol)
hba = num_hba(mol)
rot_bonds = num_rotatable_bonds(mol)
rings = num_rings(mol)
heavy_atoms = heavy_atom_count(mol)
```

### Lipinski's Rule of Five Analysis

```julia
function check_lipinski(mol::Molecule)
    """Check if molecule passes Lipinski's Rule of Five"""
    if !mol.valid
        return false, "Invalid molecule"
    end

    mw = molecular_weight(mol)
    logp_val = logp(mol)
    hbd = num_hbd(mol)
    hba = num_hba(mol)

    violations = []

    if mw > 500
        push!(violations, "MW > 500 ($(round(mw, digits=1)))")
    end
    if logp_val > 5
        push!(violations, "LogP > 5 ($(round(logp_val, digits=1)))")
    end
    if hbd > 5
        push!(violations, "HBD > 5 ($hbd)")
    end
    if hba > 10
        push!(violations, "HBA > 10 ($hba)")
    end

    passes = isempty(violations)
    message = passes ? "Passes Lipinski's Rule of Five" : "Violations: " * join(violations, ", ")

    return passes, message
end

# Test Lipinski's rule on drug-like molecules
drug_smiles = [
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("Atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2F)c(-c2ccccc2)c(-c2ccc(F)cc2)n1C[C@H](O)C[C@H](O)CC(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
]

println("\nLipinski's Rule of Five Analysis:")
for (name, smiles) in drug_smiles
    mol = mol_from_smiles(smiles)
    passes, message = check_lipinski(mol)
    status = passes ? "PASS" : "FAIL"
    println("$name: $status")
    println("$message")
    println()
end
```

## Substructure Analysis

### Finding Functional Groups

```julia
mol = mol_from_smiles("C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N") # Ibrutinib

functional_groups_to_check = [
    :alcohol, :carboxylic_acid, :amide, :amine_primary,
    :benzene, :phenol, :ketone, :aldehyde
]

println("Functional Group Analysis:")
for fg in functional_groups_to_check
    has_fg = has_functional_group(mol, fg)
    print(fg, ": " ,has_fg, "\n")
end
```

### Custom Substructure Searches

```julia
# Define custom SMARTS patterns
custom_patterns = [
    ("Hydroxyl", "[OH]"),
    ("Carbonyl", "C=O"),
    ("Aromatic Ring", "c1ccccc1"),
    ("Primary Amine", "[NH2]"),
    ("Ester", "[CX3](=O)[OX2H0]"),
    ("Benzyl", "c1ccccc1[CH2]"),
    ("Trifluoromethyl", "C(F)(F)F")
]

test_mol = mol_from_smiles("CC(=O)Oc1ccc(cc1)C(=O)O")  # Aspirin

println("Custom substructure analysis for Aspirin:")
for (name, pattern) in custom_patterns
    matches = get_substructure_matches(test_mol, pattern)
    if !isempty(matches)
        println("$name: Found $(length(matches)) match(es)")
        for (i, match) in enumerate(matches)
            println("  Match $i: atoms $(join(match, ", "))")
        end
    else
        println("$name: Not found")
    end
end
```

## Similarity and Fingerprints

### Comprehensive Similarity Analysis

```julia
# Calculate multiple similarity metrics
mol1 = mol_from_smiles("CCO")  # Ethanol
mol2 = mol_from_smiles("CCC")  # Propane
mol3 = mol_from_smiles("CCCO") # Propanol

molecules = [("Ethanol", mol1), ("Propane", mol2), ("Propanol", mol3)]

# Calculate different fingerprint types
fingerprint_types = [:morgan, :rdk, :maccs]

println("Similarity Analysis:")
for fp_type in fingerprint_types
    println("\n$fp_type Fingerprints:")

    # Calculate pairwise similarities
    for i in eachindex(molecules)
        for j in (i+1):length(molecules)
            name1, mol_i = molecules[i]
            name2, mol_j = molecules[j]

            # Calculate different similarity metrics
            tanimoto = tanimoto_similarity(mol_i, mol_j, fingerprint_type=fp_type)
            dice = dice_similarity(mol_i, mol_j, fingerprint_type=fp_type)
            cosine = cosine_similarity(mol_i, mol_j, fingerprint_type=fp_type)

            println("  $name1 vs $name2:")
            println("    Tanimoto: $(round(tanimoto, digits=3))")
            println("    Dice: $(round(dice, digits=3))")
            println("    Cosine: $(round(cosine, digits=3))")
        end
    end
end
```

### Bulk Similarity and Matrix Calculations

```julia
# Create a dataset of similar molecules
alcohol_smiles = ["CCO", "CCCO", "CC(C)O", "CCCCO", "CC(C)CO", "c1ccc(O)cc1"]
alcohol_mols = mol_from_smiles(alcohol_smiles)
valid_alcohols = filter(mol -> mol.valid, alcohol_mols)

# Query molecule
query = mol_from_smiles("CCO")

# Bulk similarity calculation
similarities = bulk_similarity(query, valid_alcohols)
println("Similarities to ethanol:")
for (i, sim) in enumerate(similarities)
    smiles = mol_to_smiles(valid_alcohols[i])
    println("  $smiles: $(round(sim, digits=3))")
end

# Similarity matrix
sim_matrix = similarity_matrix(valid_alcohols[1:4])  # First 4 molecules
println("\nSimilarity Matrix (4x4):")
for i in 1:4
    row = [round(sim_matrix[i, j], digits=3) for j in 1:4]
    println("  $(join(row, "  "))")
end
```

## Molecular Visualization

### Basic Visualization Examples

```julia
# Create various molecular images
molecules_to_draw = [
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Morphine", "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@]3([C@H]1C5)")
]

for (name, smiles) in molecules_to_draw
    mol = mol_from_smiles(smiles)
    if mol.valid
        # Generate basic image
        img = mol_to_image(mol, size=(400, 400))

        # Generate SVG
        svg = mol_to_svg(mol, size=(400, 400))

        # Save images (if desired)
        # save_molecule_image(mol, "$name.svg")

        println("Generated images for $name")
    end
end
```

### Advanced Visualization with Highlighting

```julia
# Highlight functional groups
aspirin = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O")

# Find and highlight ester group
ester_matches = get_substructure_matches(aspirin, "[CX3](=O)[OX2H0]")
if !isempty(ester_matches)
    highlighted_img = mol_to_image(aspirin,
                                  highlight_atoms=ester_matches[1],
                                  size=(500, 500))
    println("Highlighted ester group in aspirin")
end

# Create grid of molecules with legends
molecules = [mol_from_smiles(smiles) for (_, smiles) in molecules_to_draw]
names = [name for (name, _) in molecules_to_draw]

grid_img = mols_to_grid_image(molecules,
                             legends=names,
                             mols_per_row=2,
                             sub_img_size=(300, 300))
println("Created molecular grid")

# Functional group highlighting
acetaminophen = mol_from_smiles("CC(=O)Nc1ccc(O)cc1")
fg_img = draw_functional_groups(acetaminophen,
                               functional_groups=["[OH]", "C=O", "[NH]"],
                               colors=["red", "blue", "green"])
println("Highlighted functional groups in acetaminophen")
```

## Dataset Analysis

### Processing Molecular Datasets

```julia
# Simulate a small molecular dataset
dataset_smiles = [
    "CCO", "CCC", "CCCO", "CC(C)O", "CCCCO",
    "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1",
    "CC(=O)O", "CC(=O)N", "CC(=O)c1ccccc1",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
]

# Create molecules and filter valid ones
molecules = mol_from_smiles(dataset_smiles)
valid_molecules = filter(mol -> mol.valid, molecules)

println("Dataset Analysis:")
println("Total molecules: $(length(dataset_smiles))")
println("Valid molecules: $(length(valid_molecules))")

# Calculate properties for the dataset
properties = []
for mol in valid_molecules
    props = (
        mw = molecular_weight(mol),
        logp = logp(mol),
        tpsa = tpsa(mol),
        hbd = num_hbd(mol),
        hba = num_hba(mol),
        rings = num_rings(mol),
        smiles = mol_to_smiles(mol)
    )
    push!(properties, props)
end

mws = [p.mw for p in properties]
logps = [p.logp for p in properties]
```

### Clustering by Similarity

```julia
# Simple clustering based on Tanimoto similarity
function simple_clustering(molecules, threshold=0.7)
    clusters = []
    assigned = falses(length(molecules))

    for i in eachindex(molecules)
        if assigned[i]
            continue
        end

        cluster = [i]
        assigned[i] = true

        for j in (i+1):length(molecules)
            if assigned[j]
                continue
            end

            sim = tanimoto_similarity(molecules[i], molecules[j])
            if sim >= threshold
                push!(cluster, j)
                assigned[j] = true
            end
        end

        push!(clusters, cluster)
    end

    return clusters
end

# Apply clustering
clusters = simple_clustering(valid_molecules, 0.5)
println("\nClustering Results (threshold=0.5):")
for (i, cluster) in enumerate(clusters)
    println("Cluster $i ($(length(cluster)) molecules):")
    for idx in cluster
        smiles = mol_to_smiles(valid_molecules[idx])
        println("  $smiles")
    end
    println()
end
```

## Molecular Standardization

### Comprehensive Standardization Workflow

```julia
messy_molecules = [
    ("Salt mixture", "CCO.Cl"),
    ("Charged molecule", "CC(=O)[O-].[Na+]"),
    ("Tautomer", "CC(O)=CC(=O)C"),
    ("Stereochemistry", "C[C@H](O)C")
]

println("Molecular Standardization Examples:")
for (name, smiles) in messy_molecules
    mol = mol_from_smiles(smiles)
    if mol.valid
        println("\n$name: $smiles")

        # Apply different standardization steps
        stripped = strip_salts(mol)
        canonical = canonical_tautomer(mol)
        neutralized = neutralize_charges(mol)
        normalized = normalize_molecule(mol)

        # Full standardization
        standardized = standardize_molecule(mol)

        println("Original:      $(mol_to_smiles(mol))")
        println("Salt stripped: $(mol_to_smiles(stripped))")
        println("Canonical:     $(mol_to_smiles(canonical))")
        println("Neutralized:   $(mol_to_smiles(neutralized))")
        println("Normalized:    $(mol_to_smiles(normalized))")
        println("Standardized:  $(mol_to_smiles(standardized))")
    end
end

# Batch standardization
batch_smiles = ["CCO.Cl", "CC(=O)[O-].[Na+]", "CC(O)=CC(=O)C"]
batch_mols = mol_from_smiles(batch_smiles)
standardized_mols = [standardize_molecule(mol) for mol in batch_mols]

println("\nBatch Standardization:")
for (i, (orig, std)) in enumerate(zip(batch_mols, standardized_mols))
    println("  $i: $(mol_to_smiles(orig)) → $(mol_to_smiles(std))")
end
```

## 3D Conformer Generation

### Generating and Analyzing Conformers

```julia
# Generate conformers for a flexible molecule
flexible_mol = mol_from_smiles("CCCCCCCC")  # Octane

# Generate multiple conformers
conformers = generate_3d_conformers(flexible_mol, 10, optimize=true)

if !isempty(conformers)
    println("Generated $(length(conformers)) conformers for octane:")

    for (i, conf_mol) in enumerate(conformers[1:min(5, length(conformers))])  # Show first 5 or fewer
        conf_result = conf_mol.conformer_result
        energy = conf_result.energy
        converged = conf_result.converged ? "✓" : "✗"

        println("  Conformer $i: $(round(energy, digits=2)) kcal/mol (converged: $converged)")

        # Access 3D coordinates
        coords = conf_mol.molecule.props[:coordinates_3d]
        println("    Coordinates shape: $(size(coords))")
    end

    # Find lowest energy conformer
    best_conf = conformers[1]  # Already sorted by energy
    println("\nBest conformer energy: $(round(best_conf.conformer_result.energy, digits=2)) kcal/mol")
end

# 2D coordinate generation
mol_2d = mol_from_smiles("c1ccc2c(c1)ccc(=O)c2=O")  # Naphthoquinone
conformers_2d = generate_2d_conformers(mol_2d)

if !isempty(conformers_2d)
    coords_2d = conformers_2d[1].molecule.props[:coordinates_2d]
    println("\n2D coordinates for naphthoquinone: $(size(coords_2d))")
end
```

### Conformer Energy Analysis

```julia
# Analyze conformer energies and geometry
if !isempty(conformers)
    energies = [conf.conformer_result.energy for conf in conformers]

    println("\nConformer Energy Analysis:")
    println("  Number of conformers: $(length(energies))")
    println("  Energy range: $(round(minimum(energies), digits=2)) - $(round(maximum(energies), digits=2)) kcal/mol")
    println("  Energy span: $(round(maximum(energies) - minimum(energies), digits=2)) kcal/mol")

    # Count converged conformers
    converged_count = sum(conf.conformer_result.converged for conf in conformers)
    println("  Converged conformers: $converged_count/$(length(conformers))")

    # Energy distribution
    low_energy = sum(e < (minimum(energies) + 2.0) for e in energies)
    println("  Low energy conformers (< 2 kcal/mol above minimum): $low_energy")
end
```