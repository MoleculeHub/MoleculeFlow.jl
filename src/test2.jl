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


ethanol_smiles = mol_to_smiles(ethanol)
println("Ethanol SMILES: $ethanol_smiles")

# Convert multiple molecules
output_smiles = mol_to_smiles(molecules)
println("Output SMILES: $output_smiles")

# Generate InChI representation
ethanol_inchi = mol_to_inchi(ethanol)
println("Ethanol InChI: $ethanol_inchi")

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

charged_mol = mol_from_smiles("[NH3+]CCO")  # Protonated ethylamine
compute_gasteiger_charges!(charged_mol)

atoms = get_atoms(charged_mol)
println("Gasteiger partial charges:")
for (i, atom) in enumerate(atoms)
    symbol = get_symbol(atom)
    charge = get_gasteiger_charge(atom)
    println("Atom $i ($symbol): $charge")
end


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