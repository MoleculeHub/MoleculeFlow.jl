using Pkg
Pkg.activate(".")
using MoleculeFlow

# Generate conformers for a flexible molecule
flexible_mol = mol_from_smiles("C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N") # Ibrutinib

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