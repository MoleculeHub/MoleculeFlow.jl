using Test
using MoleculeFlow

@testset "Molecular Descriptors" begin
    mol = mol_from_smiles("CC(C)O")

    @test molecular_weight(mol) ≈ 60.096 atol = 0.01
    @test heavy_atom_count(mol) == 4
    @test num_hbd(mol) == 1  # OH group
    @test num_hba(mol) == 1  # O atom
    @test num_rotatable_bonds(mol) == 0  # Isopropanol has no rotatable bonds by RDKit definition

    benzene = mol_from_smiles("c1ccccc1")
    @test num_rings(benzene) == 1
    @test num_aromatic_rings(benzene) == 1

    mols = mol_from_smiles(["CCO", "CCC", "CC(C)O"])
    weights = molecular_weight(mols)
    @test length(weights) == 3
    @test all(w -> isa(w, Float64), weights)
end

@testset "Additional Molecular Descriptors" begin
    @testset "Chi Connectivity Indices" begin
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")

        # Test Chi0n series
        @test isa(chi0n(ethanol), Float64)
        @test isa(chi1n(ethanol), Float64)
        @test isa(chi2n(ethanol), Float64)
        @test isa(chi3n(ethanol), Float64)
        @test isa(chi4n(ethanol), Float64)

        @test chi0n(ethanol) > 0
        @test chi1n(ethanol) > 0
        @test chi2n(ethanol) >= 0

        # Test Chi valence series
        @test isa(chi1v(ethanol), Float64)
        @test isa(chi2v(ethanol), Float64)
        @test isa(chi3v(ethanol), Float64)
        @test isa(chi4v(ethanol), Float64)

        @test chi1v(ethanol) > 0
        @test chi2v(ethanol) >= 0

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test chi0n(invalid_mol) === missing
        @test chi1v(invalid_mol) === missing
    end

    @testset "Kappa Shape Descriptors" begin
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")

        @test isa(kappa2(ethanol), Float64)
        @test isa(kappa3(ethanol), Float64)
        @test kappa2(ethanol) > 0
        @test kappa3(ethanol) > 0

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test kappa2(invalid_mol) === missing
        @test kappa3(invalid_mol) === missing
    end

    @testset "E-State Descriptors" begin
        ethanol = mol_from_smiles("CCO")

        max_estate = max_e_state_index(ethanol)
        min_estate = min_e_state_index(ethanol)

        @test isa(max_estate, Float64)
        @test isa(min_estate, Float64)
        @test max_estate >= min_estate

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test max_e_state_index(invalid_mol) === missing
        @test min_e_state_index(invalid_mol) === missing
    end

    @testset "Atom Counts" begin
        ethanol = mol_from_smiles("CCO")  # C2H6O
        chloroform = mol_from_smiles("CCl")  # CH3Cl
        pyridine = mol_from_smiles("c1cccnc1")  # C5H5N
        dmso = mol_from_smiles("CS(=O)C")  # C2H6OS

        # Test carbon counts
        @test num_carbons(ethanol) == 2
        @test num_carbons(chloroform) == 1
        @test num_carbons(pyridine) == 5

        # Test nitrogen counts
        @test num_nitrogens(ethanol) == 0
        @test num_nitrogens(pyridine) == 1

        # Test oxygen counts
        @test num_oxygens(ethanol) == 1
        @test num_oxygens(pyridine) == 0
        @test num_oxygens(dmso) == 1

        # Test sulfur counts
        @test num_sulfurs(ethanol) == 0
        @test num_sulfurs(dmso) == 1

        # Test halogen counts
        @test num_halogens(ethanol) == 0
        @test num_halogens(chloroform) == 1

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test num_carbons(invalid_mol) === missing
        @test num_halogens(invalid_mol) === missing
    end

    @testset "Information Content (IPC)" begin
        ethanol = mol_from_smiles("CCO")
        benzene = mol_from_smiles("c1ccccc1")

        @test isa(ipc(ethanol), Float64)
        @test isa(ipc(benzene), Float64)
        @test ipc(ethanol) > 0
        @test ipc(benzene) > 0

        # Test with invalid molecule
        invalid_mol = mol_from_smiles("invalid_smiles")
        @test ipc(invalid_mol) === missing
    end

    @testset "Vectorized Operations for New Descriptors" begin
        mols = [
            mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1"), mol_from_smiles("CC(=O)N")
        ]

        # Test vectorized Chi indices
        chi0n_values = chi0n(mols)
        chi1v_values = chi1v(mols)
        @test length(chi0n_values) == 3
        @test all(isa(v, Float64) for v in chi0n_values)

        # Test vectorized atom counts
        carbon_counts = num_carbons(mols)
        nitrogen_counts = num_nitrogens(mols)
        @test carbon_counts == [2, 6, 2]  # CCO, benzene, acetamide
        @test nitrogen_counts == [0, 0, 1]  # CCO, benzene, acetamide

        # Test with missing values
        mols_with_missing = [mol_from_smiles("CCO"), missing, mol_from_smiles("c1ccccc1")]
        carbon_mixed = num_carbons(mols_with_missing)
        @test length(carbon_mixed) == 3
        @test carbon_mixed[1] == 2
        @test carbon_mixed[2] === missing
        @test carbon_mixed[3] == 6
    end
end
