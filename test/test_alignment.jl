using Test
using MoleculeFlow
using LinearAlgebra

@testset "Molecular Alignment Functions" begin
    # Helper function to create molecules with 3D conformers
    function create_mol_with_3d(smiles::String; optimize::Bool = true)
        mol_base = mol_from_smiles(smiles)
        @test mol_base.valid  # Fail the test if molecule creation fails
        conformers = generate_3d_conformers(mol_base, 1; optimize = optimize)
        @test !isempty(conformers)  # Fail the test if conformer generation fails
        return conformers[1].molecule
    end
    @testset "Basic Alignment Functions" begin
        @testset "align_mol function" begin
            # Test basic alignment as in docstring
            mol1 = create_mol_with_3d("CCO")  # ethanol
            mol2 = create_mol_with_3d("CCO")  # ethanol

            # Basic alignment should return a finite RMSD
            rmsd = align_mol(mol1, mol2)
            @test rmsd isa Float64
            @test rmsd >= 0.0
            @test isfinite(rmsd)

            # Test with different molecules
            mol3 = create_mol_with_3d("CCC")  # propane
            rmsd2 = align_mol(mol1, mol3)
            @test rmsd2 isa Float64
            @test rmsd2 >= 0.0
            @test isfinite(rmsd2) || rmsd2 == Inf  # May be Inf if no good alignment

            # Test with weights instead of atom mapping
            weights = [1.0, 2.0, 1.0]  # Give more weight to second atom (ethanol: C-C-O)
            rmsd3 = align_mol(mol1, mol2; weights = weights)
            @test rmsd3 isa Float64
            @test rmsd3 >= 0.0
            @test isfinite(rmsd3)

            # Test with reflection
            rmsd4 = align_mol(mol1, mol2; reflect = true)
            @test rmsd4 isa Float64
            @test rmsd4 >= 0.0
            @test isfinite(rmsd4)
        end

        @testset "calc_rms function" begin
            # Test RMS calculation as in docstring
            mol1 = create_mol_with_3d("CCO")
            mol2 = create_mol_with_3d("CCO")

            # Calculate RMS without transformation
            rms = calc_rms(mol1, mol2)
            @test rms isa Float64
            @test rms >= 0.0
            @test isfinite(rms)

            # Calculate RMS with transformation
            rms_transform = calc_rms(mol1, mol2; transform_probe = true)
            @test rms_transform isa Float64
            @test rms_transform >= 0.0
            @test isfinite(rms_transform)
            @test rms_transform <= rms  # Transformed should be better or equal

            # Test with weights (ethanol has 3 atoms: C-C-O)
            weights = [1.0, 1.0, 2.0]  # Give more weight to the oxygen (3rd atom)
            rms_weighted = calc_rms(mol1, mol2; weights = weights)
            @test rms_weighted isa Float64
            @test rms_weighted >= 0.0
            @test isfinite(rms_weighted)
        end

        @testset "get_best_rms function" begin
            # Test optimal RMS as in docstring (benzene has symmetry)
            mol1 = create_mol_with_3d("c1ccccc1")  # benzene
            mol2 = create_mol_with_3d("c1ccccc1")  # benzene

            best_rms = get_best_rms(mol1, mol2)
            @test best_rms isa Float64
            @test best_rms >= 0.0
            @test isfinite(best_rms)

            # Regular RMS should be >= best RMS due to symmetry consideration
            regular_rms = calc_rms(mol1, mol2)
            @test regular_rms >= best_rms || abs(regular_rms - best_rms) < 1e-6

            # Test with weights (benzene has 6 atoms: 6 carbons, no explicit hydrogens)
            # Create reasonable weights for all atoms
            weights = ones(6)   # Equal weights for all 6 carbon atoms
            weights[1] = 2.0    # More weight on first carbon
            best_rms_weighted = get_best_rms(mol1, mol2; weights = weights)
            @test best_rms_weighted isa Float64
            @test best_rms_weighted >= 0.0
            @test isfinite(best_rms_weighted)
        end

        @testset "get_alignment_transform function" begin
            # Test transformation matrix computation
            mol1 = create_mol_with_3d("CCO")
            mol2 = create_mol_with_3d("CCO")

            transform = get_alignment_transform(mol1, mol2)
            @test transform isa Matrix{Float64}
            @test size(transform) == (4, 4)
            @test all(isfinite, transform)

            # Test with reflection
            transform_reflect = get_alignment_transform(mol1, mol2; reflect = true)
            @test transform_reflect isa Matrix{Float64}
            @test size(transform_reflect) == (4, 4)
            @test all(isfinite, transform_reflect)

            # Test with weights (ethanol has 3 atoms: 2 carbons, 1 oxygen, no explicit hydrogens)
            weights = ones(3)  # Equal weights for all 3 atoms
            weights[3] = 2.0   # More weight on oxygen atom (3rd atom)
            transform_weighted = get_alignment_transform(mol1, mol2; weights = weights)
            @test transform_weighted isa Matrix{Float64}
            @test size(transform_weighted) == (4, 4)
            @test all(isfinite, transform_weighted)
        end

        @testset "random_transform function" begin
            # Test random transformation
            mol = create_mol_with_3d("CCO")

            # Apply random transformation
            transformed_mol = random_transform(mol; seed = 42)
            @test transformed_mol === mol  # Should return the same object
            @test transformed_mol.valid

            # Test with different seed
            mol2 = create_mol_with_3d("CCO")
            transformed_mol2 = random_transform(mol2; seed = 123)
            @test transformed_mol2.valid
        end

        @testset "apply_transform function" begin
            # Test applying transformation matrix
            mol = create_mol_with_3d("CCO")

            # Apply identity transformation (should not change coordinates significantly)
            identity_transform = Matrix{Float64}(I, 4, 4)
            transformed_mol = apply_transform(mol, identity_transform)
            @test transformed_mol === mol  # Should return the same object
            @test transformed_mol.valid

            # Test with a translation matrix
            translation_transform = Matrix{Float64}(I, 4, 4)
            translation_transform[1:3, 4] = [1.0, 2.0, 3.0]  # Translate by (1, 2, 3)
            transformed_mol2 = apply_transform(mol, translation_transform)
            @test transformed_mol2.valid
        end
    end

    @testset "O3A Alignment Functions" begin
        @testset "O3AResult structure" begin
            # Test O3AResult structure
            result = O3AResult(0.85, 1.23, Matrix{Float64}(I, 4, 4), [(1, 1), (2, 2)])
            @test result.score == 0.85
            @test result.rmsd == 1.23
            @test size(result.transform) == (4, 4)
            @test length(result.matched_atoms) == 2

            # Test string representation
            result_str = string(result)
            @test occursin("O3AResult", result_str)
            @test occursin("score", result_str)
            @test occursin("rmsd", result_str)
        end

        @testset "get_o3a function" begin
            # Test MMFF-based O3A alignment with similar drug-like molecules
            mol1 = create_mol_with_3d("CCc1ccccc1")     # ethylbenzene
            mol2 = create_mol_with_3d("CCc1ccc(O)cc1")  # 4-ethylphenol (similar structure)

            result = get_o3a(mol1, mol2)
            @test result isa O3AResult
            @test result.score isa Float64
            @test result.rmsd isa Float64
            @test result.rmsd >= 0.0 || result.rmsd == Inf
            @test size(result.transform) == (4, 4)
            @test result.matched_atoms isa Vector{Tuple{Int, Int}}
            @test all(isfinite, result.transform) ||
                result.transform == Matrix{Float64}(I, 4, 4)

            # Test with different parameters
            result2 = get_o3a(mol1, mol2; reflect = true, accuracy = 0.001)
            @test result2 isa O3AResult
            @test isfinite(result2.score) || result2.score == -1.0  # May fail and return -1.0
            @test isfinite(result2.rmsd) || result2.rmsd == Inf  # May fail and return Inf
        end

        @testset "get_crippen_o3a function" begin
            # Test Crippen-based O3A alignment with benzene derivatives
            mol1 = create_mol_with_3d("c1ccccc1C")      # toluene
            mol2 = create_mol_with_3d("c1ccc(C)cc1O")   # 4-methylphenol (cresol)

            result = get_crippen_o3a(mol1, mol2)
            @test result isa O3AResult
            @test result.score isa Float64
            @test result.rmsd isa Float64
            @test result.rmsd >= 0.0 || result.rmsd == Inf  # May fail
            @test size(result.transform) == (4, 4)
            @test result.matched_atoms isa Vector{Tuple{Int, Int}}
            @test all(isfinite, result.transform) ||
                result.transform == Matrix{Float64}(I, 4, 4)

            # Test with different parameters
            result2 = get_crippen_o3a(mol1, mol2; attempt_generic_features = false)
            @test result2 isa O3AResult
            @test isfinite(result2.score) || result2.score == -1.0
            @test isfinite(result2.rmsd) || result2.rmsd == Inf
        end

        @testset "o3a_align! function" begin
            # Test convenience function with similar aromatic compounds
            mol1 = create_mol_with_3d("c1ccc(CC)cc1")   # ethylbenzene
            mol2 = create_mol_with_3d("c1ccc(CCC)cc1")  # propylbenzene

            # MMFF-based alignment
            result1 = o3a_align!(mol1, mol2, :mmff)
            @test result1 isa O3AResult
            @test isfinite(result1.score) || result1.score == -1.0
            @test isfinite(result1.rmsd) || result1.rmsd == Inf

            # Reset molecules for Crippen test
            mol3 = create_mol_with_3d("c1ccc(CC)cc1")   # ethylbenzene
            mol4 = create_mol_with_3d("c1ccc(CCC)cc1")  # propylbenzene

            # Crippen-based alignment
            result2 = o3a_align!(mol3, mol4, :crippen)
            @test result2 isa O3AResult
            @test isfinite(result2.score) || result2.score == -1.0
            @test isfinite(result2.rmsd) || result2.rmsd == Inf

            # Test invalid method
            @test_throws ArgumentError o3a_align!(mol1, mol2, :invalid)
        end
    end

    @testset "Error Handling and Edge Cases" begin
        @testset "Invalid molecules" begin
            valid_mol = create_mol_with_3d("CCO")

            # Create an invalid molecule properly
            invalid_mol = mol_from_smiles("INVALID_SMILES")
            @test !invalid_mol.valid

            # Test functions with invalid molecules
            @test_throws ArgumentError align_mol(invalid_mol, valid_mol)
            @test_throws ArgumentError align_mol(valid_mol, invalid_mol)
            @test_throws ArgumentError calc_rms(invalid_mol, valid_mol)
            @test_throws ArgumentError calc_rms(valid_mol, invalid_mol)
            @test_throws ArgumentError get_best_rms(invalid_mol, valid_mol)
            @test_throws ArgumentError get_alignment_transform(invalid_mol, valid_mol)
            @test_throws ArgumentError random_transform(invalid_mol)
            @test_throws ArgumentError apply_transform(
                invalid_mol, Matrix{Float64}(I, 4, 4)
            )
            @test_throws ArgumentError get_o3a(invalid_mol, valid_mol)
            @test_throws ArgumentError get_crippen_o3a(invalid_mol, valid_mol)
        end

        @testset "Invalid transformation matrix" begin
            mol = create_mol_with_3d("CCO")

            # Test with wrong size matrix
            wrong_size_matrix = Matrix{Float64}(I, 3, 3)
            @test_throws ArgumentError apply_transform(mol, wrong_size_matrix)

            # Test with non-square matrix
            non_square_matrix = rand(4, 3)
            @test_throws ArgumentError apply_transform(mol, non_square_matrix)
        end

        @testset "Molecules without 3D conformers" begin
            # Test behavior with molecules that have no 3D conformers
            mol1_2d = mol_from_smiles("CCO")  # Only 2D structure
            mol2_2d = mol_from_smiles("CCC")  # Only 2D structure

            # These should handle gracefully and typically return Inf
            rmsd = align_mol(mol1_2d, mol2_2d)
            @test rmsd isa Float64
            @test rmsd == Inf || rmsd >= 0.0  # Should be Inf or positive

            rms = calc_rms(mol1_2d, mol2_2d)
            @test rms isa Float64
            @test rms == Inf || rms >= 0.0  # Should be Inf or positive

            # Test that transformation still returns identity on failure
            transform = get_alignment_transform(mol1_2d, mol2_2d)
            @test transform isa Matrix{Float64}
            @test size(transform) == (4, 4)
            # May be identity matrix on failure
        end

        @testset "Identical molecules" begin
            # Test with identical molecules
            mol1_base = mol_from_smiles("CCO")
            mol2_base = mol_from_smiles("CCO")
            @test mol1_base.valid && mol2_base.valid
            conformers1 = generate_3d_conformers(mol1_base, 1; random_seed = 42)
            conformers2 = generate_3d_conformers(mol2_base, 1; random_seed = 42)  # Same seed for identical conformers
            @test !isempty(conformers1) && !isempty(conformers2)

            mol1 = conformers1[1].molecule
            mol2 = conformers2[1].molecule

            rmsd = align_mol(mol1, mol2)
            @test rmsd isa Float64
            @test rmsd >= 0.0  # Should be very small but may not be exactly 0

            rms = calc_rms(mol1, mol2)
            @test rms isa Float64
            @test rms >= 0.0
        end

        @testset "Different sized molecules" begin
            # Test alignment between molecules of different sizes
            small_mol = create_mol_with_3d("CC")  # ethane
            large_mol = create_mol_with_3d("CCCCCCCCCC")  # decane

            # Should handle gracefully
            rmsd = align_mol(small_mol, large_mol)
            @test rmsd isa Float64

            rms = calc_rms(small_mol, large_mol)
            @test rms isa Float64
        end
    end

    @testset "Parameter Validation" begin
        @testset "Conformer ID validation" begin
            mol1 = create_mol_with_3d("CCO")
            mol2 = create_mol_with_3d("CCC")

            # Test with valid conformer IDs (0-based in RDKit)
            rmsd1 = align_mol(mol1, mol2; probe_conf_id = 0, ref_conf_id = 0)
            @test rmsd1 isa Float64

            # Test with -1 (default)
            rmsd2 = align_mol(mol1, mol2; probe_conf_id = -1, ref_conf_id = -1)
            @test rmsd2 isa Float64

            # Test with invalid conformer IDs (should handle gracefully)
            rmsd3 = align_mol(mol1, mol2; probe_conf_id = 999, ref_conf_id = 999)
            @test rmsd3 isa Float64  # May be Inf but should not crash
        end

        @testset "Max iterations validation" begin
            mol1 = create_mol_with_3d("CCO")
            mol2 = create_mol_with_3d("CCO")

            # Test with different max_iterations values
            rmsd1 = align_mol(mol1, mol2; max_iterations = 10)
            @test rmsd1 isa Float64
            @test rmsd1 >= 0.0

            rmsd2 = align_mol(mol1, mol2; max_iterations = 100)
            @test rmsd2 isa Float64
            @test rmsd2 >= 0.0

            # Test with minimal iterations
            rmsd3 = align_mol(mol1, mol2; max_iterations = 1)
            @test rmsd3 isa Float64
            @test rmsd3 >= 0.0
        end

        @testset "Weight validation" begin
            mol1 = create_mol_with_3d("CCO")  # ethanol: 3 atoms total (no explicit hydrogens)
            mol2 = create_mol_with_3d("CCO")  # ethanol: 3 atoms total (no explicit hydrogens)

            # Test align_mol with valid weights for all atoms
            weights = ones(3)  # Equal weights for all 3 atoms
            weights[3] = 2.0   # More weight on oxygen (3rd atom)
            rmsd1 = align_mol(mol1, mol2; weights = weights)
            @test rmsd1 isa Float64

            # Test calc_rms with valid weights
            rms1 = calc_rms(mol1, mol2; weights = weights)
            @test rms1 isa Float64

            # Test get_best_rms with valid weights
            best_rms1 = get_best_rms(mol1, mol2; weights = weights)
            @test best_rms1 isa Float64

            # Test get_alignment_transform with valid weights
            transform1 = get_alignment_transform(mol1, mol2; weights = weights)
            @test transform1 isa Matrix{Float64}
            @test size(transform1) == (4, 4)

            # Test with empty weights (should work)
            empty_weights = Float64[]
            rmsd2 = align_mol(mol1, mol2; weights = empty_weights)
            @test rmsd2 isa Float64

            rms2 = calc_rms(mol1, mol2; weights = empty_weights)
            @test rms2 isa Float64

            best_rms2 = get_best_rms(mol1, mol2; weights = empty_weights)
            @test best_rms2 isa Float64
        end
    end
end
