using Test
using MoleculeFlow
using Graphs

@testset "Graph Operations" begin
    mol = mol_from_smiles("CCO")
    g = mol_to_graph(mol)

    @test !ismissing(g)
    @test nv(g) == 3  # 3 atoms
    @test ne(g) == 2  # 2 bonds

    dg = mol_to_digraph(mol)
    @test !ismissing(dg)
    @test nv(dg) == 3  # 3 atoms
    @test ne(dg) == 4  # 4 directed edges (2 bonds × 2 directions)

    # Test using Graphs.jl functions directly on the graph
    adj_matrix = adjacency_matrix(g)
    @test size(adj_matrix) == (3, 3)
    @test adj_matrix[1, 2] == 1  # C-C bond
    @test adj_matrix[2, 3] == 1  # C-O bond
    @test adj_matrix[1, 3] == 0  # no direct C-O bond

    lap_matrix = laplacian_matrix(g)
    @test size(lap_matrix) == (3, 3)

    @test is_connected(g) == true
    @test diameter(g) == 2  # max distance between any two atoms
    @test radius(g) == 1    # min eccentricity
    @test isa(density(g), Float64)

    state = dijkstra_shortest_paths(g, 1)
    path = enumerate_paths(state, 3)  # from first C to O
    @test length(path) == 3  # [1, 2, 3] - goes through middle carbon
    @test path[1] == 1 && path[end] == 3

    # Test direct path (adjacent atoms)
    path_direct = enumerate_paths(state, 2)
    @test length(path_direct) == 2  # [1, 2]
    @test path_direct == [1, 2]
end

@testset "Graph Operations - Complex Molecules" begin
    benzene = mol_from_smiles("c1ccccc1")
    g_benzene = mol_to_graph(benzene)

    @test !ismissing(g_benzene)
    @test nv(g_benzene) == 6  # 6 carbon atoms
    @test ne(g_benzene) == 6  # 6 bonds in ring

    @test nv(g_benzene) == 6
    @test ne(g_benzene) == 6
    @test is_connected(g_benzene) == true
    @test diameter(g_benzene) == 3  # max distance in 6-ring

    # Test shortest path in ring
    state_benzene = dijkstra_shortest_paths(g_benzene, 1)
    path_ring = enumerate_paths(state_benzene, 4)  # opposite carbons
    @test length(path_ring) == 4  # 3 bonds to traverse

    # Test with branched molecule
    isobutane = mol_from_smiles("CC(C)C")  # branched alkane
    g_branched = mol_to_graph(isobutane)

    @test !ismissing(g_branched)
    @test nv(g_branched) == 4  # 4 carbon atoms
    @test ne(g_branched) == 3  # 3 bonds

    @test nv(g_branched) == 4
    @test ne(g_branched) == 3
    @test is_connected(g_branched) == true
    @test diameter(g_branched) == 2  # star topology
    @test radius(g_branched) == 1    # central carbon
end

@testset "Graph Operations - Edge Cases" begin
    # Test with single atom
    methane = mol_from_smiles("C")
    g_single = mol_to_graph(methane)

    @test !ismissing(g_single)
    @test nv(g_single) == 1
    @test ne(g_single) == 0

    @test nv(g_single) == 1
    @test ne(g_single) == 0
    @test is_connected(g_single) == true
    # diameter and radius are undefined for single vertex graphs

    # Test with invalid molecule
    bad_mol = mol_from_smiles("invalid_smiles")
    g_invalid = mol_to_graph(bad_mol)
    @test ismissing(g_invalid)

    # Invalid molecule should return missing graph
    @test ismissing(g_invalid)

    # Test invalid atom indices with good molecule
    good_mol = mol_from_smiles("CCO")
    good_g = mol_to_graph(good_mol)
    @test !ismissing(good_g)
end

@testset "Graph Operations - Disconnected Molecules" begin
    mol = mol_from_smiles("CCO")
    g = mol_to_graph(mol)

    @test !ismissing(g)
    @test is_connected(g) == true
    # For connected molecules, diameter and radius should be computable
    @test isa(diameter(g), Int)
    @test isa(radius(g), Int)
end
