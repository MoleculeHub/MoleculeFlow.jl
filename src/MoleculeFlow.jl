module MoleculeFlow

using PythonCall
using Images
using Graphs

# Basic molecule operations
export Molecule
export mol_from_smiles, mol_to_smiles, mol_to_inchi
export mol_from_molblock, mol_from_inchi, read_sdf, read_sdf_lazy
export get_address

# Drawing
export mol_to_image, mols_to_grid_image, mol_to_svg
export highlight_substructure, draw_molecule_with_atom_labels, draw_similarity_map
export save_molecule_image, draw_functional_groups

# Molecular descriptors
export molecular_weight, exact_molecular_weight, heavy_atom_count, num_heteroatoms
export num_rotatable_bonds, num_hbd, num_hba, logp, tpsa, slogp_vsa
export num_rings, num_aromatic_rings, num_saturated_rings
export bertz_ct, balaban_j, chi0v, kappa1, calc_all_descriptors

# Fingerprints
export morgan_fingerprint, rdk_fingerprint, maccs_fingerprint
export atom_pair_fingerprint,
    topological_torsion_fingerprint, fcfp_fingerprint, pattern_fingerprint

# Similarity calculations
export tanimoto_similarity, dice_similarity, cosine_similarity, sokal_similarity
export bulk_similarity, similarity_matrix

# Atom operations
export Atom, get_atoms, get_atom
export get_atomic_number, get_symbol, get_degree, get_valence, get_formal_charge
export get_hybridization, get_num_explicit_hs, get_num_implicit_hs, get_total_num_hs
export get_mass, get_isotope, is_aromatic, is_in_ring, is_in_ring_size, get_chiral_tag
export get_neighbors, get_bonds_from_atom, compute_gasteiger_charges!, get_gasteiger_charge
export get_num_radical_electrons, is_hetero, get_cip_code, is_chiral_center
export is_hydrogen_donor, is_hydrogen_acceptor, get_ring_size
export get_crippen_log_p_contribution, get_crippen_molar_refractivity_contribution
export get_tpsa_contribution, get_labute_asa_contribution, get_all_atom_properties

# Bond operations  
export Bond, get_bond_type, get_begin_atom_idx, get_end_atom_idx

# Substructure search
export has_substructure_match, get_substructure_matches, get_substructure_match
export maximum_common_substructure, has_substructure_matches, filter_by_substructure
export has_functional_group, get_functional_groups, get_ring_info, is_ring_aromatic
export FUNCTIONAL_GROUPS

# Graph operations
export mol_to_graph, mol_to_digraph
# Re-export useful Graphs.jl functions for convenience
export nv, ne, density, is_connected, diameter, radius
export adjacency_matrix, laplacian_matrix

# Progress tracking for array operations
export ProgressTracker, update_progress!, display_progress
export with_progress, map_with_progress, format_duration

# Molecular standardization
export strip_salts, enumerate_tautomers, canonical_tautomer, standardize_molecule
export neutralize_charges, normalize_molecule

# Conformer generation (2D and 3D)
export ConformerResult, ConformerMolecule, generate_3d_conformers, generate_2d_conformers

# Molecular fragmentation
export brics_decompose, recap_decompose, get_murcko_scaffold, get_generic_scaffold
export fragment_by_bonds, get_fragment_count, split_fragments, get_largest_fragment

include("./config.jl")
include("./utils.jl")
include("./molecule/molecule.jl")
include("./molecule/descriptors.jl")
include("./molecule/fingerprints.jl")
include("./molecule/similarity.jl")
include("./molecule/graph.jl")
include("./molecule/conformers.jl")
include("./atom/atom.jl")
include("./atom/descriptors.jl")
include("./bond/bond.jl")
include("./substructure/substructure.jl")
include("./io/read.jl")
include("./io/write.jl")
include("./draw/draw.jl")
include("./standardization/standardization.jl")
include("./fragmentation/fragmentation.jl")
include("./progress.jl")

function __init__()
    _rdkit_draw[] = @pyconst(pyimport("rdkit.Chem.Draw"))
end

end
