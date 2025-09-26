module MoleculeFlow

using PythonCall
using Images: RGB, N0f8
using Graphs

# Basic molecule operations
export Molecule
export mol_from_smiles, mol_to_smiles, mol_to_inchi
export mol_from_molblock, mol_from_inchi, read_sdf, read_sdf_lazy
export get_address

# Hydrogen manipulation
export add_hs, remove_hs

# Extended file I/O
export mol_from_pdb_block, mol_from_pdb_file, mol_to_pdb_block
export mol_to_inchi_key, mol_to_molblock
export mol_from_xyz_file, mol_from_xyz_block, mol_to_xyz_block
export mol_from_mol2_file, mol_from_mol2_block

# Molecular editing and manipulation
export combine_mols, delete_substructs, replace_substructs

# Stereochemistry operations
export assign_stereochemistry!, find_chiral_centers

# Ring analysis and molecular operations
export fast_find_rings!, canonical_atom_ranks, get_ring_info, find_atom_environment
export mol_fragment_to_smarts, mol_fragment_to_smiles, renumber_atoms
export remove_stereochemistry!, sanitize_mol!, compute_2d_coords!

# Pattern matching
export quick_smarts_match, mol_fragment_to_smarts

# Drawing
export mol_to_image, mols_to_grid_image, mol_to_svg
export highlight_substructure, draw_molecule_with_atom_labels, draw_similarity_map
export save_molecule_image, draw_functional_groups

# Molecular descriptors
export molecular_weight, exact_molecular_weight, heavy_atom_count, num_heteroatoms
export num_rotatable_bonds, num_hbd, num_hba, logp, tpsa, slogp_vsa
export num_rings, num_aromatic_rings, num_saturated_rings
export bertz_ct, balaban_j, chi0v, kappa1, calc_all_descriptors

# Additional molecular connectivity and shape descriptors
export chi0n, chi1n, chi2n, chi3n, chi4n, chi1v, chi2v, chi3v, chi4v
export kappa2, kappa3, max_e_state_index, min_e_state_index, ipc

# Atom counts
export num_carbons, num_nitrogens, num_oxygens, num_sulfurs, num_halogens

# Advanced drug-like and ADMET descriptors
export qed, fraction_csp3, labute_asa, molar_refractivity, synthetic_accessibility

# Advanced ring and structure counts
export num_aliphatic_carbocycles, num_aromatic_carbocycles, num_aromatic_heterocycles
export num_atom_stereo_centers,
    num_amide_bonds, num_aliphatic_heterocycles, num_saturated_heterocycles
export num_saturated_carbocycles,
    num_unspecified_atom_stereo_centers, num_spiro_atoms, num_bridgehead_atoms
export hall_kier_alpha, num_aliphatic_rings, num_heterocycles

# BCUT descriptors
export bcut2d_mwlow,
    bcut2d_mwhi,
    bcut2d_chglow,
    bcut2d_chghi,
    bcut2d_logplow,
    bcut2d_logphi,
    bcut2d_mrlow,
    bcut2d_mrhi

# VSA descriptors
export slogp_vsa2,
    slogp_vsa3,
    slogp_vsa4,
    slogp_vsa5,
    slogp_vsa6,
    slogp_vsa7,
    slogp_vsa8,
    slogp_vsa9,
    slogp_vsa10,
    slogp_vsa11,
    slogp_vsa12
export smr_vsa1,
    smr_vsa2,
    smr_vsa3,
    smr_vsa4,
    smr_vsa5,
    smr_vsa6,
    smr_vsa7,
    smr_vsa8,
    smr_vsa9,
    smr_vsa10
export peoe_vsa1,
    peoe_vsa2,
    peoe_vsa3,
    peoe_vsa4,
    peoe_vsa5,
    peoe_vsa6,
    peoe_vsa7,
    peoe_vsa8,
    peoe_vsa9,
    peoe_vsa10,
    peoe_vsa11,
    peoe_vsa12,
    peoe_vsa13,
    peoe_vsa14

# Additional E-state descriptors
export max_absolute_e_state_index, min_absolute_e_state_index

# 3D descriptors
export asphericity,
    radius_of_gyration, pmi1, pmi2, pmi3, eccentricity, inertial_shape_factor
export spherocity_index,
    getaway_descriptors, whim_descriptors, rdf_descriptors, morse_descriptors

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

# Chemical reactions
export Reaction, reaction_from_smarts, run_reaction, reaction_to_smarts
export reaction_from_rxn_file, reaction_from_rxn_block, reaction_to_rxn_block
export validate_reaction!, is_reaction_valid, sanitize_reaction!
export get_num_reactant_templates, get_num_product_templates
export get_reactant_template, get_product_template
export has_reactant_substructure_match, get_reacting_atoms
export reaction_fingerprint, reaction_structural_fingerprint, reaction_center_fingerprint
export reaction_similarity, enumerate_library, reaction_info
export reaction_complexity, reaction_type_classification, find_similar_reactions
export is_balanced, get_atom_mapping_numbers, set_atom_mapping_numbers!
export remove_unmapped_reactant_templates!, remove_unmapped_product_templates!
export preprocess_reaction!, compute_atom_mapping!, is_template_molecule_agent

# Pharmacophore features
export FeatureFactory, ChemicalFeature
export create_feature_factory, get_mol_features, pharmacophore_fingerprint
export get_feature_families, filter_features_by_family, get_pharmacophore_3d

# Molecular alignment functions
export align_mol, calc_rms, get_best_rms
export get_alignment_transform, random_transform, apply_transform
export O3AResult, get_o3a, get_crippen_o3a, o3a_align!

include("./config.jl")
include("./utils.jl")
include("./rdkit.jl")
include("./molecule/molecule.jl")
include("./molecule/operations.jl")
include("./molecule/descriptors.jl")
include("./molecule/fingerprints.jl")
include("./molecule/similarity.jl")
include("./molecule/graph.jl")
include("./molecule/conformers.jl")
include("./molecule/alignment.jl")
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
include("./reaction/reaction.jl")
include("./pharmacophore/pharmacophore.jl")

end
