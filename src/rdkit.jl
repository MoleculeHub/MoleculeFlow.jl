# RDKit Python bindings

# Basic RDKit Chem module functions
_mol_from_smiles(smiles::String) = @pyconst(pyimport("rdkit.Chem").MolFromSmiles)(smiles)
_mol_from_inchi(inchi::String) = @pyconst(pyimport("rdkit.Chem").MolFromInchi)(inchi)
function _mol_from_molfile(molfile::String)
    @pyconst(pyimport("rdkit.Chem").MolFromMolFile)(molfile)
end
function _mol_from_molblock(molblock::String)
    @pyconst(pyimport("rdkit.Chem").MolFromMolBlock)(molblock)
end
function _sdf_supplier(filename::String)
    @pyconst(pyimport("rdkit.Chem").SDMolSupplier)(filename)
end
_mol_to_smiles(mol::Py) = @pyconst(pyimport("rdkit.Chem").MolToSmiles)(mol)
_mol_to_inchi(mol::Py) = @pyconst(pyimport("rdkit.Chem").MolToInchi)(mol)
_mol_from_smarts(smarts::String) = @pyconst(pyimport("rdkit.Chem").MolFromSmarts)(smarts)
_mol_copy(mol::Py) = @pyconst(pyimport("rdkit.Chem").Mol)(mol)

# Molecular descriptors
function _calc_mol_descriptors(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").CalcMolDescriptors)(mol)
end
function _mol_wt(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MolWt)(mol)
end
function _exact_mol_wt(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").ExactMolWt)(mol)
end
function _heavy_atom_count(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").HeavyAtomCount)(mol)
end
function _num_heteroatoms(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumHeteroatoms)(mol)
end
function _num_hdonors(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumHDonors)(mol)
end
function _num_hacceptors(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumHAcceptors)(mol)
end
function _ring_count(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").RingCount)(mol)
end
function _num_aromatic_rings(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAromaticRings)(mol)
end
function _num_saturated_rings(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumSaturatedRings)(mol)
end
function _bertz_ct(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BertzCT)(mol)
end
function _balaban_j(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BalabanJ)(mol)
end
function _chi0v(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi0v)(mol)
end
function _kappa1(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Kappa1)(mol)
end
function _tpsa(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").TPSA)(mol)
end
function _slogp_vsa1(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA1)(mol)
end

# Advanced descriptors and drug-like properties
function _qed(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").qed)(mol)
end
function _fraction_csp3(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").FractionCSP3)(mol)
end
function _labute_asa(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").LabuteASA)(mol)
end
function _mol_mr(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MolMR)(mol)
end
function _hall_kier_alpha(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").HallKierAlpha)(mol)
end

# Functional group and structure counts
function _num_aliphatic_carbocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAliphaticCarbocycles)(mol)
end
function _num_aliphatic_heterocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAliphaticHeterocycles)(mol)
end
function _num_aromatic_carbocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAromaticCarbocycles)(mol)
end
function _num_aromatic_heterocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAromaticHeterocycles)(mol)
end
function _num_saturated_carbocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumSaturatedCarbocycles)(mol)
end
function _num_saturated_heterocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumSaturatedHeterocycles)(mol)
end
function _num_atom_stereo_centers(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAtomStereoCenters)(mol)
end
function _num_unspecified_atom_stereo_centers(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumUnspecifiedAtomStereoCenters)(mol)
end
function _num_amide_bonds(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAmideBonds)(mol)
end
function _num_spiro_atoms(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumSpiroAtoms)(mol)
end
function _num_bridgehead_atoms(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumBridgeheadAtoms)(mol)
end

# BCUT descriptors
function _bcut2d_mwlow(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_MWLOW)(mol)
end
function _bcut2d_mwhi(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_MWHI)(mol)
end
function _bcut2d_chglow(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_CHGLO)(mol)
end
function _bcut2d_chghi(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_CHGHI)(mol)
end
function _bcut2d_logplow(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_LOGPLOW)(mol)
end
function _bcut2d_logphi(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_LOGPHI)(mol)
end
function _bcut2d_mrlow(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_MRLOW)(mol)
end
function _bcut2d_mrhi(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BCUT2D_MRHI)(mol)
end

# 3D descriptors
function _asphericity(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").Asphericity)(mol; confId = confId)
end
function _eccentricity(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").Eccentricity)(mol; confId = confId)
end
function _radius_of_gyration(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").RadiusOfGyration)(mol; confId = confId)
end
function _pmi1(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").PMI1)(mol; confId = confId)
end
function _pmi2(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").PMI2)(mol; confId = confId)
end
function _pmi3(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").PMI3)(mol; confId = confId)
end
function _inertial_shape_factor(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").InertialShapeFactor)(mol; confId = confId)
end

# Advanced 3D descriptors
function _spherocity_index(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Descriptors3D").SpherocityIndex)(mol; confId = confId)
end

function _calc_getaway(
    mol::Py; confId::Int = -1, precision::Int = 2, custom_atom_property::String = ""
)
    if custom_atom_property == ""
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcGETAWAY)(
            mol; confId = confId, precision = precision
        )
    else
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcGETAWAY)(
            mol;
            confId = confId,
            precision = precision,
            CustomAtomProperty = custom_atom_property,
        )
    end
end

function _calc_whim(
    mol::Py; confId::Int = -1, thresh::Float64 = 0.001, custom_atom_property::String = ""
)
    if custom_atom_property == ""
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcWHIM)(
            mol; confId = confId, thresh = thresh
        )
    else
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcWHIM)(
            mol; confId = confId, thresh = thresh, CustomAtomProperty = custom_atom_property
        )
    end
end

function _calc_rdf(mol::Py; confId::Int = -1, custom_atom_property::String = "")
    if custom_atom_property == ""
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcRDF)(mol; confId = confId)
    else
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcRDF)(
            mol; confId = confId, CustomAtomProperty = custom_atom_property
        )
    end
end

function _calc_morse(mol::Py; confId::Int = -1, custom_atom_property::String = "")
    if custom_atom_property == ""
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcMORSE)(mol; confId = confId)
    else
        @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").CalcMORSE)(
            mol; confId = confId, CustomAtomProperty = custom_atom_property
        )
    end
end

# Synthetic Accessibility Score
function _sascore(mol::Py)
    @pyconst(pyimport("rdkit.Contrib.SA_Score.sascorer").calculateScore)(mol)
end

# Crippen descriptors
function _mol_logp(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Crippen").MolLogP)(mol)
end
function _get_atom_contribs(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Crippen")._GetAtomContribs)(mol)
end

# Lipinski descriptors
function _num_rotatable_bonds(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Lipinski").NumRotatableBonds)(mol)
end

# RDMolDescriptors
function _tpsa_contribs(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdMolDescriptors")._TPSAContribs)(mol)
end
function _labute_asa_contribs(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdMolDescriptors")._LabuteASAContribs)(mol)
end

# Partial charges
function _compute_gasteiger_charges(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdPartialCharges").ComputeGasteigerCharges)(mol)
end

# Conformer generation
function _etkdg_v3()
    @pyconst(pyimport("rdkit.Chem.rdDistGeom").ETKDGv3)()
end
function _embed_molecule(mol::Py, params::Py)
    @pyconst(pyimport("rdkit.Chem.rdDistGeom").EmbedMolecule)(mol, params)
end
function _embed_multiple_confs_with_kwargs(mol::Py; numConfs::Int, params::Py)
    @pyconst(pyimport("rdkit.Chem.rdDistGeom").EmbedMultipleConfs)(
        mol; numConfs = numConfs, params = params
    )
end
function _embed_multiple_confs(mol::Py, numConfs::Int, params::Py)
    @pyconst(pyimport("rdkit.Chem.rdDistGeom").EmbedMultipleConfs)(mol, numConfs, params)
end
function _mmff_get_molecule_force_field(mol::Py, props::Py; confId::Int)
    @pyconst(pyimport("rdkit.Chem.rdForceFieldHelpers").MMFFGetMoleculeForceField)(
        mol, props; confId = confId
    )
end
function _mmff_get_molecule_properties(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdForceFieldHelpers").MMFFGetMoleculeProperties)(mol)
end
function _uff_get_molecule_force_field(mol::Py; confId::Int)
    @pyconst(pyimport("rdkit.Chem.rdForceFieldHelpers").UFFGetMoleculeForceField)(
        mol; confId = confId
    )
end
function _add_hs(mol::Py)
    @pyconst(pyimport("rdkit.Chem.AllChem").AddHs)(mol)
end
function _remove_hs(mol::Py)
    @pyconst(pyimport("rdkit.Chem.AllChem").RemoveHs)(mol)
end
function _compute_2d_coords(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdDepictor").Compute2DCoords)(mol)
end

# Additional molecular creation/conversion
function _mol_from_pdb_block(pdb_block::String)
    @pyconst(pyimport("rdkit.Chem").MolFromPDBBlock)(pdb_block)
end
function _mol_from_pdb_file(filename::String)
    @pyconst(pyimport("rdkit.Chem").MolFromPDBFile)(filename)
end
function _mol_from_xyz_block(xyz_block::String)
    @pyconst(pyimport("rdkit.Chem").MolFromXYZBlock)(xyz_block)
end
function _mol_from_xyz_file(filename::String)
    @pyconst(pyimport("rdkit.Chem").MolFromXYZFile)(filename)
end
function _mol_from_mol2_block(mol2_block::String)
    @pyconst(pyimport("rdkit.Chem").MolFromMol2Block)(mol2_block)
end
function _mol_from_mol2_file(filename::String)
    @pyconst(pyimport("rdkit.Chem").MolFromMol2File)(filename)
end
function _mol_to_pdb_block(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToPDBBlock)(mol)
end
function _mol_to_pdb_file(mol::Py, filename::String)
    @pyconst(pyimport("rdkit.Chem").MolToPDBFile)(mol, filename)
end
function _mol_to_xyz_block(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToXYZBlock)(mol)
end
function _mol_to_xyz_file(mol::Py, filename::String)
    @pyconst(pyimport("rdkit.Chem").MolToXYZFile)(mol, filename)
end
function _mol_to_inchi_key(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToInchiKey)(mol)
end
function _mol_to_inchi_and_aux_info(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToInchiAndAuxInfo)(mol)
end
function _mol_to_molblock(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToMolBlock)(mol)
end
function _mol_to_v3k_molblock(mol::Py)
    @pyconst(pyimport("rdkit.Chem").MolToV3KMolBlock)(mol)
end

# Molecular editing and manipulation
function _combine_mols(mol1::Py, mol2::Py)
    @pyconst(pyimport("rdkit.Chem").CombineMols)(mol1, mol2)
end
function _delete_substructs(mol::Py, query::Py)
    @pyconst(pyimport("rdkit.Chem").DeleteSubstructs)(mol, query)
end
function _replace_substructs(mol::Py, query::Py, replacement::Py)
    @pyconst(pyimport("rdkit.Chem").ReplaceSubstructs)(mol, query, replacement)
end
function _renumber_atoms(mol::Py, new_order::Vector{Int})
    @pyconst(pyimport("rdkit.Chem").RenumberAtoms)(mol, new_order)
end
function _split_mol_by_pdb_chain_id(mol::Py)
    @pyconst(pyimport("rdkit.Chem").SplitMolByPDBChainId)(mol)
end
function _split_mol_by_pdb_residues(mol::Py)
    @pyconst(pyimport("rdkit.Chem").SplitMolByPDBResidues)(mol)
end

# Stereochemistry and 3D operations
function _assign_stereochemistry(mol::Py; cleanIt::Bool = true, force::Bool = false)
    @pyconst(pyimport("rdkit.Chem").AssignStereochemistry)(
        mol; cleanIt = cleanIt, force = force
    )
end
function _assign_stereochemistry_from_3d(mol::Py; confId::Int = -1)
    @pyconst(pyimport("rdkit.Chem").AssignStereochemistryFrom3D)(mol; confId = confId)
end
function _detect_bond_stereochemistry(mol::Py, bond_idx::Int)
    @pyconst(pyimport("rdkit.Chem").DetectBondStereochemistry)(mol, bond_idx)
end
function _find_mol_chiral_centers(mol::Py; includeUnassigned::Bool = false)
    @pyconst(pyimport("rdkit.Chem").FindMolChiralCenters)(
        mol; includeUnassigned = includeUnassigned
    )
end
function _find_potential_stereo(mol::Py)
    @pyconst(pyimport("rdkit.Chem").FindPotentialStereo)(mol)
end
function _canonicalize_enhanced_stereo(mol::Py)
    @pyconst(pyimport("rdkit.Chem").CanonicalizeEnhancedStereo)(mol)
end
function _set_bond_stereo_from_directions(mol::Py)
    @pyconst(pyimport("rdkit.Chem").SetBondStereoFromDirections)(mol)
end
function _wedge_mol_bonds(mol::Py; wedgeBonds::Bool = true)
    @pyconst(pyimport("rdkit.Chem").WedgeMolBonds)(mol; wedgeBonds = wedgeBonds)
end
function _set_terminal_atom_coords(
    mol::Py, confId::Int, terminalAtomIdx::Int, otherAtomIdx::Int, bondLength::Float64
)
    @pyconst(pyimport("rdkit.Chem").SetTerminalAtomCoords)(
        mol, confId, terminalAtomIdx, otherAtomIdx, bondLength
    )
end

# Ring analysis and aromaticity
function _fast_find_rings(mol::Py)
    @pyconst(pyimport("rdkit.Chem").FastFindRings)(mol)
end
function _find_ring_families(mol::Py)
    @pyconst(pyimport("rdkit.Chem").FindRingFamilies)(mol)
end
function _get_ring_info(mol::Py)
    mol.GetRingInfo()
end
function _canonical_rank_atoms(mol::Py)
    @pyconst(pyimport("rdkit.Chem").CanonicalRankAtoms)(mol)
end
function _canonical_rank_atoms_in_fragment(mol::Py, atomsToUse::Vector{Int})
    @pyconst(pyimport("rdkit.Chem").CanonicalRankAtomsInFragment)(mol, atomsToUse)
end
function _find_atom_environment_of_radius_n(mol::Py, radius::Int, atomIdx::Int)
    @pyconst(pyimport("rdkit.Chem").FindAtomEnvironmentOfRadiusN)(mol, radius, atomIdx)
end

# Fragmentation
function _brics_decompose(mol::Py)
    @pyconst(pyimport("rdkit.Chem.BRICS").BRICSDecompose)(mol)
end
function _recap_decompose(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Recap").Decompose)(mol)
end
function _get_scaffold_for_mol(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Scaffolds.MurckoScaffold").GetScaffoldForMol)(mol)
end
function _make_scaffold_generic(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Scaffolds.MurckoScaffold").MakeScaffoldGeneric)(mol)
end
function _fragment_on_bonds(mol::Py, bond_indices::Vector{Int})
    @pyconst(pyimport("rdkit.Chem").FragmentOnBonds)(mol, bond_indices)
end
function _get_mol_frags(mol::Py; asMols::Bool = false)
    @pyconst(pyimport("rdkit.Chem").GetMolFrags)(mol; asMols = asMols)
end

# Substructure search and pattern matching
function _find_mcs(mol_list::Vector)
    @pyconst(pyimport("rdkit.Chem.rdFMCS").FindMCS)(mol_list)
end
function _quick_smarts_match(mol::Py, smarts::String)
    @pyconst(pyimport("rdkit.Chem").QuickSmartsMatch)(mol, smarts)
end
function _mol_fragment_to_smarts(mol::Py, atomsToUse::Vector{Int})
    @pyconst(pyimport("rdkit.Chem").MolFragmentToSmarts)(mol, atomsToUse)
end
function _mol_fragment_to_smiles(mol::Py, atomsToUse::Vector{Int})
    @pyconst(pyimport("rdkit.Chem").MolFragmentToSmiles)(mol, atomsToUse)
end
function _get_most_substituted_core_match(mol::Py, core::Py)
    @pyconst(pyimport("rdkit.Chem").GetMostSubstitutedCoreMatch)(mol, core)
end

# Fingerprints
function _rdkit_fingerprint(mol::Py)
    @pyconst(pyimport("rdkit.Chem").RDKFingerprint)(mol)
end
function _get_rdk_fingerprint(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Fingerprints.FingerprintMols").GetRDKFingerprint)(mol)
end
function _morgan_fingerprint(mol::Py, radius::Int)
    @pyconst(pyimport("rdkit.Chem").GetMorganFingerprint)(mol, radius)
end
function _morgan_fingerprint_as_bit_vect(mol::Py, radius::Int, nBits::Int)
    @pyconst(pyimport("rdkit.Chem").GetMorganFingerprintAsBitVect)(mol, radius, nBits)
end
function _atom_pair_fingerprint(mol::Py)
    @pyconst(pyimport("rdkit.Chem").GetAtomPairFingerprint)(mol)
end
function _topological_torsion_fingerprint(mol::Py)
    @pyconst(pyimport("rdkit.Chem").GetTopologicalTorsionFingerprint)(mol)
end
function _maccs_keys(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").GetMACCSKeysFingerprint)(mol)
end
function _get_morgan_generator(; radius::Int, fpSize::Int)
    @pyconst(pyimport("rdkit.Chem.rdFingerprintGenerator").GetMorganGenerator)(;
        radius = radius, fpSize = fpSize
    )
end
function _get_morgan_feature_atom_inv_gen()
    @pyconst(pyimport("rdkit.Chem.rdFingerprintGenerator").GetMorganFeatureAtomInvGen)()
end
function _get_hashed_atom_pair_fingerprint_as_bit_vect(mol::Py; nBits::Int)
    @pyconst(pyimport("rdkit.Chem.rdMolDescriptors").GetHashedAtomPairFingerprintAsBitVect)(
        mol; nBits = nBits
    )
end
function _get_hashed_topological_torsion_fingerprint_as_bit_vect(mol::Py; nBits::Int)
    @pyconst(
        pyimport("rdkit.Chem.rdMolDescriptors").GetHashedTopologicalTorsionFingerprintAsBitVect
    )(
        mol; nBits = nBits
    )
end
function _pattern_fingerprint(mol::Py; fpSize::Int)
    @pyconst(pyimport("rdkit.Chem").PatternFingerprint)(mol; fpSize = fpSize)
end

# Standardization
function _cleanup(mol::Py)
    @pyconst(pyimport("rdkit.Chem.MolStandardize").Cleanup)(mol)
end
function _standardize_smiles(smiles::String)
    @pyconst(pyimport("rdkit.Chem.MolStandardize").standardize_smiles)(smiles)
end
function _charge_parent(mol::Py)
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").ChargeParent)(mol)
end
function _fragment_parent(mol::Py)
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").FragmentParent)(mol)
end
function _largest_fragment_chooser()
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").LargestFragmentChooser)()
end
function _tautomer_enumerator()
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").TautomerEnumerator)()
end
function _uncharger()
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").Uncharger)()
end
function _normalizer()
    @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize").Normalizer)()
end
function _remove_stereochemistry(mol::Py)
    @pyconst(pyimport("rdkit.Chem").RemoveStereochemistry)(mol)
end
function _sanitize_mol(mol::Py)
    @pyconst(pyimport("rdkit.Chem").SanitizeMol)(mol)
end

# Additional commonly used descriptors - low hanging fruit
function _chi0n(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi0n)(mol)
end
function _chi1n(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi1n)(mol)
end
function _chi2n(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi2n)(mol)
end
function _chi3n(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi3n)(mol)
end
function _chi4n(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi4n)(mol)
end
function _chi1v(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi1v)(mol)
end
function _chi2v(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi2v)(mol)
end
function _chi3v(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi3v)(mol)
end
function _chi4v(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Chi4v)(mol)
end

# Kappa descriptors
function _kappa2(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Kappa2)(mol)
end
function _kappa3(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Kappa3)(mol)
end

# VSA descriptors
function _slogp_vsa2(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA2)(mol)
end
function _slogp_vsa3(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA3)(mol)
end
function _slogp_vsa4(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA4)(mol)
end
function _slogp_vsa5(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA5)(mol)
end
function _slogp_vsa6(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA6)(mol)
end
function _slogp_vsa7(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA7)(mol)
end
function _slogp_vsa8(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA8)(mol)
end
function _slogp_vsa9(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA9)(mol)
end
function _slogp_vsa10(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA10)(mol)
end
function _slogp_vsa11(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA11)(mol)
end
function _slogp_vsa12(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SlogP_VSA12)(mol)
end

# SMR VSA descriptors
function _smr_vsa1(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA1)(mol)
end
function _smr_vsa2(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA2)(mol)
end
function _smr_vsa3(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA3)(mol)
end
function _smr_vsa4(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA4)(mol)
end
function _smr_vsa5(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA5)(mol)
end
function _smr_vsa6(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA6)(mol)
end
function _smr_vsa7(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA7)(mol)
end
function _smr_vsa8(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA8)(mol)
end
function _smr_vsa9(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA9)(mol)
end
function _smr_vsa10(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").SMR_VSA10)(mol)
end

# PEOE VSA descriptors
function _peoe_vsa1(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA1)(mol)
end
function _peoe_vsa2(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA2)(mol)
end
function _peoe_vsa3(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA3)(mol)
end
function _peoe_vsa4(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA4)(mol)
end
function _peoe_vsa5(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA5)(mol)
end
function _peoe_vsa6(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA6)(mol)
end
function _peoe_vsa7(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA7)(mol)
end
function _peoe_vsa8(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA8)(mol)
end
function _peoe_vsa9(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA9)(mol)
end
function _peoe_vsa10(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA10)(mol)
end
function _peoe_vsa11(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA11)(mol)
end
function _peoe_vsa12(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA12)(mol)
end
function _peoe_vsa13(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA13)(mol)
end
function _peoe_vsa14(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").PEOE_VSA14)(mol)
end

# EState descriptors
function _max_e_state_index(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MaxEStateIndex)(mol)
end
function _min_e_state_index(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MinEStateIndex)(mol)
end
function _max_absolute_e_state_index(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MaxAbsEStateIndex)(mol)
end
function _min_absolute_e_state_index(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").MinAbsEStateIndex)(mol)
end

function _num_aliphatic_rings(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumAliphaticRings)(mol)
end
function _num_heterocycles(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").NumHeterocycles)(mol)
end
function _ipc(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").Ipc)(mol)
end

# Molecular complexity
function _balabanJ(mol::Py)
    @pyconst(pyimport("rdkit.Chem.Descriptors").BalabanJ)(mol)
end

# Additional utility functions
function _get_num_atoms(mol::Py)
    @pyconst(pyimport("rdkit.Chem").Mol.GetNumAtoms)(mol)
end
function _get_num_bonds(mol::Py)
    @pyconst(pyimport("rdkit.Chem").Mol.GetNumBonds)(mol)
end
function _get_num_conformers(mol::Py)
    @pyconst(pyimport("rdkit.Chem").Mol.GetNumConformers)(mol)
end

# Simple counts that are very useful - implemented using Python list comprehension for reliability
function _num_carbons(mol::Py)
    pyeval = @pyconst(pyimport("builtins").eval)
    code = "sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNumber() == 6)"
    return pyeval(code, Dict("mol" => mol))
end
function _num_nitrogens(mol::Py)
    pyeval = @pyconst(pyimport("builtins").eval)
    code = "sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNumber() == 7)"
    return pyeval(code, Dict("mol" => mol))
end
function _num_oxygens(mol::Py)
    pyeval = @pyconst(pyimport("builtins").eval)
    code = "sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNumber() == 8)"
    return pyeval(code, Dict("mol" => mol))
end
function _num_sulfurs(mol::Py)
    pyeval = @pyconst(pyimport("builtins").eval)
    code = "sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNumber() == 16)"
    return pyeval(code, Dict("mol" => mol))
end
function _num_halogens(mol::Py)
    pyeval = @pyconst(pyimport("builtins").eval)
    code = "sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNumber() in [9, 17, 35, 53, 85])"
    return pyeval(code, Dict("mol" => mol))
end

# Drawing functions
function _mol_to_image(mol::Py; kwargs...)
    @pyconst(pyimport("rdkit.Chem.Draw").MolToImage)(mol; kwargs...)
end
function _mols_to_grid_image(mols; kwargs...)
    @pyconst(pyimport("rdkit.Chem.Draw").MolsToGridImage)(mols; kwargs...)
end

# Chemical reactions
function _reaction_from_smarts(smarts::String)
    @pyconst(pyimport("rdkit.Chem.AllChem").ReactionFromSmarts)(smarts)
end

function _reaction_from_rxn_file(filename::String)
    @pyconst(pyimport("rdkit.Chem.AllChem").ReactionFromRxnFile)(filename)
end

function _reaction_from_rxn_block(rxnblock::String)
    @pyconst(pyimport("rdkit.Chem.AllChem").ReactionFromRxnBlock)(rxnblock)
end

function _reaction_to_rxn_block(rxn::Py)
    @pyconst(pyimport("rdkit.Chem.AllChem").ReactionToRxnBlock)(rxn)
end

function _reaction_validate(rxn::Py)
    @pyconst(pyimport("rdkit.Chem.AllChem").SanitizeRxn)(rxn)
end

function _reaction_get_num_reactant_templates(rxn::Py)
    rxn.GetNumReactantTemplates()
end

function _reaction_get_num_product_templates(rxn::Py)
    rxn.GetNumProductTemplates()
end

function _reaction_get_reactant_template(rxn::Py, idx::Int)
    rxn.GetReactantTemplate(idx)
end

function _reaction_get_product_template(rxn::Py, idx::Int)
    rxn.GetProductTemplate(idx)
end

function _reaction_has_reactant_substructure_match(rxn::Py, mol::Py)
    try
        products = rxn.RunReactants(pylist([mol]))
        return length(products) > 0
    catch
        return false
    end
end

function _reaction_get_reacting_atoms(rxn::Py)
    return rxn.GetReactingAtoms()
end

function _reaction_fingerprint(rxn::Py, fp_size::Int = 2048)
    rdkit_reactions = @pyconst(pyimport("rdkit.Chem.rdChemReactions"))
    params = rdkit_reactions.ReactionFingerprintParams()
    params.fpSize = fp_size
    return rdkit_reactions.CreateDifferenceFingerprintForReaction(rxn, params)
end

function _reaction_structural_fingerprint(rxn::Py, fp_size::Int = 2048)
    rdkit_reactions = @pyconst(pyimport("rdkit.Chem.rdChemReactions"))
    params = rdkit_reactions.ReactionFingerprintParams()
    params.fpSize = fp_size
    return rdkit_reactions.CreateStructuralFingerprintForReaction(rxn, params)
end

function _compute_reaction_center_fingerprint(rxn::Py, fp_size::Int = 2048)
    rdkit_reactions = @pyconst(pyimport("rdkit.Chem.rdChemReactions"))
    params = rdkit_reactions.ReactionFingerprintParams()
    params.fpSize = fp_size
    return rdkit_reactions.CreateDifferenceFingerprintForReaction(rxn, params)
end

function _reaction_run_reactants_inline_properties(
    rxn::Py, reactants::Vector{Py}, max_products::Int = 1000
)
    rxn.RunReactants(pylist(reactants), max_products)
end

function _reaction_enumerate_library_from_reaction(
    rxn::Py, reactant_lists::Vector{Vector{Py}}
)
    @pyconst(pyimport("rdkit.Chem.AllChem").EnumerateLibraryFromReaction)(
        rxn, pylist([pylist(r) for r in reactant_lists])
    )
end

function _reaction_to_smarts(rxn::Py)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").ReactionToSmarts)(rxn)
end

function _reaction_compute_atom_mapping(rxn::Py)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").ReduceProductToSideChains)(rxn)
end

function _reaction_sanitize_reaction(rxn::Py, sanitize_ops::Int = 15)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").SanitizeRxn)(rxn, sanitize_ops)
end

function _reaction_remove_unmapped_reactant_templates(rxn::Py, mode::Int = 1)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").RemoveUnmappedReactantTemplates)(
        rxn, mode
    )
end

function _reaction_remove_unmapped_product_templates(rxn::Py, mode::Int = 1)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").RemoveUnmappedProductTemplates)(
        rxn, mode
    )
end

function _reaction_preprocess(rxn::Py)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").PreprocessReaction)(rxn)
end

function _reaction_is_template_molecule_agent(mol::Py)
    @pyconst(pyimport("rdkit.Chem.rdChemReactions").IsTemplateMoleculeAgent)(mol)
end

function _get_atom_mapping_numbers(mol::Py)
    atom_map_nums = []
    num_atoms = pyconvert(Int, mol.GetNumAtoms())
    for i in 0:(num_atoms - 1)
        atom = mol.GetAtomWithIdx(i)
        map_num = pyconvert(Int, atom.GetAtomMapNum())
        push!(atom_map_nums, map_num)
    end
    return atom_map_nums
end

function _set_atom_mapping_numbers(mol::Py, map_nums::Vector{Int})
    for (i, map_num) in enumerate(map_nums)
        atom = mol.GetAtomWithIdx(i-1)
        atom.SetAtomMapNum(map_num)
    end
end

# Pharmacophore and Chemical Features
function _build_feature_factory(filename::String)
    @pyconst(pyimport("rdkit.Chem.ChemicalFeatures").BuildFeatureFactory)(filename)
end

function _build_feature_factory_from_string(fdef_string::String)
    @pyconst(pyimport("rdkit.Chem.ChemicalFeatures").BuildFeatureFactoryFromString)(
        fdef_string
    )
end

function _get_features_for_mol(factory::Py, mol::Py; conf_id::Int = -1)
    if conf_id == -1
        factory.GetFeaturesForMol(mol)
    else
        factory.GetFeaturesForMol(mol; confId = conf_id)
    end
end

function _get_feature_families(factory::Py)
    factory.GetFeatureFamilies()
end

function _get_feature_defs(factory::Py)
    factory.GetFeatureDefs()
end

function _get_num_feature_defs(factory::Py)
    factory.GetNumFeatureDefs()
end

function _mol_from_ph4(pharmacophore::Py)
    @pyconst(pyimport("rdkit.Chem.Pharm3D.EmbedLib").EmbedMol)(pharmacophore)
end

function _get_pharmacophore_fingerprint(mol::Py, factory::Py, sig_factory::Py)
    @pyconst(pyimport("rdkit.Chem.Pharm2D.Generate").Gen2DFingerprint)(mol, sig_factory)
end

function _create_sig_factory(
    factory::Py; min_point_count::Int = 2, max_point_count::Int = 3
)
    sig_factory = @pyconst(pyimport("rdkit.Chem.Pharm2D.SigFactory").SigFactory)(
        factory, min_point_count, max_point_count
    )
    # Use non-overlapping distance bins that avoid boundary issues
    sig_factory.SetBins([(0, 2), (2, 6), (6, 12)])
    sig_factory.Init()
    return sig_factory
end

function _get_rdconfig_data_dir()
    @pyconst(pyimport("rdkit.RDConfig").RDDataDir)
end

# 3D Pharmacophore functions
function _pharmacophore_from_mol(mol::Py, feature_factory::Py; conf_id::Int = -1)
    @pyconst(pyimport("rdkit.Chem.Pharm3D.Pharmacophore").Pharmacophore)()
end

function _explicit_pharmacophore_from_mol(mol::Py, feature_factory::Py; conf_id::Int = -1)
    features = _get_features_for_mol(feature_factory, mol; conf_id = conf_id)

    # Extract feature information
    feature_list = Tuple{String, Vector{Float64}}[]
    for feature in features
        family = pyconvert(String, feature.GetFamily())
        pos_obj = feature.GetPos()
        pos = [
            pyconvert(Float64, pos_obj.x),
            pyconvert(Float64, pos_obj.y),
            pyconvert(Float64, pos_obj.z),
        ]
        push!(feature_list, (family, pos))
    end

    return feature_list
end
