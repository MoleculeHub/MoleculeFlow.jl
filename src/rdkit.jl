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

# Substructure search
function _find_mcs(mol_list::Vector)
    @pyconst(pyimport("rdkit.Chem.rdFMCS").FindMCS)(mol_list)
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

# Drawing functions
function _mol_to_image(mol::Py; kwargs...)
    @pyconst(pyimport("rdkit.Chem.Draw").MolToImage)(mol; kwargs...)
end
function _mols_to_grid_image(mols; kwargs...)
    @pyconst(pyimport("rdkit.Chem.Draw").MolsToGridImage)(mols; kwargs...)
end
