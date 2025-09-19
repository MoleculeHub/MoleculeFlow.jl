#######################################################
# Substructure search functionality
#######################################################

"""
    has_substructure_match(mol::Molecule, pattern::Molecule) -> Union{Bool,Missing}
    has_substructure_match(mol::Molecule, pattern_smarts::String) -> Union{Bool,Missing}

Check if a molecule contains a substructure pattern.

# Arguments

  - `mol::Molecule`: The molecule to search in
  - `pattern::Molecule` or `pattern_smarts::String`: The substructure pattern as a Molecule or SMARTS string

# Returns

  - `Union{Bool,Missing}`: true if pattern is found, false otherwise, missing if molecules are invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
has_oh = has_substructure_match(mol, "[OH]")
has_benzene = has_substructure_match(mol, "c1ccccc1")  # false
```

# Notes

  - SMARTS patterns are more flexible than SMILES for substructure searching
  - Case-sensitive: 'c' = aromatic carbon, 'C' = aliphatic carbon
"""
function has_substructure_match(mol::Molecule, pattern::Molecule)
    !mol.valid && return missing
    !pattern.valid && return missing

    return pyconvert(Bool, mol._rdkit_mol.HasSubstructMatch(pattern._rdkit_mol))
end

function has_substructure_match(mol::Molecule, pattern_smarts::String)
    !mol.valid && return missing

    pattern_mol = _mol_from_smarts(pattern_smarts)
    if pynot(pattern_mol)
        throw(ArgumentError("Invalid SMARTS pattern: $pattern_smarts"))
    end

    return pyconvert(Bool, mol._rdkit_mol.HasSubstructMatch(pattern_mol))
end

# Get substructure matches
"""
    get_substructure_matches(mol::Molecule, pattern::Union{Molecule,String}; unique_matches::Bool=true) -> Union{Vector{Vector{Int}},Missing}

Find all substructure matches and return atom indices for each match.

# Arguments

  - `mol::Molecule`: The molecule to search in
  - `pattern::Union{Molecule,String}`: The substructure pattern
  - `unique_matches::Bool=true`: Whether to return only unique matches

# Returns

  - `Union{Vector{Vector{Int}},Missing}`: Vector of vectors containing 1-based atom indices for each match

# Examples

```julia
mol = mol_from_smiles("c1ccc(O)cc1")  # Phenol
matches = get_substructure_matches(mol, "[OH]")
# Returns atom indices where OH pattern is found
```

# Notes

  - Returns 1-based atom indices (Julia convention)
  - Each inner vector contains atom indices for one match
  - Empty vector if no matches found
"""
function get_substructure_matches(
    mol::Molecule, pattern::Molecule; unique_matches::Bool = true
)
    !mol.valid && return missing
    !pattern.valid && return missing

    matches = mol._rdkit_mol.GetSubstructMatches(
        pattern._rdkit_mol; uniquify = unique_matches
    )
    # Convert to Julia arrays and 1-based indexing
    return [pyconvert(Vector{Int}, match) .+ 1 for match in matches]
end

function get_substructure_matches(
    mol::Molecule, pattern_smarts::String; unique_matches::Bool = true
)
    !mol.valid && return missing

    pattern_mol = _mol_from_smarts(pattern_smarts)
    if pynot(pattern_mol)
        throw(ArgumentError("Invalid SMARTS pattern: $pattern_smarts"))
    end

    matches = mol._rdkit_mol.GetSubstructMatches(pattern_mol; uniquify = unique_matches)
    # Convert to Julia arrays and 1-based indexing
    return [pyconvert(Vector{Int}, match) .+ 1 for match in matches]
end

"""
    get_substructure_match(mol::Molecule, pattern::Union{Molecule,String}) -> Union{Vector{Int}, Missing}

Get the first substructure match of a pattern in a molecule.

# Arguments

  - `mol::Molecule`: The molecule to search in
  - `pattern::Union{Molecule,String}`: The substructure pattern as a Molecule or SMARTS string

# Returns

  - `Vector{Int}`: 1-based atom indices of the first match
  - `missing`: If no match is found or molecules are invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
match = get_substructure_match(mol, "[OH]")  # Returns [3] (oxygen index)
```

# Notes

  - Returns only the first match found
  - Use `get_substructure_matches` to get all matches
"""
function get_substructure_match(mol::Molecule, pattern::Molecule)
    !mol.valid && return missing
    !pattern.valid && return missing

    match = mol._rdkit_mol.GetSubstructMatch(pattern._rdkit_mol)
    if length(match) == 0
        return missing
    end
    # Convert to Julia array and 1-based indexing
    return pyconvert(Vector{Int}, match) .+ 1
end

function get_substructure_match(mol::Molecule, pattern_smarts::String)
    !mol.valid && return missing

    pattern_mol = _mol_from_smarts(pattern_smarts)
    if pynot(pattern_mol)
        throw(ArgumentError("Invalid SMARTS pattern: $pattern_smarts"))
    end

    match = mol._rdkit_mol.GetSubstructMatch(pattern_mol)
    if length(match) == 0
        return missing
    end
    # Convert to Julia array and 1-based indexing
    return pyconvert(Vector{Int}, match) .+ 1
end

"""
    maximum_common_substructure(mol1::Molecule, mol2::Molecule) -> Union{Molecule, Missing}

Find the maximum common substructure (MCS) between two molecules.

# Arguments

  - `mol1::Molecule`: First molecule
  - `mol2::Molecule`: Second molecule

# Returns

  - `Molecule`: The maximum common substructure as a new Molecule object
  - `missing`: If no common substructure is found or molecules are invalid

# Examples

```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
mcs = maximum_common_substructure(mol1, mol2)  # Common C-C substructure
```

# Notes

  - Uses RDKit's FindMCS algorithm
  - Returns the largest substructure common to both molecules
  - Useful for scaffold analysis and lead optimization
"""
function maximum_common_substructure(mol1::Molecule, mol2::Molecule)
    !mol1.valid && return missing
    !mol2.valid && return missing

    # Find MCS
    mcs = _find_mcs([mol1._rdkit_mol, mol2._rdkit_mol])

    if pyconvert(Int, mcs.numAtoms) == 0
        return missing
    end

    mcs_mol = _mol_from_smarts(pyconvert(String, mcs.smartsString))
    if pynot(mcs_mol)
        return missing
    end

    return Molecule(;
        _rdkit_mol = mcs_mol, valid = true, source = pyconvert(String, mcs.smartsString)
    )
end

"""
    has_substructure_matches(mols::Vector{Union{Molecule,Missing}}, pattern::Union{Molecule,String}) -> Vector{Union{Bool,Missing}}

Check for substructure matches across a vector of molecules.

# Arguments

  - `mols::Vector{Union{Molecule,Missing}}`: Vector of molecules to search
  - `pattern::Union{Molecule,String}`: The substructure pattern to search for

# Returns

  - `Vector{Union{Bool,Missing}}`: Boolean vector indicating matches for each molecule

# Examples

```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("CCC"), mol_from_smiles("CC(O)C")]
matches = has_substructure_matches(mols, "[OH]")  # [true, false, true]
```

# Notes

  - Vectorized version of `has_substructure_match`
  - Useful for filtering large molecular datasets
"""
function has_substructure_matches(
    mols::Vector{Union{Molecule, Missing}}, pattern::Union{Molecule, String}
)
    return [has_substructure_match(mol, pattern) for mol in mols]
end

"""
    filter_by_substructure(mols::Vector{Union{Molecule,Missing}}, pattern::Union{Molecule,String}) -> Vector{Union{Molecule,Missing}}

Filter molecules that contain a specific substructure pattern.

# Arguments

  - `mols::Vector{Union{Molecule,Missing}}`: Vector of molecules to filter
  - `pattern::Union{Molecule,String}`: The substructure pattern to filter by

# Returns

  - `Vector{Union{Molecule,Missing}}`: Filtered vector containing only molecules with the pattern

# Examples

```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("CCC"), mol_from_smiles("CC(O)C")]
alcohols = filter_by_substructure(mols, "[OH]")  # Returns [CCO, CC(O)C]
```

# Notes

  - Convenient wrapper around `has_substructure_matches`
  - Useful for creating focused molecular libraries
"""
function filter_by_substructure(
    mols::Vector{Union{Molecule, Missing}}, pattern::Union{Molecule, String}
)
    mask = has_substructure_matches(mols, pattern)
    return mols[mask .== true]
end

"""
    FUNCTIONAL_GROUPS

A comprehensive dictionary containing 100+ predefined functional group SMARTS patterns organized by category.

# Basic Functional Groups

  - `:alcohol`, `:phenol`, `:carboxylic_acid`, `:ester`, `:ether`, `:aldehyde`, `:ketone`
  - `:amine_primary`, `:amine_secondary`, `:amine_tertiary`, `:amide`, `:nitrile`

# Sulfur-Containing Groups

  - `:thiol`, `:sulfide`, `:disulfide`, `:sulfoxide`, `:sulfone`
  - `:sulfonamide`, `:sulfonate`, `:sulfonic_acid`, `:sulfonyl_chloride`

# Phosphorus-Containing Groups

  - `:phosphate`, `:phosphonate`, `:phosphine`, `:phosphine_oxide`

# Halogen-Containing Groups

  - `:fluoride`, `:chloride`, `:bromide`, `:iodide`, `:trifluoromethyl`, `:trichloromethyl`

# Advanced Nitrogen Groups

  - `:nitro`, `:nitroso`, `:azide`, `:diazo`, `:hydrazine`, `:hydroxylamine`
  - `:imine`, `:oxime`, `:enamine`, `:guanidine`, `:urea`, `:carbamate`
  - `:isocyanate`, `:isothiocyanate`

# Advanced Oxygen Groups

  - `:peroxide`, `:acetal`, `:ketal`, `:hemiacetal`, `:hemiketal`
  - `:anhydride`, `:carbonate`, `:carbamate_ester`

# Carbon-Carbon Multiple Bonds

  - `:alkene`, `:alkyne`, `:allene`, `:conjugated_diene`

# 5-Membered Aromatic Heterocycles

  - `:furan`, `:thiophene`, `:pyrrole`, `:imidazole`, `:pyrazole`
  - `:oxazole`, `:isoxazole`, `:thiazole`, `:isothiazole`
  - `:triazole_1_2_3`, `:triazole_1_2_4`, `:tetrazole`

# 6-Membered Aromatic Heterocycles

  - `:benzene`, `:pyridine`, `:pyrimidine`, `:pyrazine`, `:pyridazine`, `:triazine`

# Fused Aromatic Systems

  - `:naphthalene`, `:anthracene`, `:phenanthrene`, `:quinoline`, `:isoquinoline`
  - `:indole`, `:benzofuran`, `:benzothiophene`, `:purine`

# Saturated Heterocycles

  - `:tetrahydrofuran`, `:tetrahydropyran`, `:pyrrolidine`, `:piperidine`
  - `:morpholine`, `:piperazine`, `:azetidine`, `:oxetane`, `:thietane`

# Biomolecule Patterns

  - `:glucose_like`, `:anomeric_carbon`, `:glycosidic_bond` (sugars)
  - `:peptide_bond`, `:alpha_amino_acid`, `:proline_like` (proteins)
  - `:fatty_acid`, `:fatty_acid_long`, `:triglyceride_like` (lipids)

# Pharmaceutical Patterns

  - `:benzodiazepine_core`, `:beta_lactam`, `:sulfonamide_drug`, `:barbiturate_core`

# Natural Product Patterns

  - `:steroid_core`, `:flavonoid_core`, `:coumarin`, `:chromone`

# Reactive Groups and Electrophiles

  - `:epoxide`, `:aziridine`, `:cyclopropane`, `:michael_acceptor`
  - `:alpha_beta_unsaturated_carbonyl`

# Protecting Groups (Synthetic Chemistry)

  - `:tert_butyl`, `:benzyl`, `:acetyl`, `:benzoyl`, `:tosyl`
  - `:boc`, `:cbz`, `:fmoc`

# Examples

```julia
mol = mol_from_smiles("CCO")
has_functional_group(mol, :alcohol)  # true
# Basic functional groups
drug = mol_from_smiles("CC(=O)Nc1ccc(O)cc1")  # Acetaminophen
has_functional_group(drug, :phenol)     # true
has_functional_group(drug, :amide)      # true
has_functional_group(drug, :benzene)    # true

# Get all functional groups at once
all_groups = get_functional_groups(drug)
```

# Notes

  - Contains 100+ functional group patterns covering basic to advanced organic chemistry
  - SMARTS patterns are carefully designed for specificity and broad applicability
  - Useful for drug discovery, natural product analysis, and chemical library filtering    # Sulfur-containing groups
  - Patterns organized by chemical similarity and complexity
"""
const FUNCTIONAL_GROUPS = Dict{Symbol, String}(
    # Basic functional groups
    :alcohol => "[OH1]",
    :phenol => "[OH1][c]",
    :carboxylic_acid => "[CX3](=O)[OX2H1]",
    :ester => "[#6][CX3](=O)[OX2H0][#6]",
    :ether => "[OD2]([#6])[#6]",
    :aldehyde => "[CX3H1](=O)[#6]",
    :ketone => "[#6][CX3](=O)[#6]",
    :amine_primary => "[NX3;H2;!\$(NC=O)]",
    :amine_secondary => "[NX3;H1;!\$(NC=O)]",
    :amine_tertiary => "[NX3;H0;!\$(NC=O)]",
    :amide => "[NX3][CX3](=[OX1])[#6]",
    :nitrile => "[NX1]#[CX2]",

    # Sulfur-containing groups
    :thiol => "[SH1]",
    :sulfide => "[SD2]([#6])[#6]",
    :disulfide => "[#6][SD2][SD2][#6]",
    :sulfoxide => "[#6][SX3](=O)[#6]",
    :sulfone => "[#6][SX4](=O)(=O)[#6]",
    :sulfonamide => "[SX4](=O)(=O)[NX3]",
    :sulfonate => "[SX4](=O)(=O)[O-,OH]",
    :sulfonic_acid => "[SX4](=O)(=O)[OH]",
    :sulfonyl_chloride => "[SX4](=O)(=O)[Cl]",

    # Phosphorus-containing groups
    :phosphate => "[PX4](=O)([OH,O-])([OH,O-])[OH,O-]",
    :phosphonate => "[PX4](=O)([OH,O-])[#6]",
    :phosphine => "[PX3]([#6])([#6])[#6]",
    :phosphine_oxide => "[PX4](=O)([#6])([#6])[#6]",

    # Halogen-containing groups
    :fluoride => "[F]",
    :chloride => "[Cl]",
    :bromide => "[Br]",
    :iodide => "[I]",
    :trifluoromethyl => "[CX4]([F])([F])[F]",
    :trichloromethyl => "[CX4]([Cl])([Cl])[Cl]",

    # Nitrogen-containing advanced groups
    :nitro => "[NX3+](=O)[O-]",
    :nitroso => "[NX2](=O)",
    :azide => "[NX2-][NX2+]#[NX1]",
    :diazo => "[#6][NX2+]#[NX1-]",
    :hydrazine => "[NX3][NX3]",
    :hydroxylamine => "[NX3][OH]",
    :imine => "[CX3]=[NX2]",
    :oxime => "[CX3]=[NX2][OH]",
    :enamine => "[NX3][CX3]=[CX3]",
    :guanidine => "[NX3][CX3](=[NX3+])[NX3]",
    :urea => "[NX3][CX3](=O)[NX3]",
    :carbamate => "[NX3][CX3](=O)[OX2]",
    :isocyanate => "[NX2]=[CX2]=[OX1]",
    :isothiocyanate => "[NX2]=[CX2]=[SX1]",

    # Oxygen-containing advanced groups
    :peroxide => "[OX2][OX2]",
    :acetal => "[CX4]([OX2])([OX2])[#6]",
    :ketal => "[CX4]([OX2])([OX2])([#6])[#6]",
    :hemiacetal => "[CX4]([OX2H])([OX2])[#6]",
    :hemiketal => "[CX4]([OX2H])([OX2])([#6])[#6]",
    :anhydride => "[CX3](=O)[OX2][CX3](=O)",
    :carbonate => "[OX2][CX3](=O)[OX2]",
    :carbamate_ester => "[OX2][CX3](=O)[NX3]",

    # Carbon-carbon multiple bonds
    :alkene => "[CX3]=[CX3]",
    :alkyne => "[CX2]#[CX2]",
    :allene => "[CX3]=[CX2]=[CX3]",
    :conjugated_diene => "[CX3]=[CX3][CX3]=[CX3]",

    # Aromatic heterocycles (5-membered)
    :furan => "c1ccoc1",
    :thiophene => "c1ccsc1",
    :pyrrole => "c1cc[nH]c1",
    :imidazole => "c1cnc[nH]1",
    :pyrazole => "c1cn[nH]c1",
    :oxazole => "c1cnoc1",
    :isoxazole => "c1noc[cH]1",
    :thiazole => "c1cnsc1",
    :isothiazole => "c1nsc[cH]1",
    :triazole_1_2_3 => "c1nn[nH]c1",
    :triazole_1_2_4 => "c1ncn[nH]1",
    :tetrazole => "c1nnn[nH]1",

    # Aromatic heterocycles (6-membered)
    :benzene => "c1ccccc1",
    :pyridine => "c1ccncc1",
    :pyrimidine => "c1cncnc1",
    :pyrazine => "c1cnccn1",
    :pyridazine => "c1ccnnc1",
    :triazine => "c1ncncn1",

    # Fused aromatic systems
    :naphthalene => "c1ccc2ccccc2c1",
    :anthracene => "c1ccc2cc3ccccc3cc2c1",
    :phenanthrene => "c1ccc2c(c1)ccc3ccccc32",
    :quinoline => "c1ccc2ncccc2c1",
    :isoquinoline => "c1cnc2ccccc2c1",
    :indole => "c1ccc2[nH]ccc2c1",
    :benzofuran => "c1ccc2occc2c1",
    :benzothiophene => "c1ccc2sccc2c1",
    :purine => "c1nc2c([nH]1)ncn2",
    :pyrimidine_fused => "c1nc2nccnc2[nH]1",

    # Saturated heterocycles
    :tetrahydrofuran => "C1CCOC1",
    :tetrahydropyran => "C1CCOCC1",
    :pyrrolidine => "C1CCN[CH2]1",
    :piperidine => "C1CCNCC1",
    :morpholine => "C1COCCN1",
    :piperazine => "C1CNCCN1",
    :azetidine => "C1CNC1",
    :oxetane => "C1COC1",
    :thietane => "C1CSC1",

    # Sugar and carbohydrate patterns
    :glucose_like => "[CH1]([OH])[CH1]([OH])[CH1]([OH])",
    :anomeric_carbon => "[CH1]([OH])[OX2]",
    :glycosidic_bond => "[CH1][OX2][CH1]",

    # Peptide and protein patterns
    :peptide_bond => "[NX3][CX3](=O)",
    :alpha_amino_acid => "[NX3][CH1]([CX3](=O))[*]",
    :proline_like => "N1[CH2][CH2][CH2][CH1]1[CX3](=O)",

    # Lipid patterns
    :fatty_acid => "[CH3][CH2][CH2][CX3](=O)[OH]",
    :fatty_acid_long => "[CH3][CH2][CH2][CH2][CH2][CH2][CX3](=O)[OH]",
    :triglyceride_like => "[CH2][OX2][CX3](=O)[CH2]",

    # Pharmaceutical patterns
    :benzodiazepine_core => "c1ccc2c(c1)[nH]c(=O)c3ccccc3n2",
    :beta_lactam => "[CH2][CX4][NX3][CX3]=O",
    :sulfonamide_drug => "c1ccc(cc1)[SX4](=O)(=O)[NX3]",
    :barbiturate_core => "[NX3][CX3](=O)[NX3][CX3](=O)",

    # Natural product patterns
    :steroid_core => "C1CC2CCC3C(CCC4CCCCC34)C2CC1",
    :flavonoid_core => "c1ccc(cc1)[CX3]2[OX2]c3ccccc3[CX3](=O)[CH2]2",
    :coumarin => "c1ccc2oc(=O)ccc2c1",
    :chromone => "c1ccc2c(c1)oc(=O)cc2",

    # Reactive intermediates and electrophiles
    :epoxide => "C1OC1",
    :aziridine => "C1NC1",
    :cyclopropane => "C1CC1",
    :michael_acceptor => "[CX3]=[CX3][CX3]=O",
    :alpha_beta_unsaturated_carbonyl => "[CX3]=[CX3][CX3]=O",

    # Protecting groups (common in synthesis)
    :tert_butyl => "[CH3]C([CH3])([CH3])",
    :benzyl => "[CH2]c1ccccc1",
    :acetyl => "[CH3][CX3]=O",
    :benzoyl => "c1ccccc1[CX3]=O",
    :tosyl => "c1ccc(cc1)[CH3][SX4](=O)(=O)",
    :boc => "[CH3]C([CH3])([CH3])[OX2][CX3](=O)[NX3]",
    :cbz => "c1ccccc1[CH2][OX2][CX3](=O)[NX3]",
    :fmoc => "c1ccc2c(c1)c(c3ccccc3c2)[CH2][OX2][CX3](=O)[NX3]",
)

"""
    has_functional_group(mol::Molecule, group::Symbol) -> Union{Bool, Missing}

Check if a molecule contains a specific functional group.

# Arguments

  - `mol::Molecule`: The molecule to analyze
  - `group::Symbol`: The functional group to search for (see `FUNCTIONAL_GROUPS` for available groups)

# Returns

  - `Union{Bool, Missing}`: true if the functional group is present, false otherwise, missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
has_functional_group(mol, :alcohol)  # true
has_functional_group(mol, :ketone)   # false
```

# Available Functional Groups

See `FUNCTIONAL_GROUPS` constant for the complete list of 100+ available functional groups including:

  - `:alcohol`, `:carboxylic_acid`, `:ester`, `:ether`, `:aldehyde`, `:ketone`
  - `:amine_primary`, `:amine_secondary`, `:amine_tertiary`, `:amide`, `:nitrile`
  - `:benzene`, `:pyridine`, `:furan`, `:thiophene`, `:imidazole`, `:pyrrole`
  - `:thiol`, `:sulfide`, `:disulfide`, `:sulfoxide`, `:sulfone`, `:sulfonamide`
  - `:phosphate`, `:phosphonate`, `:phosphine`, `:phosphine_oxide`
  - `:fluoride`, `:chloride`, `:bromide`, `:iodide`, `:trifluoromethyl`
  - `:nitro`, `:nitroso`, `:azide`, `:epoxide`, `:steroid_core`, `:beta_lactam`
  - `:morpholine`, `:piperidine`, `:pyrrolidine`, `:tetrahydrofuran`
  - `:peptide_bond`, `:alpha_amino_acid`, `:glucose_like`, `:fatty_acid`
  - `:benzodiazepine_core`, `:barbiturate_core`, `:sulfonamide_drug`
"""
function has_functional_group(mol::Molecule, group::Symbol)
    if haskey(FUNCTIONAL_GROUPS, group)
        return has_substructure_match(mol, FUNCTIONAL_GROUPS[group])
    else
        throw(
            ArgumentError(
                "Unknown functional group: $group. Available groups: $(keys(FUNCTIONAL_GROUPS))",
            ),
        )
    end
end

"""
    get_functional_groups(mol::Molecule) -> Union{Dict{Symbol, Bool}, Missing}

Get a dictionary of all functional groups present in a molecule.

# Arguments

  - `mol::Molecule`: The molecule to analyze

# Returns

  - `Dict{Symbol, Bool}`: Dictionary mapping functional group symbols to their presence (true/false)
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CC(=O)O")  # Acetic acid
groups = get_functional_groups(mol)
# Returns Dict(:carboxylic_acid => true, :alcohol => false, ...)
```

# Notes

  - Checks for all functional groups defined in `FUNCTIONAL_GROUPS`
  - Useful for quick functional group profiling
"""
function get_functional_groups(mol::Molecule)
    results = Dict{Symbol, Bool}()
    for (group, pattern) in FUNCTIONAL_GROUPS
        results[group] = has_substructure_match(mol, pattern)
    end
    return results
end

"""
    get_ring_info(mol::Molecule) -> Union{Dict, Missing}

Get detailed information about rings in a molecule.

# Arguments

  - `mol::Molecule`: The molecule to analyze

# Returns

  - `Dict`: Dictionary containing ring information with keys:

      + `:num_rings`: Total number of rings
      + `:atom_rings`: Vector of vectors containing 1-based atom indices for each ring
      + `:bond_rings`: Vector of vectors containing 1-based bond indices for each ring

  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
ring_info = get_ring_info(mol)
# Returns Dict(:num_rings => 1, :atom_rings => [[1,2,3,4,5,6]], :bond_rings => [[1,2,3,4,5,6]])
```

# Notes

  - Provides comprehensive ring analysis
  - Useful for understanding molecular topology
"""
function get_ring_info(mol::Molecule)
    !mol.valid && return missing

    ring_info = mol._rdkit_mol.GetRingInfo()

    return Dict(
        :num_rings => pyconvert(Int, ring_info.NumRings()),
        :atom_rings =>
            [pyconvert(Vector{Int}, ring) .+ 1 for ring in ring_info.AtomRings()],
        :bond_rings =>
            [pyconvert(Vector{Int}, ring) .+ 1 for ring in ring_info.BondRings()],
    )
end

"""
    is_ring_aromatic(mol::Molecule, ring_atoms::Vector{Int}) -> Union{Bool, Missing}

Check if a specific ring in a molecule is aromatic.

# Arguments

  - `mol::Molecule`: The molecule containing the ring
  - `ring_atoms::Vector{Int}`: 1-based indices of atoms forming the ring

# Returns

  - `Bool`: true if the ring is aromatic, false otherwise
  - `missing`: If molecule is invalid

# Examples

```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
ring_info = get_ring_info(mol)
is_aromatic = is_ring_aromatic(mol, ring_info[:atom_rings][1])  # true
```

# Notes

  - Checks if all atoms in the ring have aromatic character
  - Useful for distinguishing aromatic from aliphatic rings
"""
function is_ring_aromatic(mol::Molecule, ring_atoms::Vector{Int})
    !mol.valid && return missing

    # Convert to 0-based indexing
    ring_atoms_0based = ring_atoms .- 1

    # Check if all atoms in the ring are aromatic
    all_aromatic = true
    for atom_idx in ring_atoms_0based
        atom = mol._rdkit_mol.GetAtomWithIdx(atom_idx)
        if !pyconvert(Bool, atom.GetIsAromatic())
            all_aromatic = false
            break
        end
    end

    return all_aromatic
end
