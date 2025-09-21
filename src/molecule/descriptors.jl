#######################################################
# Molecular descriptors
#######################################################

# Basic molecular properties
"""
    molecular_weight(mol::Molecule) -> Union{Float64,Missing}

Calculate the molecular weight of a molecule in Daltons (g/mol).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Molecular weight in g/mol, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
mw = molecular_weight(mol)    # ≈ 46.07 g/mol
```
"""
function molecular_weight(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _mol_wt(mol._rdkit_mol))
end

"""
    exact_molecular_weight(mol::Molecule) -> Union{Float64,Missing}

Calculate the exact molecular weight using isotopic masses.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Exact molecular weight, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
exact_mw = exact_molecular_weight(mol)  # More precise than molecular_weight
```

# Notes

  - Uses exact isotopic masses rather than average atomic weights
  - More precise than `molecular_weight` for mass spectrometry applications
"""
function exact_molecular_weight(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _exact_mol_wt(mol._rdkit_mol))
end

"""
    heavy_atom_count(mol::Molecule) -> Union{Int,Missing}

Count the number of heavy atoms (non-hydrogen atoms) in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of heavy atoms, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
count = heavy_atom_count(mol)  # 3 (two carbons and one oxygen)
```

# Notes

  - Heavy atoms include all atoms except hydrogen
  - Useful for drug-like property calculations
"""
function heavy_atom_count(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _heavy_atom_count(mol._rdkit_mol))
end

"""
    num_heteroatoms(mol::Molecule) -> Union{Int,Missing}

Count the number of heteroatoms (non-carbon heavy atoms) in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of heteroatoms, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
heteroatoms = num_heteroatoms(mol)  # 1 (oxygen)
```

# Notes

  - Heteroatoms include N, O, S, P, halogens, etc. (everything except C and H)
  - Important for drug-like property calculations
"""
function num_heteroatoms(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_heteroatoms(mol._rdkit_mol))
end

"""
    num_rotatable_bonds(mol::Molecule) -> Union{Int,Missing}

Count the number of rotatable bonds in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of rotatable bonds, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCCC")  # Butane
rot_bonds = num_rotatable_bonds(mol)  # 3
```

# Notes

  - Rotatable bonds are single bonds that can freely rotate
  - Excludes bonds in rings and bonds to terminal atoms
  - Important for molecular flexibility and drug-like properties
"""
function num_rotatable_bonds(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_rotatable_bonds(mol._rdkit_mol))
end

"""
    num_hbd(mol::Molecule) -> Union{Int,Missing}

Count the number of hydrogen bond donors in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of hydrogen bond donors, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
hbd = num_hbd(mol)  # 1 (the OH group)
```

# Notes

# Hydrogen bond donors

  - Used in Lipinski's Rule of Five (≤5 donors)
  - Important for drug-like properties  # Hydrogen bond donors
"""
function num_hbd(mol::Molecule)  # Hydrogen bond donors
    !mol.valid && return missing
    return pyconvert(Int, _num_hdonors(mol._rdkit_mol))
end

"""
    num_hba(mol::Molecule) -> Union{Int,Missing}

Count the number of hydrogen bond acceptors in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of hydrogen bond acceptors, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
hba = num_hba(mol)  # 1 (the oxygen)
```

# Notes

# Hydrogen bond acceptors

  - Used in Lipinski's Rule of Five (≤10 acceptors)
  - Important for drug-like properties  # Hydrogen bond acceptors
"""
function num_hba(mol::Molecule)  # Hydrogen bond acceptors
    !mol.valid && return missing
    return pyconvert(Int, _num_hacceptors(mol._rdkit_mol))
end

# Lipinski descriptors
"""
    logp(mol::Molecule) -> Union{Float64,Missing}

Calculate the octanol-water partition coefficient (LogP) using Crippen's method.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: LogP value, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
lp = logp(mol)  # Approximately -0.31
```

# Notes

  - LogP measures lipophilicity (fat-loving vs water-loving)
  - Positive values indicate lipophilic molecules
  - Negative values indicate hydrophilic molecules
  - Important for drug absorption and distribution
"""
function logp(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _mol_logp(mol._rdkit_mol))
end

"""
    tpsa(mol::Molecule) -> Union{Float64,Missing}

Calculate the Topological Polar Surface Area (TPSA) in Ų.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: TPSA in Ų, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
area = tpsa(mol)  # Approximately 20.23 Ų
```

# Notes

  - TPSA is the polar surface area based on fragment contributions
  - Important for predicting drug permeability  # Topological polar surface area
  - Values < 60 Ų typically indicate good oral bioavailability
  - Used in Lipinski's Rule of Five
"""
function tpsa(mol::Molecule)  # Topological polar surface area
    !mol.valid && return missing
    return pyconvert(Float64, _tpsa(mol._rdkit_mol))
end

"""
    slogp_vsa(mol::Molecule) -> Union{Float64,Missing}

Calculate the SlogP_VSA1 descriptor (MOE-type descriptor).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: SlogP_VSA1 value, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
vsa = slogp_vsa(mol)
```

# Notes

  - Part of the MOE-type descriptors family
  - Combines surface area and LogP information
  - Used in QSAR modeling
"""
function slogp_vsa(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _slogp_vsa1(mol._rdkit_mol))
end

"""
    num_rings(mol::Molecule) -> Union{Int,Missing}

Count the total number of rings in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of rings, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("c1ccc2ccccc2c1")  # Naphthalene
rings = num_rings(mol)  # 2
```

# Notes

  - Counts all ring systems (aromatic and aliphatic)
  - Important for drug-like property calculations
"""
function num_rings(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _ring_count(mol._rdkit_mol))
end

"""
    num_aromatic_rings(mol::Molecule) -> Union{Int,Missing}

Count the number of aromatic rings in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of aromatic rings, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
aromatic_rings = num_aromatic_rings(mol)  # 1
```

# Notes

  - Only counts rings with aromatic character
  - Important for drug design and π-π interactions
"""
function num_aromatic_rings(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_aromatic_rings(mol._rdkit_mol))
end

"""
    num_saturated_rings(mol::Molecule) -> Union{Int,Missing}

Count the number of saturated (aliphatic) rings in a molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of saturated rings, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("C1CCCCC1")  # Cyclohexane
saturated_rings = num_saturated_rings(mol)  # 1
```

# Notes

  - Only counts rings without aromatic character
  - Includes cycloalkanes and saturated heterocycles
"""
function num_saturated_rings(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_saturated_rings(mol._rdkit_mol))
end

"""
    bertz_ct(mol::Molecule) -> Union{Float64,Missing}

Calculate the BertzCT molecular complexity descriptor.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: BertzCT complexity score, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
complexity = bertz_ct(mol)
```

# Notes

  - Measures molecular complexity based on graph theory  # BertzCT complexity
  - Higher values indicate more complex molecular structures
  - Useful for drug design and synthesis planning
"""
function bertz_ct(mol::Molecule)  # BertzCT complexity
    !mol.valid && return missing
    return pyconvert(Float64, _bertz_ct(mol._rdkit_mol))
end

"""
    balaban_j(mol::Molecule) -> Union{Float64,Missing}

Calculate the Balaban J topological index.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Balaban J index, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCCC")
balaban = balaban_j(mol)
```

# Notes

  - Balaban J index is a topological descriptor
  - Measures molecular branching and connectivity
  - Used in QSAR studies
"""
function balaban_j(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _balaban_j(mol._rdkit_mol))
end

"""
    chi0v(mol::Molecule) -> Union{Float64,Missing}

Calculate the Chi0v valence connectivity index.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Chi0v connectivity index, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
chi = chi0v(mol)
```

# Notes

  - Chi0v is a valence connectivity index  # Connectivity index
  - Describes molecular connectivity considering valence electrons
  - Useful for predicting molecular properties
"""
function chi0v(mol::Molecule)  # Connectivity index
    !mol.valid && return missing
    return pyconvert(Float64, _chi0v(mol._rdkit_mol))
end

"""
    kappa1(mol::Molecule) -> Union{Float64,Missing}

Calculate the Kappa1 molecular shape index.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Kappa1 shape index, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCCC")  # Linear molecule
kappa = kappa1(mol)  # Higher for linear molecules
```

# Notes

  - Kappa1 describes molecular shape and branching  # Kappa shape index
  - Higher values indicate more linear structures
  - Part of the Kier and Hall molecular shape indices
"""
function kappa1(mol::Molecule)  # Kappa shape index
    !mol.valid && return missing
    return pyconvert(Float64, _kappa1(mol._rdkit_mol))
end

#######################################################
# Advanced Drug-like and ADMET Properties
#######################################################

"""
    qed(mol::Molecule) -> Union{Float64,Missing}

Calculate the Quantitative Estimate of Drug-likeness (QED) score.

QED combines multiple molecular properties into a single drug-likeness score
ranging from 0 (non-drug-like) to 1 (highly drug-like).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: QED score (0-1), or missing if molecule is invalid

# Examples

```julia
aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
drug_score = qed(aspirin)  # ≈ 0.73 (fairly drug-like)
```

# Notes

  - Based on desirability functions for MW, LogP, HBD, HBA, PSA, rotatable bonds, aromatic rings, and alerts
  - Higher scores indicate more drug-like properties
  - Commonly used threshold: QED > 0.5 for drug-like compounds
"""
function qed(mol::Molecule)
    !mol.valid && return missing
    try
        return pyconvert(Float64, _qed(mol._rdkit_mol))
    catch
        return missing
    end
end

"""
    fraction_csp3(mol::Molecule) -> Union{Float64,Missing}

Calculate the fraction of sp3 hybridized carbons.

Higher sp3 fraction correlates with increased drug-likeness and 3D character.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Fraction of sp3 carbons (0-1), or missing if molecule is invalid

# Examples

```julia
cyclohexane = mol_from_smiles("C1CCCCC1")
fsp3 = fraction_csp3(cyclohexane)  # 1.0 (all carbons are sp3)

benzene = mol_from_smiles("c1ccccc1")
fsp3_benzene = fraction_csp3(benzene)  # 0.0 (all carbons are sp2)
```

# Notes

  - Values range from 0 (fully aromatic/planar) to 1 (fully saturated)
  - Drug-like compounds typically have Fsp3 > 0.25
  - Important for assessing molecular complexity and 3D character
"""
function fraction_csp3(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _fraction_csp3(mol._rdkit_mol))
end

"""
    labute_asa(mol::Molecule) -> Union{Float64,Missing}

Calculate the Labute Accessible Surface Area.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Accessible surface area in Ų, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
asa = labute_asa(mol)  # Accessible surface area
```

# Notes

  - Estimates the solvent-accessible surface area
  - Important for understanding molecular size and shape
  - Correlates with solubility and membrane permeability
"""
function labute_asa(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _labute_asa(mol._rdkit_mol))
end

"""
    molar_refractivity(mol::Molecule) -> Union{Float64,Missing}

Calculate the molar refractivity.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: Molar refractivity, or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
mr = molar_refractivity(mol)  # Molar refractivity value
```

# Notes

  - Related to polarizability and molecular volume
  - Important for QSAR modeling
  - Correlates with London dispersion forces
"""
function molar_refractivity(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Float64, _mol_mr(mol._rdkit_mol))
end

"""
    num_aliphatic_carbocycles(mol::Molecule) -> Union{Int,Missing}

Count the number of aliphatic carbocycles (saturated carbon-only rings).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of aliphatic carbocycles, or missing if molecule is invalid

# Examples

```julia
cyclohexane = mol_from_smiles("C1CCCCC1")
count = num_aliphatic_carbocycles(cyclohexane)  # 1
```
"""
function num_aliphatic_carbocycles(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_aliphatic_carbocycles(mol._rdkit_mol))
end

"""
    num_aromatic_carbocycles(mol::Molecule) -> Union{Int,Missing}

Count the number of aromatic carbocycles (aromatic carbon-only rings).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of aromatic carbocycles, or missing if molecule is invalid

# Examples

```julia
benzene = mol_from_smiles("c1ccccc1")
count = num_aromatic_carbocycles(benzene)  # 1
```
"""
function num_aromatic_carbocycles(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_aromatic_carbocycles(mol._rdkit_mol))
end

"""
    num_aromatic_heterocycles(mol::Molecule) -> Union{Int,Missing}

Count the number of aromatic heterocycles (aromatic rings containing heteroatoms).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of aromatic heterocycles, or missing if molecule is invalid

# Examples

```julia
pyridine = mol_from_smiles("c1cccnc1")
count = num_aromatic_heterocycles(pyridine)  # 1
```
"""
function num_aromatic_heterocycles(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_aromatic_heterocycles(mol._rdkit_mol))
end

"""
    num_atom_stereo_centers(mol::Molecule) -> Union{Int,Missing}

Count the number of defined atom stereocenters.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of defined stereocenters, or missing if molecule is invalid

# Examples

```julia
chiral_mol = mol_from_smiles("C[C@H](O)C")
count = num_atom_stereo_centers(chiral_mol)  # 1
```
"""
function num_atom_stereo_centers(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_atom_stereo_centers(mol._rdkit_mol))
end

"""
    num_amide_bonds(mol::Molecule) -> Union{Int,Missing}

Count the number of amide bonds in the molecule.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Int,Missing}`: Number of amide bonds, or missing if molecule is invalid

# Examples

```julia
acetamide = mol_from_smiles("CC(=O)N")
count = num_amide_bonds(acetamide)  # 1
```
"""
function num_amide_bonds(mol::Molecule)
    !mol.valid && return missing
    return pyconvert(Int, _num_amide_bonds(mol._rdkit_mol))
end

"""
    asphericity(mol::Molecule; conf_id::Int=-1) -> Union{Float64,Missing}

Calculate the asphericity of a molecule from its 3D coordinates.

Asphericity describes how much a molecule deviates from a spherical shape.

# Arguments

  - `mol::Molecule`: Input molecule (must have 3D coordinates)
  - `conf_id::Int`: Conformer ID to use (-1 for default)

# Returns

  - `Union{Float64,Missing}`: Asphericity value, or missing if molecule is invalid or lacks 3D coordinates

# Examples

```julia
mol = mol_from_smiles("CCO")
conformers = generate_3d_conformers(mol, 1)
if !isempty(conformers)
    mol_3d = conformers[1].molecule
    asp = asphericity(mol_3d)
end
```

# Notes

  - Requires 3D coordinates to be present
  - Values range from 0 (perfect sphere) to 1 (linear molecule)
  - Useful for describing molecular shape and compactness
"""
function asphericity(mol::Molecule; conf_id::Int = -1)
    !mol.valid && return missing
    try
        return pyconvert(Float64, _asphericity(mol._rdkit_mol; confId = conf_id))
    catch
        return missing
    end
end

"""
    radius_of_gyration(mol::Molecule; conf_id::Int=-1) -> Union{Float64,Missing}

Calculate the radius of gyration from 3D coordinates.

# Arguments

  - `mol::Molecule`: Input molecule (must have 3D coordinates)
  - `conf_id::Int`: Conformer ID to use (-1 for default)

# Returns

  - `Union{Float64,Missing}`: Radius of gyration, or missing if molecule is invalid or lacks 3D coordinates

# Notes

  - Requires 3D coordinates to be present
  - Measures molecular compactness
  - Useful for comparing molecular sizes and shapes
"""
function radius_of_gyration(mol::Molecule; conf_id::Int = -1)
    !mol.valid && return missing
    try
        return pyconvert(Float64, _radius_of_gyration(mol._rdkit_mol; confId = conf_id))
    catch
        return missing
    end
end

"""
    synthetic_accessibility(mol::Molecule) -> Union{Float64,Missing}

Calculate the Synthetic Accessibility Score (SAscore).

SAscore estimates how difficult a compound would be to synthesize, ranging from
1 (very easy) to 10 (very difficult).

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Float64,Missing}`: SAscore (1-10 scale), or missing if molecule is invalid

# Examples

```julia
# Simple molecules are easy to synthesize
ethanol = mol_from_smiles("CCO")
sa_score = synthetic_accessibility(ethanol)  # ≈ 1.98 (easy)

# Complex natural products are difficult
paclitaxel = mol_from_smiles(
    "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
)
if paclitaxel.valid
    complex_score = synthetic_accessibility(paclitaxel)  # ≈ 8+ (very difficult)
end
```

# Notes

  - Based on fragment contributions and structural complexity
  - Scores: 1-3 (easy), 4-6 (moderate), 7-10 (difficult)
  - Trained on known synthetic compounds vs. non-synthesizable structures
  - Essential for virtual screening and drug design
  - Helps prioritize synthesizable compounds in large libraries
"""
function synthetic_accessibility(mol::Molecule)
    !mol.valid && return missing
    try
        return pyconvert(Float64, _sascore(mol._rdkit_mol))
    catch
        return missing
    end
end

"""
    calc_all_descriptors(mol::Molecule) -> Union{Dict{Symbol,Any},Missing}

Calculate all available molecular descriptors for a molecule using RDKit's CalcMolDescriptors.

# Arguments

  - `mol::Molecule`: Input molecule

# Returns

  - `Union{Dict{Symbol,Any},Missing}`: Dictionary with descriptor names as keys and values (Float64, Int, String, etc.), or missing if molecule is invalid

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
descriptors = calc_all_descriptors(mol)
println(descriptors[:MolWt])       # Molecular weight
println(descriptors[:NumHDonors])  # Hydrogen bond donors
println(descriptors[:TPSA])        # Topological polar surface area
```

# Notes

  - Returns all available descriptors in a single call
  - More efficient than calling individual descriptor functions
  - Descriptor names are converted to symbols for Julia compatibility
  - Contains 200+ descriptors including all Lipinski, Crippen, and topological descriptors
"""
function calc_all_descriptors(mol::Molecule)
    !mol.valid && return missing

    # Get all descriptors as a Python dictionary
    py_descriptors = _calc_mol_descriptors(mol._rdkit_mol)

    # Convert to Julia Dict with Symbol keys
    result = Dict{Symbol, Any}()
    for item in py_descriptors.items()
        key = item[0]
        value = item[1]

        # Convert key to Symbol
        key_str = pyconvert(String, key)

        # Try to convert value to appropriate Julia type
        try
            # First try as Float64 (most common)
            val_converted = pyconvert(Float64, value)
            result[Symbol(key_str)] = val_converted
        catch
            try
                # Try as Int if Float64 fails
                val_converted = pyconvert(Int, value)
                result[Symbol(key_str)] = val_converted
            catch
                try
                    # Try as String if numeric types fail
                    val_converted = pyconvert(String, value)
                    result[Symbol(key_str)] = val_converted
                catch
                    # If all else fails, keep as Python object
                    result[Symbol(key_str)] = value
                end
            end
        end
    end

    return result
end

# Vectorized versions for multiple molecules
for func in [
    :molecular_weight,
    :exact_molecular_weight,
    :heavy_atom_count,
    :num_heteroatoms,
    :num_rotatable_bonds,
    :num_hbd,
    :num_hba,
    :logp,
    :tpsa,
    :slogp_vsa,
    :num_rings,
    :num_aromatic_rings,
    :num_saturated_rings,
    :bertz_ct,
    :balaban_j,
    :chi0v,
    :kappa1,
    :calc_all_descriptors,
    # Advanced drug-like and ADMET descriptors
    :qed,
    :synthetic_accessibility,
    :fraction_csp3,
    :labute_asa,
    :molar_refractivity,
    # Advanced ring and structure counts
    :num_aliphatic_carbocycles,
    :num_aromatic_carbocycles,
    :num_aromatic_heterocycles,
    :num_atom_stereo_centers,
    :num_amide_bonds,
    # 3D descriptors (note: these may fail for molecules without 3D coordinates)
    :asphericity,
    :radius_of_gyration,
]
    @eval function $(func)(mols::Vector{Union{Molecule, Missing}})
        return [mol === missing ? missing : $(func)(mol) for mol in mols]
    end
    @eval function $(func)(mols::Vector{Molecule})
        return [$(func)(mol) for mol in mols]
    end
end

# Address extraction (moved from old descriptors.jl)
"""
    get_address(mol::Union{Molecule, Missing}) -> Union{String, Missing}
    get_address(mol_list::Vector{Union{Molecule, Missing}}) -> Vector{Union{String, Missing}}

Extract the memory address of the underlying RDKit molecule object.

# Arguments

  - `mol`: A Molecule object or Missing
  - `mol_list`: Vector of Molecule objects or Missing values

# Returns

  - `String`: Hexadecimal memory address (e.g., "0x1a2b3c4d")
  - `missing`: If molecule is missing or invalid

# Examples

```julia
mol = mol_from_smiles("CCO")
address = get_address(mol)  # Returns something like "0x1a2b3c4d"
```

# Notes

  - Useful for debugging and tracking molecule objects
  - Each molecule instance has a unique memory address
  - Addresses change between different program runs
"""
function get_address(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    return extract_address(pyconvert(String, mol._rdkit_mol.__str__()))
end

function get_address(mol_list::Vector{Union{Molecule, Missing}})
    results = Vector{Py}(undef, length(mol_list))
    @inbounds for i in eachindex(mol_list)
        if ismissing(mol_list[i])
            results[i] = missing
        else
            results[i] = mol_list[i]._rdkit_mol.__str__()
        end
    end
    strings = pyconvert(Vector{Union{String, Missing}}, results)
    return map(extract_address, strings)
end

# Additional Chi indices - molecular connectivity
"""
    chi0n(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi0n molecular connectivity index.
"""
function chi0n(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi0n(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi0n: $e"
        return missing
    end
end

"""
    chi1n(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi1n molecular connectivity index.
"""
function chi1n(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi1n(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi1n: $e"
        return missing
    end
end

"""
    chi2n(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi2n molecular connectivity index.
"""
function chi2n(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi2n(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi2n: $e"
        return missing
    end
end

"""
    chi3n(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi3n molecular connectivity index.
"""
function chi3n(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi3n(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi3n: $e"
        return missing
    end
end

"""
    chi4n(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi4n molecular connectivity index.
"""
function chi4n(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi4n(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi4n: $e"
        return missing
    end
end

"""
    chi1v(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi1v valence molecular connectivity index.
"""
function chi1v(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi1v(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi1v: $e"
        return missing
    end
end

"""
    chi2v(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi2v valence molecular connectivity index.
"""
function chi2v(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi2v(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi2v: $e"
        return missing
    end
end

"""
    chi3v(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi3v valence molecular connectivity index.
"""
function chi3v(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi3v(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi3v: $e"
        return missing
    end
end

"""
    chi4v(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Chi4v valence molecular connectivity index.
"""
function chi4v(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _chi4v(mol._rdkit_mol))
    catch e
        @warn "Error calculating Chi4v: $e"
        return missing
    end
end

# Additional Kappa descriptors
"""
    kappa2(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Kappa2 shape index.
"""
function kappa2(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _kappa2(mol._rdkit_mol))
    catch e
        @warn "Error calculating Kappa2: $e"
        return missing
    end
end

"""
    kappa3(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Kappa3 shape index.
"""
function kappa3(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _kappa3(mol._rdkit_mol))
    catch e
        @warn "Error calculating Kappa3: $e"
        return missing
    end
end

# EState descriptors
"""
    max_e_state_index(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the maximum E-state index.
"""
function max_e_state_index(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _max_e_state_index(mol._rdkit_mol))
    catch e
        @warn "Error calculating MaxEStateIndex: $e"
        return missing
    end
end

"""
    min_e_state_index(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the minimum E-state index.
"""
function min_e_state_index(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _min_e_state_index(mol._rdkit_mol))
    catch e
        @warn "Error calculating MinEStateIndex: $e"
        return missing
    end
end

# Simple atom counts
"""
    num_carbons(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of carbon atoms.
"""
function num_carbons(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        atoms = get_atoms(mol)
        return count(atom -> get_atomic_number(atom) == 6, atoms)
    catch e
        @warn "Error counting carbons: $e"
        return missing
    end
end

"""
    num_nitrogens(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of nitrogen atoms.
"""
function num_nitrogens(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        atoms = get_atoms(mol)
        return count(atom -> get_atomic_number(atom) == 7, atoms)
    catch e
        @warn "Error counting nitrogens: $e"
        return missing
    end
end

"""
    num_oxygens(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of oxygen atoms.
"""
function num_oxygens(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        atoms = get_atoms(mol)
        return count(atom -> get_atomic_number(atom) == 8, atoms)
    catch e
        @warn "Error counting oxygens: $e"
        return missing
    end
end

"""
    num_sulfurs(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of sulfur atoms.
"""
function num_sulfurs(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        atoms = get_atoms(mol)
        return count(atom -> get_atomic_number(atom) == 16, atoms)
    catch e
        @warn "Error counting sulfurs: $e"
        return missing
    end
end

"""
    num_halogens(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of halogen atoms (F, Cl, Br, I, At).
"""
function num_halogens(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        atoms = get_atoms(mol)
        halogen_nums = [9, 17, 35, 53, 85]  # F, Cl, Br, I, At
        return count(atom -> get_atomic_number(atom) in halogen_nums, atoms)
    catch e
        @warn "Error counting halogens: $e"
        return missing
    end
end

"""
    ipc(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Information Content of the distance degree sequence (IPC).
"""
function ipc(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _ipc(mol._rdkit_mol))
    catch e
        @warn "Error calculating IPC: $e"
        return missing
    end
end

# Additional missing descriptors that need high-level wrappers

"""
    num_aliphatic_heterocycles(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of aliphatic heterocycles in the molecule.
"""
function num_aliphatic_heterocycles(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_aliphatic_heterocycles(mol._rdkit_mol))
    catch e
        @warn "Error calculating aliphatic heterocycles: $e"
        return missing
    end
end

"""
    num_saturated_heterocycles(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of saturated heterocycles in the molecule.
"""
function num_saturated_heterocycles(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_saturated_heterocycles(mol._rdkit_mol))
    catch e
        @warn "Error calculating saturated heterocycles: $e"
        return missing
    end
end

"""
    num_saturated_carbocycles(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of saturated carbocycles in the molecule.
"""
function num_saturated_carbocycles(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_saturated_carbocycles(mol._rdkit_mol))
    catch e
        @warn "Error calculating saturated carbocycles: $e"
        return missing
    end
end

"""
    num_unspecified_atom_stereo_centers(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of unspecified atom stereo centers.
"""
function num_unspecified_atom_stereo_centers(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_unspecified_atom_stereo_centers(mol._rdkit_mol))
    catch e
        @warn "Error calculating unspecified stereo centers: $e"
        return missing
    end
end

"""
    num_spiro_atoms(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of spiro atoms in the molecule.
"""
function num_spiro_atoms(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_spiro_atoms(mol._rdkit_mol))
    catch e
        @warn "Error calculating spiro atoms: $e"
        return missing
    end
end

"""
    num_bridgehead_atoms(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of bridgehead atoms in the molecule.
"""
function num_bridgehead_atoms(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_bridgehead_atoms(mol._rdkit_mol))
    catch e
        @warn "Error calculating bridgehead atoms: $e"
        return missing
    end
end

"""
    hall_kier_alpha(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the Hall-Kier alpha descriptor.
"""
function hall_kier_alpha(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _hall_kier_alpha(mol._rdkit_mol))
    catch e
        @warn "Error calculating Hall-Kier alpha: $e"
        return missing
    end
end

"""
    eccentricity(mol::Union{Molecule, Missing}; conf_id::Int = -1) -> Union{Float64, Missing}

Calculate the eccentricity of the molecule's 3D structure.
"""
function eccentricity(mol::Union{Molecule, Missing}; conf_id::Int = -1)
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _eccentricity(mol._rdkit_mol; confId = conf_id))
    catch e
        @warn "Error calculating eccentricity: $e"
        return missing
    end
end

"""
    inertial_shape_factor(mol::Union{Molecule, Missing}; conf_id::Int = -1) -> Union{Float64, Missing}

Calculate the inertial shape factor of the molecule.
"""
function inertial_shape_factor(mol::Union{Molecule, Missing}; conf_id::Int = -1)
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _inertial_shape_factor(mol._rdkit_mol; confId = conf_id))
    catch e
        @warn "Error calculating inertial shape factor: $e"
        return missing
    end
end

# BCUT descriptors
"""
    bcut2d_mwlow(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_MWLOW descriptor.
"""
function bcut2d_mwlow(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_mwlow(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_MWLOW: $e"
        return missing
    end
end

"""
    bcut2d_mwhi(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_MWHI descriptor.
"""
function bcut2d_mwhi(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_mwhi(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_MWHI: $e"
        return missing
    end
end

"""
    bcut2d_chglow(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_CHGLO descriptor.
"""
function bcut2d_chglow(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_chglow(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_CHGLO: $e"
        return missing
    end
end

"""
    bcut2d_chghi(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_CHGHI descriptor.
"""
function bcut2d_chghi(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_chghi(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_CHGHI: $e"
        return missing
    end
end

"""
    bcut2d_logplow(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_LOGPLOW descriptor.
"""
function bcut2d_logplow(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_logplow(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_LOGPLOW: $e"
        return missing
    end
end

"""
    bcut2d_logphi(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_LOGPHI descriptor.
"""
function bcut2d_logphi(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_logphi(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_LOGPHI: $e"
        return missing
    end
end

"""
    bcut2d_mrlow(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_MRLOW descriptor.
"""
function bcut2d_mrlow(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_mrlow(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_MRLOW: $e"
        return missing
    end
end

"""
    bcut2d_mrhi(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the BCUT2D_MRHI descriptor.
"""
function bcut2d_mrhi(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _bcut2d_mrhi(mol._rdkit_mol))
    catch e
        @warn "Error calculating BCUT2D_MRHI: $e"
        return missing
    end
end

"""
    max_absolute_e_state_index(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the maximum absolute E-state index.
"""
function max_absolute_e_state_index(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _max_absolute_e_state_index(mol._rdkit_mol))
    catch e
        @warn "Error calculating max absolute E-state index: $e"
        return missing
    end
end

"""
    min_absolute_e_state_index(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the minimum absolute E-state index.
"""
function min_absolute_e_state_index(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _min_absolute_e_state_index(mol._rdkit_mol))
    catch e
        @warn "Error calculating min absolute E-state index: $e"
        return missing
    end
end

# Removed num_sp3_heavy_atoms as the underlying RDKit function doesn't exist

"""
    num_aliphatic_rings(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of aliphatic rings in the molecule.
"""
function num_aliphatic_rings(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_aliphatic_rings(mol._rdkit_mol))
    catch e
        @warn "Error calculating aliphatic rings: $e"
        return missing
    end
end

"""
    num_heterocycles(mol::Union{Molecule, Missing}) -> Union{Int, Missing}

Count the number of heterocycles in the molecule.
"""
function num_heterocycles(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Int, _num_heterocycles(mol._rdkit_mol))
    catch e
        @warn "Error calculating heterocycles: $e"
        return missing
    end
end

# Vectorized functions for new descriptors (moved to bottom for organization)
const NEW_DESCRIPTOR_FUNCTIONS = [
    # Chi connectivity indices
    :chi0n,
    :chi1n,
    :chi2n,
    :chi3n,
    :chi4n,
    :chi1v,
    :chi2v,
    :chi3v,
    :chi4v,
    # Kappa shape indices
    :kappa2,
    :kappa3,
    # E-state descriptors
    :max_e_state_index,
    :min_e_state_index,
    :max_absolute_e_state_index,
    :min_absolute_e_state_index,
    # Atom counts
    :num_carbons,
    :num_nitrogens,
    :num_oxygens,
    :num_sulfurs,
    :num_halogens,
    # Additional ring counts
    :num_aliphatic_heterocycles,
    :num_saturated_heterocycles,
    :num_saturated_carbocycles,
    :num_aliphatic_rings,
    :num_heterocycles,
    # Stereo and structure counts
    :num_unspecified_atom_stereo_centers,
    :num_spiro_atoms,
    :num_bridgehead_atoms,
    # 3D descriptors
    :eccentricity,
    :inertial_shape_factor,
    # BCUT descriptors
    :bcut2d_mwlow,
    :bcut2d_mwhi,
    :bcut2d_chglow,
    :bcut2d_chghi,
    :bcut2d_logplow,
    :bcut2d_logphi,
    :bcut2d_mrlow,
    :bcut2d_mrhi,
    # Other descriptors
    :hall_kier_alpha,
    # Complexity measures
    :ipc,
]

for func in NEW_DESCRIPTOR_FUNCTIONS
    @eval function $(func)(mols::Vector{Union{Molecule, Missing}})
        return [mol === missing ? missing : $(func)(mol) for mol in mols]
    end
    @eval function $(func)(mols::Vector{Molecule})
        return [$(func)(mol) for mol in mols]
    end
end

# Additional 3D descriptors
"""
    pmi1(mol::Union{Molecule, Missing}; conf_id::Int = -1) -> Union{Float64, Missing}

Calculate the first principal moment of inertia (PMI1).
"""
function pmi1(mol::Union{Molecule, Missing}; conf_id::Int = -1)
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _pmi1(mol._rdkit_mol; confId = conf_id))
    catch e
        @warn "Error calculating PMI1: $e"
        return missing
    end
end

"""
    pmi2(mol::Union{Molecule, Missing}; conf_id::Int = -1) -> Union{Float64, Missing}

Calculate the second principal moment of inertia (PMI2).
"""
function pmi2(mol::Union{Molecule, Missing}; conf_id::Int = -1)
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _pmi2(mol._rdkit_mol; confId = conf_id))
    catch e
        @warn "Error calculating PMI2: $e"
        return missing
    end
end

"""
    pmi3(mol::Union{Molecule, Missing}; conf_id::Int = -1) -> Union{Float64, Missing}

Calculate the third principal moment of inertia (PMI3).
"""
function pmi3(mol::Union{Molecule, Missing}; conf_id::Int = -1)
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _pmi3(mol._rdkit_mol; confId = conf_id))
    catch e
        @warn "Error calculating PMI3: $e"
        return missing
    end
end

# VSA descriptors (SlogP_VSA series)
"""
    slogp_vsa2(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA2 descriptor.
"""
function slogp_vsa2(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa2(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA2: $e"
        return missing
    end
end

"""
    slogp_vsa3(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA3 descriptor.
"""
function slogp_vsa3(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa3(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA3: $e"
        return missing
    end
end

"""
    slogp_vsa4(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA4 descriptor.
"""
function slogp_vsa4(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa4(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA4: $e"
        return missing
    end
end

"""
    slogp_vsa5(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA5 descriptor.
"""
function slogp_vsa5(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa5(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA5: $e"
        return missing
    end
end

"""
    slogp_vsa6(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA6 descriptor.
"""
function slogp_vsa6(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa6(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA6: $e"
        return missing
    end
end

"""
    slogp_vsa7(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA7 descriptor.
"""
function slogp_vsa7(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa7(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA7: $e"
        return missing
    end
end

"""
    slogp_vsa8(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA8 descriptor.
"""
function slogp_vsa8(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa8(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA8: $e"
        return missing
    end
end

"""
    slogp_vsa9(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA9 descriptor.
"""
function slogp_vsa9(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa9(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA9: $e"
        return missing
    end
end

"""
    slogp_vsa10(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA10 descriptor.
"""
function slogp_vsa10(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa10(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA10: $e"
        return missing
    end
end

"""
    slogp_vsa11(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA11 descriptor.
"""
function slogp_vsa11(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa11(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA11: $e"
        return missing
    end
end

"""
    slogp_vsa12(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SlogP_VSA12 descriptor.
"""
function slogp_vsa12(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _slogp_vsa12(mol._rdkit_mol))
    catch e
        @warn "Error calculating SlogP_VSA12: $e"
        return missing
    end
end

# SMR_VSA descriptors
"""
    smr_vsa1(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA1 descriptor.
"""
function smr_vsa1(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa1(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA1: $e"
        return missing
    end
end

"""
    smr_vsa2(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA2 descriptor.
"""
function smr_vsa2(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa2(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA2: $e"
        return missing
    end
end

"""
    smr_vsa3(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA3 descriptor.
"""
function smr_vsa3(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa3(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA3: $e"
        return missing
    end
end

"""
    smr_vsa4(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA4 descriptor.
"""
function smr_vsa4(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa4(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA4: $e"
        return missing
    end
end

"""
    smr_vsa5(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA5 descriptor.
"""
function smr_vsa5(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa5(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA5: $e"
        return missing
    end
end

"""
    smr_vsa6(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA6 descriptor.
"""
function smr_vsa6(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa6(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA6: $e"
        return missing
    end
end

"""
    smr_vsa7(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA7 descriptor.
"""
function smr_vsa7(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa7(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA7: $e"
        return missing
    end
end

"""
    smr_vsa8(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA8 descriptor.
"""
function smr_vsa8(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa8(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA8: $e"
        return missing
    end
end

"""
    smr_vsa9(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA9 descriptor.
"""
function smr_vsa9(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa9(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA9: $e"
        return missing
    end
end

"""
    smr_vsa10(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the SMR_VSA10 descriptor.
"""
function smr_vsa10(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _smr_vsa10(mol._rdkit_mol))
    catch e
        @warn "Error calculating SMR_VSA10: $e"
        return missing
    end
end

# PEOE_VSA descriptors
"""
    peoe_vsa1(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA1 descriptor.
"""
function peoe_vsa1(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa1(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA1: $e"
        return missing
    end
end

"""
    peoe_vsa2(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA2 descriptor.
"""
function peoe_vsa2(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa2(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA2: $e"
        return missing
    end
end

"""
    peoe_vsa3(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA3 descriptor.
"""
function peoe_vsa3(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa3(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA3: $e"
        return missing
    end
end

"""
    peoe_vsa4(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA4 descriptor.
"""
function peoe_vsa4(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa4(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA4: $e"
        return missing
    end
end

"""
    peoe_vsa5(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA5 descriptor.
"""
function peoe_vsa5(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa5(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA5: $e"
        return missing
    end
end

"""
    peoe_vsa6(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA6 descriptor.
"""
function peoe_vsa6(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa6(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA6: $e"
        return missing
    end
end

"""
    peoe_vsa7(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA7 descriptor.
"""
function peoe_vsa7(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa7(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA7: $e"
        return missing
    end
end

"""
    peoe_vsa8(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA8 descriptor.
"""
function peoe_vsa8(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa8(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA8: $e"
        return missing
    end
end

"""
    peoe_vsa9(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA9 descriptor.
"""
function peoe_vsa9(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa9(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA9: $e"
        return missing
    end
end

"""
    peoe_vsa10(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA10 descriptor.
"""
function peoe_vsa10(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa10(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA10: $e"
        return missing
    end
end

"""
    peoe_vsa11(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA11 descriptor.
"""
function peoe_vsa11(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa11(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA11: $e"
        return missing
    end
end

"""
    peoe_vsa12(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA12 descriptor.
"""
function peoe_vsa12(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa12(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA12: $e"
        return missing
    end
end

"""
    peoe_vsa13(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA13 descriptor.
"""
function peoe_vsa13(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa13(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA13: $e"
        return missing
    end
end

"""
    peoe_vsa14(mol::Union{Molecule, Missing}) -> Union{Float64, Missing}

Calculate the PEOE_VSA14 descriptor.
"""
function peoe_vsa14(mol::Union{Molecule, Missing})
    isa(mol, Missing) && return missing
    !mol.valid && return missing
    try
        return pyconvert(Float64, _peoe_vsa14(mol._rdkit_mol))
    catch e
        @warn "Error calculating PEOE_VSA14: $e"
        return missing
    end
end

# Update vectorized function list with the new descriptors
const ADDITIONAL_VSA_DESCRIPTOR_FUNCTIONS = [
    # 3D descriptors
    :pmi1,
    :pmi2,
    :pmi3,
    # SlogP_VSA series (2-12)
    :slogp_vsa2,
    :slogp_vsa3,
    :slogp_vsa4,
    :slogp_vsa5,
    :slogp_vsa6,
    :slogp_vsa7,
    :slogp_vsa8,
    :slogp_vsa9,
    :slogp_vsa10,
    :slogp_vsa11,
    :slogp_vsa12,
    # SMR_VSA series (1-10)
    :smr_vsa1,
    :smr_vsa2,
    :smr_vsa3,
    :smr_vsa4,
    :smr_vsa5,
    :smr_vsa6,
    :smr_vsa7,
    :smr_vsa8,
    :smr_vsa9,
    :smr_vsa10,
    # PEOE_VSA series (1-14)
    :peoe_vsa1,
    :peoe_vsa2,
    :peoe_vsa3,
    :peoe_vsa4,
    :peoe_vsa5,
    :peoe_vsa6,
    :peoe_vsa7,
    :peoe_vsa8,
    :peoe_vsa9,
    :peoe_vsa10,
    :peoe_vsa11,
    :peoe_vsa12,
    :peoe_vsa13,
    :peoe_vsa14,
]

for func in ADDITIONAL_VSA_DESCRIPTOR_FUNCTIONS
    @eval function $(func)(mols::Vector{Union{Molecule, Missing}})
        return [mol === missing ? missing : $(func)(mol) for mol in mols]
    end
    @eval function $(func)(mols::Vector{Molecule})
        return [$(func)(mol) for mol in mols]
    end
end
