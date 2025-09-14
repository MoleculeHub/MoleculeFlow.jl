#######################################################
# Molecular descriptors
#######################################################

# Import RDKit descriptor modules
function _descriptors()
    return @pyconst(pyimport("rdkit.Chem.Descriptors"))
end

function _crippen()
    return @pyconst(pyimport("rdkit.Chem.Crippen"))
end

function _lipinski()
    return @pyconst(pyimport("rdkit.Chem.Lipinski"))
end

function _rdmoldescs()
    return @pyconst(pyimport("rdkit.Chem.rdMolDescriptors"))
end

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
    return pyconvert(Float64, _descriptors().MolWt(mol._rdkit_mol))
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
    return pyconvert(Float64, _descriptors().ExactMolWt(mol._rdkit_mol))
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
    return pyconvert(Int, _descriptors().HeavyAtomCount(mol._rdkit_mol))
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
    return pyconvert(Int, _descriptors().NumHeteroatoms(mol._rdkit_mol))
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
    return pyconvert(Int, _lipinski().NumRotatableBonds(mol._rdkit_mol))
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
- Used in Lipinski's Rule of Five (≤5 donors)
- Important for drug-like properties
"""
function num_hbd(mol::Molecule)  # Hydrogen bond donors
    !mol.valid && return missing
    return pyconvert(Int, _descriptors().NumHDonors(mol._rdkit_mol))
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
- Used in Lipinski's Rule of Five (≤10 acceptors)
- Important for drug-like properties
"""
function num_hba(mol::Molecule)  # Hydrogen bond acceptors
    !mol.valid && return missing
    return pyconvert(Int, _descriptors().NumHAcceptors(mol._rdkit_mol))
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
    return pyconvert(Float64, _crippen().MolLogP(mol._rdkit_mol))
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
- Important for predicting drug permeability
- Values < 60 Ų typically indicate good oral bioavailability
- Used in Lipinski's Rule of Five
"""
function tpsa(mol::Molecule)  # Topological polar surface area
    !mol.valid && return missing
    return pyconvert(Float64, _descriptors().TPSA(mol._rdkit_mol))
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
    return pyconvert(Float64, _descriptors().SlogP_VSA1(mol._rdkit_mol))
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
    return pyconvert(Int, _descriptors().RingCount(mol._rdkit_mol))
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
    return pyconvert(Int, _descriptors().NumAromaticRings(mol._rdkit_mol))
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
    return pyconvert(Int, _descriptors().NumSaturatedRings(mol._rdkit_mol))
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
- Measures molecular complexity based on graph theory
- Higher values indicate more complex molecular structures
- Useful for drug design and synthesis planning
"""
function bertz_ct(mol::Molecule)  # BertzCT complexity
    !mol.valid && return missing
    return pyconvert(Float64, _descriptors().BertzCT(mol._rdkit_mol))
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
    return pyconvert(Float64, _descriptors().BalabanJ(mol._rdkit_mol))
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
- Chi0v is a valence connectivity index
- Describes molecular connectivity considering valence electrons
- Useful for predicting molecular properties
"""
function chi0v(mol::Molecule)  # Connectivity index
    !mol.valid && return missing
    return pyconvert(Float64, _descriptors().Chi0v(mol._rdkit_mol))
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
- Kappa1 describes molecular shape and branching
- Higher values indicate more linear structures
- Part of the Kier and Hall molecular shape indices
"""
function kappa1(mol::Molecule)  # Kappa shape index
    !mol.valid && return missing
    return pyconvert(Float64, _descriptors().Kappa1(mol._rdkit_mol))
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
    py_descriptors = _descriptors().CalcMolDescriptors(mol._rdkit_mol)

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
]
    @eval function $(func)(mols::Vector{Union{Molecule, Missing}})
        return [$(func)(mol) for mol in mols]
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
