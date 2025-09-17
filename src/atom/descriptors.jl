#######################################################
# Atom descriptors
#######################################################

"""
    get_atomic_number(atom::Atom) -> Int

Get the atomic number of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The atomic number (number of protons)

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 1)
atomic_num = get_atomic_number(atom)  # Returns 6 for carbon
```
"""
function get_atomic_number(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetAtomicNum())
end

"""
    get_symbol(atom::Atom) -> String

Get the atomic symbol of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `String`: The atomic symbol (e.g., "C", "N", "O")

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 3)
symbol = get_symbol(atom)  # Returns "O" for oxygen
```
"""
function get_symbol(atom::Atom)
    return pyconvert(String, atom._rdkit_atom.GetSymbol())
end

"""
    get_degree(atom::Atom) -> Int

Get the degree (number of bonds) of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The number of bonds connected to this atom

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 1)  # First carbon
degree = get_degree(atom)  # Returns 1 (connected to second carbon)
```
"""
function get_degree(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetDegree())
end

"""
    get_valence(atom::Atom) -> Int

Get the total valence of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The total valence (sum of bond orders)

# Example
```julia
mol = mol_from_smiles("C=O")
atom = get_atom(mol, 1)  # Carbon
valence = get_valence(atom)  # Returns 4
```
"""
function get_valence(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetTotalValence())
end

function get_implicit_valence(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetImplicitValence())
end

function get_explicit_valence(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetExplicitValence())
end

"""
    get_formal_charge(atom::Atom) -> Int

Get the formal charge of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The formal charge (-1, 0, +1, etc.)

# Example
```julia
mol = mol_from_smiles("[NH3+]")
atom = get_atom(mol, 1)
charge = get_formal_charge(atom)  # Returns +1
```
"""
function get_formal_charge(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetFormalCharge())
end

"""
    get_hybridization(atom::Atom) -> String

Get the hybridization state of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `String`: The hybridization state ("SP3", "SP2", "SP", "UNSPECIFIED", etc.)

# Example
```julia
mol = mol_from_smiles("C=C")
atom = get_atom(mol, 1)  # Carbon with double bond
hybridization = get_hybridization(atom)  # Returns "SP2"
```
"""
function get_hybridization(atom::Atom)
    return pyconvert(String, atom._rdkit_atom.GetHybridization().__str__())
end

"""
    get_num_explicit_hs(atom::Atom) -> Int

Get the number of explicit hydrogens attached to an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The number of explicitly specified hydrogen atoms

# Example
```julia
mol = mol_from_smiles("C[H]")
atom = get_atom(mol, 1)  # Carbon
explicit_hs = get_num_explicit_hs(atom)  # Returns 1
```
"""
function get_num_explicit_hs(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetNumExplicitHs())
end

"""
    get_num_implicit_hs(atom::Atom) -> Int

Get the number of implicit hydrogens attached to an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The number of implicitly specified hydrogen atoms

# Example
```julia
mol = mol_from_smiles("C")
atom = get_atom(mol, 1)  # Carbon
implicit_hs = get_num_implicit_hs(atom)  # Returns 4 (CH4)
```
"""
function get_num_implicit_hs(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetNumImplicitHs())
end

"""
    get_total_num_hs(atom::Atom) -> Int

Get the total number of hydrogens (explicit + implicit) attached to an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The total number of hydrogen atoms attached

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 1)  # First carbon
total_hs = get_total_num_hs(atom)  # Returns 3 (CH3)
```
"""
function get_total_num_hs(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetTotalNumHs())
end

"""
    get_mass(atom::Atom) -> Float64

Get the atomic mass of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Float64`: The atomic mass in atomic mass units (amu)

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 1)  # Carbon
mass = get_mass(atom)  # Returns ~12.011
```
"""
function get_mass(atom::Atom)
    return pyconvert(Float64, atom._rdkit_atom.GetMass())
end

"""
    get_isotope(atom::Atom) -> Int

Get the isotope number of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The isotope number (0 for most common isotope)

# Example
```julia
mol = mol_from_smiles("[13C]")
atom = get_atom(mol, 1)  # Carbon-13
isotope = get_isotope(atom)  # Returns 13
```
"""
function get_isotope(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetIsotope())
end

"""
    is_aromatic(atom::Atom) -> Bool

Check if an atom is aromatic.

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is aromatic, false otherwise

# Example
```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
atom = get_atom(mol, 1)  # Aromatic carbon
is_aromatic(atom)  # Returns true
```
"""
function is_aromatic(atom::Atom)
    return pyconvert(Bool, atom._rdkit_atom.GetIsAromatic())
end

"""
    is_in_ring(atom::Atom) -> Bool

Check if an atom is part of a ring.

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is in a ring, false otherwise

# Example
```julia
mol = mol_from_smiles("c1ccccc1")
atom = get_atom(mol, 1)
is_in_ring(atom)  # Returns true
```
"""
function is_in_ring(atom::Atom)
    return pyconvert(Bool, atom._rdkit_atom.IsInRing())
end

"""
    is_in_ring_size(atom::Atom, size::Int) -> Bool

Check if an atom is part of a ring of specific size.

# Arguments
- `atom`: An Atom object
- `size`: Ring size to check for

# Returns
- `Bool`: true if the atom is in a ring of the specified size

# Example
```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene (6-membered ring)
atom = get_atom(mol, 1)
is_in_ring_size(atom, 6)  # Returns true
is_in_ring_size(atom, 5)  # Returns false
```
"""
function is_in_ring_size(atom::Atom, size::Int)
    return pyconvert(Bool, atom._rdkit_atom.IsInRingSize(size))
end

"""
    get_chiral_tag(atom::Atom) -> String

Get the chirality tag of an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `String`: The chirality tag ("CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW", etc.)

# Example
```julia
mol = mol_from_smiles("C[C@H](O)N")  # Chiral carbon
atom = get_atom(mol, 2)
chiral_tag = get_chiral_tag(atom)  # Returns chirality information
```
"""
function get_chiral_tag(atom::Atom)
    return pyconvert(String, atom._rdkit_atom.GetChiralTag().__str__())
end

"""
    get_atoms(mol::Molecule) -> Union{Vector{Atom}, Missing}

Get all atoms from a molecule as a vector of Atom objects.

# Arguments
- `mol`: A Molecule object

# Returns
- `Vector{Atom}`: Vector of all atoms in the molecule
- `missing`: If the molecule is invalid

# Example
```julia
mol = mol_from_smiles("CCO")
atoms = get_atoms(mol)  # Returns vector of 3 atoms
length(atoms)  # Returns 3
```
"""
function get_atoms(mol::Molecule)
    !mol.valid && return missing
    rdkit_atoms = mol._rdkit_mol.GetAtoms()
    return [Atom(; _rdkit_atom=atom) for atom in rdkit_atoms]
end

"""
    get_atom(mol::Molecule, idx::Int) -> Union{Atom, Missing}

Get a specific atom from a molecule by its index (1-based).

# Arguments
- `mol`: A Molecule object
- `idx`: 1-based index of the atom to retrieve

# Returns
- `Atom`: The atom at the specified index
- `missing`: If the molecule is invalid

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 1)  # Returns first carbon atom
symbol = get_symbol(atom)  # Returns "C"
```
"""
function get_atom(mol::Molecule, idx::Int)
    !mol.valid && return missing
    try
        # Convert to 0-based indexing for RDKit
        rdkit_atom = mol._rdkit_mol.GetAtomWithIdx(idx - 1)
        return Atom(; _rdkit_atom=rdkit_atom)
    catch
        return missing
    end
end

# Atom environment descriptors
"""
    get_neighbors(mol::Molecule, atom_idx::Int) -> Union{Vector{Int}, Missing}

Get the indices of neighboring atoms for a given atom.

# Arguments
- `mol::Molecule`: The molecule containing the atom
- `atom_idx::Int`: 1-based index of the atom

# Returns
- `Vector{Int}`: Vector of 1-based atom indices that are neighbors to the specified atom
- `missing`: If the molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
neighbors = get_neighbors(mol, 1)  # Neighbors of first carbon: [2]
neighbors = get_neighbors(mol, 2)  # Neighbors of second carbon: [1, 3]
```
"""
function get_neighbors(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing
    atom = mol._rdkit_mol.GetAtomWithIdx(atom_idx - 1)  # Convert to 0-based
    neighbors = atom.GetNeighbors()
    return [pyconvert(Int, neighbor.GetIdx()) + 1 for neighbor in neighbors]  # Convert back to 1-based
end

"""
    get_bonds_from_atom(mol::Molecule, atom_idx::Int) -> Union{Vector{Bond}, Missing}

Get all bonds connected to a specific atom.

# Arguments
- `mol::Molecule`: The molecule containing the atom
- `atom_idx::Int`: 1-based index of the atom

# Returns
- `Vector{Bond}`: Vector of Bond objects connected to the specified atom
- `missing`: If the molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
bonds = get_bonds_from_atom(mol, 2)  # Bonds from second carbon
for bond in bonds
    bond_type = get_bond_type(bond)
    println("Bond type: ", bond_type)
end
```
"""
function get_bonds_from_atom(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing
    atom = mol._rdkit_mol.GetAtomWithIdx(atom_idx - 1)  # Convert to 0-based
    bonds = atom.GetBonds()
    return [Bond(; _rdkit_bond=bond) for bond in bonds]
end

# Gasteiger partial charges (if available)
"""
    compute_gasteiger_charges!(mol::Molecule) -> Nothing

Compute Gasteiger partial charges for all atoms in a molecule.

# Arguments
- `mol::Molecule`: The molecule to compute charges for (modified in-place)

# Returns
- `Nothing`

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
compute_gasteiger_charges!(mol)

# Now can access charges for individual atoms
atoms = get_atoms(mol)
charge = get_gasteiger_charge(atoms[1])  # Charge on first carbon
```
"""
function compute_gasteiger_charges!(mol::Molecule)
    !mol.valid && return nothing
    @pyconst(pyimport("rdkit.Chem.rdPartialCharges").ComputeGasteigerCharges)(
        mol._rdkit_mol
    )
    return nothing
end

"""
    get_gasteiger_charge(atom::Atom) -> Float64

Get the Gasteiger partial charge for a specific atom.

# Arguments
- `atom::Atom`: The atom to get the charge for

# Returns
- `Float64`: The Gasteiger partial charge

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
compute_gasteiger_charges!(mol)  # Must compute charges first

atoms = get_atoms(mol)
charge_carbon = get_gasteiger_charge(atoms[1])
charge_oxygen = get_gasteiger_charge(atoms[3])
```
"""
function get_gasteiger_charge(atom::Atom)
    return pyconvert(Float64, atom._rdkit_atom.GetDoubleProp("_GasteigerCharge"))
end

# Atom type information
function get_atom_type(atom::Atom)
    try
        return pyconvert(String, atom._rdkit_atom.GetAtomType())
    catch
        return missing
    end
end

# Vectorized functions for multiple atoms
function get_atomic_numbers(atoms::Vector{Atom})
    return [get_atomic_number(atom) for atom in atoms]
end

function get_symbols(atoms::Vector{Atom})
    return [get_symbol(atom) for atom in atoms]
end

function get_degrees(atoms::Vector{Atom})
    return [get_degree(atom) for atom in atoms]
end

function get_formal_charges(atoms::Vector{Atom})
    return [get_formal_charge(atom) for atom in atoms]
end

"""
    get_num_radical_electrons(atom::Atom) -> Int

Get the number of radical electrons on an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Int`: The number of radical electrons

# Example
```julia
mol = mol_from_smiles("C[CH2]")  # Radical carbon
atom = get_atom(mol, 2)
radicals = get_num_radical_electrons(atom)  # Returns number of unpaired electrons
```
"""
function get_num_radical_electrons(atom::Atom)
    return pyconvert(Int, atom._rdkit_atom.GetNumRadicalElectrons())
end

"""
    is_hetero(atom::Atom) -> Bool

Check if an atom is a heteroatom (not carbon or hydrogen).

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is a heteroatom, false otherwise

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 3)  # Oxygen
is_hetero(atom)  # Returns true
```
"""
function is_hetero(atom::Atom)
    symbol = get_symbol(atom)
    return !(symbol in ["C", "H"])
end

"""
    get_cip_code(atom::Atom) -> Union{String,Missing}

Get the CIP (Cahn-Ingold-Prelog) stereochemistry code for an atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Union{String,Missing}`: The CIP code ("R", "S", or missing if not assigned)

# Example
```julia
mol = mol_from_smiles("C[C@H](O)N")  # Chiral carbon
atom = get_atom(mol, 2)
cip = get_cip_code(atom)  # Returns "R" or "S"
```
"""
function get_cip_code(atom::Atom)
    try
        cip_prop = atom._rdkit_atom.GetProp("_CIPCode")
        return pyconvert(String, cip_prop)
    catch
        return missing
    end
end

"""
    is_chiral_center(atom::Atom) -> Bool

Check if an atom is a chiral center.

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is a chiral center, false otherwise

# Example
```julia
mol = mol_from_smiles("C[C@H](O)N")  # Chiral carbon
atom = get_atom(mol, 2)
is_chiral_center(atom)  # Returns true
```
"""
function is_chiral_center(atom::Atom)
    chiral_tag = get_chiral_tag(atom)
    return !(chiral_tag in ["CHI_UNSPECIFIED", "UNSPECIFIED"])
end

"""
    is_hydrogen_donor(atom::Atom) -> Bool

Check if an atom can act as a hydrogen bond donor.

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is a hydrogen bond donor, false otherwise

# Example
```julia
mol = mol_from_smiles("CCO")
atom = get_atom(mol, 3)  # Oxygen
is_hydrogen_donor(atom)  # Returns true (OH can donate hydrogen)
```
"""
function is_hydrogen_donor(atom::Atom)
    # An atom is a hydrogen donor if it's O, N, or S with hydrogens attached
    symbol = get_symbol(atom)
    if symbol in ["O", "N", "S"]
        return get_total_num_hs(atom) > 0
    end
    return false
end

"""
    is_hydrogen_acceptor(atom::Atom) -> Bool

Check if an atom can act as a hydrogen bond acceptor.

# Arguments
- `atom`: An Atom object

# Returns
- `Bool`: true if the atom is a hydrogen bond acceptor, false otherwise

# Example
```julia
mol = mol_from_smiles("C=O")
atom = get_atom(mol, 2)  # Oxygen
is_hydrogen_acceptor(atom)  # Returns true (C=O oxygen can accept hydrogen bonds)
```
"""
function is_hydrogen_acceptor(atom::Atom)
    # An atom is a hydrogen acceptor if it's O, N, or F with lone pairs
    symbol = get_symbol(atom)
    if symbol in ["O", "N", "F"]
        return true
    end
    # Some other atoms can also be acceptors (S, Cl, etc.)
    return symbol in ["S", "Cl", "Br", "I"]
end

"""
    get_ring_size(atom::Atom) -> Union{Int,Missing}

Get the size of the smallest ring containing the atom.

# Arguments
- `atom`: An Atom object

# Returns
- `Union{Int,Missing}`: The size of the smallest ring, or missing if not in a ring

# Example
```julia
mol = mol_from_smiles("c1ccccc1")  # Benzene
atom = get_atom(mol, 1)
ring_size = get_ring_size(atom)  # Returns 6
```
"""
function get_ring_size(atom::Atom)
    if !is_in_ring(atom)
        return missing
    end

    # Check ring sizes from 3 to 20 (reasonable range)
    for size in 3:20
        if is_in_ring_size(atom, size)
            return size
        end
    end

    return missing  # If no ring size found in reasonable range
end

"""
    get_crippen_log_p_contribution(mol::Molecule, atom_idx::Int) -> Union{Float64,Missing}

Get the Crippen LogP contribution for a specific atom.

# Arguments
- `mol`: A Molecule object
- `atom_idx`: 1-based index of the atom

# Returns
- `Union{Float64,Missing}`: The LogP contribution, or missing if not available

# Example
```julia
mol = mol_from_smiles("CCO")
logp_contrib = get_crippen_log_p_contribution(mol, 1)
```
"""
function get_crippen_log_p_contribution(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing

    try
        crippen = @pyconst(pyimport("rdkit.Chem.Crippen"))
        contribs = crippen._GetAtomContribs(mol._rdkit_mol)
        # Convert to 0-based indexing for RDKit
        return pyconvert(Float64, contribs[atom_idx - 1][0])  # LogP contribution
    catch e
        return missing
    end
end

"""
    get_crippen_molar_refractivity_contribution(mol::Molecule, atom_idx::Int) -> Union{Float64,Missing}

Get the Crippen molar refractivity contribution for a specific atom.

# Arguments
- `mol`: A Molecule object
- `atom_idx`: 1-based index of the atom

# Returns
- `Union{Float64,Missing}`: The molar refractivity contribution, or missing if not available

# Example
```julia
mol = mol_from_smiles("CCO")
mr_contrib = get_crippen_molar_refractivity_contribution(mol, 1)
```
"""
function get_crippen_molar_refractivity_contribution(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing

    try
        crippen = @pyconst(pyimport("rdkit.Chem.Crippen"))
        contribs = crippen._GetAtomContribs(mol._rdkit_mol)
        # Convert to 0-based indexing for RDKit
        return pyconvert(Float64, contribs[atom_idx - 1][1])  # MR contribution
    catch e
        return missing
    end
end

"""
    get_tpsa_contribution(mol::Molecule, atom_idx::Int) -> Union{Float64,Missing}

Get the TPSA (Topological Polar Surface Area) contribution for a specific atom.

# Arguments
- `mol`: A Molecule object
- `atom_idx`: 1-based index of the atom

# Returns
- `Union{Float64,Missing}`: The TPSA contribution, or missing if not available

# Example
```julia
mol = mol_from_smiles("CCO")
tpsa_contrib = get_tpsa_contribution(mol, 3)  # Oxygen contribution
```
"""
function get_tpsa_contribution(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing

    try
        rdmoldescriptors = @pyconst(pyimport("rdkit.Chem.rdMolDescriptors"))
        contribs = rdmoldescriptors._TPSAContribs(mol._rdkit_mol)
        # Convert to 0-based indexing for RDKit
        return pyconvert(Float64, contribs[atom_idx - 1])
    catch e
        return missing
    end
end

"""
    get_labute_asa_contribution(mol::Molecule, atom_idx::Int) -> Union{Float64,Missing}

Get the Labute ASA (Accessible Surface Area) contribution for a specific atom.

# Arguments
- `mol`: A Molecule object
- `atom_idx`: 1-based index of the atom

# Returns
- `Union{Float64,Missing}`: The ASA contribution, or missing if not available

# Example
```julia
mol = mol_from_smiles("CCO")
asa_contrib = get_labute_asa_contribution(mol, 1)
```
"""
function get_labute_asa_contribution(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing

    try
        rdmoldescriptors = @pyconst(pyimport("rdkit.Chem.rdMolDescriptors"))
        contribs = rdmoldescriptors._LabuteASAContribs(mol._rdkit_mol)
        # Convert to 0-based indexing for RDKit
        return pyconvert(Float64, contribs[atom_idx - 1])
    catch e
        return missing
    end
end

"""
    get_all_atom_properties(mol::Molecule, atom_idx::Int) -> Union{Dict{Symbol,Any},Missing}

Get all available properties for a specific atom as a dictionary.

# Arguments
- `mol`: A Molecule object
- `atom_idx`: 1-based index of the atom

# Returns
- `Union{Dict{Symbol,Any},Missing}`: Dictionary with all atom properties, or missing if invalid

# Example
```julia
mol = mol_from_smiles("CCO")
props = get_all_atom_properties(mol, 1)
println(props[:symbol])         # "C"
println(props[:hybridization])  # "SP3"
println(props[:formal_charge])  # 0
```
"""
function get_all_atom_properties(mol::Molecule, atom_idx::Int)
    !mol.valid && return missing

    atom = get_atom(mol, atom_idx)
    atom === missing && return missing

    props = Dict{Symbol, Any}()

    # Basic properties
    props[:symbol] = get_symbol(atom)
    props[:hybridization] = get_hybridization(atom)
    props[:formal_charge] = get_formal_charge(atom)
    props[:total_num_hs] = get_total_num_hs(atom)
    props[:total_valence] = get_valence(atom)
    props[:num_radical_electrons] = get_num_radical_electrons(atom)
    props[:degree] = get_degree(atom)
    props[:aromatic] = is_aromatic(atom)
    props[:hetero] = is_hetero(atom)
    props[:hydrogen_donor] = is_hydrogen_donor(atom)
    props[:hydrogen_acceptor] = is_hydrogen_acceptor(atom)
    props[:ring] = is_in_ring(atom)
    props[:ring_size] = get_ring_size(atom)
    props[:chiral_center] = is_chiral_center(atom)
    props[:cip_code] = get_cip_code(atom)

    # Contribution-based properties (require molecule context)
    props[:crippen_log_p_contribution] = get_crippen_log_p_contribution(mol, atom_idx)
    props[:crippen_molar_refractivity_contribution] = get_crippen_molar_refractivity_contribution(mol, atom_idx)
    props[:tpsa_contribution] = get_tpsa_contribution(mol, atom_idx)
    props[:labute_asa_contribution] = get_labute_asa_contribution(mol, atom_idx)

    # Gasteiger charge (requires computation)
    try
        props[:gasteiger_charge] = get_gasteiger_charge(atom)
    catch
        # Might not be computed yet
        props[:gasteiger_charge] = missing
    end

    return props
end
