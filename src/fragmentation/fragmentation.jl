#######################################################
# Molecular Fragmentation
#######################################################

# Import RDKit fragmentation modules
function _brics()
    return @pyconst(pyimport("rdkit.Chem.BRICS"))
end

function _recap()
    return @pyconst(pyimport("rdkit.Chem.Recap"))
end

function _fragmentmatchers()
    return @pyconst(pyimport("rdkit.Chem.FragmentMatcher"))
end

function _rdchem()
    return @pyconst(pyimport("rdkit.Chem"))
end

"""
    brics_decompose(mol::Molecule; min_fragment_size::Int=1, max_fragment_size::Union{Int,Nothing}=nothing) -> Union{Vector{String},Missing}

Decompose a molecule using BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures) rules.

# Arguments
- `mol::Molecule`: Input molecule to decompose
- `min_fragment_size::Int=1`: Minimum number of atoms in returned fragments
- `max_fragment_size::Union{Int,Nothing}=nothing`: Maximum number of atoms in returned fragments (nothing for no limit)

# Returns
- `Union{Vector{String},Missing}`: Vector of SMILES strings representing fragments, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCCOCCc1ccccc1")
fragments = brics_decompose(mol)
# Returns fragments like ["CCCO", "CCc1ccccc1", etc.]

# With size constraints
fragments = brics_decompose(mol, min_fragment_size=3, max_fragment_size=10)
```

# Notes
- BRICS uses retrosynthetic rules to identify breakable bonds
- Returns canonical SMILES for each fragment
- Useful for fragment-based drug design and retrosynthetic analysis
"""
function brics_decompose(mol::Molecule; min_fragment_size::Int=1, max_fragment_size::Union{Int,Nothing}=nothing)
    !mol.valid && return missing

    try
        # Get BRICS fragments
        brics_frags = _brics().BRICSDecompose(mol._rdkit_mol; minFragmentSize=min_fragment_size)

        # Convert to Julia vector of strings
        fragments = Vector{String}()
        for frag in brics_frags
            frag_smiles = pyconvert(String, frag)

            # Apply max fragment size filter if specified
            if max_fragment_size !== nothing
                # Create temporary molecule to count atoms
                temp_mol = mol_from_smiles(frag_smiles)
                if temp_mol.valid
                    atom_count = heavy_atom_count(temp_mol)
                    if atom_count !== missing && atom_count <= max_fragment_size
                        push!(fragments, frag_smiles)
                    end
                end
            else
                push!(fragments, frag_smiles)
            end
        end

        return fragments
    catch e
        @warn "Error in BRICS decomposition: $e"
        return String[]
    end
end

"""
    recap_decompose(mol::Molecule; min_fragment_size::Int=1) -> Union{Vector{String},Missing}

Decompose a molecule using RECAP (Retrosynthetic Combinatorial Analysis Procedure) rules.

# Arguments
- `mol::Molecule`: Input molecule to decompose
- `min_fragment_size::Int=1`: Minimum number of atoms in returned fragments

# Returns
- `Union{Vector{String},Missing}`: Vector of SMILES strings representing fragments, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCCOCCc1ccccc1")
fragments = recap_decompose(mol)
# Returns RECAP fragments

fragments = recap_decompose(mol, min_fragment_size=5)
```

# Notes
- RECAP uses different fragmentation rules compared to BRICS
- Focuses on bonds that are commonly formed in synthetic chemistry
- Returns canonical SMILES for each fragment
"""
function recap_decompose(mol::Molecule; min_fragment_size::Int=1)
    !mol.valid && return missing

    try
        # Create RECAP decomposition tree
        recap_tree = _recap().RecapDecompose(mol._rdkit_mol)

        # Extract leaf nodes (final fragments)
        fragments = Vector{String}()
        leaves = recap_tree.GetLeaves()

        for leaf in leaves.values()
            frag_mol = leaf.mol
            if !pynot(frag_mol)
                # Check fragment size
                temp_molecule = Molecule(; _rdkit_mol=frag_mol, valid=true, source="recap_fragment")
                atom_count = heavy_atom_count(temp_molecule)

                if atom_count !== missing && atom_count >= min_fragment_size
                    frag_smiles = mol_to_smiles(temp_molecule)
                    if !ismissing(frag_smiles)
                        push!(fragments, frag_smiles)
                    end
                end
            end
        end

        return unique(fragments)  # Remove duplicates
    catch e
        @warn "Error in RECAP decomposition: $e"
        return String[]
    end
end

"""
    get_murcko_scaffold(mol::Molecule) -> Union{String,Missing}

Extract the Murcko scaffold (framework) from a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{String,Missing}`: SMILES string of the Murcko scaffold, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCCOCCc1ccc(CC)cc1")
scaffold = get_murcko_scaffold(mol)
# Returns the core ring system: "c1ccccc1"
```

# Notes
- Murcko scaffolds represent the core ring systems of molecules
- Side chains and functional groups are removed
- Useful for scaffold-based drug design and classification
"""
function get_murcko_scaffold(mol::Molecule)
    !mol.valid && return missing

    try
        scaffolds = @pyconst(pyimport("rdkit.Chem.Scaffolds.MurckoScaffold"))
        scaffold_mol = scaffolds.GetScaffoldForMol(mol._rdkit_mol)

        if pynot(scaffold_mol)
            return missing
        end

        scaffold_molecule = Molecule(; _rdkit_mol=scaffold_mol, valid=true, source="murcko_scaffold")
        return mol_to_smiles(scaffold_molecule)
    catch e
        @warn "Error extracting Murcko scaffold: $e"
        return missing
    end
end

"""
    get_generic_scaffold(mol::Molecule) -> Union{String,Missing}

Extract the generic Murcko scaffold from a molecule (heteroatoms converted to carbon).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{String,Missing}`: SMILES string of the generic scaffold, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCCOCCc1ccc(N)cc1")
generic_scaffold = get_generic_scaffold(mol)
# Returns generic scaffold with all atoms as carbon
```

# Notes
- Generic scaffolds normalize heteroatoms to carbon
- Useful for finding structurally similar scaffolds regardless of heteroatom identity
- Better for broad scaffold-based clustering
"""
function get_generic_scaffold(mol::Molecule)
    !mol.valid && return missing

    try
        scaffolds = @pyconst(pyimport("rdkit.Chem.Scaffolds.MurckoScaffold"))
        generic_mol = scaffolds.MakeScaffoldGeneric(scaffolds.GetScaffoldForMol(mol._rdkit_mol))

        if pynot(generic_mol)
            return missing
        end

        generic_molecule = Molecule(; _rdkit_mol=generic_mol, valid=true, source="generic_scaffold")
        return mol_to_smiles(generic_molecule)
    catch e
        @warn "Error extracting generic scaffold: $e"
        return missing
    end
end

"""
    fragment_by_bonds(mol::Molecule, bond_indices::Vector{Int}) -> Union{Vector{String},Missing}

Fragment a molecule by breaking specific bonds.

# Arguments
- `mol::Molecule`: Input molecule
- `bond_indices::Vector{Int}`: Indices of bonds to break (0-based indexing)

# Returns
- `Union{Vector{String},Missing}`: Vector of SMILES strings representing fragments, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCCOCCc1ccccc1")
# Break bonds at indices 3 and 7
fragments = fragment_by_bonds(mol, [3, 7])
```

# Notes
- Bond indices use 0-based indexing (RDKit convention)
- Fragments are returned as canonical SMILES
- Useful for custom fragmentation strategies
"""
function fragment_by_bonds(mol::Molecule, bond_indices::Vector{Int})
    !mol.valid && return missing

    try
        # Fragment the molecule
        fragmented = _rdchem().FragmentOnBonds(mol._rdkit_mol, bond_indices)

        # Get individual fragments
        frags = _rdchem().GetMolFrags(fragmented; asMols=true)

        fragments = Vector{String}()
        for frag in frags
            if !pynot(frag)
                frag_molecule = Molecule(; _rdkit_mol=frag, valid=true, source="bond_fragment")
                frag_smiles = mol_to_smiles(frag_molecule)
                if !ismissing(frag_smiles)
                    push!(fragments, frag_smiles)
                end
            end
        end

        return fragments
    catch e
        @warn "Error in bond-based fragmentation: $e"
        return String[]
    end
end

"""
    get_fragment_count(mol::Molecule) -> Union{Int,Missing}

Count the number of disconnected fragments in a molecule.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{Int,Missing}`: Number of disconnected fragments, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO.CCC")  # Two disconnected fragments
count = get_fragment_count(mol)   # Returns 2
```

# Notes
- Counts disconnected components in the molecular graph
- Useful for identifying salts, mixtures, or reaction products
"""
function get_fragment_count(mol::Molecule)
    !mol.valid && return missing

    try
        # Get fragment indices and count them
        frags = _rdchem().GetMolFrags(mol._rdkit_mol; asMols=false)
        return length(frags)
    catch e
        @warn "Error counting fragments: $e"
        return missing
    end
end

"""
    split_fragments(mol::Molecule) -> Union{Vector{Molecule},Missing}

Split a molecule into its disconnected fragments.

# Arguments
- `mol::Molecule`: Input molecule (potentially with disconnected fragments)

# Returns
- `Union{Vector{Molecule},Missing}`: Vector of individual fragment molecules, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO.CCC")  # Two disconnected fragments
fragments = split_fragments(mol)
# Returns [Molecule(CCO), Molecule(CCC)]
```

# Notes
- Separates disconnected molecular components
- Each fragment is returned as a separate Molecule object
- Useful for processing mixtures or salt forms
"""
function split_fragments(mol::Molecule)
    !mol.valid && return missing

    try
        # Get fragments as molecule objects
        frags = _rdchem().GetMolFrags(mol._rdkit_mol; asMols=true)

        fragments = Vector{Molecule}()
        for (i, frag) in enumerate(frags)
            if !pynot(frag)
                frag_molecule = Molecule(; _rdkit_mol=frag, valid=true, source="fragment_$i")
                push!(fragments, frag_molecule)
            end
        end

        return fragments
    catch e
        @warn "Error splitting fragments: $e"
        return Molecule[]
    end
end

"""
    get_largest_fragment(mol::Molecule) -> Union{Molecule,Missing}

Get the largest fragment from a molecule (by heavy atom count).

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Union{Molecule,Missing}`: The largest fragment as a Molecule object, or missing if molecule is invalid

# Examples
```julia
mol = mol_from_smiles("CCO.CCCCCCCc1ccccc1")  # Small and large fragments
largest = get_largest_fragment(mol)
# Returns the larger fragment (CCCCCCCc1ccccc1)
```

# Notes
- Identifies the fragment with the most heavy atoms
- Useful for removing salts or small contaminants
- Returns the original molecule if already connected
"""
function get_largest_fragment(mol::Molecule)
    !mol.valid && return missing

    fragments = split_fragments(mol)
    if fragments === missing || isempty(fragments)
        return missing
    end

    # If only one fragment, return it
    if length(fragments) == 1
        return fragments[1]
    end

    # Find the fragment with the most heavy atoms
    max_atoms = 0
    largest_frag = fragments[1]

    for frag in fragments
        atom_count = heavy_atom_count(frag)
        if atom_count !== missing && atom_count > max_atoms
            max_atoms = atom_count
            largest_frag = frag
        end
    end

    return largest_frag
end

# Vectorized versions for multiple molecules
for func in [
    :brics_decompose,
    :recap_decompose,
    :get_murcko_scaffold,
    :get_generic_scaffold,
    :get_fragment_count,
    :split_fragments,
    :get_largest_fragment,
]
    @eval function $(func)(mols::Vector{Union{Molecule, Missing}}, args...; kwargs...)
        return [$(func)(mol, args...; kwargs...) for mol in mols]
    end
    @eval function $(func)(mols::Vector{Molecule}, args...; kwargs...)
        return [$(func)(mol, args...; kwargs...) for mol in mols]
    end
end