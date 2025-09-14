# Molecular standardization functionality
# Includes salt stripping, tautomer enumeration, and other standardization operations

# Import RDKit modules for standardization
function _standardize_module()
    return @pyconst(pyimport("rdkit.Chem.MolStandardize"))
end

function _tautomer_module()
    return @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize"))
end

function _fragment_module()
    return @pyconst(pyimport("rdkit.Chem.MolStandardize.rdMolStandardize"))
end

function _chem_standardize()
    return @pyconst(pyimport("rdkit.Chem"))
end

"""
    strip_salts(mol::Molecule; largest_fragment_only::Bool=true) -> Molecule

Remove salts and other small fragments from a molecule, keeping only the largest fragment.

# Arguments
- `mol::Molecule`: Input molecule
- `largest_fragment_only::Bool=true`: If true, return only the largest fragment

# Returns
- `Molecule`: Molecule with salts removed

# Example
```julia
mol = mol_from_smiles("CCO.Cl")  # Ethanol with chloride salt
clean_mol = strip_salts(mol)     # Returns just the ethanol
```
"""
function strip_salts(mol::Molecule; largest_fragment_only::Bool=true)
    if !mol.valid
        return mol
    end

    try
        # Get the fragment remover
        fragment_remover = _fragment_module().LargestFragmentChooser()

        # Remove fragments/salts
        clean_mol = fragment_remover.choose(mol._rdkit_mol)

        if pynot(clean_mol)
            return Molecule(;
                _rdkit_mol=clean_mol, valid=false, source="stripped_salts_$(mol.source)"
            )
        end

        return Molecule(;
            _rdkit_mol=clean_mol,
            valid=true,
            source="stripped_salts_$(mol.source)",
            props=copy(mol.props),
        )
    catch e
        @warn "Failed to strip salts from molecule: $e"
        return Molecule(;
            _rdkit_mol=mol._rdkit_mol, valid=false, source="error_$(mol.source)"
        )
    end
end

"""
    enumerate_tautomers(mol::Molecule; max_tautomers::Int=1000) -> Vector{Molecule}

Enumerate all reasonable tautomers of a molecule.

# Arguments
- `mol::Molecule`: Input molecule
- `max_tautomers::Int=1000`: Maximum number of tautomers to generate

# Returns
- `Vector{Molecule}`: Vector of tautomer molecules

# Example
```julia
mol = mol_from_smiles("CC(=O)CC(=O)C")  # Acetylacetone
tautomers = enumerate_tautomers(mol)     # Generates keto-enol tautomers
```
"""
function enumerate_tautomers(mol::Molecule; max_tautomers::Int=1000)
    if !mol.valid
        return [mol]
    end

    try
        # Get the tautomer enumerator
        enumerator = _tautomer_module().TautomerEnumerator()

        # Enumerate tautomers
        tautomers = enumerator.Enumerate(mol._rdkit_mol)

        # Convert to Julia molecules
        result = Molecule[]
        count = 0

        for taut in tautomers
            if count >= max_tautomers
                break
            end

            # All tautomers from RDKit are valid molecules
            taut_smiles = mol_to_smiles(
                Molecule(; _rdkit_mol=taut, valid=true, source="tautomer")
            )
            push!(
                result,
                Molecule(;
                    _rdkit_mol=taut,
                    valid=true,
                    source="tautomer_$(taut_smiles)",
                    props=copy(mol.props),
                ),
            )
            count += 1
        end

        return result
    catch e
        @warn "Failed to enumerate tautomers: $e"
        return [mol]
    end
end

"""
    canonical_tautomer(mol::Molecule) -> Molecule

Get the canonical tautomer of a molecule according to RDKit's tautomer rules.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Canonical tautomer

# Example
```julia
mol = mol_from_smiles("CC(O)=CC(=O)C")   # Enol form
canonical = canonical_tautomer(mol)      # Returns keto form
```
"""
function canonical_tautomer(mol::Molecule)
    if !mol.valid
        return mol
    end

    try
        # Get the tautomer enumerator
        enumerator = _tautomer_module().TautomerEnumerator()

        # Get canonical tautomer
        canonical = enumerator.Canonicalize(mol._rdkit_mol)

        if pynot(canonical)
            return Molecule(;
                _rdkit_mol=canonical, valid=false, source="canonical_$(mol.source)"
            )
        end

        return Molecule(;
            _rdkit_mol=canonical,
            valid=true,
            source="canonical_$(mol.source)",
            props=copy(mol.props),
        )
    catch e
        @warn "Failed to get canonical tautomer: $e"
        return mol
    end
end

"""
    standardize_molecule(mol::Molecule; 
                        strip_salts_flag::Bool=true,
                        canonical_tautomer_flag::Bool=true,
                        sanitize::Bool=true,
                        remove_stereochemistry::Bool=false) -> Molecule

Perform comprehensive molecular standardization including salt stripping, 
tautomer canonicalization, and sanitization.

# Arguments
- `mol::Molecule`: Input molecule
- `strip_salts_flag::Bool=true`: Whether to strip salts and keep largest fragment
- `canonical_tautomer_flag::Bool=true`: Whether to canonicalize tautomers
- `sanitize::Bool=true`: Whether to sanitize the molecule
- `remove_stereochemistry::Bool=false`: Whether to remove stereochemistry information

# Returns
- `Molecule`: Standardized molecule

# Example
```julia
mol = mol_from_smiles("CC(O)=CC(=O)C.Na")  # Enol form with sodium salt
std_mol = standardize_molecule(mol)         # Returns clean, canonical form
```
"""
function standardize_molecule(
    mol::Molecule;
    strip_salts_flag::Bool=true,
    canonical_tautomer_flag::Bool=true,
    sanitize::Bool=true,
    remove_stereochemistry::Bool=false,
)
    if !mol.valid
        return mol
    end

    result_mol = mol

    try
        # Step 1: Strip salts if requested
        if strip_salts_flag
            result_mol = strip_salts(result_mol)
            if !result_mol.valid
                return result_mol
            end
        end

        # Step 2: Get canonical tautomer if requested
        if canonical_tautomer_flag
            result_mol = canonical_tautomer(result_mol)
            if !result_mol.valid
                return result_mol
            end
        end

        # Step 3: Remove stereochemistry if requested
        if remove_stereochemistry
            try
                # RemoveStereochemistry modifies the molecule in-place and returns None
                # So we need to make a copy first to avoid modifying the original
                mol_copy = _chem_standardize().Mol(result_mol._rdkit_mol)
                _chem_standardize().RemoveStereochemistry(mol_copy)
                result_mol = Molecule(;
                    _rdkit_mol=mol_copy,
                    valid=true,
                    source="no_stereo_$(result_mol.source)",
                    props=copy(result_mol.props),
                )
            catch e
                @warn "RemoveStereochemistry failed: $e"
                return Molecule(;
                    _rdkit_mol=result_mol._rdkit_mol,
                    valid=false,
                    source="stereo_removal_failed_$(result_mol.source)",
                )
            end
        end

        # Step 4: Sanitize if requested
        if sanitize
            # Only sanitize if we have a valid RDKit molecule
            if result_mol.valid && result_mol._rdkit_mol !== nothing
                try
                    _chem_standardize().SanitizeMol(result_mol._rdkit_mol)
                catch e
                    @warn "Sanitization failed: $e"
                    return Molecule(;
                        _rdkit_mol=result_mol._rdkit_mol,
                        valid=false,
                        source="sanitize_failed_$(result_mol.source)",
                    )
                end
            end
        end

        # Update source to reflect standardization
        new_source = "standardized_$(mol.source)"
        return Molecule(;
            _rdkit_mol=result_mol._rdkit_mol,
            valid=true,
            source=new_source,
            props=copy(mol.props),
        )

    catch e
        @warn "Standardization failed: $e"
        return Molecule(;
            _rdkit_mol=mol._rdkit_mol, valid=false, source="error_$(mol.source)"
        )
    end
end

"""
    neutralize_charges(mol::Molecule) -> Molecule

Neutralize charges in a molecule by adding/removing protons.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Neutralized molecule

# Example
```julia
mol = mol_from_smiles("CC(=O)[O-]")     # Acetate anion
neutral = neutralize_charges(mol)       # Returns acetic acid
```
"""
function neutralize_charges(mol::Molecule)
    if !mol.valid
        return mol
    end

    try
        # Get the charge parent
        uncharger = _standardize_module().rdMolStandardize.Uncharger()
        neutral_mol = uncharger.uncharge(mol._rdkit_mol)

        if pynot(neutral_mol)
            return Molecule(;
                _rdkit_mol=neutral_mol, valid=false, source="neutralized_$(mol.source)"
            )
        end

        return Molecule(;
            _rdkit_mol=neutral_mol,
            valid=true,
            source="neutralized_$(mol.source)",
            props=copy(mol.props),
        )
    catch e
        @warn "Failed to neutralize charges: $e"
        return mol
    end
end

"""
    normalize_molecule(mol::Molecule) -> Molecule

Apply RDKit's normalizer to fix common issues in molecular structures.

# Arguments
- `mol::Molecule`: Input molecule

# Returns
- `Molecule`: Normalized molecule

# Example
```julia
mol = mol_from_smiles("c1ccc2c(c1)ccc(=O)c2=O")  # Quinone
normalized = normalize_molecule(mol)               # Applies normalization rules
```
"""
function normalize_molecule(mol::Molecule)
    if !mol.valid
        return mol
    end

    try
        # Apply normalization
        normalizer = _standardize_module().rdMolStandardize.Normalizer()
        normalized_mol = normalizer.normalize(mol._rdkit_mol)

        if pynot(normalized_mol)
            return Molecule(;
                _rdkit_mol=normalized_mol, valid=false, source="normalized_$(mol.source)"
            )
        end

        return Molecule(;
            _rdkit_mol=normalized_mol,
            valid=true,
            source="normalized_$(mol.source)",
            props=copy(mol.props),
        )
    catch e
        @warn "Failed to normalize molecule: $e"
        return mol
    end
end
