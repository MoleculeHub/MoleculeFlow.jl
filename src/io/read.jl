#######################################################
# Reading Mols
#######################################################

"""
    mol_from_smiles(smiles::String) -> Molecule

Convert a SMILES string to a Molecule object.

# Arguments
- `smiles::String`: A valid SMILES (Simplified Molecular Input Line Entry System) string

# Returns
- `Molecule`: A Molecule object. Check `mol.valid` to see if parsing was successful.

# Examples
```julia
mol = mol_from_smiles("CCO")  # Ethanol
if mol.valid
    println("Successfully parsed molecule")
end
```

# Notes
- Invalid SMILES strings will return a Molecule with `valid=false`
- **Important**: For SMILES containing backslashes,
  use `raw"..."` strings to avoid Julia's escape sequence parsing, or double the backslashes `\\`
"""
@inline function mol_from_smiles(smiles::String)
    # Handle escaped backslashes in SMILES strings (e.g., \\C=C\\ -> \C=C\)
    cleaned_smiles = replace(smiles, "\\\\" => "\\")
    rdkit_mol = _mol_from_smiles(cleaned_smiles)
    if pynot(rdkit_mol)
        return Molecule(; _rdkit_mol=rdkit_mol, valid=false, source=smiles)
    end
    return Molecule(; _rdkit_mol=rdkit_mol, valid=true, source=smiles)
end

"""
    mol_from_smiles(smiles_list::Vector{String}) -> Vector{Molecule}

Convert a vector of SMILES strings to a vector of Molecule objects.

# Arguments
- `smiles_list::Vector{String}`: Vector of SMILES strings

# Returns
- `Vector{Molecule}`: Vector of Molecule objects

# Examples
```julia
smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
mols = mol_from_smiles(smiles)
valid_mols = filter(m -> m.valid, mols)
```
"""
function mol_from_smiles(smiles_list::Vector{String})
    results = Vector{Molecule}(undef, length(smiles_list))

    @inbounds for i in eachindex(smiles_list)
        results[i] = mol_from_smiles(smiles_list[i])
    end

    return results
end

"""
    mol_from_molblock(molblock::String) -> Molecule

Convert a MOL block string (MDL format) to a Molecule object.

# Arguments
- `molblock::String`: A valid MOL block string in MDL format

# Returns
- `Molecule`: A Molecule object. Check `mol.valid` to see if parsing was successful.

# Examples
```julia
molblock = \"\"\"

  RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
\"\"\"
mol = mol_from_molblock(molblock)
```

# Notes
- Invalid MOL blocks will return a Molecule with `valid=false`
- Commonly used for reading individual molecular structures from files
"""
function mol_from_molblock(molblock::String)
    rdkit_mol = _mol_from_molblock(molblock)
    if pynot(rdkit_mol)
        return Molecule(; _rdkit_mol=rdkit_mol, valid=false, source="molblock")
    end
    return Molecule(; _rdkit_mol=rdkit_mol, valid=true, source="molblock")
end

"""
    read_sdf(filename::String; max_mols::Union{Int,Nothing}=nothing,
            include_props::Bool=true, sanitize::Bool=true) -> Vector{Molecule}

Read molecules from an SDF (Structure Data File) format file.

# Arguments
- `filename::String`: Path to the SDF file
- `max_mols::Union{Int,Nothing}=nothing`: Maximum number of molecules to read (nothing for all)
- `include_props::Bool=true`: Whether to include SDF properties in molecule props
- `sanitize::Bool=true`: Whether to sanitize molecules during reading

# Returns
- `Vector{Molecule}`: Vector of Molecule objects read from the SDF file

# Examples
```julia
# Read all molecules from an SDF file
molecules = read_sdf("compounds.sdf")
valid_mols = filter(mol -> mol.valid, molecules)

# Read only first 100 molecules
molecules = read_sdf("large_database.sdf", max_mols=100)

# Read without properties for faster parsing
molecules = read_sdf("structures_only.sdf", include_props=false)
```

# Notes
- Invalid molecules in the SDF will be included with `valid=false`
- SDF properties are stored in the `props` field of each molecule
- Large SDF files are processed efficiently using RDKit's supplier mechanism
"""
function read_sdf(filename::String;
                  max_mols::Union{Int,Nothing}=nothing,
                  include_props::Bool=true,
                  sanitize::Bool=true)
    if !isfile(filename)
        throw(ArgumentError("File not found: $filename"))
    end

    try
        # Create SDF supplier
        supplier = _sdf_supplier(filename)

        molecules = Molecule[]
        mol_count = 0

        for rdkit_mol in supplier
            # Check if we've reached the maximum number of molecules
            if max_mols !== nothing && mol_count >= max_mols
                break
            end

            mol_count += 1

            if pynot(rdkit_mol)
                # Invalid molecule - still add it with valid=false
                push!(molecules, Molecule(; _rdkit_mol=rdkit_mol, valid=false, source="sdf_mol_$mol_count"))
                continue
            end

            # Create molecule properties dictionary
            mol_props = Dict{Symbol, Any}()

            if include_props
                # Extract SDF properties
                try
                    prop_names = rdkit_mol.GetPropNames()
                    for prop_name in prop_names
                        prop_value = rdkit_mol.GetProp(pyconvert(String, prop_name))
                        # Convert property name to symbol and store value
                        mol_props[Symbol(pyconvert(String, prop_name))] = pyconvert(String, prop_value)
                    end
                catch e
                    @warn "Error reading properties for molecule $mol_count: $e"
                end
            end

            # Add molecule index as property
            mol_props[:sdf_index] = mol_count
            mol_props[:source_file] = filename

            # Create Molecule object
            mol = Molecule(; _rdkit_mol=rdkit_mol, valid=true, source="sdf_mol_$mol_count", props=mol_props)
            push!(molecules, mol)
        end

        return molecules

    catch e
        throw(ArgumentError("Error reading SDF file '$filename': $e"))
    end
end

"""
    read_sdf_lazy(filename::String; sanitize::Bool=true) -> Function

Create a lazy iterator for reading molecules from an SDF file one at a time.

# Arguments
- `filename::String`: Path to the SDF file
- `sanitize::Bool=true`: Whether to sanitize molecules during reading

# Returns
- `Function`: A function that returns the next molecule when called, or `nothing` when done

# Examples
```julia
# Process large SDF file without loading all molecules into memory
next_mol = read_sdf_lazy("huge_database.sdf")

while true
    mol = next_mol()
    if mol === nothing
        break  # End of file
    end

    if mol.valid
        # Process the molecule
        println("Processing: \$(mol_to_smiles(mol))")
    end
end
```

# Notes
- Memory efficient for large SDF files
- Returns `nothing` when no more molecules are available
- Each molecule includes SDF properties in its `props` field
"""
function read_sdf_lazy(filename::String; sanitize::Bool=true)
    if !isfile(filename)
        throw(ArgumentError("File not found: $filename"))
    end

    supplier = _sdf_supplier(filename)
    supplier_iter = pyiter(supplier)
    mol_count = Ref(0)  # Use Ref to make it mutable in closure

    function next_molecule()
        try
            rdkit_mol = pynext(supplier_iter)
            mol_count[] += 1

            if pynot(rdkit_mol)
                # Invalid molecule
                return Molecule(; _rdkit_mol=rdkit_mol, valid=false, source="sdf_mol_$(mol_count[])")
            end

            # Extract SDF properties
            mol_props = Dict{Symbol, Any}()
            try
                prop_names = rdkit_mol.GetPropNames()
                for prop_name in prop_names
                    prop_value = rdkit_mol.GetProp(pyconvert(String, prop_name))
                    mol_props[Symbol(pyconvert(String, prop_name))] = pyconvert(String, prop_value)
                end
            catch e
                @warn "Error reading properties for molecule $(mol_count[]): $e"
            end

            mol_props[:sdf_index] = mol_count[]
            mol_props[:source_file] = filename

            return Molecule(; _rdkit_mol=rdkit_mol, valid=true, source="sdf_mol_$(mol_count[])", props=mol_props)

        catch e
            if pyisinstance(e, pybuiltins.StopIteration)
                return nothing
            else
                @warn "Error reading molecule $(mol_count[]): $e"
                return nothing
            end
        end
    end

    return next_molecule
end

"""
    mol_from_inchi(inchi::String) -> Molecule

Convert an InChI string to a Molecule object.

# Arguments
- `inchi::String`: A valid InChI (International Chemical Identifier) string

# Returns
- `Molecule`: A Molecule object. Check `mol.valid` to see if parsing was successful.

# Examples
```julia
mol = mol_from_inchi("InChI=1S/C2H6/c1-2/h1-2H3")  # Ethane
if mol.valid
    println("Successfully parsed molecule")
end

mol2 = mol_from_inchi("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")  # Benzene
```

# Notes
- Invalid InChI strings will return a Molecule with `valid=false`
- The original InChI string is stored in the molecule's props as `:InChI`
"""
function mol_from_inchi(inchi::String)
    rdkit_mol = _mol_from_inchi(inchi)

    if pynot(rdkit_mol)
        return Molecule(; _rdkit_mol=rdkit_mol, valid=false, source=inchi, props=Dict(:InChI => inchi))
    end

    return Molecule(; _rdkit_mol=rdkit_mol, valid=true, source=inchi, props=Dict(:InChI => inchi))
end

"""
    mol_from_inchi(inchi_list::Vector{String}) -> Vector{Molecule}

Convert a vector of InChI strings to a vector of Molecule objects.

# Arguments
- `inchi_list::Vector{String}`: Vector of InChI strings

# Returns
- `Vector{Molecule}`: Vector of Molecule objects

# Examples
```julia
inchis = [
    "InChI=1S/C2H6/c1-2/h1-2H3",     # Ethane
    "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3",  # Propane
    "invalid_inchi",                  # Invalid
    "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"  # Benzene
]
mols = mol_from_inchi(inchis)
valid_mols = filter(m -> m.valid, mols)
```
"""
function mol_from_inchi(inchi_list::Vector{String})
    results = Vector{Molecule}(undef, length(inchi_list))

    @inbounds for i in eachindex(inchi_list)
        inchi = inchi_list[i]
        results[i] = mol_from_inchi(inchi)
    end

    return results
end
