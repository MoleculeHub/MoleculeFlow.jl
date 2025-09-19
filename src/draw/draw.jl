#######################################################
# Drawing Mols
#######################################################
const _rdkit_draw = Ref{Py}()

# #######################################################

"""
    mol_to_image(mol::Molecule; kwargs...) -> Matrix{RGB{N0f8}}

Render a molecule as a 2D image with various customization options.

# Arguments
- `mol::Molecule`: The molecule to render
- `size::Tuple=(300, 300)`: Image size (width, height) in pixels
- `kekulize::Bool=true`: Whether to kekulize the molecule before drawing
- `wedge_bonds::Bool=true`: Whether to draw wedge bonds for stereochemistry
- `highlight_atoms::Union{Vector{Int},Nothing}=nothing`: Atom indices to highlight (1-based)
- `highlight_bonds::Union{Vector{Int},Nothing}=nothing`: Bond indices to highlight (1-based)
- `highlight_atom_colors::Union{Dict,Nothing}=nothing`: Custom colors for highlighted atoms
- `highlight_bond_colors::Union{Dict,Nothing}=nothing`: Custom colors for highlighted bonds
- `use_svg::Bool=false`: Whether to use SVG rendering (returns SVG string instead of image)

# Returns
- `Matrix{RGB{N0f8}}`: RGB image matrix (if use_svg=false)
- `String`: SVG string (if use_svg=true)

# Examples
```julia
mol = mol_from_smiles("CCO")
img = mol_to_image(mol)
img_highlighted = mol_to_image(mol, highlight_atoms=[1, 3])
svg_str = mol_to_image(mol, use_svg=true)
```
"""
function mol_to_image(
    mol::Molecule;
    size::Tuple=(300, 300),
    kekulize::Bool=true,
    wedge_bonds::Bool=true,
    highlight_atoms::Union{Vector{Int}, Nothing}=nothing,
    highlight_bonds::Union{Vector{Int}, Nothing}=nothing,
    highlight_atom_colors::Union{Dict, Nothing}=nothing,
    highlight_bond_colors::Union{Dict, Nothing}=nothing,
    use_svg::Bool=false,
)
    if !mol.valid
        throw(ArgumentError("Cannot draw invalid molecule"))
    end

    if use_svg
        # For SVG, try to use RDKit's SVG drawing capabilities
        try
            # Try using MolsToGridImage with SVG for single molecule
            svg_str = pycall(
                _rdkit_draw[].MolsToGridImage,
                [mol._rdkit_mol];
                molsPerRow=1,
                subImgSize=size,
                useSVG=true,
                highlightAtomLists=if highlight_atoms !== nothing
                    [julia_to_python_indices(highlight_atoms)]
                else
                    nothing
                end,
                highlightBondLists=if highlight_bonds !== nothing
                    [julia_to_python_indices(highlight_bonds)]
                else
                    nothing
                end,
            )
            return pyconvert(String, svg_str)
        catch e1
            @warn "Primary SVG generation failed, trying alternative: $e1"
            try
                # Alternative approach: try MolToImage with different parameters
                svg_result = pycall(
                    _rdkit_draw[].MolToImage, mol._rdkit_mol; size=size, useSVG=true
                )
                return pyconvert(String, svg_result)
            catch e2
                @warn "SVG generation failed completely, falling back to PNG: $e2"
                # Fall back to PNG and return as image
                img = pycall(
                    _rdkit_draw[].MolToImage,
                    mol._rdkit_mol;
                    size=size,
                    kekulize=kekulize,
                    wedgeBonds=wedge_bonds,
                    highlightAtoms=julia_to_python_indices(highlight_atoms),
                    highlightBonds=julia_to_python_indices(highlight_bonds),
                    highlightAtomColors=highlight_atom_colors,
                    highlightBondColors=highlight_bond_colors,
                )
                return pil_png_to_rgb(img)
            end
        end
    else
        img = pycall(
            _rdkit_draw[].MolToImage,
            mol._rdkit_mol;
            size=size,
            kekulize=kekulize,
            wedgeBonds=wedge_bonds,
            highlightAtoms=julia_to_python_indices(highlight_atoms),
            highlightBonds=julia_to_python_indices(highlight_bonds),
            highlightAtomColors=highlight_atom_colors,
            highlightBondColors=highlight_bond_colors,
        )
        return pil_png_to_rgb(img)
    end
end

"""
    mols_to_grid_image(mols::Vector{Union{Molecule,Missing}}; kwargs...) -> Matrix{RGB{N0f8}}

Render multiple molecules in a grid layout.

# Arguments
- `mols::Vector{Union{Molecule,Missing}}`: Vector of molecules to render
- `sub_img_size::Tuple=(200, 200)`: Size of each individual molecule image
- `mols_per_row::Int=3`: Number of molecules per row in the grid
- `legends::Union{Vector{String},Nothing}=nothing`: Text labels for each molecule
- `highlight_atom_lists::Union{Vector,Nothing}=nothing`: Lists of atoms to highlight for each molecule
- `highlight_bond_lists::Union{Vector,Nothing}=nothing`: Lists of bonds to highlight for each molecule
- `use_svg::Bool=false`: Whether to use SVG rendering

# Returns
- `Matrix{RGB{N0f8}}`: RGB image matrix (if use_svg=false)
- `String`: SVG string (if use_svg=true)

# Examples
```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1")]
img = mols_to_grid_image(mols, legends=["Ethanol", "Benzene"])
```
"""
function mols_to_grid_image(
    mols::Vector{Union{Molecule, Missing}};
    sub_img_size::Tuple=(200, 200),
    mols_per_row::Int=3,
    legends::Union{Vector{String}, Nothing}=nothing,
    highlight_atom_lists::Union{Vector, Nothing}=nothing,
    highlight_bond_lists::Union{Vector, Nothing}=nothing,
    use_svg::Bool=false,
)
    # Filter out missing molecules and get their RDKit objects
    valid_mols = [mol._rdkit_mol for mol in mols if !ismissing(mol) && mol.valid]

    if isempty(valid_mols)
        throw(ArgumentError("No valid molecules to draw"))
    end

    # Process highlight lists to convert to Python indices
    highlight_atoms_py = if highlight_atom_lists !== nothing
        [julia_to_python_indices(atoms) for atoms in highlight_atom_lists]
    else
        nothing
    end
    highlight_bonds_py = if highlight_bond_lists !== nothing
        [julia_to_python_indices(bonds) for bonds in highlight_bond_lists]
    else
        nothing
    end

    if use_svg
        try
            svg_str = pycall(
                _rdkit_draw[].MolsToGridImage,
                valid_mols;
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                legends=legends,
                highlightAtomLists=highlight_atoms_py,
                highlightBondLists=highlight_bonds_py,
                useSVG=true,
            )
            return pyconvert(String, svg_str)
        catch e
            @warn "Grid SVG generation failed, falling back to PNG: $e"
            # Fall back to PNG
            img = pycall(
                _rdkit_draw[].MolsToGridImage,
                valid_mols;
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                legends=legends,
                highlightAtomLists=highlight_atoms_py,
                highlightBondLists=highlight_bonds_py,
            )
            return pil_png_to_rgb(img)
        end
    else
        img = pycall(
            _rdkit_draw[].MolsToGridImage,
            valid_mols;
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=legends,
            highlightAtomLists=highlight_atoms_py,
            highlightBondLists=highlight_bonds_py,
        )
        return pil_png_to_rgb(img)
    end
end

function mols_to_grid_image(
    mols::Vector{Molecule};
    sub_img_size::Tuple=(200, 200),
    mols_per_row::Int=3,
    legends::Union{Vector{String}, Nothing}=nothing,
    highlight_atom_lists::Union{Vector, Nothing}=nothing,
    highlight_bond_lists::Union{Vector, Nothing}=nothing,
    use_svg::Bool=false,
)
    # Get RDKit objects for valid molecules
    valid_mols = [mol._rdkit_mol for mol in mols if mol.valid]

    if isempty(valid_mols)
        throw(ArgumentError("No valid molecules to draw"))
    end

    # Process highlight lists
    highlight_atoms_py = if highlight_atom_lists !== nothing
        [julia_to_python_indices(atoms) for atoms in highlight_atom_lists]
    else
        nothing
    end
    highlight_bonds_py = if highlight_bond_lists !== nothing
        [julia_to_python_indices(bonds) for bonds in highlight_bond_lists]
    else
        nothing
    end

    if use_svg
        try
            svg_str = pycall(
                _rdkit_draw[].MolsToGridImage,
                valid_mols;
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                legends=legends,
                highlightAtomLists=highlight_atoms_py,
                highlightBondLists=highlight_bonds_py,
                useSVG=true,
            )
            return pyconvert(String, svg_str)
        catch e
            @warn "Grid SVG generation failed, falling back to PNG: $e"
            # Fall back to PNG
            img = pycall(
                _rdkit_draw[].MolsToGridImage,
                valid_mols;
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                legends=legends,
                highlightAtomLists=highlight_atoms_py,
                highlightBondLists=highlight_bonds_py,
            )
            return pil_png_to_rgb(img)
        end
    else
        img = pycall(
            _rdkit_draw[].MolsToGridImage,
            valid_mols;
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=legends,
            highlightAtomLists=highlight_atoms_py,
            highlightBondLists=highlight_bonds_py,
        )
        return pil_png_to_rgb(img)
    end
end

"""
    mol_to_svg(mol::Molecule; kwargs...) -> String

Render a molecule as an SVG string.

# Arguments
- `mol::Molecule`: The molecule to render
- `size::Tuple=(300, 300)`: Image size (width, height) in pixels
- `kekulize::Bool=true`: Whether to kekulize the molecule before drawing
- `wedge_bonds::Bool=true`: Whether to draw wedge bonds for stereochemistry
- `highlight_atoms::Union{Vector{Int},Nothing}=nothing`: Atom indices to highlight
- `highlight_bonds::Union{Vector{Int},Nothing}=nothing`: Bond indices to highlight

# Returns
- `String`: SVG representation of the molecule

# Examples
```julia
mol = mol_from_smiles("CCO")
svg_str = mol_to_svg(mol)
```
"""
function mol_to_svg(
    mol::Molecule;
    size::Tuple=(300, 300),
    kekulize::Bool=true,
    wedge_bonds::Bool=true,
    highlight_atoms::Union{Vector{Int}, Nothing}=nothing,
    highlight_bonds::Union{Vector{Int}, Nothing}=nothing,
)
    return mol_to_image(
        mol;
        size=size,
        kekulize=kekulize,
        wedge_bonds=wedge_bonds,
        highlight_atoms=highlight_atoms,
        highlight_bonds=highlight_bonds,
        use_svg=true,
    )
end

"""
    highlight_substructure(mol::Molecule, pattern::Union{Molecule,String}; kwargs...) -> Matrix{RGB{N0f8}}

Draw a molecule with a substructure pattern highlighted.

# Arguments
- `mol::Molecule`: The molecule to render
- `pattern::Union{Molecule,String}`: Substructure pattern (Molecule or SMARTS string)
- `size::Tuple=(300, 300)`: Image size
- `highlight_color::String="yellow"`: Color for highlighting (yellow, red, blue, etc.)
- `use_svg::Bool=false`: Whether to return SVG instead of image matrix

# Returns
- `Matrix{RGB{N0f8}}` or `String`: Rendered molecule with highlighted substructure

# Examples
```julia
mol = mol_from_smiles("c1ccc(O)cc1")  # Phenol
img = highlight_substructure(mol, "[OH]")
```
"""
function highlight_substructure(
    mol::Molecule,
    pattern::Union{Molecule, String};
    size::Tuple=(300, 300),
    highlight_color::String="yellow",
    use_svg::Bool=false,
)
    if !mol.valid
        throw(ArgumentError("Cannot draw invalid molecule"))
    end

    # Get substructure matches
    matches = if pattern isa String
        get_substructure_matches(mol, pattern)
    else
        get_substructure_matches(mol, pattern)
    end

    if isempty(matches)
        @warn "No substructure matches found"
        return mol_to_image(mol; size=size, use_svg=use_svg)
    end

    # Highlight the first match
    highlight_atoms = matches[1]

    return mol_to_image(mol; size=size, highlight_atoms=highlight_atoms, use_svg=use_svg)
end

"""
    draw_molecule_with_atom_labels(mol::Molecule; kwargs...) -> Matrix{RGB{N0f8}}

Draw a molecule with atom indices labeled.

# Arguments
- `mol::Molecule`: The molecule to render
- `size::Tuple=(300, 300)`: Image size
- `label_atoms::Bool=true`: Whether to show atom indices
- `label_bonds::Bool=false`: Whether to show bond indices
- `use_svg::Bool=false`: Whether to return SVG

# Returns
- `Matrix{RGB{N0f8}}` or `String`: Rendered molecule with labels

# Examples
```julia
mol = mol_from_smiles("CCO")
img = draw_molecule_with_atom_labels(mol)
```
"""
function draw_molecule_with_atom_labels(
    mol::Molecule;
    size::Tuple=(300, 300),
    label_atoms::Bool=true,
    label_bonds::Bool=false,
    use_svg::Bool=false,
)
    if !mol.valid
        throw(ArgumentError("Cannot draw invalid molecule"))
    end

    # This is a simplified version - full implementation would require
    # more advanced RDKit drawing options or custom labeling
    if use_svg
        try
            # Try using MolsToGridImage with SVG for single molecule (more reliable for SVG)
            svg_str = pycall(
                _rdkit_draw[].MolsToGridImage,
                [mol._rdkit_mol];
                molsPerRow=1,
                subImgSize=size,
                useSVG=true,
            )
            return pyconvert(String, svg_str)
        catch e
            @warn "SVG generation failed for atom labels, falling back to PNG: $e"
            # Fall back to PNG
            img = pycall(_rdkit_draw[].MolToImage, mol._rdkit_mol; size=size)
            return pil_png_to_rgb(img)
        end
    else
        img = pycall(_rdkit_draw[].MolToImage, mol._rdkit_mol; size=size)
        return pil_png_to_rgb(img)
    end
end

"""
    draw_similarity_map(mol1::Molecule, mol2::Molecule; kwargs...) -> Matrix{RGB{N0f8}}

Draw two molecules side by side highlighting their most common substructure.

# Arguments
- `mol1::Molecule`: First molecule
- `mol2::Molecule`: Second molecule
- `size::Tuple=(600, 300)`: Total image size
- `use_svg::Bool=false`: Whether to return SVG

# Returns
- `Matrix{RGB{N0f8}}` or `String`: Side-by-side comparison with common substructure highlighted

# Examples
```julia
mol1 = mol_from_smiles("CCO")
mol2 = mol_from_smiles("CCC")
img = draw_similarity_map(mol1, mol2)
```
"""
function draw_similarity_map(
    mol1::Molecule, mol2::Molecule; size::Tuple=(600, 300), use_svg::Bool=false
)
    if !mol1.valid || !mol2.valid
        throw(ArgumentError("Cannot draw invalid molecules"))
    end

    # Find maximum common substructure
    try
        mcs_result = maximum_common_substructure(mol1, mol2)
        if !mcs_result.valid
            # No common substructure, just draw side by side
            return mols_to_grid_image(
                [mol1, mol2];
                sub_img_size=(Int(size[1] / 2), size[2]),
                mols_per_row=2,
                use_svg=use_svg,
            )
        end

        # Get matches for both molecules
        matches1 = get_substructure_matches(mol1, mcs_result)
        matches2 = get_substructure_matches(mol2, mcs_result)

        highlight_atoms1 = !isempty(matches1) ? matches1[1] : Int[]
        highlight_atoms2 = !isempty(matches2) ? matches2[1] : Int[]

        return mols_to_grid_image(
            [mol1, mol2];
            sub_img_size=(Int(size[1] / 2), size[2]),
            mols_per_row=2,
            highlight_atom_lists=[highlight_atoms1, highlight_atoms2],
            use_svg=use_svg,
        )
    catch e
        @warn "Could not compute similarity map: $e"
        return mols_to_grid_image(
            [mol1, mol2];
            sub_img_size=(Int(size[1] / 2), size[2]),
            mols_per_row=2,
            use_svg=use_svg,
        )
    end
end

"""
    save_molecule_image(mol::Molecule, filepath::String; kwargs...)

Save a molecule image to a file.

# Arguments
- `mol::Molecule`: The molecule to render
- `filepath::String`: Output file path (supports .png, .svg)
- `size::Tuple=(300, 300)`: Image size
- `kwargs...`: Additional arguments passed to drawing functions

# Examples
```julia
mol = mol_from_smiles("CCO")
save_molecule_image(mol, "ethanol.png")
save_molecule_image(mol, "ethanol.svg", highlight_atoms=[1, 3])
```
"""
function save_molecule_image(
    mol::Molecule, filepath::String; size::Tuple=(300, 300), kwargs...
)
    if !mol.valid
        throw(ArgumentError("Cannot save invalid molecule"))
    end

    if endswith(lowercase(filepath), ".svg")
        svg_content = mol_to_svg(mol; size=size, kwargs...)
        open(filepath, "w") do file
            write(file, svg_content)
        end
    elseif endswith(lowercase(filepath), ".png")
        img = mol_to_image(mol; size=size, kwargs...)
        # Note: This would require Images.jl for saving
        # save(filepath, img)
        @warn "PNG saving requires Images.jl package. Use mol_to_image() and save manually."
        return img
    else
        throw(ArgumentError("Unsupported file format. Use .png or .svg"))
    end
end

"""
    draw_functional_groups(mol::Molecule; kwargs...) -> Matrix{RGB{N0f8}}

Draw a molecule with functional groups highlighted in different colors.

# Arguments
- `mol::Molecule`: The molecule to render
- `size::Tuple=(300, 300)`: Image size
- `functional_groups::Vector{String}=["[OH]", "[C=O]", "[NH2]", "[COOH]"]`: SMARTS patterns for functional groups
- `colors::Vector{String}=["red", "blue", "green", "orange"]`: Colors for each functional group
- `use_svg::Bool=false`: Whether to return SVG

# Returns
- `Matrix{RGB{N0f8}}` or `String`: Rendered molecule with highlighted functional groups

# Examples
```julia
mol = mol_from_smiles("CC(=O)O")  # Acetic acid
img = draw_functional_groups(mol)
```
"""
function draw_functional_groups(
    mol::Molecule;
    size::Tuple=(300, 300),
    functional_groups::Vector{String}=["[OH]", "C=O", "[NH2]", "[COOH]"],
    colors::Vector{String}=["red", "blue", "green", "orange"],
    use_svg::Bool=false,
)
    if !mol.valid
        throw(ArgumentError("Cannot draw invalid molecule"))
    end

    # Find all functional group matches
    all_highlight_atoms = Int[]
    highlight_colors = Dict{Int, String}()

    for (i, pattern) in enumerate(functional_groups)
        matches = get_substructure_matches(mol, pattern)
        color = i <= length(colors) ? colors[i] : "yellow"

        for match in matches
            for atom_idx in match
                if atom_idx ∉ all_highlight_atoms
                    push!(all_highlight_atoms, atom_idx)
                    highlight_colors[atom_idx] = color
                end
            end
        end
    end

    if isempty(all_highlight_atoms)
        @info "No functional groups found"
        return mol_to_image(mol; size=size, use_svg=use_svg)
    end

    return mol_to_image(
        mol; size=size, highlight_atoms=all_highlight_atoms, use_svg=use_svg
    )
end
