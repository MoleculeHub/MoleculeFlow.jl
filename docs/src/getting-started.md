# Getting Started

## Basic Usage

```julia
using MoleculeFlow

# Create a molecule from SMILES
mol = mol_from_smiles("CCO")  # Ethanol

# Check if it's valid
mol.valid  # true

# Calculate molecular properties
mw = molecular_weight(mol)     # 46.07 g/mol
logp_val = logp(mol)          # -0.31
hbd = num_hbd(mol)            # 1

# Generate fingerprints
fp = morgan_fingerprint(mol)   # ECFP4 fingerprint
maccs = maccs_fingerprint(mol) # MACCS keys

# Substructure search
has_oh = has_substructure_match(mol, "[OH]")  # true

# Convert back to SMILES
smiles = mol_to_smiles(mol)    # "CCO"
```

## Working with SMILES Containing Backslashes

When working with SMILES strings that contain backslashes (e.g., stereochemistry markers), you need to handle Julia's string parsing carefully:

```julia
# This will cause a ParseError:
# mol = mol_from_smiles("C\C=C\C")

# Use raw strings (recommended):
# mol = mol_from_smiles(raw"C\C=C\C")

# Or escape the backslashes:
# mol = mol_from_smiles("C\\C=C\\C")
```

**Why this happens**: Julia treats `\C` as an invalid escape sequence. Raw strings (`raw"..."`) tell Julia not to process escape sequences, while doubling backslashes (`\\`) creates a literal backslash.

## Working with Chemical Reactions

MoleculeFlow.jl provides comprehensive support for chemical reaction processing:

```julia
# Create a reaction from SMARTS
rxn = reaction_from_smarts("[C:1](=O)[O:2][C:3]>>[C:1](=O)[O-].[C:3][O+]")

# Apply reaction to molecules
reactant = mol_from_smiles("CC(=O)OCC")  # Ethyl acetate
if has_reactant_substructure_match(rxn, reactant)
    products = run_reaction(rxn, [reactant])
    for product_set in products
        for product in product_set
            println(mol_to_smiles(product))
        end
    end
end

# Analyze reaction properties
complexity = reaction_complexity(rxn)
classification = reaction_type_classification(rxn)
fingerprint = reaction_fingerprint(rxn)

# Compare reactions
rxn2 = reaction_from_smarts("[C:1][OH:2]>>[C:1][O-]")
similarity = reaction_similarity(rxn, rxn2)
```

!!! tip "Type Safety with Reactions"
    When building reaction databases, use typed arrays to avoid method dispatch errors:
    ```julia
    # Correct
    reaction_db = Reaction[]

    # Incorrect
    reaction_db = []  # Creates Vector{Any}
    ```