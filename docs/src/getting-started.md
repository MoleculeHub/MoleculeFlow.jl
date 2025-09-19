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
mol = mol_from_smiles(raw"C\C=C\C")

# Or escape the backslashes:
mol = mol_from_smiles("C\\C=C\\C")

# Complex example with stereochemistry:
mol = mol_from_smiles(raw"CN(C)C\C=C\C(=O)Nc3cc1c(Nc(cc2Cl)ccc2F)ncnc1cc3OC4COCC4")
```

**Why this happens**: Julia treats `\C` as an invalid escape sequence. Raw strings (`raw"..."`) tell Julia not to process escape sequences, while doubling backslashes (`\\`) creates a literal backslash.