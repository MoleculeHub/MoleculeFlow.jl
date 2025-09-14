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

More detailed tutorials coming soon!