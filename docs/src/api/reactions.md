```@meta
CurrentModule = MoleculeFlow
```

# Chemical Reactions

Functions for creating, analyzing, and manipulating chemical reactions.

## Reaction Type

```@docs
Reaction
```

## Creating Reactions

### From SMARTS

```@docs
reaction_from_smarts
```

### From RXN Files

```@docs
reaction_from_rxn_file
```

### From RXN Blocks

```@docs
reaction_from_rxn_block
```

## Converting Reactions

### To SMARTS

```@docs
reaction_to_smarts
```

### To RXN Block

```@docs
reaction_to_rxn_block
```

## Running Reactions

```@docs
run_reaction
```

## Reaction Validation

```@docs
validate_reaction!
```

```@docs
is_reaction_valid
```

```@docs
sanitize_reaction!
```

## Template Analysis

### Template Counts

```@docs
get_num_reactant_templates
```

```@docs
get_num_product_templates
```

### Template Access

```@docs
get_reactant_template
```

```@docs
get_product_template
```

### Substructure Matching

```@docs
has_reactant_substructure_match
```

```@docs
get_reacting_atoms
```

## Reaction Fingerprinting

### Difference Fingerprint

```@docs
reaction_fingerprint
```

### Structural Fingerprint

```@docs
reaction_structural_fingerprint
```

### Reaction Center Fingerprint

```@docs
reaction_center_fingerprint
```

### Similarity Analysis

```@docs
reaction_similarity
```

## Atom Mapping

### Get Mapping Numbers

```@docs
get_atom_mapping_numbers
```

### Set Mapping Numbers

```@docs
set_atom_mapping_numbers!
```

### Template Management

```@docs
remove_unmapped_reactant_templates!
```

```@docs
remove_unmapped_product_templates!
```

### Reaction Preprocessing

```@docs
preprocess_reaction!
```

```@docs
compute_atom_mapping!
```

### Agent Detection

```@docs
is_template_molecule_agent
```

## Reaction Analysis

### Information

```@docs
reaction_info
```

### Classification

```@docs
reaction_type_classification
```

### Complexity

```@docs
reaction_complexity
```

### Balance Check

```@docs
is_balanced
```

## Library Enumeration

```@docs
enumerate_library
```

## Database Operations

```@docs
find_similar_reactions
```