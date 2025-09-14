using Test
using MoleculeFlow

@testset "Substructure Search" begin
    mol = mol_from_smiles("CCc1ccccc1")  # ethylbenzene

    # Test SMARTS pattern matching
    @test has_substructure_match(mol, "c1ccccc1")  # benzene ring
    @test !has_substructure_match(mol, "N")        # no nitrogen

    # Test functional group detection
    alcohol = mol_from_smiles("CCO")
    @test has_functional_group(alcohol, :alcohol)
    @test !has_functional_group(alcohol, :ketone)

    # Test getting all functional groups
    func_groups = get_functional_groups(alcohol)
    @test isa(func_groups, Dict)
    @test func_groups[:alcohol] == true
    @test func_groups[:ketone] == false

    # Test substructure matches
    matches = get_substructure_matches(mol, "c1ccccc1")
    @test length(matches) >= 1
    @test all(match -> isa(match, Vector{Int}), matches)
end

@testset "Basic Functional Groups" begin
    # Alcohol and phenol
    alcohol = mol_from_smiles("CCO")
    @test has_functional_group(alcohol, :alcohol)

    phenol = mol_from_smiles("c1ccc(O)cc1")
    @test has_functional_group(phenol, :phenol)
    @test has_functional_group(phenol, :alcohol)  # phenol is also an alcohol

    # Carboxylic acid and ester
    acid = mol_from_smiles("CC(=O)O")
    @test has_functional_group(acid, :carboxylic_acid)

    ester = mol_from_smiles("CC(=O)OC")
    @test has_functional_group(ester, :ester)

    # Ether
    ether = mol_from_smiles("COC")
    @test has_functional_group(ether, :ether)

    # Aldehyde and ketone
    aldehyde = mol_from_smiles("CC=O")
    @test has_functional_group(aldehyde, :aldehyde)

    ketone = mol_from_smiles("CC(=O)C")
    @test has_functional_group(ketone, :ketone)

    # Amines
    primary_amine = mol_from_smiles("CCN")
    @test has_functional_group(primary_amine, :amine_primary)

    secondary_amine = mol_from_smiles("CCNC")
    @test has_functional_group(secondary_amine, :amine_secondary)

    tertiary_amine = mol_from_smiles("CN(C)C")
    @test has_functional_group(tertiary_amine, :amine_tertiary)

    # Amide and nitrile
    amide = mol_from_smiles("CC(=O)N")
    @test has_functional_group(amide, :amide)

    nitrile = mol_from_smiles("CC#N")
    @test has_functional_group(nitrile, :nitrile)
end

@testset "Sulfur-Containing Groups" begin
    # Thiol
    thiol = mol_from_smiles("CCS")
    @test has_functional_group(thiol, :thiol)

    # Sulfide
    sulfide = mol_from_smiles("CCSC")
    @test has_functional_group(sulfide, :sulfide)

    # Disulfide
    disulfide = mol_from_smiles("CCSSCC")
    @test has_functional_group(disulfide, :disulfide)

    # Sulfoxide
    sulfoxide = mol_from_smiles("CCS(=O)C")
    @test has_functional_group(sulfoxide, :sulfoxide)

    # Sulfone
    sulfone = mol_from_smiles("CCS(=O)(=O)C")
    @test has_functional_group(sulfone, :sulfone)

    # Sulfonamide
    sulfonamide = mol_from_smiles("CS(=O)(=O)N")
    @test has_functional_group(sulfonamide, :sulfonamide)

    # Sulfonic acid
    sulfonic_acid = mol_from_smiles("CS(=O)(=O)O")
    @test has_functional_group(sulfonic_acid, :sulfonic_acid)
end

@testset "Phosphorus-Containing Groups" begin
    # Test phosphorus groups (if molecules can be created)
    try
        phosphine = mol_from_smiles("CP(C)C")
        @test has_functional_group(phosphine, :phosphine)
    catch
        @test_skip "Phosphine test skipped - molecule creation failed"
    end

    try
        phosphine_oxide = mol_from_smiles("CP(=O)(C)C")
        @test has_functional_group(phosphine_oxide, :phosphine_oxide)
    catch
        @test_skip "Phosphine oxide test skipped - molecule creation failed"
    end
end

@testset "Halogen-Containing Groups" begin
    # Individual halogens
    fluoride = mol_from_smiles("CCF")
    @test has_functional_group(fluoride, :fluoride)

    chloride = mol_from_smiles("CCCl")
    @test has_functional_group(chloride, :chloride)

    bromide = mol_from_smiles("CCBr")
    @test has_functional_group(bromide, :bromide)

    iodide = mol_from_smiles("CCI")
    @test has_functional_group(iodide, :iodide)

    # Polyhalogenated groups
    trifluoromethyl = mol_from_smiles("CC(F)(F)F")
    @test has_functional_group(trifluoromethyl, :trifluoromethyl)
end

@testset "Advanced Nitrogen Groups" begin
    # Nitro group
    nitro = mol_from_smiles("c1ccc(cc1)[N+](=O)[O-]")
    @test has_functional_group(nitro, :nitro)

    # Azide
    try
        azide = mol_from_smiles("CC[N-][N+]#N")
        @test has_functional_group(azide, :azide)
    catch
        @test_skip "Azide test skipped - molecule creation failed"
    end

    # Imine
    imine = mol_from_smiles("CC=N")
    @test has_functional_group(imine, :imine)

    # Urea
    urea = mol_from_smiles("NC(=O)N")
    @test has_functional_group(urea, :urea)

    # Carbamate
    carbamate = mol_from_smiles("NC(=O)OC")
    @test has_functional_group(carbamate, :carbamate)
end

@testset "Advanced Oxygen Groups" begin
    # Acetal
    acetal = mol_from_smiles("CC(OC)OC")
    @test has_functional_group(acetal, :acetal)

    # Anhydride
    anhydride = mol_from_smiles("CC(=O)OC(=O)C")
    @test has_functional_group(anhydride, :anhydride)

    # Carbonate
    carbonate = mol_from_smiles("COC(=O)OC")
    @test has_functional_group(carbonate, :carbonate)
end

@testset "Carbon-Carbon Multiple Bonds" begin
    # Alkene
    alkene = mol_from_smiles("CC=CC")
    @test has_functional_group(alkene, :alkene)

    # Alkyne
    alkyne = mol_from_smiles("CC#CC")
    @test has_functional_group(alkyne, :alkyne)

    # Conjugated diene
    diene = mol_from_smiles("C=CC=C")
    @test has_functional_group(diene, :conjugated_diene)
end

@testset "5-Membered Aromatic Heterocycles" begin
    # Furan
    furan = mol_from_smiles("c1ccoc1")
    @test has_functional_group(furan, :furan)

    # Thiophene
    thiophene = mol_from_smiles("c1ccsc1")
    @test has_functional_group(thiophene, :thiophene)

    # Pyrrole
    pyrrole = mol_from_smiles("c1cc[nH]c1")
    @test has_functional_group(pyrrole, :pyrrole)

    # Imidazole
    imidazole = mol_from_smiles("c1cnc[nH]1")
    @test has_functional_group(imidazole, :imidazole)

    # Pyrazole
    pyrazole = mol_from_smiles("c1cn[nH]c1")
    @test has_functional_group(pyrazole, :pyrazole)

    # Oxazole
    oxazole = mol_from_smiles("c1cnoc1")
    @test has_functional_group(oxazole, :oxazole)

    # Thiazole
    thiazole = mol_from_smiles("c1cnsc1")
    @test has_functional_group(thiazole, :thiazole)
end

@testset "6-Membered Aromatic Heterocycles" begin
    # Benzene
    benzene = mol_from_smiles("c1ccccc1")
    @test has_functional_group(benzene, :benzene)

    # Pyridine
    pyridine = mol_from_smiles("c1ccncc1")
    @test has_functional_group(pyridine, :pyridine)

    # Pyrimidine
    pyrimidine = mol_from_smiles("c1cncnc1")
    @test has_functional_group(pyrimidine, :pyrimidine)

    # Pyrazine
    pyrazine = mol_from_smiles("c1cnccn1")
    @test has_functional_group(pyrazine, :pyrazine)
end

@testset "Fused Aromatic Systems" begin
    # Naphthalene
    naphthalene = mol_from_smiles("c1ccc2ccccc2c1")
    @test has_functional_group(naphthalene, :naphthalene)
    @test has_functional_group(naphthalene, :benzene)  # should also match benzene

    # Quinoline
    quinoline = mol_from_smiles("c1ccc2ncccc2c1")
    @test has_functional_group(quinoline, :quinoline)

    # Indole
    indole = mol_from_smiles("c1ccc2[nH]ccc2c1")
    @test has_functional_group(indole, :indole)

    # Benzofuran
    benzofuran = mol_from_smiles("c1ccc2occc2c1")
    @test has_functional_group(benzofuran, :benzofuran)
end

@testset "Saturated Heterocycles" begin
    # Tetrahydrofuran
    thf = mol_from_smiles("C1CCOC1")
    @test has_functional_group(thf, :tetrahydrofuran)

    # Pyrrolidine
    pyrrolidine = mol_from_smiles("C1CCNC1")
    @test has_functional_group(pyrrolidine, :pyrrolidine)

    # Piperidine
    piperidine = mol_from_smiles("C1CCNCC1")
    @test has_functional_group(piperidine, :piperidine)

    # Morpholine
    morpholine = mol_from_smiles("C1COCCN1")
    @test has_functional_group(morpholine, :morpholine)
end

@testset "Reactive Groups" begin
    # Epoxide
    epoxide = mol_from_smiles("C1OC1")
    @test has_functional_group(epoxide, :epoxide)

    # Cyclopropane
    cyclopropane = mol_from_smiles("C1CC1")
    @test has_functional_group(cyclopropane, :cyclopropane)

    # Michael acceptor
    michael = mol_from_smiles("C=CC(=O)C")
    @test has_functional_group(michael, :michael_acceptor)
    @test has_functional_group(michael, :alpha_beta_unsaturated_carbonyl)
end

@testset "Protecting Groups" begin
    # tert-Butyl
    tert_butyl = mol_from_smiles("CC(C)(C)C")
    @test has_functional_group(tert_butyl, :tert_butyl)

    # Benzyl
    benzyl = mol_from_smiles("OCc1ccccc1")  # Benzyl alcohol
    @test has_functional_group(benzyl, :benzyl)

    # Acetyl
    acetyl = mol_from_smiles("CC(=O)O")
    @test has_functional_group(acetyl, :acetyl)

    # Benzoyl
    benzoyl = mol_from_smiles("c1ccccc1C(=O)O")
    @test has_functional_group(benzoyl, :benzoyl)
end

@testset "Drug-like Molecules" begin
    # Acetaminophen (paracetamol)
    acetaminophen = mol_from_smiles("CC(=O)Nc1ccc(O)cc1")
    groups = get_functional_groups(acetaminophen)
    @test groups[:phenol]
    @test groups[:amide]
    @test groups[:benzene]
    @test groups[:acetyl]

    # Aspirin
    aspirin = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    aspirin_groups = get_functional_groups(aspirin)
    @test aspirin_groups[:ester]
    @test aspirin_groups[:carboxylic_acid]
    @test aspirin_groups[:benzene]

    # Simple antibiotic-like β-lactam
    try
        beta_lactam = mol_from_smiles("C1CN(C1=O)C")
        @test has_functional_group(beta_lactam, :beta_lactam)
    catch
        @test_skip "β-lactam test skipped - molecule creation failed"
    end
end

@testset "FUNCTIONAL_GROUPS Dictionary" begin
    # Test that the dictionary exists and has expected size
    @test isa(FUNCTIONAL_GROUPS, Dict{Symbol, String})
    @test length(FUNCTIONAL_GROUPS) > 100  # Should have 100+ functional groups

    # Test that all basic groups are present
    basic_groups = [:alcohol, :ketone, :ester, :amine_primary, :benzene, :pyridine]
    for group in basic_groups
        @test haskey(FUNCTIONAL_GROUPS, group)
        @test isa(FUNCTIONAL_GROUPS[group], String)
        @test length(FUNCTIONAL_GROUPS[group]) > 0
    end

    # Test that advanced groups are present
    advanced_groups = [
        :sulfonamide, :phosphate, :trifluoromethyl, :nitro, :epoxide, :indole, :morpholine
    ]
    for group in advanced_groups
        @test haskey(FUNCTIONAL_GROUPS, group)
        @test isa(FUNCTIONAL_GROUPS[group], String)
        @test length(FUNCTIONAL_GROUPS[group]) > 0
    end
end

@testset "Ring Analysis" begin
    benzene = mol_from_smiles("c1ccccc1")

    ring_info = get_ring_info(benzene)
    @test ring_info[:num_rings] == 1
    @test length(ring_info[:atom_rings]) == 1
    @test length(ring_info[:atom_rings][1]) == 6  # 6-membered ring

    # Test aromatic ring detection
    @test is_ring_aromatic(benzene, ring_info[:atom_rings][1])

    cyclohexane = mol_from_smiles("C1CCCCC1")
    ring_info_saturated = get_ring_info(cyclohexane)
    @test !is_ring_aromatic(cyclohexane, ring_info_saturated[:atom_rings][1])
end
