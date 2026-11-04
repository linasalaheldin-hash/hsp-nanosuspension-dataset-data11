"""
hansen_group_contribution_v4.py

v4 CHANGE FROM v3 (the version already used for the published Table 1/Table 2 numbers):
  Adds a group-contribution term for the 1,4-dihydropyridine (DHP) ring N-H,
  which v3 did not match at all. In v3, this atom is silently uncounted for
  every DHP-class drug (Nifedipine, Felodipine, Nitrendipine, Nicardipine,
  Lercanidipine) because:
    - it is a ring N-H, so it fails the acyclic "amine_NH" pattern (!R), and
    - it is not bonded to a carbonyl carbon, so it fails "amide_or_lactam_NH", and
    - it is not aromatic (the DHP ring is non-aromatic due to the sp3 C4
      "dihydro" carbon), so it fails the aromatic "[nX3H1]" pattern too.
  This is exactly the gap disclosed in the manuscript's Limitations section.

WHAT WAS ADDED:
  1. A new SMARTS pattern "dihydropyridine_NH" matching the specific six-atom
     ring pattern N(H)-C(sp2)=C-C(sp3)-C=C(sp2), i.e. the ring nitrogen
     flanked on both sides by a ring C=C (the two "vinylogous"/enamine arms
     characteristic of 1,4-DHPs), with the sp3 "dihydro" carbon (bearing the
     aryl substituent) at the far side of the ring.
  2. A new GROUPS entry "dihydropyridine_NH" using Hansen's plain secondary-
     amine values, (4.5, 1150, 100, 750) -- IMPORTANT: Hansen's table has no
     dedicated entry for this group (it is not a standard amine, amide, or
     aromatic N-H), so this is a deliberate, disclosed approximation: the DHP
     nitrogen is fundamentally a secondary amine (R2N-H) whose two ring
     substituents happen to be enamine-conjugated carbons rather than plain
     alkyl/aryl groups. Reusing the existing, Hansen-tabulated secondary-amine
     value is the most conservative choice available -- it does not invent a
     new number, it reuses one already in Hansen's own table for the closest
     precedented functional class. This should be stated explicitly whenever
     these corrected values are reported.

VERIFICATION BUILT INTO THIS SCRIPT (same discipline as the v2->v3 fix):
  The __main__ block below runs on ALL drugs in the DRUGS list (not just the
  5 DHPs), and explicitly prints which drugs matched the new pattern. Before
  trusting the corrected numbers, confirm:
    (a) exactly Nifedipine, Felodipine, Nitrendipine, Nicardipine, and
        Lercanidipine show a nonzero "dihydropyridine_NH" count (1 each --
        there is only one ring N per DHP core), and
    (b) every other drug's dD/dP/dH is IDENTICAL to its v3 output (i.e. this
        change is additive and does not alter anything for non-DHP molecules).
  If either check fails, do not use these numbers -- flag it and we debug
  the SMARTS together, the same way the v2->v3 carboxylic-acid fix was
  verified before being reported in the manuscript.

HOW TO RUN:
  1. pip install rdkit (if not already installed)
  2. python3 hansen_group_contribution_v4.py
  3. Paste the FULL output back (including the verification section) so the
     new dD/dP/dH/Ra values for the 5 DHP compounds can be checked into
     Table 1, Table 2, and the descriptor CSV, and the Limitations bullet
     can be updated from "disclosed" to "corrected."

Everything else in this file (all other SMARTS, all other GROUPS entries,
all 56 drug SMILES, the v2->v3 carboxylic-acid fix) is UNCHANGED from v3.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import math

# ============================================================
# Hansen's Table 1.1 group contributions (aliphatic values, cal/mol)
# Format: (delta_V, delta_V_dD2, delta_V_dP2, delta_V_dH2)
# ============================================================
GROUPS = {
    "CH3":          (33.5, 1125, 0, 0),
    "CH2":          (16.1, 1180, 0, 0),
    "CH":           (-1.0, 820, 0, 0),
    "C_quat":       (-19.2, 350, 0, 0),
    "CH2_olefin":   (28.5, 850, 25, 180),
    "CH_olefin":    (13.5, 875, 18, 180),
    "C_olefin":     (-5.5, 800, 60, 180),
    "phenyl":       (71.4, 7530, 50, 50),      # monosubstituted benzene ring (C6H5-)
    "ring6":        (16.0, 250, 0, 0),         # saturated 6-membered ring closure
    "ring5":        (16.0, 250, 0, 0),         # saturated 5-membered ring closure
    "F":            (18.0, 0, 1000, 0),
    "Cl":           (24.0, 1400, 1250, 100),
    "Br":           (30.0, 1950, 1250, 500),
    "I":            (31.5, 2350, 1250, 1000),
    "ether_O":      (3.8, 0, 500, 450),
    "ketone":       (10.8, 2350, 1000, 800),
    "aldehyde":     (23.2, 950, 2100, 1000),
    "ester":        (18.0, 640, 1350, 1250),
    "carboxylic_acid": (28.5, 3350, 500, 2750),
    "hydroxyl":     (10.0, 1770, 700, 4650),
    "hydroxyl_twin":(26.0, 0, 1500, 9000),
    "nitrile":      (24.0, 1600, 4000, 500),
    "nitro":        (24.0, 3000, 3600, 400),
    "amine_NH2":    (19.2, 1050, 600, 1350),
    "amine_NH":     (4.5, 1150, 100, 750),
    "amide_NH2":    (6.7, 0, 0, 2700),
    "amide_or_lactam_NH": (10.0, 1150, 400, 2700),
    # NEW (v4): 1,4-dihydropyridine ring N-H. No dedicated Hansen table entry
    # exists for this group; the plain secondary-amine values are reused as
    # the most conservative, precedented approximation (see module docstring).
    "dihydropyridine_NH": (4.5, 1150, 100, 750),
}

# ============================================================
# SMARTS patterns to identify each group
# ============================================================
SMARTS = {
    "CH3":          "[CX4H3]",
    "CH2":          "[CX4H2;!R]",
    "CH":           "[CX4H1;!R]",
    "C_quat":       "[CX4H0;!R]",
    "CH2_olefin":   "[CX3H2]=[CX3]",
    "CH_olefin":    "[CX3H1]=[CX3]",
    "phenyl":       "c1ccccc1",
    "ring6_sat":    "[R6;!a]",
    "ring5_sat":    "[R5;!a]",
    "F":            "[F]",
    "Cl":           "[Cl]",
    "Br":           "[Br]",
    "I":            "[I]",
    "ether_O":      "[OX2;!$(O=*);!$([OX2H]);!$([OX2][C]=O)]([#6])[#6]",
    "ketone":       "[#6][CX3](=O)[#6]",
    "aldehyde":     "[CX3H1](=O)[#6]",
    "ester":        "[#6][CX3](=O)[OX2][#6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "hydroxyl":     "[OX2H1;!$(O=*);!$([OX2H1][CX3](=[OX1]))]",
    "nitrile":      "[CX2]#[NX1]",
    "nitro":        "[NX3](=O)=O",
    "amine_NH2":    "[NX3H2;!$(NC=O)]",
    "amine_NH":     "[NX3H1;!$(NC=O);!R]",
    # NEW (v4): 1,4-dihydropyridine ring N-H. Matches N(H) in a ring, flanked
    # by two ring C=C double bonds with an sp3 carbon (the DHP "C4") between
    # them on the far side of the ring -- i.e. the six-membered
    # N-C(sp2)=C-C(sp3)-C=C(sp2)- ring characteristic of Hantzsch DHPs.
    "dihydropyridine_NH": "[NX3H1;R1]1[CX3]=[CX3][CX4][CX3]=[CX3]1",
}

AMIDE_OR_AROMATIC_NH_PATTERNS = [
    "[#7X3H1]C=O",
    "[nX3H1]",
]

def count_groups(mol):
    """Count matches for each SMARTS pattern, returning a dict of group->count."""
    counts = {}
    for name, smarts in SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            print(f"  WARNING: could not parse SMARTS for {name}: {smarts}")
            continue
        matches = mol.GetSubstructMatches(patt)
        counts[name] = len(matches)

    nh_atom_indices = set()
    for smarts in AMIDE_OR_AROMATIC_NH_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            print(f"  WARNING: could not parse SMARTS: {smarts}")
            continue
        for match in mol.GetSubstructMatches(patt):
            nh_atom_indices.add(match[0])
    counts["amide_or_lactam_NH"] = len(nh_atom_indices)

    return counts

def calculate_hsp(smiles, drug_name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"!!! ERROR: could not parse SMILES for {drug_name}: {smiles}")
        return None

    mol = Chem.AddHs(mol)
    counts = count_groups(mol)

    total_V = 0.0
    total_dD2 = 0.0
    total_dP2 = 0.0
    total_dH2 = 0.0

    mapping = {
        "CH3": "CH3", "CH2": "CH2", "CH": "CH", "C_quat": "C_quat",
        "CH2_olefin": "CH2_olefin", "CH_olefin": "CH_olefin",
        "phenyl": "phenyl", "ring6_sat": "ring6", "ring5_sat": "ring5",
        "F": "F", "Cl": "Cl", "Br": "Br", "I": "I",
        "ether_O": "ether_O", "ketone": "ketone", "aldehyde": "aldehyde",
        "ester": "ester", "carboxylic_acid": "carboxylic_acid",
        "hydroxyl": "hydroxyl", "nitrile": "nitrile", "nitro": "nitro",
        "amine_NH2": "amine_NH2", "amine_NH": "amine_NH",
        "amide_or_lactam_NH": "amide_or_lactam_NH",
        "dihydropyridine_NH": "dihydropyridine_NH",  # NEW (v4)
    }

    detail_lines = []
    for smarts_key, n in counts.items():
        if n == 0:
            continue
        group_key = mapping.get(smarts_key)
        if group_key not in GROUPS:
            continue
        dV, dD2, dP2, dH2 = GROUPS[group_key]
        total_V += n * dV
        total_dD2 += n * dD2
        total_dP2 += n * dP2
        total_dH2 += n * dH2
        detail_lines.append(f"    {smarts_key}: {n} x (V={dV}, dD2={dD2}, dP2={dP2}, dH2={dH2})")

    if total_V <= 0:
        print(f"!!! ERROR: total molar volume <=0 for {drug_name} -- group counting likely failed")
        return None

    CAL_TO_MPA = 2.0455
    dD = math.sqrt(total_dD2 / total_V) * CAL_TO_MPA
    dP = math.sqrt(total_dP2 / total_V) * CAL_TO_MPA
    dH = math.sqrt(total_dH2 / total_V) * CAL_TO_MPA

    dhp_count = counts.get("dihydropyridine_NH", 0)

    print(f"\n--- {drug_name} ---")
    print(f"  SMILES: {smiles}")
    print(f"  Matched groups:")
    for line in detail_lines:
        print(line)
    print(f"  Total molar volume (V): {total_V:.1f} cm3/mol")
    print(f"  Hansen dD: {dD:.2f} MPa^0.5")
    print(f"  Hansen dP: {dP:.2f} MPa^0.5")
    print(f"  Hansen dH: {dH:.2f} MPa^0.5")
    if dhp_count > 0:
        print(f"  >>> dihydropyridine_NH matched {dhp_count}x on this molecule <<<")

    return (dD, dP, dH, dhp_count)


# ============================================================
# All drugs from v3, UNCHANGED (56 SMILES)
# ============================================================
DRUGS = [
    ("Naproxen", "COc1ccc2cc(ccc2c1)[C@@H](C)C(=O)O"),
    ("Indomethacin", "CC1=C(C2=CC(=CC=C2N1C(=O)C3=CC=C(C=C3)Cl)OC)CC(=O)O"),
    ("Itraconazole", "CCC(C)N1N=CN(N1)C2=CC=C(C=C2)N3CCN(CC3)C4=CC=C(C=C4)OCC5COC(O5)(CN6C=NC=N6)C7=CC=C(C=C7)Cl"),
    ("Cinnarizine", "C1CN(CCN1CC=CC2=CC=CC=C2)C(C3=CC=CC=C3)C4=CC=CC=C4"),
    ("Fenofibrate", "CC(C)OC(=O)C(C)(C)OC1=CC=C(C=C1)C(=O)C2=CC=C(C=C2)Cl"),
    ("Fisetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(O2)C=C(C=C3)O)O)O)O"),
    ("Celecoxib", "CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F"),
    ("Danazol", "CC12CCC3C(C1CCC2(C#C)O)CCC4=CC(=O)C=CC34C"),
    ("Ketoconazole", "CC(=O)N1CCN(CC1)C1=CC=C(C=C1)OCC1COC(CN2C=CN=C2)(O1)C1=CC=C(Cl)C=C1Cl"),
    ("Nifedipine", "CC1=C(C(C(=C(N1)C)C(=O)OC)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC"),
    ("Piroxicam", "CN1C(=C(O)c2ccccc2S1(=O)=O)C(=O)Nc1ccccn1"),
    ("Diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("Flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2F)cc1"),
    ("Baicalein", "O=C\\1c3c(O)c(O)c(O)cc3O/C(=C/1)c2ccccc2"),
    ("Spironolactone", "CC(=O)SC1CCC2(C)C3CCC4=CC(=O)CCC4(C)C3CCC12C"),
    ("Probucol", "CC(C)(C)c1cc(SC(C)(C)c2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)c(O)c(C(C)(C)C)c1"),
    ("Felodipine", "CCOC(=O)C1=C(C)NC(=C(C1c1ccccc1Cl)C(=O)OC)C"),
    ("Nimodipine", "CC(C)OC(=O)C1=C(C)NC(=C(C1c1cccc([N+](=O)[O-])c1)C(=O)OCCOC)C"),
    ("Aprepitant", "O=C1N(C(c2ccccc2F)C(F)(F)F)N=C(N1)N1CCN(CC1)Cc1nnnn1C"),
    ("Meloxicam", "CC1=CN=C(S1)NC(=O)C1=C(O)c2ccccc2S(=O)(=O)N1C"),
    ("Ezetimibe", "O=C1N(c2ccc(F)cc2)C(c2ccc(O)cc2)C1CC(O)c1ccc(F)cc1"),
    ("Glibenclamide", "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"),
    ("Cilostazol", "O=C1CCCCN1CCCCOc1ccc2[nH]c(nc2c1)C1CCCCC1"),
    ("Aripiprazole", "Clc1cccc(c1N1CCN(CCCCOc2ccc3c(c2)CCC(=O)N3)CC1)Cl"),
    ("Ritonavir", "CC(C)C1=NC(=CS1)CN(C)C(=O)NC(C(C)C)C(=O)NC(CC2=CC=CC=C2)CC(C(CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O"),
    ("Simvastatin", "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C"),
    ("Lovastatin", "CCC(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C"),
    ("Atovaquone", "OC=2C(=O)c1ccccc1C(=O)C=2C3CCC(CC3)c4ccc(Cl)cc4"),
    ("Paclitaxel", "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"),
    ("Loratadine", "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1"),
    ("Repaglinide", "CCOC1=C(C=CC(=C1)CC(=O)NC(CC(C)C)C2=CC=CC=C2N3CCCCC3)C(=O)O"),
    ("Nitrendipine", "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1cccc([N+](=O)[O-])c1"),
    ("Rifampicin", "CN1CCN(CC1)N=Cc2c(O)c3c5C(=O)C4(C)OC=CC(OC)C(C)C(OC(C)=O)C(C)C(O)C(C)C(O)C(C)C=CC=C(C)C(=O)Nc2c(O)c3c(O)c(C)c5O4"),
    ("Clotrimazole", "C1=CC=C(C=C1)C(C2=CC=CC=C2)(C3=CC=CC=C3Cl)N4C=CN=C4"),
    ("Miconazole", "C1=CC(=C(C=C1Cl)Cl)COC(CN2C=CN=C2)C3=C(C=C(C=C3)Cl)Cl"),
    ("Albendazole", "CCCSC1=CC2=C(C=C1)N=C(N2)NC(=O)OC"),
    ("Mebendazole", "COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CC=C3"),
    ("Haloperidol", "C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F"),
    ("Domperidone", "C1CN(CCC1N2C3=C(C=C(C=C3)Cl)NC2=O)CCCN4C5=CC=CC=C5NC4=O"),
    ("Raloxifene", "C1CCN(CC1)CCOC2=CC=C(C=C2)C(=O)C3=C(SC4=C3C=CC(=C4)O)C5=CC=C(C=C5)O"),
    ("Tamoxifen", "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3"),
    ("Nicardipine", "CC1=C(C(C(=C(N1)C)C(=O)OCCN(C)CC2=CC=CC=C2)C3=CC(=CC=C3)[N+](=O)[O-])C(=O)OC"),
    ("Diosmin", "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3=CC(=C4C(=C3)OC(=CC4=O)C5=CC(=C(C=C5)OC)O)O)O)O)O)O)O)O"),
    ("Ciprofloxacin", "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"),
    ("Rutin", "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4"),
    ("Candesartan_Cilexetil", "CCOC1=NC2=CC=CC(=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NNN=N5)C(=O)OC(C)OC(=O)OC6CCCCC6"),
    ("Hydrochlorothiazide", "C1NC2=CC(=C(C=C2S(=O)(=O)N1)S(=O)(=O)N)Cl"),
    ("Tadalafil", "CN1CC(=O)N2C(C1=O)CC3=C(C2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36"),
    ("Quercetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"),
    ("Rosuvastatin", "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C"),
    ("Curcumin", "COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O"),
    ("Silibinin", "COC1=C(C=CC(=C1)C2C(OC3=C(O2)C=C(C=C3)C4C(C(=O)C5=C(C=C(C=C5O4)O)O)O)CO)O"),
    ("Lercanidipine", "CC1=C(C(C(=C(N1)C)C(=O)OC(C)(C)CN(C)CCC(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC(=CC=C4)[N+](=O)[O-])C(=O)OC"),
    ("Luteolin", "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O"),
    ("Furosemide", "C1=COC(=C1)CNC2=CC(=C(C=C2C(=O)O)S(=O)(=O)N)Cl"),
]

EXPECTED_DHP_DRUGS = {"Nifedipine", "Felodipine", "Nitrendipine", "Nicardipine", "Lercanidipine"}

if __name__ == "__main__":
    results = {}
    dhp_matches = {}
    for name, smi in DRUGS:
        r = calculate_hsp(smi, name)
        if r:
            dD, dP, dH, dhp_count = r
            results[name] = (dD, dP, dH)
            if dhp_count > 0:
                dhp_matches[name] = dhp_count

    print("\n\n========== SUMMARY TABLE (all 56 drugs, v4) ==========")
    print(f"{'Drug':<25} {'dD':>8} {'dP':>8} {'dH':>8}")
    for name, (dD, dP, dH) in results.items():
        print(f"{name:<25} {dD:>8.2f} {dP:>8.2f} {dH:>8.2f}")

    print("\n\n========== VERIFICATION: dihydropyridine_NH matches ==========")
    print(f"Expected to match (exactly these 5, count=1 each): {sorted(EXPECTED_DHP_DRUGS)}")
    print(f"Actually matched: {dhp_matches}")
    matched_names = set(dhp_matches.keys())
    if matched_names == EXPECTED_DHP_DRUGS and all(v == 1 for v in dhp_matches.values()):
        print("PASS: exactly the 5 expected DHP drugs matched, once each.")
    else:
        missing = EXPECTED_DHP_DRUGS - matched_names
        unexpected = matched_names - EXPECTED_DHP_DRUGS
        print("*** MISMATCH -- DO NOT USE THESE NUMBERS YET ***")
        if missing:
            print(f"  Missing (expected but did not match): {sorted(missing)}")
        if unexpected:
            print(f"  Unexpected (matched but should not have): {sorted(unexpected)}")
        wrong_count = {k: v for k, v in dhp_matches.items() if v != 1}
        if wrong_count:
            print(f"  Wrong match count (expected 1): {wrong_count}")

    print("\n\n========== VERIFICATION: non-DHP drugs unchanged from v3 ==========")
    print("Manually spot-check a few of these dD/dP/dH values above against the")
    print("v3 output you already have on file (e.g. Naproxen, Indomethacin,")
    print("Diclofenac, Ibuprofen) -- they should be IDENTICAL, since v4 only")
    print("adds a new pattern and changes nothing else.")
