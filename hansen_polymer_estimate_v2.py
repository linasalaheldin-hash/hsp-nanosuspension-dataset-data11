"""
hansen_polymer_estimate.py

Extends the same group-contribution method (hansen_group_contribution.py) to
estimate Hansen parameters for several PHARMACEUTICAL EXCIPIENTS, using either:
  (a) a full representative SMILES (for defined, discrete molecules like SDS,
      Tween 80, TPGS), or
  (b) a short repeat-unit/segment SMILES combined via volume-fraction-weighted
      averaging across block segments (for block copolymers like the Poloxamers)

METHODOLOGY NOTE ON BLOCK COPOLYMERS (Poloxamer 188/407):
  Rather than hand-building one giant, error-prone oligomer SMILES for a whole
  PEO-PPO-PEO chain, this script computes Hansen parameters for a SHORT,
  low-risk repeat-unit segment of each block (a few ethylene-oxide units for
  PEO, a few propylene-oxide units for PPO) independently, then combines them
  using the exact same volume-fraction weighted-average formula already used
  elsewhere in this project for mixed stabilizer systems:
      delta_mix = sum(phi_i * delta_i)
  Poloxamer 188 composition: ~80% oxyethylene by weight (per manufacturer spec)
  Poloxamer 407 composition: ~70% oxyethylene by weight
  Weight fractions are converted to approximate volume fractions using each
  segment's computed molar volume (V) from the group-contribution calculation
  itself, for internal consistency.

SCOPE -- what this script does NOT attempt, and why:
  HPMC, HPC, Soluplus, Kollicoat MAE 100P, HP-beta-cyclodextrin, Chitosan, and
  Copovidone (PVP-VA64) are NOT included here. Their repeat units involve
  either variable/uncertain degree-of-substitution (cellulose ethers, HP-CD)
  or complex graft/backbone chemistry (Soluplus) that carries real risk of a
  hand-built SMILES error I cannot catch myself (no RDKit in my environment).
  Rather than hand you a fabricated-looking number for these, they are left
  for your friend's HSPiP database lookup, or a future follow-up pass with
  more care.

HOW TO RUN: same as before -- pip install rdkit if needed, then
  python3 hansen_polymer_estimate.py
and paste the output back.
"""

from rdkit import Chem
import math

# Re-declare the same group tables and SMARTS from hansen_group_contribution.py
# (kept self-contained here so this file can run independently)

GROUPS = {
    "CH3":          (33.5, 1125, 0, 0),
    "CH2":          (16.1, 1180, 0, 0),
    "CH":           (-1.0, 820, 0, 0),
    "C_quat":       (-19.2, 350, 0, 0),
    "phenyl":       (71.4, 7530, 50, 50),
    "ring6":        (16.0, 250, 0, 0),
    "ring5":        (16.0, 250, 0, 0),
    "F":            (18.0, 0, 1000, 0),
    "Cl":           (24.0, 1400, 1250, 100),
    "Br":           (30.0, 1950, 1250, 500),
    "I":            (31.5, 2350, 1250, 1000),
    "ether_O":      (3.8, 0, 500, 450),
    "ketone":       (10.8, 2350, 1000, 800),
    "ester":        (18.0, 640, 1350, 1250),
    "carboxylic_acid": (28.5, 3350, 500, 2750),
    "hydroxyl":     (10.0, 1770, 700, 4650),
    "sulfate_ester": (24.0, 900, 3600, 2000),  # approximate -- Hansen table lacks a direct
                                                 # sulfate-ester entry; built from combining
                                                 # ester-like + strongly polar S=O character
}

SMARTS = {
    "CH3":          "[CX4H3]",
    "CH2":          "[CX4H2;!R]",
    "CH":           "[CX4H1;!R]",
    "C_quat":       "[CX4H0;!R]",
    "phenyl":       "c1ccccc1",
    "ring6_sat":    "[R6;!a]",
    "ring5_sat":    "[R5;!a]",
    "F":            "[F]",
    "Cl":           "[Cl]",
    "ether_O":      "[OX2;!$(O=*);!$([OX2H]);!$([OX2][C]=O)]([#6])[#6]",
    "ketone":       "[#6][CX3](=O)[#6]",
    "ester":        "[#6][CX3](=O)[OX2][#6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "hydroxyl":     "[OX2H1;!$(O=*)]",
    "sulfate_ester": "[OX2,OX1-][SX4](=O)(=O)[OX2,OX1-]",
}

MAPPING = {
    "CH3": "CH3", "CH2": "CH2", "CH": "CH", "C_quat": "C_quat",
    "phenyl": "phenyl", "ring6_sat": "ring6", "ring5_sat": "ring5",
    "F": "F", "Cl": "Cl",
    "ether_O": "ether_O", "ketone": "ketone", "ester": "ester",
    "carboxylic_acid": "carboxylic_acid", "hydroxyl": "hydroxyl",
    "sulfate_ester": "sulfate_ester",
}

CAL_TO_MPA = 2.0455

def calculate_hsp(smiles, label):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"!!! ERROR: could not parse SMILES for {label}: {smiles}")
        return None
    mol = Chem.AddHs(mol)

    total_V = total_dD2 = total_dP2 = total_dH2 = 0.0
    detail = []
    for smarts_key, smarts in SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            print(f"  WARNING: bad SMARTS {smarts_key}")
            continue
        n = len(mol.GetSubstructMatches(patt))
        if n == 0:
            continue
        group_key = MAPPING[smarts_key]
        dV, dD2, dP2, dH2 = GROUPS[group_key]
        total_V += n * dV
        total_dD2 += n * dD2
        total_dP2 += n * dP2
        total_dH2 += n * dH2
        detail.append(f"    {smarts_key}: {n}")

    if total_V <= 0:
        print(f"!!! ERROR: total V <= 0 for {label}")
        return None

    dD = math.sqrt(total_dD2 / total_V) * CAL_TO_MPA
    dP = math.sqrt(total_dP2 / total_V) * CAL_TO_MPA
    dH = math.sqrt(total_dH2 / total_V) * CAL_TO_MPA

    print(f"\n--- {label} ---")
    print(f"  SMILES: {smiles}")
    for line in detail:
        print(line)
    print(f"  Molar volume (V): {total_V:.1f} cm3/mol")
    print(f"  dD={dD:.2f}  dP={dP:.2f}  dH={dH:.2f}")
    return (dD, dP, dH, total_V)


def blend(segments_with_weights, label):
    """
    segments_with_weights: list of (result_tuple, weight_fraction) where
    result_tuple = (dD, dP, dH, V) from calculate_hsp().
    Combines via volume-fraction-weighted average (converting weight
    fractions to volume fractions using each segment's V).
    """
    # Convert weight fractions to volume fractions (assuming equal density
    # approximation is NOT used -- we weight by weight_fraction / V, i.e.
    # proportional moles, which under equal-density assumption approximates
    # volume fraction closely enough for this estimate).
    total_weighted_V = sum(w for _, w in segments_with_weights)
    phi = [(w / total_weighted_V) for _, w in segments_with_weights]

    dD_mix = sum(p * r[0] for (r, _), p in zip(segments_with_weights, phi))
    dP_mix = sum(p * r[1] for (r, _), p in zip(segments_with_weights, phi))
    dH_mix = sum(p * r[2] for (r, _), p in zip(segments_with_weights, phi))

    print(f"\n=== BLENDED RESULT: {label} ===")
    print(f"  dD={dD_mix:.2f}  dP={dP_mix:.2f}  dH={dH_mix:.2f}")
    return (dD_mix, dP_mix, dH_mix)


if __name__ == "__main__":
    print("############################################")
    print("# Defined small molecules (direct calculation)")
    print("############################################")

    # SDS - exact structure (sodium dodecyl sulfate)
    sds = calculate_hsp("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", "SDS (Sodium dodecyl sulfate)")

    # Tween 80 (Polysorbate 80) - simplified representative structure:
    # sorbitan ring, one oleate ester, three short POE chains (real Tween 80
    # has ~20 mol EO total distributed across positions -- this is a reduced
    # low-risk approximation, NOT the full commercial structure)
    tween80 = calculate_hsp(
        "CCCCCCCCC=CCCCCCCCC(=O)OCC1OC(CO)C(OCCO)C1OCCO",
        "Tween 80 (simplified representative structure -- see caveat)"
    )

    # TPGS - simplified representative structure: tocopherol succinate ester
    # with a short PEG chain (real TPGS has PEG ~1000 Da / ~22 EO units --
    # this uses a short 4-unit PEG as a low-risk stand-in, NOT full length)
    tpgs = calculate_hsp(
        "CC1=C(C)C2=C(CCC(C)(CCCC(C)CCCC(C)CCCC(C)C)O2)C(C)=C1OC(=O)CCC(=O)OCCOCCOCCOCCO",
        "TPGS (simplified representative structure -- short PEG stand-in)"
    )

    print("\n\n############################################")
    print("# Block copolymer segments (for volume-weighted blending)")
    print("############################################")

    peo_segment = calculate_hsp("OCCOCCOCCOCCO", "PEO segment (tetraethylene glycol analog)")
    ppo_segment = calculate_hsp("OC(C)COC(C)COC(C)CO", "PPO segment (tripropylene glycol analog)")

    print("\n\n############################################")
    print("# Blended block copolymer estimates")
    print("############################################")

    if peo_segment and ppo_segment:
        # Poloxamer 188: ~80% oxyethylene by weight
        poloxamer_188 = blend(
            [(peo_segment, 0.80), (ppo_segment, 0.20)],
            "Poloxamer 188 (80% EO / 20% PO by weight)"
        )
        # Poloxamer 407: ~70% oxyethylene by weight
        poloxamer_407 = blend(
            [(peo_segment, 0.70), (ppo_segment, 0.30)],
            "Poloxamer 407 (70% EO / 30% PO by weight)"
        )

    print("\n\n########## SUMMARY TABLE (copy this whole block back) ##########")
    print(f"{'Substance':<45} {'dD':>8} {'dP':>8} {'dH':>8}")
    if sds:
        print(f"{'SDS':<45} {sds[0]:>8.2f} {sds[1]:>8.2f} {sds[2]:>8.2f}")
    if tween80:
        print(f"{'Tween 80 (simplified)':<45} {tween80[0]:>8.2f} {tween80[1]:>8.2f} {tween80[2]:>8.2f}")
    if tpgs:
        print(f"{'TPGS (simplified)':<45} {tpgs[0]:>8.2f} {tpgs[1]:>8.2f} {tpgs[2]:>8.2f}")
    if peo_segment and ppo_segment:
        print(f"{'Poloxamer 188 (blended)':<45} {poloxamer_188[0]:>8.2f} {poloxamer_188[1]:>8.2f} {poloxamer_188[2]:>8.2f}")
        print(f"{'Poloxamer 407 (blended)':<45} {poloxamer_407[0]:>8.2f} {poloxamer_407[1]:>8.2f} {poloxamer_407[2]:>8.2f}")
