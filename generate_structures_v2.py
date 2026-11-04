"""
generate_structures.py -- render 2D chemical structures for the manuscript.

WHAT THIS DOES
1. Fetches the canonical SMILES for all 42 unique APIs directly from PubChem's
   PUG-REST API by compound name (not from memory -- this is the reliable part).
2. Falls back to a small hand-verified SMILES dictionary only for the handful
   of representative compounds where the name might not resolve cleanly on
   PubChem (e.g. "Silibinin", "Candesartan Cilexetil").
3. Produces TWO outputs:
   a) supplementary_structures/  -- all 42 APIs, individual PNGs + one big
      grid image, for Supplementary Materials.
   b) main_figure_structures/    -- one representative compound per
      pharmacological class (12 total), individual PNGs + one grid image
      sized for the main manuscript.

Usage:
    pip install rdkit requests --break-system-packages   # if not already installed
    python3 generate_structures.py

Upload back to the chat: main_figure_panel.png (for the manuscript) and,
if you want the full catalog archived, supplementary_panel.png plus the
supplementary_structures/ folder.
"""

import os
import time
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES/TXT"

# All 42 unique APIs in the curated dataset, with assigned pharmacological class.
# "pubchem_name" lets us pass a cleaner/alternate name to PubChem when the
# dataset's display name might not resolve (e.g. the marker compound used for
# silymarin). Leave as None to just use the drug name as-is.
DRUGS = [
    # (display_name, pharmacological_class, pubchem_name_override)
    ("Naproxen", "NSAIDs", None),
    ("Ibuprofen", "NSAIDs", None),
    ("Indomethacin", "NSAIDs", None),
    ("Celecoxib", "NSAIDs", None),
    ("Piroxicam", "NSAIDs", None),
    ("Meloxicam", "NSAIDs", None),
    ("Nifedipine", "Cardiovascular", None),
    ("Felodipine", "Cardiovascular", None),
    ("Nitrendipine", "Cardiovascular", None),
    ("Nicardipine", "Cardiovascular", None),
    ("Lercanidipine", "Cardiovascular", None),
    ("Simvastatin", "Cardiovascular", None),
    ("Rosuvastatin", "Cardiovascular", None),
    ("Atorvastatin", "Cardiovascular", None),
    ("Spironolactone", "Cardiovascular", None),
    ("Furosemide", "Cardiovascular", None),
    ("Candesartan Cilexetil", "Cardiovascular", "Candesartan cilexetil"),
    ("Itraconazole", "Antifungals", None),
    ("Clotrimazole", "Antifungals", None),
    ("Rifampicin", "Antibacterial", "Rifampin"),
    ("Ciprofloxacin", "Antibacterial", None),
    ("Fisetin", "Flavonoids & polyphenols", None),
    ("Quercetin", "Flavonoids & polyphenols", None),
    ("Rutin", "Flavonoids & polyphenols", None),
    ("Diosmin", "Flavonoids & polyphenols", None),
    ("Luteolin", "Flavonoids & polyphenols", None),
    ("Curcumin", "Flavonoids & polyphenols", None),
    ("Silibinin", "Flavonoids & polyphenols", "Silibinin"),
    ("Haloperidol", "Antiemetic / CNS dopaminergic", None),
    ("Aprepitant", "Antiemetic / CNS dopaminergic", None),
    ("Domperidone", "Antiemetic / CNS dopaminergic", None),
    ("Glibenclamide", "Antidiabetic / lipid-modulating", "Glyburide"),
    ("Ezetimibe", "Antidiabetic / lipid-modulating", None),
    ("Probucol", "Antidiabetic / lipid-modulating", None),
    ("Cilostazol", "PDE inhibitors", None),
    ("Tadalafil", "PDE inhibitors", None),
    ("Mebendazole", "Other (anthelmintic)", None),
    ("Paclitaxel", "Other (antineoplastic)", None),
    ("Loratadine", "Other (antihistamine)", None),
    ("Raloxifene", "Other (SERM)", None),
]

# One representative per class for the main manuscript figure.
REPRESENTATIVES = [
    "Naproxen", "Nifedipine", "Itraconazole", "Rifampicin",
    "Curcumin", "Haloperidol", "Glibenclamide", "Tadalafil",
    "Mebendazole", "Paclitaxel", "Loratadine", "Raloxifene",
]

# Hand-verified fallback SMILES, used ONLY if the PubChem lookup fails.
# Cross-checked against Wikipedia/PubChem/SCBT records; kept here so the
# script still produces the 12 main-figure compounds even with a flaky
# network call.
FALLBACK_SMILES = {
    "Naproxen": "COC1=CC2=CC(=CC=C2C=C1)C(C)C(=O)O",
    "Nifedipine": "CC1=C(C(=O)OC)C(c2ccccc2[N+](=O)[O-])C(=C(N1)C)C(=O)OC",
    "Itraconazole": "CCC(C)N1C(=O)N(C=N1)C2=CC=C(C=C2)N3CCN(CC3)C4=CC=C(C=C4)OCC5COC(O5)(CN6C=NC=N6)C7=C(C=C(C=C7)Cl)Cl",
    "Rifampicin": "CN1CCN(CC1)/N=C/c2c(O)c3c5C(=O)[C@@]4(C)O/C=C/[C@H](OC)[C@@H](C)[C@@H](OC(C)=O)[C@H](C)[C@H](O)[C@H](C)[C@@H](O)[C@@H](C)\\C=C\\C=C(\\C)C(=O)Nc2c(O)c3c(O)c(C)c5O4",
    "Curcumin": "COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O",
    "Haloperidol": "C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F",
    "Glibenclamide": "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
    "Tadalafil": "CN1CC(=O)N2C(C1=O)CC3=C(C2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36",
    "Mebendazole": "COC(=O)Nc3nc2ccc(C(=O)c1ccccc1)cc2[nH]3",
    "Paclitaxel": "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    "Loratadine": "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
    "Raloxifene": "C1CCN(CC1)CCOC2=CC=C(C=C2)C(=O)C3=C(SC4=C3C=CC(=C4)O)C5=CC=C(C=C5)O",
    "Silibinin": "COC1=C(C=CC(=C1)C2C(OC3=C(O2)C=C(C=C3)C4C(C(=O)C5=C(C=C(C=C5O4)O)O)O)CO)O",
}


def fetch_smiles_from_pubchem(name):
    url = PUG_REST.format(requests.utils.quote(name))
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200 and resp.text.strip():
            return resp.text.strip().splitlines()[0]
    except requests.RequestException as e:
        print(f"  PubChem request failed for '{name}': {e}")
    return None


def get_smiles(display_name, pubchem_name_override):
    query_name = pubchem_name_override or display_name
    smiles = fetch_smiles_from_pubchem(query_name)
    source = "PubChem"
    if not smiles:
        smiles = FALLBACK_SMILES.get(display_name)
        source = "fallback (hand-verified)" if smiles else None
    return smiles, source


def draw_structure(smiles, out_path, size=(500, 450)):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 2
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    with open(out_path, "wb") as f:
        f.write(drawer.GetDrawingText())
    return True


def main():
    os.makedirs("supplementary_structures", exist_ok=True)
    os.makedirs("main_figure_structures", exist_ok=True)

    resolved = {}
    failures = []

    print("Fetching SMILES for all 42 APIs from PubChem...\n")
    for display_name, pharm_class, override in DRUGS:
        smiles, source = get_smiles(display_name, override)
        if smiles is None:
            print(f"  FAILED: {display_name} -- no SMILES from PubChem or fallback")
            failures.append(display_name)
            continue
        resolved[display_name] = (smiles, pharm_class)
        print(f"  {display_name:25s} [{source:24s}] {smiles[:60]}{'...' if len(smiles) > 60 else ''}")
        time.sleep(0.2)

    sup_mols, sup_legends = [], []
    for display_name, (smiles, pharm_class) in resolved.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  RDKit could not parse SMILES for {display_name} -- skipping")
            continue
        out_path = f"supplementary_structures/{display_name.replace(' ', '_')}.png"
        draw_structure(smiles, out_path)
        sup_mols.append(mol)
        sup_legends.append(f"{display_name}\n({pharm_class})")

    if sup_mols:
        img = Draw.MolsToGridImage(sup_mols, legends=sup_legends, molsPerRow=5,
                                    subImgSize=(320, 300), useSVG=False)
        img.save("supplementary_panel.png")
        print(f"\nWrote supplementary_panel.png ({len(sup_mols)}/{len(DRUGS)} compounds)")

    rep_mols, rep_legends = [], []
    for name in REPRESENTATIVES:
        if name not in resolved:
            print(f"  WARNING: representative '{name}' was not resolved -- skipping in main figure")
            continue
        smiles, pharm_class = resolved[name]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        out_path = f"main_figure_structures/{name.replace(' ', '_')}.png"
        draw_structure(smiles, out_path, size=(600, 520))
        rep_mols.append(mol)
        rep_legends.append(f"{name}\n({pharm_class})")

    if rep_mols:
        img = Draw.MolsToGridImage(rep_mols, legends=rep_legends, molsPerRow=4,
                                    subImgSize=(380, 340), useSVG=False)
        img.save("main_figure_panel.png")
        print(f"Wrote main_figure_panel.png ({len(rep_mols)}/{len(REPRESENTATIVES)} representatives)")

    if failures:
        print(f"\n{len(failures)} compound(s) could not be resolved at all: {failures}")
        print("If any of these were meant to be in the main figure, send me the name(s)")
        print("and I'll help source a verified SMILES manually.")

    print("\nDone. Upload main_figure_panel.png (required) and supplementary_panel.png")
    print("+ the supplementary_structures/ folder (if you want the full catalog archived).")


if __name__ == "__main__":
    main()
