# HSP nanosuspension dataset: code and data

Companion code for "Hansen Solubility Parameters for Drug-Excipient Compatibility and Nanosuspension Stability: A Source-Traced Dataset and Grouped Validation Across Drugs".

## Contents
| File | Purpose |
|---|---|
| Supplementary_Dataset_S1_v2.csv | Curated dataset, 54 formulations (43 Stable, 11 Unstable), 38 columns with audit fields (including `Stabilizer_System`, the harmonized stabilizer system used for the 21 and 27 system counts) |
| Route_check_16_entries.csv | Literature check of the 16 entries with no recorded preparation route |
| perm_group2.txt, perm_row2.txt | Null AUCs from the drug-preserving (400) and row-level (500) permutation runs; each line is the AUC and the number of Unstable rows |
| hansen_group_contribution_v4.py | Drug-side Hansen parameters (group contribution, RDKit SMARTS) and the Tetronic 1107 fragment |
| hansen_polymer_estimate_v2.py | Stabilizer-side group-contribution values for SDS, Tween 80 and TPGS (simplified representative structures) and the PEO and PPO segments used for the poloxamer estimates that were later replaced by manufacturer values |
| pipe.py | Shared pipeline: 15 features, in-fold median imputation, RF and GB settings, metrics |
| main_final.py | Table 5, feature-set comparison and sensitivity analyses (Supplementary Tables S4 to S6) |
| add_system.py | Writes the `Stabilizer_System` column (components only, ignoring grade, concentration and ratio) |
| fig2.py | Draws the curation flow chart (Figure 2) |
| table5_check.py | Quick regeneration of Table 5 only |
| perm_final.py | Permutation tests (row-level and drug-preserving) |
| generate_structures_v2.py | Structure panel for all 42 APIs (Supplementary Figure S1) |

## Run
    pip install -r requirements.txt
    python perm_final.py obs 0 0          # observed leave-one-drug-out RF ROC-AUC, expected 0.6047
    python table5_check.py                # Table 5 only, about 30 seconds, writes table5_check.csv
    python main_final.py                  # writes table5_final.csv, table_feature_hierarchy.csv, table_sensitivity.csv (about 30 minutes)
    python perm_final.py group2 101 300   # drug-preserving permutations: arguments are mode, random seed, number of runs; exactly 11 Unstable rows on 4 drugs (stored file: seed 101 x 300, then seed 303 x 100)
    python perm_final.py row 11 250       # row-level permutations (stored file: seed 11 x 250, then seed 22 x 250)

The two permutation files hold the null AUCs of several batches run with different seeds (the `perm_final.py` arguments above show the form of each call), so rerunning one command reproduces a batch, not the whole file.

Empirical p = (1 + number of null AUCs at or above the observed AUC) / (permutations + 1).

## Notes
- `table5_check.csv` is the output of `table5_check.py` (Table 5 of the manuscript).
- The dihydropyridine N-H assignment in the group-contribution method is approximate (sensitivity analysis in Supplementary Table S4).
- The drug-side descriptor columns (MW, LogP, HBD, HBA, TPSA) were calculated with RDKit from canonical SMILES. The script that wrote those columns is not included. The SMILES are in the `DRUGS` list of `hansen_group_contribution_v4.py`, so the columns can be recomputed with RDKit (Crippen LogP, Lipinski donors and acceptors, topological polar surface area).
- `ID` runs 1 to 75. IDs 1 to 62 are the original compilation (the gaps are the 21 rows excluded, listed in Supplementary Table S3); IDs 63 to 75 are the 13 failure-specific search rows added later (manuscript Section 2.5). IDs are read as text (`dtype={'ID': str}` in `pipe.py`).
- `Endpoint_Category` and `Formulation_Type` state where the evidence for a label or a preparation route was not recorded; `Route_check_route` and `Route_check_status` give the later literature check of those entries.
- `Stabilizer_HSP_Source_Detail` gives, for every row, the source of each stabilizer component. The values for SDS, Tween 80 and TPGS come from `hansen_polymer_estimate_v2.py`, which uses simplified representative structures (Tween 80: one oleate ester with short polyoxyethylene chains; TPGS: a four-unit PEG chain), so they approximate the commercial materials. The values for PVA, PEG 400, hydroxypropyl-beta-cyclodextrin, ethyl cellulose, HPC and HPMC are literature-consensus values or values whose derivation was not recorded.
- `Source_Location_in_Publication` is filled in only for the four indomethacin rows from the Taylor et al. screening study; every row gives its citation and DOI in `Primary_Citation`.
- `generate_structures_v2.py` downloads SMILES from PubChem and needs an internet connection.

## Licence
Code: MIT (see LICENSE). Dataset (`Supplementary_Dataset_S1_v2.csv` and `Route_check_16_entries.csv`): CC BY 4.0. Cite the archived release with its Zenodo DOI (listed on the repository's release page) and the metadata in `CITATION.cff`.

## Software versions
- Descriptor and Hansen steps (`hansen_group_contribution_v4.py`, `hansen_polymer_estimate_v2.py`; author's MacBook): Python 3.14.7, RDKit 2026.03.5 (descriptor step first run in IDLE on 6 to 7 September 2026; entries added afterwards were processed in the same way).
- All model analyses (`pipe.py`, `table5_check.py`, `main_final.py`, `perm_final.py`): Python 3.13.17, scikit-learn 1.9.1, pandas 3.0.5, NumPy 2.5.3, SciPy 1.18.1; the figure scripts use Matplotlib 3.11.2 and `generate_structures_v2.py` uses Requests 2.34.2. `requirements.txt` pins these versions (RDKit is needed only for the descriptor, Hansen and structure scripts).

## Reproducing Table 5
`pipe.py` holds the settings (Random Forest 300 trees, depth 5, class-balanced; Gradient Boosting 200 estimators, depth 3, learning rate 0.05; 15 features; in-fold median imputation; seed 42) and `table5_check.py` runs them.
