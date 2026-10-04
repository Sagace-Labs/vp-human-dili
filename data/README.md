# Human DILI data

`human_dili.parquet` is an explicitly fetched training table and is ignored by
Git. `example/human_dili_example.parquet` is a synthetic, stratified fixture;
its labels have no biological interpretation. `validation/` holds the frozen
FDA cohort, source assertions, structure-resolution exclusions and split keys.

## Training source

[Seal et al. source](https://github.com/srijitseal/DILI/tree/98b7a1b7fba4547cefe936b9890b44baab10e490)
provides 1,111 DILIst-derived entries with `smiles_r` and binary `TOXICITY`.
Normalize structures, retain the largest fragment, neutralize, reject invalid
structures, group by InChIKey connectivity, exclude conflicting labels, retain
the first consistent identity and sort by full InChIKey. This yields 1,096
compounds: 708 positive and 388 negative. Source rank, published predictions,
dose and exposure columns do not enter the feature matrix.

    python -m vp_human_dili.data fetch --verify

Fetching verifies the pinned raw snapshot, writes `sources/gold.csv` and a
retrieval ledger, and checks the parsed `inchikey`/`smiles`/`label` hash.
The full table and raw snapshot are not redistributed; see `../LICENSE-DATA`.

## FDA validation cohort

The [FDA DILIrank 2.0 table](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-20-dataset)
contains 1,336 entries: Most 217, Less 351, No 414 and Ambiguous 354.
`validation/sources/fda_web_snapshot.txt` retains the official table used for
this frozen cohort. Of 169 New entries with definite labels, PubChem resolves 155;
14 lack usable structures. Most/Less are positive and No negative. Ambiguous
labels are excluded. Connectivity conflicts and duplicate identities are
excluded; all training identities and scaffolds, including the empty scaffold,
are removed. The resulting `fda_new.parquet` has 119 drugs, 35 positive and
84 negative, across 115 scaffolds.

`resolutions.csv`, `exclusions.csv` and `sources/pubchem_assertions.json` retain
per-drug name-resolution evidence. `manifest.json` pins files and the standardized
cohort hash. `split_keys.csv` freezes every primary training/calibration key and
external key for seed 101; FDA labels do not select hyperparameters or thresholds.
Run `python -m vp_human_dili.validation` to rebuild this cohort offline from
the retained FDA and PubChem assertions after fetching the primary table.

Automatic parent normalization does not verify every formulation or biologic.
Large peptide structures remain; 100 test drugs have nearest-training binary
Morgan Tanimoto below 0.4. This related-source transfer experiment does not
provide independent clinical adjudication. The FDA test cohort is distinct
from the primary benchmark; scores are never compared across them.
