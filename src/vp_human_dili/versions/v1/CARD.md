# human-dili v1

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-10-04 · signature 1

**Why this version.** Initial counted-Morgan forest for human drug-level DILI labels with held-out calibration and unseen FDA drug validation

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `human_dili` | float32 | 0.0–1.0 | Structure-based probability of a positive drug-level DILIst-derived label, calibrated on held-out training-source drugs. It is not patient-specific DILI incidence or a clinical safety decision. |

Missing values: NaN when a nonempty molecular structure cannot be parsed and standardized

## Performance — `scaffold-balanced-5seed@1`

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-10-04 on n_train=822, n_val=110, n_test=164.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.705 | 0.030 | 0.694, 0.706, 0.703, 0.759, 0.665 |
| auprc | 0.807 | 0.015 | 0.822, 0.828, 0.791, 0.801, 0.794 |
| mcc | 0.139 | 0.078 | 0.219, 0.204, 0.142, 0.000, 0.130 |
| brier | 0.212 | 0.009 | 0.205, 0.200, 0.210, 0.222, 0.223 |

Model-family selection used previously examined test splits; these scaffold scores remain exploratory. Standard deviations summarize seeds, not confidence intervals.

Validation-selected threshold means: mcc 0.230, sensitivity 0.797, specificity 0.399. Fixed-0.5 MCC is in the metric table.

## Performance — `human-dili-fda-new-unseen-scaffold@1`

Protocol `human-dili-fda-new-unseen-scaffold@1` — Frozen 1,096-compound DILIst-derived primary cohort, with a seed-101 15% scaffold calibration carve, and 119 FDA DILIrank 2.0 New drugs after primary identity and scaffold exclusion. Most/Less concern is positive, No concern negative. One saved deployment model; related-source transfer, not independent clinical adjudication. Split identities are pinned by hash.

Evaluated 2026-10-04 on n_train=932, n_val=164, n_test=119.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.617 | 0.000 | 0.617 |
| auprc | 0.352 | 0.000 | 0.352 |
| mcc | 0.149 | 0.000 | 0.149 |
| brier | 0.339 | 0.000 | 0.339 |

This fixed FDA New cohort contains 35 positive and 84 negative drugs after primary identity/scaffold exclusion. Its related label sources provide source-transfer evidence, not independently adjudicated clinical validation. The single-seed standard deviation of zero is not an uncertainty estimate.

At validation-selected threshold 0.692414: TP=22, FP=34, FN=13, TN=50; sensitivity 0.629, specificity 0.595, precision 0.393.

Raw-score Brier 0.272; calibrated Brier 0.339. Calibration was fitted using primary-source labels only. This diagnostic does not select a replacement pipeline.

> Comparable only with metrics carrying the same protocol id.

## Data

DILIst-derived human drug-level gold-standard labels — Seal et al. pinned gold-standard table, standardized to connectivity identities; conflicting labels excluded. Retrieved 2026-10-04, licensed source-specific; not redistributed, not redistributable.

`1096` compounds, positive rate `0.646`, table SHA-256 `9892ba3810d6122a…`

Regenerate and check for upstream drift with `python -m vp_human_dili.data fetch --verify`.

## Model

balanced-random-forest on `morgan2_count_2048` features. Shipped weights: 300 trees fit on 932 compounds; 164 compounds in a 15% scaffold carve calibrate scores and select the classification threshold only

`weights.joblib` SHA-256 `4af0a6f8c1303e3c…`

## Intended use

Research screening of molecular structures for association with human drug-level DILI labels. Inference requires SMILES only.

The score is not patient-specific risk, dose-response or proof of clinical safety. Model-family selection used previously examined scaffold test splits. The FDA validation cohort is an unseen-drug source-transfer check using related label sources, not prospective or independently adjudicated clinical validation. Calibration transfers poorly; low structure similarity and unresolved formulations limit coverage.

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0, sklearn 1.9.0, numpy 2.4.6, pandas 3.0.5, joblib 1.5.3, platform Windows, device cpu.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
