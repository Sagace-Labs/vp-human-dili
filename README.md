# vp-human-dili

Predicts association with human drug-level DILI labels from a SMILES string
using a counted-Morgan random forest. Intended for research screening.

## Install

    pip install git+https://github.com/Sagace-Labs/vp-core.git@v1.8.0
    pip install git+https://github.com/Sagace-Labs/vp-human-dili.git@v1

## Use

    from vp_human_dili import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])  # -> DataFrame[human_dili]

Returns a training-source-calibrated positive-label probability, not
patient-specific incidence or a clinical safety decision. Invalid or empty
structures return NaN. Pin with `predict(smiles, version="v1")`; list releases
with `versions()`. Dose, exposure and assay measurements are not inputs.

## Current version

**v1**, signature 1, evaluated under `scaffold-balanced-5seed@1` and a separate
frozen FDA New drug protocol. See
[`src/vp_human_dili/versions/v1/CARD.md`](src/vp_human_dili/versions/v1/CARD.md)
for metrics, operating points, calibration limits and provenance.

## Data

The pinned Seal et al. DILIst-derived table yields 1,096 standardized compounds.
Training labels are fetched explicitly and are not redistributed. A separate
119-drug FDA DILIrank 2.0 validation cohort excludes training identities and
scaffolds; its related source labels are not independent clinical adjudication.
See [`data/README.md`](data/README.md) for sources, processing and coverage.

## Retrain

    python -m vp_human_dili.data fetch --verify
    python -m vp_human_dili.train --version v2 --reason "why this version exists"
    python -m vp_human_dili.evaluate --version v2
    python -m vp_human_dili.evaluate --version v2 --protocol human-dili-fda-new-unseen-scaffold@1

Training selects the forest on training scaffolds and reserves 15% of compounds
for calibration. Evaluation refits every primary seed; the FDA protocol scores
the saved deployment model. Further tuning after seeing FDA results requires
fresh validation. Recorded hashes and environments pin reproducibility.

## Licence

Code and weights are Apache-2.0 ([`LICENSE`](LICENSE)). Source-specific data
terms and the synthetic fixture are described in [`LICENSE-DATA`](LICENSE-DATA).

## Cite

See [`CITATION.cff`](CITATION.cff).
