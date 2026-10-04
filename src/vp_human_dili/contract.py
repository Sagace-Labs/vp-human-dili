"""Drug-level label scores; the probability is not calibrated patient risk."""

SIGNATURE_VERSION = 1
INPUTS = ["smiles"]
PRIMARY = "human_dili"
INTENDED_USE = {
    "scope": "Research screening of molecular structures for association with human drug-level DILI labels. Inference requires SMILES only.",
    "limitations": "The score is not patient-specific risk, dose-response or proof of clinical safety. Model-family selection used previously examined scaffold test splits. The FDA validation cohort is an unseen-drug source-transfer check using related label sources, not prospective or independently adjudicated clinical validation. Calibration transfers poorly; low structure similarity and unresolved formulations limit coverage.",
}
OUTPUTS = [{
    "name": PRIMARY,
    "dtype": "float32",
    "range": [0.0, 1.0],
    "semantics": (
        "Structure-based probability of a positive drug-level DILIst-derived label, "
        "calibrated on held-out training-source drugs. It is not patient-specific "
        "DILI incidence or a clinical safety decision."
    ),
    "missing": "NaN when a nonempty molecular structure cannot be parsed and standardized",
}]


def column_names():
    return [output["name"] for output in OUTPUTS]


def as_manifest_table():
    return {"version": SIGNATURE_VERSION, "inputs": list(INPUTS),
            "outputs": [dict(output) for output in OUTPUTS]}
