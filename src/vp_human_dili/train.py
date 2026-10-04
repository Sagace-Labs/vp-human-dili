"""Write a new forest artifact with an independent scaffold calibration carve."""

from __future__ import annotations

import argparse
import platform
from datetime import date

import joblib
import numpy as np
import pandas as pd
import rdkit
import sklearn
import xgboost

import vp_core
from vp_core import dataset, fingerprints, hashing, manifest, protocols
from vp_core.registry import parse_version
from vp_core.splits import murcko_scaffold, scaffold_train_val
from vp_human_dili import VERSIONS_DIR, contract, data, model

DEFAULT_PROTOCOL = "scaffold-balanced-5seed@1"
DEPLOYMENT_SEED = 101


def provenance():
    return {"python": platform.python_version(), "rdkit": rdkit.__version__,
            "xgboost": xgboost.__version__, "sklearn": sklearn.__version__,
            "numpy": np.__version__, "pandas": pd.__version__,
            "joblib": joblib.__version__, "platform": platform.system(), "device": "cpu"}


def build_version(version, *, reason, protocol=DEFAULT_PROTOCOL, supersedes=None, seed=DEPLOYMENT_SEED):
    parse_version(version)
    protocols.get(protocol)
    directory = VERSIONS_DIR / version
    if directory.exists():
        raise FileExistsError("Released versions are immutable; use a new version name")
    table = data.load()
    smiles = table.smiles.tolist()
    features = fingerprints.featurize(smiles, model.FEATURES)
    labels = table.label.to_numpy(dtype=int)
    train, valid = scaffold_train_val(smiles, val_frac=0.15, seed=seed)
    groups = np.array([murcko_scaffold(value) for value in smiles])
    bundle = model.fit(features[train], labels[train], features[valid], labels[valid], groups[train], seed=seed)
    directory.mkdir(parents=True)
    weights = directory / "weights.joblib"
    joblib.dump(bundle, weights, compress=3)
    record = {
        "schema": manifest.SCHEMA_VERSION, "pathway": "human-dili", "version": version,
        "released": date.today().isoformat(), "reason": reason, "supersedes": supersedes,
        "signature": contract.as_manifest_table(),
        "dataset": {"name": "DILIst-derived human drug-level gold-standard labels",
            "source": "Seal et al. pinned gold-standard table, standardized to connectivity identities; conflicting labels excluded",
            "url": data.SOURCE, "retrieved": "2026-10-04", "licence": "source-specific; not redistributed",
            "redistributable": False, "sha256": dataset.dataset_hash(table), "n_rows": len(table),
            "base_rate": float(table.label.mean()), "labels": ["label"],
            "source_sha256": data.SOURCE_SHA256, "fetch": "python -m vp_human_dili.data fetch --verify"},
        "model": {"family": "balanced-random-forest", "features": model.FEATURES,
            "fit": f"300 trees fit on {len(train)} compounds; {len(valid)} compounds in a 15% scaffold carve calibrate scores and select the classification threshold only",
            "weights": weights.name, "sha256": hashing.sha256_file(weights), "seed": seed,
            "n_estimators": model.N_ESTIMATORS, "class_weight": "balanced",
            "selected_min_samples_leaf": bundle["chosen_params"]["min_samples_leaf"],
            "selected_max_features": bundle["chosen_params"]["max_features"],
            "classification_threshold": bundle["threshold"], "calibration_fraction": 0.15},
        "protocol": {"id": protocol, "provider": "vp-core", "core_version": vp_core.__version__},
        "provenance": provenance(),
        "intended_use": contract.INTENDED_USE,
    }
    problems = manifest.validate(record, version_dir=directory)
    if problems:
        raise ValueError("; ".join(problems))
    manifest.write(directory / "manifest.toml", record)
    print(f"Wrote {directory}; training={len(train)}, calibration={len(valid)}; selected={bundle['chosen_params']}")
    return directory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--protocol", default=DEFAULT_PROTOCOL)
    parser.add_argument("--supersedes")
    parser.add_argument("--seed", type=int, default=DEPLOYMENT_SEED)
    args = parser.parse_args()
    build_version(args.version, reason=args.reason, protocol=args.protocol,
                  supersedes=args.supersedes, seed=args.seed)


if __name__ == "__main__":
    main()
