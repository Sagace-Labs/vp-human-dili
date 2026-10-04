"""Refit every primary seed and test frozen weights on the separate FDA cohort."""

from __future__ import annotations

import argparse
from datetime import date

import numpy as np
import pandas as pd

from vp_core import card, dataset, fingerprints, forest, metrics, metrics_store, protocols
from vp_core.splits import murcko_scaffold, scaffold_train_val
from vp_human_dili import current_version, data, get, model
from vp_human_dili.train import provenance


def evaluate_version(version, *, protocol_id=None, use_example=False, write=True):
    resolved = get(version)
    protocol = protocols.get(protocol_id or resolved.protocol_id)
    protocols.require_comparable(protocol.id, protocol_id or resolved.protocol_id)
    if use_example and write:
        raise ValueError("Fixture scores cannot be recorded as release metrics")
    table = data.example() if use_example else data.load()
    if not use_example and dataset.dataset_hash(table) != resolved.manifest["dataset"]["sha256"]:
        raise ValueError("Release training dataset hash mismatch")
    bundle = resolved.load_model()
    smiles = table.smiles.tolist()
    features = fingerprints.featurize(smiles, resolved.features)
    labels = table.label.to_numpy(dtype=int)
    groups = np.array([murcko_scaffold(value) for value in smiles])
    folds, calibrated, raw_scores, diagnostics = [], [], [], []
    notes = []
    external_predictions = None
    if protocol.id == data.EXTERNAL_PROTOCOL:
        if use_example:
            raise ValueError("FDA validation requires the pinned primary cohort")
        external = data.external()
        if set(external.connectivity) & set(table.inchikey.str.split("-").str[0]):
            raise ValueError("External identity overlap")
        if set(external.smiles.map(murcko_scaffold)) & set(groups):
            raise ValueError("External scaffold overlap")
        split = pd.read_csv(data.DATA / "validation/split_keys.csv")
        expected = set(zip(split.cohort, split.inchikey, strict=True))
        actual = {("primary", key) for key in table.inchikey} | {("external", key) for key in external.inchikey}
        if expected != actual:
            raise ValueError("Frozen split identities differ from evaluation tables")
        if bundle["seed"] != 101:
            raise ValueError("The registered FDA protocol fixes the deployment seed at 101")
        train, valid = scaffold_train_val(smiles, val_frac=0.15, seed=101)
        for role, indices in (("train", train), ("calibration", valid)):
            if set(split.loc[split.role == role, "inchikey"]) != set(table.inchikey.iloc[indices]):
                raise ValueError("Frozen calibration carve changed")
        xe = fingerprints.featurize(external.smiles.tolist(), resolved.features)
        y = external.label.to_numpy(dtype=int)
        probability = forest.predict_proba(bundle, xe)
        raw = forest.predict_proba(bundle, xe, raw=True)
        test_table = external
        scored_bundle = bundle
        units = [(101, train, valid, np.arange(len(external)), y, probability, raw, scored_bundle)]
        notes.append("This fixed FDA New cohort contains 35 positive and 84 negative drugs after primary identity/scaffold exclusion. Its related label sources provide source-transfer evidence, not independently adjudicated clinical validation. The single-seed standard deviation of zero is not an uncertainty estimate.")
        external_predictions = pd.DataFrame({"inchikey": external.inchikey, "label": y,
            "probability": probability, "raw_probability": raw,
            "threshold": bundle["threshold"]})
    else:
        units = []
        test_table = table
        for seed in protocol.seeds:
            train, valid, test = protocol.split_indices(smiles, seed)
            if any(len(np.unique(labels[ix])) != 2 for ix in (train, valid, test)):
                raise ValueError(f"Unsupported class counts for seed {seed}")
            scored_bundle = model.fit(features[train], labels[train], features[valid], labels[valid],
                groups[train], seed=seed, recipe=bundle["recipe"])
            probability = forest.predict_proba(scored_bundle, features[test])
            raw = forest.predict_proba(scored_bundle, features[test], raw=True)
            units.append((seed, train, valid, test, labels[test], probability, raw, scored_bundle))
        notes.append("Model-family selection used previously examined test splits; these scaffold scores remain exploratory. Standard deviations summarize seeds, not confidence intervals.")
    for seed, train, valid, test, y, probability, raw, fitted in units:
        if len(np.unique(y)) != 2:
            raise ValueError("Unsupported single-class test fold")
        calibrated.append(metrics.binary_metrics(y, probability))
        raw_scores.append(metrics.binary_metrics(y, raw))
        folds.append({"seed": seed, "train": len(train), "val": len(valid), "test": len(test),
                      "train_positive": int(labels[train].sum()), "val_positive": int(labels[valid].sum()),
                      "test_positive": int(y.sum())})
        diagnostics.append({"seed": seed, "selection": fitted["selection"],
            "chosen_params": fitted["chosen_params"],
            "operating": forest.operating(y, probability, fitted["threshold"]),
            "validation_operating": fitted["validation_operating"]})
        print(f"Seed {seed}: AUC={calibrated[-1]['auc_roc']:.6f}, Brier={calibrated[-1]['brier']:.6f}", flush=True)
    operating = diagnostics[0]["operating"] if len(diagnostics) == 1 else None
    if operating:
        notes.append(f"At validation-selected threshold {operating['threshold']:.6f}: TP={operating['tp']}, FP={operating['fp']}, FN={operating['fn']}, TN={operating['tn']}; sensitivity {operating['sensitivity']:.3f}, specificity {operating['specificity']:.3f}, precision {operating['precision']:.3f}.")
        notes.append(f"Raw-score Brier {raw_scores[0]['brier']:.3f}; calibrated Brier {calibrated[0]['brier']:.3f}. Calibration was fitted using primary-source labels only. This diagnostic does not select a replacement pipeline.")
    else:
        notes.append("Validation-selected threshold means: " + ", ".join(
            f"{key} {np.mean([item['operating'][key] for item in diagnostics]):.3f}"
            for key in ("mcc", "sensitivity", "specificity")) + ". Fixed-0.5 MCC is in the metric table.")
    record = {"pathway": "human-dili", "version": version, "protocol_id": protocol.id,
        "dataset_sha256": dataset.dataset_hash(test_table), "dataset": "example fixture" if use_example else "full",
        "evaluated": date.today().isoformat(),
        "n": {"total": len(table), "train": folds[0]["train"], "val": folds[0]["val"], "test": folds[0]["test"]},
        "folds": folds, "test": metrics.aggregate(calibrated, protocol.metrics), "notes": notes,
        "diagnostics": {"raw_test": metrics.aggregate(raw_scores, protocol.metrics),
                        "folds": diagnostics, "environment": provenance()}}
    if any(item["n_finite"] != len(protocol.seeds) for item in record["test"].values()):
        raise ValueError("A missing seed would change the estimator")
    if write:
        metrics_store.write(resolved.directory, record)
        if external_predictions is not None:
            external_predictions.to_csv(resolved.directory / "fda_predictions.csv", index=False)
        card.write_card(resolved.directory)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default=None)
    parser.add_argument("--protocol", default=None)
    parser.add_argument("--example", action="store_true")
    args = parser.parse_args()
    evaluate_version(args.version or current_version(), protocol_id=args.protocol,
                     use_example=args.example, write=not args.example)


if __name__ == "__main__":
    main()
