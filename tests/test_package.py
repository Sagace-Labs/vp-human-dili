"""Offline contract, source preparation and retained validation safeguards."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import vp_human_dili
from vp_human_dili import contract, data

from vp_core import dataset, hashing, manifest, protocols
from vp_core.docs import check_pathway_readme


def test_released_artifact_and_frame():
    version = vp_human_dili.get("v1")
    assert manifest.validate(version.manifest, version_dir=version.directory) == []
    assert version.manifest["intended_use"] == contract.INTENDED_USE
    frame = vp_human_dili.predict(["CCO", "OCC", "not a molecule", "", None])
    assert list(frame.columns) == contract.column_names()
    assert frame.dtypes.iloc[0] == np.dtype("float32")
    assert len(frame) == 5
    assert frame.iloc[0, 0] == frame.iloc[1, 0]
    assert frame.iloc[:2, 0].between(0, 1).all()
    assert frame.iloc[2:, 0].isna().all()
    assert vp_human_dili.predict([]).empty


def test_dataset_load_never_fetches(monkeypatch, tmp_path):
    monkeypatch.setattr(data, "TABLE", tmp_path / "missing.parquet")
    monkeypatch.setattr(data.requests, "get", lambda *_a, **_k: pytest.fail("implicit network"))
    with pytest.raises(FileNotFoundError, match="fetch --verify"):
        data.load()


def test_connectivity_conflicts_are_excluded():
    raw = pd.DataFrame({"smiles_r": ["CCO", "OCC", "CCN", "NCC", "CCCl"],
                        "TOXICITY": [0, 1, 0, 0, 1], "Source": ["synthetic"] * 5})
    table = data.prepare(raw)
    assert len(table) == 2
    assert sorted(table.label) == [0, 1]
    assert table.connectivity.is_unique
    assert "CCO" not in table.smiles.tolist()


def test_synthetic_fixture_and_readme():
    assert dataset.validate_table(data.example()) == []
    readme = Path(__file__).resolve().parents[1] / "README.md"
    assert check_pathway_readme(readme.read_text(encoding="utf-8")) == []


def test_frozen_validation_sources_and_protocol():
    external = data.external()
    assert len(external) == 119 and external.label.sum() == 35
    root = data.DATA / "validation"
    ledger = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    assert dataset.dataset_hash(external) == ledger["dataset_sha256"]
    for name, expected in ledger["source_hashes"].items():
        assert hashing.sha256_file(root / name) == expected
    assertions = json.loads((root / "sources/pubchem_assertions.json").read_text(encoding="utf-8"))
    for row in pd.read_csv(root / "resolutions.csv").itertuples():
        assert hashlib.sha256(assertions[row.ltkb_id].encode()).hexdigest() == row.assertion_sha256
    with pytest.raises(ValueError, match="frozen source-transfer"):
        protocols.get(data.EXTERNAL_PROTOCOL).split_indices(external.smiles.tolist(), 101)
    with pytest.raises(ValueError, match="not comparable"):
        protocols.require_comparable(data.EXTERNAL_PROTOCOL, "scaffold-balanced-5seed@1")


def test_train_refuses_overwriting():
    from vp_human_dili.train import build_version

    with pytest.raises(FileExistsError, match="immutable"):
        build_version("v1", reason="must not replace a release")
