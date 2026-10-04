"""Explicit source fetching and connectivity deduplication with label-conflict exclusion."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem

from vp_core import dataset
from vp_core.standardise import standardise_many

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
if not (ROOT / "pyproject.toml").exists():
    DATA = Path(__file__).resolve().parent / "data"
TABLE = DATA / "human_dili.parquet"
LABELS = ("label",)
COMMIT = "98b7a1b7fba4547cefe936b9890b44baab10e490"
SOURCE = f"https://raw.githubusercontent.com/srijitseal/DILI/{COMMIT}/data/DILI_Goldstandard_1111.csv"
SOURCE_SHA256 = "a86fd3a71c6832aac5e89bcdc35bdd27c000d6449260e4467a92982ae8e89be6"
EXPECTED_HASH = "9892ba3810d6122accf23a56f1ed5c5a8922f14017cc05f76004184bfdcbc33f"
EXTERNAL_PROTOCOL = "human-dili-fda-new-unseen-scaffold@1"


def prepare(raw):
    smiles, keys = standardise_many(raw.smiles_r.astype(str).tolist())
    table = pd.DataFrame({"inchikey": keys, "smiles": smiles, "label": raw.TOXICITY,
                          "source_row": raw.index, "source": raw.Source})
    table = table.dropna(subset=["inchikey", "smiles", "label"])
    table = table.loc[table.smiles.map(lambda value: bool(value) and Chem.MolFromSmiles(value) is not None)]
    table["connectivity"] = table.inchikey.str.split("-").str[0]
    counts = table.groupby("connectivity").label.nunique()
    table = table.loc[~table.connectivity.isin(counts[counts > 1].index)]
    table = table.drop_duplicates("connectivity").sort_values("inchikey").reset_index(drop=True)
    table.label = table.label.astype(int)
    problems = dataset.validate_table(table)
    if problems:
        raise ValueError("; ".join(problems))
    return table


def load():
    if not TABLE.exists():
        raise FileNotFoundError("Training table absent. Run: python -m vp_human_dili.data fetch --verify")
    table = dataset.read_table(TABLE)
    if dataset.dataset_hash(table) != EXPECTED_HASH:
        raise ValueError("Training table does not match the pinned standardized source")
    return table


def example():
    return dataset.read_table(DATA / "example/human_dili_example.parquet")


def external():
    from vp_core import hashing, protocols

    path = DATA / "validation/fda_new.parquet"
    table = dataset.read_table(path)
    ledger = json.loads((DATA / "validation/manifest.json").read_text(encoding="utf-8"))
    if dataset.dataset_hash(table) != ledger["dataset_sha256"]:
        raise ValueError("Validation dataset hash mismatch")
    if hashing.sha256_file(DATA / "validation/split_keys.csv") != protocols.get(EXTERNAL_PROTOCOL).split_key_sha256:
        raise ValueError("Validation split-key hash mismatch")
    return table


def fetch(*, verify=False):
    response = requests.get(SOURCE, timeout=60)
    response.raise_for_status()
    content = response.content
    if hashlib.sha256(content).hexdigest() != SOURCE_SHA256:
        raise ValueError("Pinned source snapshot hash mismatch")
    from io import BytesIO

    table = prepare(pd.read_csv(BytesIO(content)))
    actual = dataset.dataset_hash(table)
    if verify and actual != EXPECTED_HASH:
        raise ValueError(f"Standardized table hash mismatch: {actual}")
    source_dir = DATA / "sources"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "gold.csv").write_bytes(content)
    dataset.write_table(table, TABLE)
    (source_dir / "retrieval.json").write_text(json.dumps({
        "url": SOURCE, "source_sha256": SOURCE_SHA256,
        "dataset_sha256": actual, "retrieved": datetime.now(UTC).isoformat(),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared {len(table)} compounds; standardized SHA-256 {actual}")
    return table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["fetch"])
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    fetch(verify=args.verify)


if __name__ == "__main__":
    main()
