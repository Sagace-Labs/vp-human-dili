"""Reconstruct the frozen FDA cohort offline from official and PubChem assertions."""

from __future__ import annotations

import hashlib
import json
import re

import pandas as pd
from rdkit import Chem

from vp_core import dataset, hashing
from vp_core.splits import murcko_scaffold
from vp_core.standardise import standardise_many
from vp_human_dili import data


def rebuild(primary):
    root = data.DATA / "validation"
    metadata = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in metadata["source_hashes"].items():
        if hashing.sha256_file(root / name) != expected:
            raise ValueError(f"Source hash mismatch: {name}")
    rows = {}
    for line in (root / "sources/fda_web_snapshot.txt").read_text(encoding="utf-8").splitlines():
        match = re.match(r"L\d+: (LT\d{5}.*)$", line)
        if match:
            cells = re.split(r"\s*\|\s*", match[1])
            if len(cells) != 6:
                raise ValueError("Changed FDA table shape")
            rows[cells[0]] = cells
    table = pd.DataFrame(rows.values(), columns=["ltkb_id", "name", "severity", "label_section", "concern", "change"])
    table["concern"] = table.concern.str.replace("^{v}", "", regex=False).str.lower()
    if table.concern.value_counts().to_dict() != {"no-dili-concern": 414,
        "ambiguous-dili-concern": 354, "less-dili-concern": 351, "most-dili-concern": 217}:
        raise ValueError("FDA category counts changed")
    selected = table.loc[(table.change == "New") & ~table.concern.str.startswith("ambiguous")]
    assertions = json.loads((root / "sources/pubchem_assertions.json").read_text(encoding="utf-8"))
    records = []
    for row in selected.to_dict(orient="records"):
        assertion = assertions[row["ltkb_id"]]
        source = json.loads(assertion)
        properties = source["body"].get("PropertyTable", {}).get("Properties", [])
        if not properties:
            continue
        smiles, keys = standardise_many([p.get("SMILES", p.get("CanonicalSMILES")) for p in properties])
        valid = [(s, k, p["CID"]) for s, k, p in zip(smiles, keys, properties, strict=True)
                 if s and k and Chem.MolFromSmiles(s) is not None]
        if len({key.split("-")[0] for _, key, _ in valid}) != 1:
            continue
        smiles, key, cid = valid[0]
        records.append({**row, "smiles": smiles, "inchikey": key, "cid": cid,
                        "assertion_sha256": hashlib.sha256(assertion.encode("utf-8")).hexdigest()})
    result = pd.DataFrame(records)
    result["label"] = result.concern.map({"most-dili-concern": 1, "less-dili-concern": 1, "no-dili-concern": 0})
    result["connectivity"] = result.inchikey.str.split("-").str[0]
    conflicts = result.groupby("connectivity").label.nunique()
    result = result.loc[~result.connectivity.isin(conflicts[conflicts > 1].index)]
    result = result.drop_duplicates("connectivity").sort_values("inchikey").reset_index(drop=True)
    keys = set(primary.inchikey.str.split("-").str[0])
    scaffolds = set(primary.smiles.map(murcko_scaffold))
    result = result.loc[~result.connectivity.isin(keys) & ~result.smiles.map(murcko_scaffold).isin(scaffolds)].reset_index(drop=True)
    if dataset.dataset_hash(result) != metadata["dataset_sha256"]:
        raise ValueError("Rebuilt FDA cohort differs from the frozen standardized hash")
    retained = pd.read_csv(root / "resolutions.csv")
    for row in retained.itertuples():
        actual = hashlib.sha256(assertions[row.ltkb_id].encode("utf-8")).hexdigest()
        if actual != row.assertion_sha256:
            raise ValueError(f"PubChem assertion mismatch: {row.ltkb_id}")
    return result


def main():
    table = rebuild(data.load())
    print(f"Verified offline FDA cohort: {len(table)} compounds, {int(table.label.sum())} positives")


if __name__ == "__main__":
    main()
