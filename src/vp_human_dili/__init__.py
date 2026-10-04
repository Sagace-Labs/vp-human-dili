"""Structure-based human drug-level DILI label predictions."""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from vp_core import fingerprints, forest
from vp_core.registry import Version, VersionedPathway
from vp_core.standardise import standardise_many

__version__ = "1.0.0"
PATHWAY = "human-dili"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"


def _predict_values(bundle, smiles, version):
    RDLogger.DisableLog("rdApp.*")
    values = np.full((len(smiles), len(version.output_names)), np.nan, dtype=np.float32)
    candidates = [i for i, value in enumerate(smiles)
                  if isinstance(value, str) and value.strip() and Chem.MolFromSmiles(value) is not None]
    standardized, _ = standardise_many([smiles[i] for i in candidates])
    valid = [(i, value) for i, value in zip(candidates, standardized, strict=True)
             if value and Chem.MolFromSmiles(value) is not None]
    if valid:
        indices, structures = zip(*valid, strict=True)
        features = fingerprints.featurize(list(structures), version.features)
        values[list(indices), 0] = forest.predict_proba(bundle, features)
    return values


_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR, predict_fn=_predict_values)
__all__ = ["PATHWAY", "VERSIONS_DIR", "Version", "__version__", "current_version",
           "get", "predict", "signature", "versions"]


def predict(smiles: list[str], *, version: str | None = None) -> pd.DataFrame:
    return _pathway.predict(smiles, version=version)


def versions():
    return _pathway.versions()


def current_version():
    return _pathway.current()


def get(version=None):
    return _pathway.get(version)


def signature(version=None):
    return get(version).manifest["signature"]
