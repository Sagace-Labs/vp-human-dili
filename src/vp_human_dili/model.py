"""Counted Morgan forest recipe with training-only scaffold model selection."""

from vp_core import forest

FEATURES = "morgan2_count_2048"
FEATURES_BY_OUTPUT = {}
GRID = [{"min_samples_leaf": 2, "max_features": "sqrt"},
        {"min_samples_leaf": 5, "max_features": 0.25}]
N_ESTIMATORS = 300


def features_for(output):
    return FEATURES_BY_OUTPUT.get(output, FEATURES)


def fit(xt, yt, xv, yv, groups, *, seed, recipe=None):
    recipe = recipe or {"grid": GRID, "n_estimators": N_ESTIMATORS, "n_jobs": 6}
    return forest.fit(xt, yt, xv, yv, groups, seed=seed, **recipe)
