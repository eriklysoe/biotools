"""Rarefaction curve calculation."""

import numpy as np


def calc_rarefaction_curves(abundance_df, n_steps=10, n_reps=10, seed=42):
    """
    Generate rarefaction curves for all samples.

    At each depth step, subsample and count observed features (mean of n_reps).

    Returns:
        dict with 'depths' (list) and 'curves' (dict of sample -> list of values)
    """
    rng = np.random.RandomState(seed)

    sample_depths = abundance_df.sum(axis=0)
    max_depth = int(sample_depths.max())

    # Create depth steps
    depths = [int(d) for d in np.linspace(0, max_depth, n_steps + 1)[1:]]
    # Ensure we include some lower depths
    if depths[0] > 100 and max_depth > 100:
        depths = [int(d) for d in np.linspace(10, max_depth, n_steps)]

    curves = {}

    for sample in abundance_df.columns:
        counts = abundance_df[sample].values.astype(int)
        total = counts.sum()
        sample_curve = []

        # Build pool once
        pool = np.repeat(np.arange(len(counts)), counts)

        for depth in depths:
            if depth > total:
                # Use full sample observed features
                sample_curve.append(int(np.sum(counts > 0)))
            else:
                obs_list = []
                for _ in range(n_reps):
                    chosen = rng.choice(pool, size=depth, replace=False)
                    obs = len(np.unique(chosen))
                    obs_list.append(obs)
                sample_curve.append(round(float(np.mean(obs_list)), 1))

        curves[str(sample)] = sample_curve

    return {"depths": depths, "curves": curves}
